from random import shuffle
from argparse import ArgumentParser
from os.path import join, exists
from os import mkdir
import contextlib
import time
import numpy as np

from myopic_mces.filter_MCES import ComputationMode

def _load_mces_matrix(dmatrix_file):
    """Load a precomputed MCES distance matrix HDF5 lookup once, for reuse across batches."""
    from scipy.spatial.distance import squareform
    import h5py
    with h5py.File(dmatrix_file, 'r') as hf:
        mces = hf['mces'][:]
        all_smiles = [s.decode() for s in hf['mces_smiles_order'][:]] # Since here the input is csv we need to decode the lookup
    mces = squareform(mces)
    smiles_index = {smiles: i for i, smiles in enumerate(all_smiles)}
    return mces, smiles_index


def filter_inputs(smiles_list, index, dmatrix_file, threshold=None, _matrix=None, _smiles_index=None):
    print('filtering for precomputed mces')
    filtered_inputs = []
    precomputed_mces = []
    precomputed_index = []
    if _matrix is None or _smiles_index is None:
        mces, smiles_index = _load_mces_matrix(dmatrix_file)
    else:
        mces, smiles_index = _matrix, _smiles_index
    for i, s1, s2 in index:
        idx1 = smiles_index.get(smiles_list[s1])
        idx2 = smiles_index.get(smiles_list[s2])
        
        if idx1 is not None and idx2 is not None:
            val = mces[idx1][idx2]
            if val != -1:
                if threshold is None or val < threshold:
                    precomputed_mces.append((i, val,-1,-1))
                    precomputed_index.append([i, s1, s2])
                    continue
        
        filtered_inputs.append([i, s1, s2])
    return np.array(filtered_inputs),np.array(precomputed_index) ,precomputed_mces


def _build_inchikey_cache(smiles_list):
    """Pre-convert unique SMILES to first block inchikeys, once, for reuse across batches."""
    from rdkit import Chem
    smiles_to_inchikey = {}
    for s in smiles_list:
        if isinstance(s, bytes):
            s = s.decode()
        if s in smiles_to_inchikey:
            continue
        mol = Chem.MolFromSmiles(s)
        smiles_to_inchikey[s] = Chem.MolToInchiKey(mol)[:14] if mol is not None else None
    return smiles_to_inchikey


def db_filter_inputs(smiles_list, index, client, threshold, always_stronger_bound, batch_size=50000, smiles_to_inchikey=None):
    """Filter precomputed MCES distances from the MCES database (mcesdb).

    Analogous to filter_inputs but queries a mcesdb PostgreSQL backend instead of
    a local HDF5 lookup matrix.

    An exact distance strictly below the threshold is reused
    directly, and a stored lower bound strictly above the threshold is reused as a
    bound. The choice between mces_match2_lower (always_stronger_bound=True, which
    corresponds to filter2 / STRONGEST_BOUND) and best_lower_bound
    (always_stronger_bound=False / dynamic mode, UNKNOWN) mirrors which bound
    apply_filter would actually have computed.

    Pairs in the DB are stored symmetrically in both orientations using first block inchikeys.

    TODO: we will have to think up something when we get multiple billion DB hits
    """
    print('filtering for precomputed mces (DB)')
    t0 = time.time()

    if smiles_to_inchikey is None:
        smiles_to_inchikey = _build_inchikey_cache(smiles_list)

    filtered_inputs = []
    precomputed_mces = []
    precomputed_index = []
    # rows whose InChIKey pair we can identify, deferred for a single batched DB query
    pending = []   # entries: (i, s1, s2, ik_a, ik_b) - canonical (min,max) InChIKey pair

    for i, s1, s2 in index:
        smi1 = smiles_list[s1]
        smi2 = smiles_list[s2]
        if isinstance(smi1, bytes):
            smi1 = smi1.decode()
        if isinstance(smi2, bytes):
            smi2 = smi2.decode()
        ik1 = smiles_to_inchikey.get(smi1)
        ik2 = smiles_to_inchikey.get(smi2)
        if ik1 is None or ik2 is None:
            # can't identify the structure - has to be (re)computed
            filtered_inputs.append([i, s1, s2])
            continue
        # canonical ordering for symmetric-row lookup
        aik, bik = (ik1, ik2) if ik1 <= ik2 else (ik2, ik1)
        pending.append((i, s1, s2, aik, bik))

    # batch the DB query so a single staging-table round-trip stays bounded
    bound_map = {}
    for chunk_start in range(0, len(pending), batch_size):
        chunk = pending[chunk_start: chunk_start + batch_size]
        resp = client.get_pair_bounds([(r[3], r[4]) for r in chunk])
        for r in resp:
            # get_pair_bounds returns the canonical (min, max) pair back; key by it
            aik = r['inchikey_a']
            bik = r['inchikey_b']
            bound_map[(aik, bik)] = r

    for i, s1, s2, aik, bik in pending:
        row = bound_map.get((aik, bik))
        if row is None:
            # pair not in the DB at all - recompute
            filtered_inputs.append([i, s1, s2])
            continue
        exact = row.get('mces_exact')
        best_lower = row.get('best_lower_bound')
        match2_lower = row.get('mces_match2_lower')
        # 1) exact distance already proven and below threshold -> reuse
        #    (uses strict < to stay consistent with filter_inputs' matrix path)
        if exact is not None and exact < threshold:
            precomputed_mces.append((i, float(exact), -1, ComputationMode.EXACT.value))
            precomputed_index.append([i, s1, s2])
            continue
        # 2) stronger bound short-circuit (always_stronger_bound=True corresponds
        #    to apply_filter always computing filter2 -> STRONGEST_BOUND)
        if always_stronger_bound and match2_lower is not None and match2_lower > threshold:
            precomputed_mces.append((i, float(match2_lower), -1, ComputationMode.STRONGEST_BOUND.value))
            precomputed_index.append([i, s1, s2])
            continue
        # 3) dynamic bound short-circuit (always_stronger_bound=False). Reuse the
        #    stored best_lower_bound; the runtime filter would have picked the same
        #    bound dynamically, so we mark the mode UNKNOWN (the specific filter
        #    that produced this bound is not recoverable from best_lower_bound)
        if (not always_stronger_bound) and best_lower is not None and best_lower > threshold:
            precomputed_mces.append((i, float(best_lower), -1, ComputationMode.UNKNOWN.value))
            precomputed_index.append([i, s1, s2])
            continue
        # 4) no usable bound known - recompute
        filtered_inputs.append([i, s1, s2])

    print(f'done, took {(time.time() - t0) / 60:.1f}min and found {len(precomputed_mces)} in DB, '
          f'{len(filtered_inputs)} to go')
    return np.array(filtered_inputs), np.array(precomputed_index), precomputed_mces


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('input_file', help='readily formatedd myopic_mces input OR list of smiles for `h5`-mode')
    parser.add_argument('out_folder')
    parser.add_argument('--hdf5_mode', action='store_true')
    parser.add_argument('--hdf5_extra_input_file', help='in HDF5 mode, with this argument an extra input file (see above) can be provided: '
                        'pairs from `input_file`-SMILES and `hdf5_extra_input_file`-SMILES will then be generated', default=None)
    parser.add_argument('--batch_size', type=int, default=224960 * 100)
    parser.add_argument('--batch_start', type=int, default=0)
    parser.add_argument('--no_shuffle', action='store_true')
    parser.add_argument('--use_matrix_lookup', help='(experimental) Use with the '
                        'path to a HDF5 file with precomputed MCES distances. Computation for these instances will be '
                        'skipped, using the provided values. HDF5 has to contain distances (key `mces`) and SMILES '
                        '(`mces_smiles_order`), like the HDF5 files produced by this script. ')
    parser.add_argument('--lookup_threshold', help='(experimental) Use with `--use_matrix_lookup`: '
                        'Precomputed values equal or greater than the threshold will be ignored; these '
                        'instances will be recomputed', default=None, type=float)
    parser.add_argument('--use_db_lookup', help='(experimental) Use the mcesdb backend as a '
                        'precomputed-MCES lookup. Connection parameters come from the standard PostgreSQL '
                        'environment variables. Pairwise matches are done by first block inchikey.',
                        action='store_true')
    parser.add_argument('--threshold', help='threshold for the distance, used by --use_db_lookup to decide '
                        'which stored bounds are reusable. Should match the --threshold of the upcoming '
                        'compute step.', default=10., type=float)
    parser.add_argument('--choose_bound_dynamically', action='store_true',
                        help='if this is set (in combination with `--use_db_lookup`), a stored best_lower_bound '
                        'above the threshold is reused as the distance. Otherwise (default), only a stored '
                        'mces_match2_lower (the strongest bound, corresponding to filter2) above the threshold '
                        'is reused. Should match the flag of the same name on the compute side.')
    args = parser.parse_args()

    if (args.hdf5_mode):
        import h5py
        from scipy.special import binom
        print('hdf5-mode...')
        # TWO input files (e.g., dataset x bio)
        if (args.hdf5_extra_input_file is not None):
            smiles_input1 = [l.strip() for l in open(args.input_file).readlines()]
            smiles_input2 = [l.strip() for l in open(args.hdf5_extra_input_file).readlines()]
            smiles_input = smiles_input1 + smiles_input2
            n1, n2 = len(smiles_input1), len(smiles_input2)
            ninstances = n1 * n2
            if ninstances >= np.iinfo(np.int64).max:
                # can probably not happen anyways?
                raise Exception('too many instances:', ninstances)
            nbatches = int(np.ceil(ninstances/args.batch_size))
            print(f'read {n1:_}*{n2:_} SMILES as input -> {ninstances:_} instances '
                  f'-> {nbatches} batches of {args.batch_size:_}')

            def unrank(k):
                return k // n2, n1 + (k % n2)
        # one input file
        else:
            smiles_input = [l.strip() for l in open(args.input_file).readlines()]
            n = len(smiles_input)
            ninstances = int(binom(n, 2))
            if ninstances >= np.iinfo(np.int64).max:
                # can probably not happen anyways?
                raise Exception('too many instances:', ninstances)
            nbatches = int(np.ceil(ninstances/args.batch_size))
            print(f'read {n:_} SMILES as input -> {ninstances:_} instances -> {nbatches} batches of {args.batch_size:_}')
            # cumulative pair-count offset of each row in triu(n, k=1) order, used to
            # unrank a linear pair ordinal into (row, col) without materializing all pairs
            row_idx = np.arange(n, dtype=np.int64)
            row_starts = row_idx*(n-1) - row_idx*(row_idx-1)//2

            def unrank(k):
                i = np.searchsorted(row_starts, k, side='right') - 1
                j = i + 1 + (k - row_starts[i])
                return i, j

        if not args.no_shuffle:
            print('NOTE: --hdf5_mode generates batches in a memory-bounded streaming '
                  'fashion and does not support shuffling; batches are contiguous '
                  'ranges of the pair ordering regardless of --no_shuffle.')

        # shared per-run lookup state, built once and reused for every batch below
        inchikey_cache = None
        matrix_cache = matrix_smiles_index_cache = None
        db_ctx = contextlib.nullcontext()
        if args.use_db_lookup:
            from mcesdb import McesDbClient
            inchikey_cache = _build_inchikey_cache(smiles_input)
            # empty conninfo -> libpq falls back to PGHOST/PGDATABASE/... env vars
            db_ctx = McesDbClient("")
        elif args.use_matrix_lookup:
            matrix_cache, matrix_smiles_index_cache = _load_mces_matrix(args.use_matrix_lookup)

        if (not exists(args.out_folder)):
            mkdir(args.out_folder)
        print('creating batches')
        # when any precomputed pairs turn up, they take slot 0 (batch0_precomputed),
        # and the to-compute batches start at 1 so they sort after it
        i_offset = 1 if (args.use_matrix_lookup or args.use_db_lookup) else 0
        precomputed_index_all = []
        precomputed_mces_all = []
        written_batches = 0

        with db_ctx as client:
            for b in range(nbatches):
                k = np.arange(b*args.batch_size, min((b+1)*args.batch_size, ninstances), dtype=np.int64)
                i_arr, j_arr = unrank(k)
                batch_index = np.stack([k, i_arr, j_arr], axis=1)

                if args.use_db_lookup:
                    batch_index, pre_idx, pre_mces = db_filter_inputs(
                        smiles_list=smiles_input, index=batch_index, client=client,
                        threshold=args.threshold, always_stronger_bound=not args.choose_bound_dynamically,
                        smiles_to_inchikey=inchikey_cache)
                    precomputed_index_all.extend(pre_idx.tolist())
                    precomputed_mces_all.extend(pre_mces)
                elif args.use_matrix_lookup:
                    batch_index, pre_idx, pre_mces = filter_inputs(
                        smiles_list=smiles_input, index=batch_index, dmatrix_file=args.use_matrix_lookup,
                        threshold=args.lookup_threshold, _matrix=matrix_cache, _smiles_index=matrix_smiles_index_cache)
                    precomputed_index_all.extend(pre_idx.tolist())
                    precomputed_mces_all.extend(pre_mces)

                if batch_index.shape[0] > 0:
                    file_path = join(args.out_folder, f'batch{written_batches + i_offset}.hdf5')
                    print(f'creating batch {written_batches + i_offset} at {file_path} '
                          f'({batch_index.shape[0]:_} instances)...', end=' ')
                    with h5py.File(file_path, 'w') as f:
                        f.create_dataset('computation_indices', data=batch_index, dtype='int64', compression='gzip')
                        f.create_dataset('smiles', data=smiles_input)
                        f.attrs['original_path'] = file_path
                    print('done')
                    written_batches += 1

        if args.use_matrix_lookup or args.use_db_lookup:
            print(f'Found {len(precomputed_index_all):_} precomputed MCES -> '
                  f'{ninstances - len(precomputed_index_all):_} instances left to compute')
            if len(precomputed_index_all) > 0:
                file_path = join(args.out_folder, 'batch0_precomputed.hdf5')
                print(f'creating batch0_precomputed at {file_path} with precomputed mces...', end=' ')
                with h5py.File(file_path, 'w') as f:
                    f.create_dataset('computation_indices', data=np.array(precomputed_index_all), dtype='int64', compression='gzip')
                    f.create_dataset('mces', data=[row[1] for row in precomputed_mces_all], compression='gzip')
                    f.create_dataset('smiles', data=smiles_input)
                    f.create_dataset('computation_times', data=[row[2] if len(row) > 2 else -1 for row in precomputed_mces_all], compression='gzip')
                    f.create_dataset('computation_modes', data=[row[3] if len(row) > 3 else -1 for row in precomputed_mces_all], compression='gzip')
                print('done')

        if written_batches == 0:
            print("Nothing to compute, everything cached")
    else:
        pairs = [line.strip() for line in open(args.input_file).readlines() if not line.startswith('i,smiles')]
        print('pairs read')
        if (not args.no_shuffle):
            shuffle(pairs)
        # make folders, if necessary
        for folder in [args.out_folder, join(args.out_folder, 'data'), join(args.out_folder, 'out')]:
            if (not exists(folder)):
                mkdir(folder)
                print('created', folder)
        pairs_counter = 0
        batch_counter = 0
        while (pairs_counter < len(pairs)):
            batch_pairs = pairs[pairs_counter:(pairs_counter+args.batch_size)]
            csv = join(args.out_folder, 'data', f'batch{args.batch_start + batch_counter}.csv')
            with open(csv, 'w') as out:
                out.write('\n'.join(batch_pairs) + '\n')
            print('wrote', csv)
            pairs_counter += len(batch_pairs)
            batch_counter += 1
        assert pairs_counter == len(pairs)
