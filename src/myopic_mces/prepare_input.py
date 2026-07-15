from random import shuffle
from argparse import ArgumentParser
from os.path import join, exists
from os import mkdir
import time
import numpy as np

from myopic_mces.filter_MCES import ComputationMode

def filter_inputs(smiles_list,index,dmatrix_file,threshold=None):
    from scipy.spatial.distance import squareform
    import h5py
    print('filtering for precomputed mces')
    filtered_inputs = []
    precomputed_mces = []
    precomputed_index = []
    with h5py.File(dmatrix_file,'r') as hf:
            mces = hf["mces"][:]
            all_smiles =[s.decode() for s in hf['mces_smiles_order'][:]] # Since here the input is csv we need to decode the lookup
    mces = squareform(mces)
    smiles_index = {}
    for i, smiles in enumerate(all_smiles):
        smiles_index[smiles] = i
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


def db_filter_inputs(smiles_list, index, client, threshold, always_stronger_bound, batch_size=50000):
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
    """
    print('filtering for precomputed mces (DB)')
    t0 = time.time()
    from rdkit import Chem

    # pre-convert unique SMILES to first block inchikeys
    smiles_to_inchikey = {}
    for s in smiles_list:
        if isinstance(s, bytes):
            s = s.decode()
        if s in smiles_to_inchikey:
            continue
        mol = Chem.MolFromSmiles(s)
        smiles_to_inchikey[s] = Chem.MolToInchiKey(mol)[:14] if mol is not None else None

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
            ninstances = len(smiles_input1) * len(smiles_input2)
            if ninstances >= np.iinfo(np.int64).max:
                # can probably not happen anyways?
                raise Exception('too many instances:', ninstances)
            nbatches = int(np.ceil(ninstances/args.batch_size))
            print(f'read {len(smiles_input1):_}*{len(smiles_input2):_} SMILES as input -> {ninstances:_} instances '
                  f'-> {nbatches} batches of {args.batch_size:_}')
            print('creating full index array')
            index_array_full = np.zeros((ninstances, 3), dtype='int64')
            index_array_full[:, 0] = np.arange(ninstances)
            index_array_full[:, 1:] = np.stack(np.meshgrid(range(0, len(smiles_input1)),
                                                           range(len(smiles_input1), len(smiles_input))),
                                               axis=-1).reshape(-1, 2)
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
            print('creating full index array')
            index_array_full = np.zeros((ninstances, 3), dtype='int64')
            index_array_full[:, 0] = np.arange(ninstances)
            index_array_full[:, 1:] = np.stack(np.triu_indices(n, k=1), axis=-1)
        if (not args.no_shuffle):
            print('shuffling full index array')
            np.random.shuffle(index_array_full)
        if args.use_db_lookup:
            from mcesdb import McesDbClient
            # empty conninfo -> libpq falls back to PGHOST/PGDATABASE/... env vars
            with McesDbClient("") as client:
                index_array_full, precomputed_index, precomputed_mces = db_filter_inputs(
                    smiles_list=smiles_input, index=index_array_full, client=client,
                    threshold=args.threshold, always_stronger_bound=not args.choose_bound_dynamically)
            print(f'Found {len(precomputed_index):_} precomputed MCES -> {len(index_array_full):_} instances left to compute'
                  f'-> {int(np.ceil(index_array_full.shape[0]/args.batch_size))} batches of {args.batch_size:_} and 1 precomputed batch')
        elif args.use_matrix_lookup:
            index_array_full, precomputed_index,precomputed_mces = filter_inputs(smiles_list=smiles_input,index=index_array_full, dmatrix_file=args.use_matrix_lookup,threshold=args.lookup_threshold)
            print(f'Found {len(precomputed_index):_} precomputed MCES -> {len(index_array_full):_} instances left to compute'
                  f'-> {int(np.ceil(index_array_full.shape[0]/args.batch_size))} batches of {args.batch_size:_} and 1 precomputed batch')
        # splitting
        if (not exists(args.out_folder)):
            mkdir(args.out_folder)
        print('creating batches')
        has_precomputed = (args.use_matrix_lookup or args.use_db_lookup) and precomputed_index.shape[0] > 0
        # when a precomputed batch exists it takes slot 0 (batch0_precomputed),
        # and the to-compute batches start at 1 so they sort after it
        i_offset = 1 if has_precomputed else 0
        if has_precomputed:
            file_path = join(args.out_folder, 'batch0_precomputed.hdf5')
            print(f'creating batch0_precomputed at {file_path} with precomputed mces...', end=' ')
            with h5py.File(file_path, 'w') as f:
                f.create_dataset('computation_indices', data=precomputed_index, dtype='int64', compression='gzip')
                f.create_dataset('mces', data=[row[1] for row in precomputed_mces], compression='gzip')
                f.create_dataset('smiles', data=smiles_input)
                f.create_dataset('computation_times', data=[row[2] if len(row) > 2 else -1 for row in precomputed_mces], compression='gzip')
                f.create_dataset('computation_modes', data=[row[3] if len(row) > 3 else -1 for row in precomputed_mces], compression='gzip')
            print('done')
        if index_array_full.shape[0] >0:
            for i, batch in enumerate(np.array_split(index_array_full, (index_array_full.shape[0]/args.batch_size)+1)):
                file_path = join(args.out_folder, f'batch{i + i_offset}.hdf5')
                print(f'creating batch {i + i_offset} at {file_path}...', end=' ')
                with h5py.File(file_path, 'w') as f:
                    f.create_dataset('computation_indices', data=batch, dtype='int64', compression='gzip')
                    f.create_dataset('smiles', data=smiles_input)
                    f.attrs['original_path'] = file_path
                print('done')
        else:
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
