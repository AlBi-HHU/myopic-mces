from random import shuffle
from argparse import ArgumentParser
from os.path import join, exists
from os import mkdir
import numpy as np

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
        if args.use_matrix_lookup:
            index_array_full, precomputed_index,precomputed_mces = filter_inputs(smiles_list=smiles_input,index=index_array_full, dmatrix_file=args.use_matrix_lookup,threshold=args.lookup_threshold)
            print(f'Found {len(precomputed_index):_} precomputed MCES -> {len(index_array_full):_} instances left to compute'
                  f'-> {int(np.ceil(index_array_full.shape[0]/args.batch_size))} batches of {args.batch_size:_} and 1 precomputed batch')
        # splitting
        if (not exists(args.out_folder)):
            mkdir(args.out_folder)
        print('creating batches')
        if index_array_full.shape[0] >0:
            for i, batch in enumerate(np.array_split(index_array_full, (index_array_full.shape[0]/args.batch_size)+1)):
                file_path = join(args.out_folder, f'batch{i}.hdf5')
                print(f'creating batch {i} at {file_path}...', end=' ')
                with h5py.File(file_path, 'w') as f:
                    f.create_dataset('computation_indices', data=batch, dtype='int64', compression='gzip')
                    f.create_dataset('smiles', data=smiles_input)
                    f.attrs['original_path'] = file_path
                print('done')
        else:
            print("Nothing to compute, everything cached")
            i = -1 # Hacky but for now needed if all smiles are found in the lookup
        if args.use_matrix_lookup:
            if precomputed_index.shape[0]>0:
                file_path = join(args.out_folder, f'batch{i+1}.hdf5')
                print(f'creating batch {i+1} at {file_path} with precomputed mces...', end=' ')
                with h5py.File(file_path, 'w') as f:
                    f.create_dataset('computation_indices', data=precomputed_index, dtype='int64', compression='gzip')
                    f.create_dataset('mces', data=[row[1] for row in precomputed_mces], compression='gzip')
                    f.create_dataset('smiles', data=smiles_input)
                    f.create_dataset('computation_times', data=[row[2] if len(row) > 2 else -1 for row in precomputed_mces], compression='gzip')
                    f.create_dataset('computation_modes', data=[row[3] if len(row) > 3 else -1 for row in precomputed_mces], compression='gzip')         
                print('done')
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
