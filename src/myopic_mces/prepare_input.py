from random import shuffle
from argparse import ArgumentParser
from os.path import join, exists
from os import mkdir
import numpy as np

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('input_file', help='readily formatedd myopic_mces input OR list of smiles for `h5`-mode')
    parser.add_argument('out_folder')
    parser.add_argument('--hdf5_mode', action='store_true')
    parser.add_argument('--batch_size', type=int, default=224960 * 100)
    parser.add_argument('--batch_start', type=int, default=0)
    parser.add_argument('--no_shuffle', action='store_true')
    args = parser.parse_args()

    if (args.hdf5_mode):
        import h5py
        from scipy.special import binom
        print('hdf5-mode...')
        # TODO: also make mode for TWO input files (e.g., dataset x bio)
        smiles_input = [l.strip() for l in open(args.input_file).readlines()]
        n = len(smiles_input)
        ninstances = int(binom(n, 2))
        if ninstances >= np.iinfo(np.int64).max:
            # can probably not happen anyways?
            raise Exception('too many instances:', ninstances)
        nbatches = int(np.ceil(ninstances/args.batch_size))
        print(f'read {n:_} SMILES as input -> {ninstances:_} instances -> {nbatches} batches of {args.batch_size:_}')
        settings = dict(threshold=args.threshold, choose_bound_dynamically=args.choose_bound_dynamically)
        print('creating full index array')
        index_array_full = np.zeros((ninstances, 3), dtype='int64')
        index_array_full[:, 0] = np.arange(ninstances)
        index_array_full[:, 1:] = np.stack(np.triu_indices(n, k=1), axis=-1)
        if (not args.no_shuffle):
            print('shuffling full index array')
            np.random.shuffle(index_array_full)
        # splitting
        if (not exists(args.out_folder)):
            mkdir(args.out_folder)
        print('creating batches')
        for i, batch in enumerate(np.array_split(index_array_full, (index_array_full.shape[0]/args.batch_size)+1)):
            file_path = join(args.out_folder, f'batch{i}.hdf5')
            print(f'creating batch {i} at {file_path}...', end=' ')
            with h5py.File(file_path, 'w') as f:
                f.create_dataset('computation_indices', data=batch, dtype='int64', compression='gzip')
                f.create_dataset('smiles', data=smiles_input)
                settings_group = f.create_group('settings')
                for k, v in settings.items():
                    settings_group[k] = v
                f.attrs['original_path'] = file_path
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
