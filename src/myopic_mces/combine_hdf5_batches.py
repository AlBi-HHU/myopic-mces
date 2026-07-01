"""Two modes: 1. one dataset, all compounds against all
              2. two datasets, all compounds from dataset 1 against all compounds from dataset 2
                 Use --two_datasets_shape!
"""


import h5py
from scipy.special import binom
import numpy as np
from argparse import ArgumentParser
from time import time

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('input_batches', help='bachted HDF5 output files from MCES computations', nargs='+')
    parser.add_argument('--out', help='HDF5 file for combined results', required=True)
    parser.add_argument('--two_datasets_shape', help='Use this if the instances are from two different datasets (e.g., biomolecules versus target dataset), '
                        'specify their sizes: e.g., 19994 2000 if the target dataset contains 2000 compounds', nargs='+', type=int, default=None)
    args = parser.parse_args()

    t0 = time()

    with h5py.File(args.input_batches[0], 'r') as f:
        smiles = list(f['smiles'])
        ninstances = int(binom(len(smiles), 2))

    if (args.two_datasets_shape is not None):
        if (len(args.two_datasets_shape) != 2):
            raise Exception('for two datasets, specify two numbers', args.two_datasets_shape)
        ninstances = args.two_datasets_shape[0] * args.two_datasets_shape[1]

    all_indices = np.empty((ninstances, 3), dtype='int64')
    all_mces = np.empty((ninstances,), dtype='float32')
        
    i = 0
    for batch in args.input_batches:
        print('adding computations from', batch)
        with h5py.File(batch, 'r') as f:
            n = len(f['computation_indices'])
            f['mces'].read_direct(all_mces[i:(i+n)])
            f['computation_indices'].read_direct(all_indices[i:(i+n)])
            i += n

    print('restoring correct order')
    order = all_indices[:, 0].argsort()
    all_mces_order = all_mces[order]
    all_indices_order = all_indices[order]

    if (args.two_datasets_shape is not None):
        print(f'reshaping computation matrix from {all_mces_order.shape} to {tuple(args.two_datasets_shape)}')
        all_mces_order = all_mces_order.reshape(tuple(args.two_datasets_shape))
        

    # TODO: sanity check, repeat and test some computations
        
    print('writing output to', args.out)
    with h5py.File(args.out, 'w') as f:
        f.create_dataset('mces', data=all_mces_order)
        f.create_dataset('mces_smiles_order', data=smiles)

    print(f'combined {len(all_indices_order)} MCES computations in {(time() - t0) / 60:.1f}min')
