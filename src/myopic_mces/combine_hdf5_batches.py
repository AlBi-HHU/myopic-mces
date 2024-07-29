import h5py
from scipy.special import binom
import numpy as np
from argparse import ArgumentParser
from time import time

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('input_batches', help='bachted HDF5 output files from MCES computations')
    parser.add_argument('--out', help='HDF5 file for combined results', required=True)
    args = parser.parse_args()

    t0 = time.time()

    with h5py.File(args.input_batches[0], 'r') as f:
        smiles = list(f['smiles'])
        ninstances = int(binom(len(smiles), 2))

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

    order = all_indices[:, 0].argsort()
    all_mces_order = all_mces[order]
    all_indices_order = all_indices[order]
    with h5py.File(args.out, 'w') as f:
        f.create_dataset('mces', data=all_mces_order)
        f.create_dataset('mces_smiles_order', data=smiles)

    print(f'combined {len(all_indices_order)} MCES computations in {(time.time() - t0) / 60:.1f}min')
