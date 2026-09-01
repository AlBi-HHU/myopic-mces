"""Two modes: 1. one dataset, all compounds against all
              2. two datasets, all compounds from dataset 1 against all compounds from dataset 2
                 Use --two_datasets_shape!
"""


import h5py
from scipy.special import binom
import numpy as np
from argparse import ArgumentParser
from time import time
from myopic_mces.myopic_mces import MCES

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('input_batches', help='bachted HDF5 output files from MCES computations', nargs='+')
    parser.add_argument('--out', help='HDF5 file for combined results', required=True)
    parser.add_argument('--two_datasets_shape', help='Use this if the instances are from two different datasets (e.g., biomolecules versus target dataset), '
                        'specify their sizes: e.g., 19994 2000 if the target dataset contains 2000 compounds', nargs='+', type=int, default=None)
    parser.add_argument('--sanity', type=int, default=None,
                        help='if set, recompute this many random pairs and compare against the combined results as a sanity check')
    args = parser.parse_args()

    t0 = time()

    with h5py.File(args.input_batches[0], 'r') as f:
        smiles = list(f['smiles'])
        computation_args = {k: f['computation_args'][k][()] for k in f['computation_args']}
        ninstances = int(binom(len(smiles), 2))

    if (args.two_datasets_shape is not None):
        if (len(args.two_datasets_shape) != 2):
            raise Exception('for two datasets, specify two numbers', args.two_datasets_shape)
        ninstances = args.two_datasets_shape[0] * args.two_datasets_shape[1]

    all_indices = np.empty((ninstances, 3), dtype='int64')
    all_mces = np.empty((ninstances,), dtype='float32')
    all_modes = np.empty((ninstances,), dtype='uint8')
        
    i = 0
    for batch in args.input_batches:
        print('adding computations from', batch)
        with h5py.File(batch, 'r') as f:
            n = len(f['computation_indices'])
            f['mces'].read_direct(all_mces[i:(i+n)])
            f['computation_modes'].read_direct(all_modes[i:(i+n)])
            f['computation_indices'].read_direct(all_indices[i:(i+n)])
            i += n

    print('restoring correct order')
    order = all_indices[:, 0].argsort()
    all_mces_order = all_mces[order]
    all_modes_order = all_modes[order]
    all_indices_order = all_indices[order]

    if (args.two_datasets_shape is not None):
        print(f'reshaping computation matrix from {all_mces_order.shape} to {tuple(args.two_datasets_shape)}')
        all_mces_order = all_mces_order.reshape(tuple(args.two_datasets_shape))
        all_modes_order = all_modes_order.reshape(tuple(args.two_datasets_shape))
        

    if (args.sanity is not None):
        print(f'sanity check: recomputing {args.sanity} random pairs and comparing to combined results')
        solver = computation_args.get('solver', 'COIN_CMD')
        if isinstance(solver, bytes):
            solver = solver.decode()
        threshold = float(computation_args.get('threshold', 10.))
        mces_flat = all_mces_order.reshape(-1)
        modes_flat = all_modes_order.reshape(-1)
        check_indices = np.random.choice(len(all_indices_order), size=min(args.sanity, len(all_indices_order)), replace=False)
        nmismatches = 0
        for idx in check_indices:
            _, s1_i, s2_i = all_indices_order[idx]
            smiles1 = smiles[s1_i].decode() if isinstance(smiles[s1_i], bytes) else smiles[s1_i]
            smiles2 = smiles[s2_i].decode() if isinstance(smiles[s2_i], bytes) else smiles[s2_i]
            _, distance, _, mode = MCES(smiles1, smiles2, threshold=threshold, solver=solver)
            if not np.isclose(distance, mces_flat[idx]) or mode != modes_flat[idx]:
                nmismatches += 1
                print(f'sanity check MISMATCH at index {idx}: recomputed distance {distance} (mode {mode}), '
                      f'stored distance {mces_flat[idx]} (mode {modes_flat[idx]})')
        print(f'sanity check done: {nmismatches}/{len(check_indices)} mismatches')

    print('writing output to', args.out)
    with h5py.File(args.out, 'w') as f:
        g = f.create_group('computation_args')
        for k, v in computation_args.items():
            g.create_dataset(k, data=v)
        f.create_dataset('mces', data=all_mces_order)
        f.create_dataset('computation_modes_order', data=all_modes_order)
        if (args.two_datasets_shape is not None):
            n1, n2 = args.two_datasets_shape
            f.create_dataset('smiles_dim1', data=smiles[n2:n1+n2])
            f.create_dataset('smiles_dim2', data=smiles[:n2])
        else:
            f.create_dataset('mces_smiles_order', data=smiles)

    print(f'combined {len(all_indices_order)} MCES computations in {(time() - t0) / 60:.1f}min')
