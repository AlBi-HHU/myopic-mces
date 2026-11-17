"""Many small one-dataset HDF5 files -> few compute batches -> results back into the small files.

Every input file holds its own SMILES list (key 'smiles' or 'mces_smiles_order'); all pairs
*within* a file are computed (no cross-file pairs). Files that already contain 'mces' are skipped.

    pack:   groups whole input files into batches of ~batch_size pairs, in the usual batch format
            ('smiles', 'computation_indices') that `myopic_mces --hdf5_mode` consumes. Each batch
            additionally remembers its source files ('source_files', 'source_pair_offsets').
    unpack: writes the results of one computed batch back into its source files, in the same
            format as combine_hdf5_batches ('mces', 'computation_modes_order', 'computation_args').
"""

import os
import sys
import h5py
import numpy as np
from argparse import ArgumentParser

SMILES_KEYS = ('smiles', 'mces_smiles_order')
STR = h5py.string_dtype()


def _read_smiles(f):
    for k in SMILES_KEYS:
        if k in f:
            return [s.decode() if isinstance(s, bytes) else s for s in f[k][:]]
    raise KeyError(f'{f.filename}: none of {SMILES_KEYS} found')


def pack(input_dir, out_dir, batch_size):
    os.makedirs(out_dir, exist_ok=True)
    files = sorted(os.path.join(os.path.abspath(input_dir), p) for p in os.listdir(input_dir)
                   if p.endswith(('.hdf5', '.h5')))
    nbatches = 0
    cur = dict(smiles=[], indices=[], sources=[], offsets=[], npairs=0)

    def flush():
        nonlocal nbatches, cur
        if cur['npairs'] == 0:
            return
        idx = np.concatenate(cur['indices'])
        with h5py.File(os.path.join(out_dir, f'batch{nbatches}.hdf5'), 'w') as f:
            f.create_dataset('computation_indices', data=idx, dtype='int64', compression='gzip')
            f.create_dataset('smiles', data=cur['smiles'], dtype=STR)
            f.create_dataset('source_files', data=cur['sources'], dtype=STR)
            f.create_dataset('source_pair_offsets', data=cur['offsets'] + [cur['npairs']], dtype='int64')
        nbatches += 1
        cur = dict(smiles=[], indices=[], sources=[], offsets=[], npairs=0)

    skipped = 0
    for path in files:
        with h5py.File(path, 'r') as f:
            if 'mces' in f:
                skipped += 1
                continue
            smiles = _read_smiles(f)
        n = len(smiles)
        i, j = np.triu_indices(n, k=1)   # same pair order as squareform / combine_hdf5_batches
        m = len(i)
        smiles_off = len(cur['smiles'])
        cur['indices'].append(np.stack([cur['npairs'] + np.arange(m), smiles_off + i, smiles_off + j], axis=1))
        cur['smiles'].extend(smiles)
        cur['sources'].append(path)
        cur['offsets'].append(cur['npairs'])
        cur['npairs'] += m
        if cur['npairs'] >= batch_size:
            flush()
    flush()
    print(f'{len(files)} files ({skipped} already computed) -> {nbatches} batches')


def unpack(batch):
    with h5py.File(batch, 'r') as f:
        mces = f['mces'][:]
        modes = f['computation_modes'][:]
        sources = [s.decode() if isinstance(s, bytes) else s for s in f['source_files'][:]]
        offsets = f['source_pair_offsets'][:]
        comp_args = {k: f['computation_args'][k][()] for k in f['computation_args']}
        assert (f['computation_indices'][:, 0] == np.arange(len(mces))).all()
    for path, start, end in zip(sources, offsets[:-1], offsets[1:]):
        with h5py.File(path, 'a') as f:
            for k in ('mces', 'computation_modes_order', 'computation_args'):
                if k in f:
                    del f[k]
            f.create_dataset('mces', data=mces[start:end])
            f.create_dataset('computation_modes_order', data=modes[start:end])
            g = f.create_group('computation_args')
            for k, v in comp_args.items():
                g.create_dataset(k, data=v)
        print(path)   # stdout: one written file per line (used by the workflow's sanity check)
    print(f'{batch}: wrote results into {len(sources)} files', file=sys.stderr)


if __name__ == '__main__':
    parser = ArgumentParser()
    sub = parser.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('pack')
    p.add_argument('input_dir')
    p.add_argument('out_dir')
    p.add_argument('--batch_size', type=int, required=True)
    u = sub.add_parser('unpack')
    u.add_argument('batches', nargs='+')
    args = parser.parse_args()
    if args.cmd == 'pack':
        pack(args.input_dir, args.out_dir, args.batch_size)
    else:
        for b in args.batches:
            unpack(b)
