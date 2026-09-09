import subprocess
import sys

import h5py
import numpy as np

from test_prepare_input import SMILES, _prepare

THRESHOLD = 10.0
SOLVER = 'PULP_CBC_CMD'


def _compute_batches(batch_files):
    for f in batch_files:
        subprocess.run(
            [sys.executable, '-m', 'myopic_mces.myopic_mces', str(f), '--hdf5_mode',
             '--threshold', str(THRESHOLD), '--solver', SOLVER],
            check=True, capture_output=True, text=True,
        )


def _combine(tmp_path, batch_files, shape=None):
    out_file = tmp_path / 'combined.hdf5'
    args = [sys.executable, '-m', 'myopic_mces.combine_hdf5_batches', *[str(f) for f in batch_files],
            '--out', str(out_file)]
    if shape is not None:
        args += ['--two_datasets_shape', *[str(s) for s in shape]]
    subprocess.run(args, check=True, capture_output=True, text=True)
    return out_file


def _naive_lookup(tmp_path, pairs):
    # keyed by (smiles1, smiles2) rather than position, so it can check the combined
    # output's actual (row, col) -> pair mapping instead of assuming an order
    input_csv = tmp_path / 'naive_input.csv'
    input_csv.write_text('\n'.join(f'{i},{s1},{s2}' for i, (s1, s2) in enumerate(pairs)) + '\n')
    output_csv = tmp_path / 'naive_output.csv'
    subprocess.run(
        [sys.executable, '-m', 'myopic_mces.myopic_mces', str(input_csv), str(output_csv),
         '--threshold', str(THRESHOLD), '--solver', SOLVER],
        check=True, capture_output=True, text=True,
    )
    rows = [line.split(',') for line in output_csv.read_text().splitlines()]
    return {pairs[int(r[0])]: (float(r[1]), int(r[3])) for r in rows}


def test_combine_matches_naive_self_mode(tmp_path):
    batch_files = _prepare(tmp_path, SMILES, batch_size=13)
    _compute_batches(batch_files)
    combined = _combine(tmp_path, batch_files)
    with h5py.File(combined, 'r') as f:
        smiles = [s.decode() for s in f['mces_smiles_order'][:]]
        mces = f['mces'][:]
        modes = f['computation_modes_order'][:]

    i_arr, j_arr = np.triu_indices(len(smiles), k=1)
    pairs = [(smiles[i], smiles[j]) for i, j in zip(i_arr, j_arr)]
    expected = _naive_lookup(tmp_path, pairs)

    for k, pair in enumerate(pairs):
        exp_dist, exp_mode = expected[pair]
        assert mces[k] == exp_dist
        assert modes[k] == exp_mode


def test_combine_matches_naive_two_file_mode(tmp_path):
    smiles_a, smiles_b = SMILES[:7], SMILES[7:]
    batch_files = _prepare(tmp_path, smiles_a, smiles_b, batch_size=13)
    _compute_batches(batch_files)
    # shape is (extra, main), matching prepare_input's id ordering -- see test_prepare_input.py
    combined = _combine(tmp_path, batch_files, shape=(len(smiles_b), len(smiles_a)))
    with h5py.File(combined, 'r') as f:
        dim1 = [s.decode() for s in f['smiles_dim1'][:]]
        dim2 = [s.decode() for s in f['smiles_dim2'][:]]
        mces = f['mces'][:]
        modes = f['computation_modes_order'][:]

    pairs = [(s1, s2) for s1 in dim1 for s2 in dim2]
    expected = _naive_lookup(tmp_path, pairs)

    for i, s1 in enumerate(dim1):
        for j, s2 in enumerate(dim2):
            exp_dist, exp_mode = expected[(s1, s2)]
            assert mces[i, j] == exp_dist
            assert modes[i, j] == exp_mode
