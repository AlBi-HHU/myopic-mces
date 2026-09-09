import math
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

SMILES = [line.split(',')[1] for line in Path('example/example_data.csv').read_text().splitlines()[:12]]


def _prepare(tmp_path, smiles_a, smiles_b=None, batch_size=1000):
    a_file = tmp_path / 'a.txt'
    a_file.write_text('\n'.join(smiles_a) + '\n')
    args = [sys.executable, '-m', 'myopic_mces.prepare_input', str(a_file), str(tmp_path / 'out'),
            '--hdf5_mode', '--batch_size', str(batch_size), '--no_shuffle']
    if smiles_b is not None:
        b_file = tmp_path / 'b.txt'
        b_file.write_text('\n'.join(smiles_b) + '\n')
        args += ['--hdf5_extra_input_file', str(b_file)]
    subprocess.run(args, check=True, capture_output=True, text=True)
    out_dir = tmp_path / 'out'
    return sorted(out_dir.glob('batch*.hdf5'), key=lambda p: int(p.stem.removeprefix('batch')))


def _all_indices(batch_files):
    return np.concatenate([h5py.File(f, 'r')['computation_indices'][:] for f in batch_files])


def test_self_mode_order(tmp_path):
    batch_files = _prepare(tmp_path, SMILES, batch_size=13)
    indices = _all_indices(batch_files)
    i_arr, j_arr = np.triu_indices(len(SMILES), k=1)
    expected = np.stack([i_arr, j_arr], axis=1)
    assert np.array_equal(indices[:, 0], np.arange(len(expected)))
    assert np.array_equal(indices[:, 1:], expected)


def test_two_file_mode_order(tmp_path):
    smiles_a, smiles_b = SMILES[:7], SMILES[7:]
    batch_files = _prepare(tmp_path, smiles_a, smiles_b, batch_size=13)
    indices = _all_indices(batch_files)
    n_a = len(smiles_a)
    i_arr = np.tile(np.arange(n_a), len(smiles_b))
    j_arr = np.repeat(np.arange(n_a, n_a + len(smiles_b)), n_a)
    expected = np.stack([i_arr, j_arr], axis=1)
    assert np.array_equal(indices[:, 0], np.arange(len(expected)))
    assert np.array_equal(indices[:, 1:], expected)


@pytest.mark.parametrize('split', [None, 7])
def test_batch_count(tmp_path, split):
    smiles_a = SMILES if split is None else SMILES[:split]
    smiles_b = None if split is None else SMILES[split:]
    batch_files = _prepare(tmp_path, smiles_a, smiles_b, batch_size=25)
    if smiles_b is None:
        ninstances = len(smiles_a) * (len(smiles_a) - 1) // 2
    else:
        ninstances = len(smiles_a) * len(smiles_b)
    assert len(batch_files) == math.ceil(ninstances / 25)
    assert sum(len(h5py.File(f, 'r')['computation_indices']) for f in batch_files) == ninstances


@pytest.mark.parametrize('split', [None, 7])
def test_batch_smiles_match_input(tmp_path, split):
    smiles_a = SMILES if split is None else SMILES[:split]
    smiles_b = None if split is None else SMILES[split:]
    batch_files = _prepare(tmp_path, smiles_a, smiles_b, batch_size=13)
    expected = smiles_a + (smiles_b or [])
    for f in batch_files:
        stored = [s.decode() for s in h5py.File(f, 'r')['smiles'][:]]
        assert stored == expected
