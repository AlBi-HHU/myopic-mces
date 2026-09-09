import shutil
import subprocess
from pathlib import Path

import h5py
import pytest

from test_prepare_input import SMILES, _prepare
from test_combine_batches import _compute_batches, _combine, _naive_lookup

REPO_ROOT = Path(__file__).resolve().parent.parent

pytestmark = pytest.mark.skipif(shutil.which('nextflow') is None, reason='nextflow not installed')


def _run_nextflow(tmp_path, smiles_a, smiles_b):
    a_file = tmp_path / 'nf_a.txt'
    b_file = tmp_path / 'nf_b.txt'
    a_file.write_text('\n'.join(smiles_a) + '\n')
    b_file.write_text('\n'.join(smiles_b) + '\n')
    out_dir = tmp_path / 'nf_out'
    subprocess.run(
        ['nextflow', '-c', str(REPO_ROOT / 'nextflow.config'),
         'run', str(REPO_ROOT / 'nextflow' / 'two_datasets.nf'),
         '--smiles_a', str(a_file), '--smiles_b', str(b_file),
         '--batch_size', '13', '--cpus', '2', '--out', str(out_dir)],
        cwd=tmp_path, check=True,
    )
    return out_dir / 'combined.hdf5'


def test_csv_matches_manual_hdf5_matches_nextflow(tmp_path):
    smiles_a, smiles_b = SMILES[:7], SMILES[7:]

    batch_files = _prepare(tmp_path, smiles_a, smiles_b, batch_size=13)
    _compute_batches(batch_files)
    manual_combined = _combine(tmp_path, batch_files, shape=(len(smiles_b), len(smiles_a)))
    with h5py.File(manual_combined, 'r') as f:
        manual_dim1 = [s.decode() for s in f['smiles_dim1'][:]]
        manual_dim2 = [s.decode() for s in f['smiles_dim2'][:]]
        manual_mces = f['mces'][:]
        manual_modes = f['computation_modes_order'][:]

    nf_combined = _run_nextflow(tmp_path, smiles_a, smiles_b)
    with h5py.File(nf_combined, 'r') as f:
        nf_dim1 = [s.decode() for s in f['smiles_dim1'][:]]
        nf_dim2 = [s.decode() for s in f['smiles_dim2'][:]]
        nf_mces = f['mces'][:]
        nf_modes = f['computation_modes_order'][:]

    assert manual_dim1 == nf_dim1 and manual_dim2 == nf_dim2
    pairs = [(s1, s2) for s1 in manual_dim1 for s2 in manual_dim2]
    expected = _naive_lookup(tmp_path, pairs)

    for i, s1 in enumerate(manual_dim1):
        for j, s2 in enumerate(manual_dim2):
            exp_dist, exp_mode = expected[(s1, s2)]
            assert manual_mces[i, j] == exp_dist
            assert manual_modes[i, j] == exp_mode

    for i, s1 in enumerate(nf_dim1):
        for j, s2 in enumerate(nf_dim2):
            exp_dist, exp_mode = expected[(s1, s2)]
            assert nf_mces[i, j] == exp_dist
            assert nf_modes[i, j] == exp_mode
