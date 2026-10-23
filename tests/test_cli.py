import subprocess
import sys

import pytest


def _run_cli(*args):
    return subprocess.run(
        [sys.executable, '-m', 'myopic_mces.myopic_mces', *args],
        capture_output=True,
        text=True,
    )


def test_output_required_without_hdf5_mode():
    proc = _run_cli('input.csv')
    assert proc.returncode == 2
    assert 'the following argument is required: output' in proc.stderr


def test_output_optional_with_hdf5_mode():
    proc = _run_cli('input.hdf5', '--hdf5_mode')
    assert proc.returncode != 0
    assert 'the following argument is required: output' not in proc.stderr


@pytest.mark.parametrize('stereo, expected', [(False, [0., 0., 1.]), (True, [0.5, 0.5, 1.5])])
def test_stereo(tmp_path, stereo, expected):
    import pulp
    solver = 'CPLEX_PY' if 'CPLEX_PY' in pulp.listSolvers(onlyAvailable=True) else 'COIN_CMD'
    input_file, output_file = tmp_path / 'input.csv', tmp_path / 'output.csv'
    input_file.write_text('0,C/C=C/C,C/C=C\\C\n'
                          '1,C[C@H](N)O,C[C@@H](N)O\n'
                          '2,N[C@@H](C)C(=O)O,N[C@H](CS)C(=O)O\n')
    proc = _run_cli(str(input_file), str(output_file), '--threshold', '-1', '--solver', solver, '--solver_no_msg',
                    *(['--stereo'] if stereo else []))
    assert proc.returncode == 0, proc.stderr
    distances = [float(line.split(',')[1]) for line in output_file.read_text().splitlines()]
    assert distances == pytest.approx(expected)
