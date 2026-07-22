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
