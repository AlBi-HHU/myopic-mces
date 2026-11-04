import subprocess
import sys
import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

SMILES = [
    "CC(=O)Oc1ccccc1C(=O)O",      # 0: Aspirin
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # 1: Ibuprofen
    "CC(=O)Nc1ccc(O)cc1",          # 2: Paracetamol
]


def _make_lookup(path, smiles, distances):
    with h5py.File(path, "w") as f:
        f.create_dataset("mces", data=np.array(distances, dtype="float32"))
        f.create_dataset("mces_smiles_order", data=[s.encode() for s in smiles])


@pytest.mark.parametrize("threshold, exp_pre, exp_comp", [
    (None, 1, 2),  # (0,1)=2.0 found; (0,2) and (1,2) to compute
    (1.0,  0, 3),  # 2.0 >= 1.0 threshold -> not reused, all 3 to compute
])
def test_hdf5_matrix_lookup(tmp_path, threshold, exp_pre, exp_comp):
    lookup = tmp_path / "lookup.hdf5"
    _make_lookup(lookup, SMILES[:2], [2.0])  # Precompute pair (0, 1)

    input_file = tmp_path / "input.txt"
    input_file.write_text("\n".join(SMILES) + "\n")
    out_dir = tmp_path / "out"

    cmd = [
        sys.executable, "-m", "myopic_mces.prepare_input",
        str(input_file), str(out_dir),
        "--hdf5_mode", "--use_matrix_lookup", str(lookup), "--no_shuffle",
    ]
    if threshold is not None:
        cmd += ["--lookup_threshold", str(threshold)]

    subprocess.run(cmd, check=True, capture_output=True, text=True)

    if exp_pre > 0:
        with h5py.File(out_dir / "batch0_precomputed.hdf5", "r") as f:
            assert len(f["computation_indices"]) == exp_pre
            assert f["mces"][0] == 2.0
    else:
        assert not (out_dir / "batch0_precomputed.hdf5").exists()

    with h5py.File(out_dir / "batch1.hdf5", "r") as f:
        assert len(f["computation_indices"]) == exp_comp
