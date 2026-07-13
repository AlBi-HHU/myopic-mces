"""Sanity check for MCES HDF5 results.

Recomputes a random sample of pairs and compares to the combined HDF5.
Auto-detects one-dataset vs two-dataset mode from the HDF5 keys.

Usage:
    python3 sanity_check.py combined.hdf5 <threshold> <always_stronger> <n_samples>
"""

import sys
import h5py
import numpy as np
from random import randint
from scipy.spatial.distance import squareform
from myopic_mces.myopic_mces import MCES


combined_file = sys.argv[1]
threshold = float(sys.argv[2])
always_stronger = sys.argv[3] == "True"
n_samples = int(sys.argv[4])


def get_mces(s1, s2):
    return MCES(s1, s2, always_stronger_bound=always_stronger, threshold=threshold)[1]


result = False
with h5py.File(combined_file, "r") as f:
    if "mces_smiles_order" in f:
        mces = squareform(f["mces"][:])
        smiles = [s.decode() if isinstance(s, bytes) else s for s in f["mces_smiles_order"][:]]
        n = len(smiles)
        tests = []
        for _ in range(n_samples):
            i1 = randint(0, n - 1)
            i2 = randint(0, n - 1)
            if i1 == i2:
                continue
            expected = mces[i1, i2]
            actual = get_mces(smiles[i1], smiles[i2])
            tests.append(np.isclose(actual, expected))
        result = all(tests) if tests else False
        print(f"one-dataset mode: {n} smiles, {len(tests)} checks, all match: {result}")
    elif "smiles_dim1" in f:
        mces = f["mces"][:]
        smiles_a = [s.decode() if isinstance(s, bytes) else s for s in f["smiles_dim1"][:]]
        smiles_b = [s.decode() if isinstance(s, bytes) else s for s in f["smiles_dim2"][:]]
        n_a, n_b = mces.shape
        tests = []
        for _ in range(n_samples):
            i1 = randint(0, n_a - 1)
            i2 = randint(0, n_b - 1)
            expected = mces[i1, i2]
            actual = get_mces(smiles_a[i1], smiles_b[i2])
            tests.append(np.isclose(actual, expected))
        result = all(tests) if tests else False
        print(f"two-dataset mode: {n_a} x {n_b}, {len(tests)} checks, all match: {result}")
    else:
        print("ERROR: unknown HDF5 format")

with open("sanity_report.txt", "w") as out:
    out.write(str(result) + "\n")
