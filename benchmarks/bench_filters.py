"""Benchmark filter0/1/2: C++ (mmces_filters) vs. pure Python (networkx).

See README.md. Requires the C++ extension (myopic_mces_cpp) to be installed.
"""
import argparse
import contextlib
import datetime
import io
import subprocess
import sys
import time
from pathlib import Path

from rdkit import Chem, RDLogger  # rdkit must be imported before mmces_filters
import mmces_filters as mf

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "src"))
with contextlib.redirect_stdout(io.StringIO()):  # filter_MCES prints loaded .so files on import
    from myopic_mces import filter_MCES as pf
from myopic_mces.graph import construct_graph

RDLogger.DisableLog("rdApp.*")

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("pairs", help="CSV with header smiles1,smiles2 or headerless id,smiles1,smiles2,...")
parser.add_argument("-n", type=int, default=None, help="only use the first n pairs")
parser.add_argument("--pairs-out", help="write per-pair values and times to this CSV")
parser.add_argument("--summary", default=str(HERE / "results.csv"),
                    help="append a summary row to this CSV (default: %(default)s)")
args = parser.parse_args()

pairs = []
with open(args.pairs) as f:
    for i, line in enumerate(f):
        parts = line.rstrip("\n").split(",")
        if i == 0 and parts[:2] == ["smiles1", "smiles2"]:
            two_cols = True
            continue
        if i == 0:
            two_cols = False
        pairs.append((str(len(pairs)), parts[0], parts[1]) if two_cols else (parts[0], parts[1], parts[2]))
        if args.n and len(pairs) >= args.n:
            break

# preprocessing (not part of filter timings)
t0 = time.perf_counter()
mols = [(Chem.MolFromSmiles(s1), Chem.MolFromSmiles(s2)) for _, s1, s2 in pairs]
t_mol = time.perf_counter() - t0
t0 = time.perf_counter()
graphs = [(construct_graph(m1), construct_graph(m2)) for m1, m2 in mols]
t_graph = time.perf_counter() - t0

FILTERS = ["filter0", "filter1", "filter2"]
impls = {
    "py": {f: getattr(pf, f) for f in FILTERS},
    "cpp": {f: getattr(mf, f) for f in FILTERS},
}
inputs = {"py": graphs, "cpp": mols}

totals = {(impl, f): 0.0 for impl in impls for f in FILTERS}
mismatches = {f: 0 for f in FILTERS}
rows = []
perf = time.perf_counter
for i, (pid, _, _) in enumerate(pairs):
    row = [pid]
    for f in FILTERS:
        vals = {}
        for impl in impls:
            fn = impls[impl][f]
            a, b = inputs[impl][i]
            t0 = perf()
            vals[impl] = fn(a, b)
            dt = perf() - t0
            totals[(impl, f)] += dt
            row += [vals[impl], dt]
        if abs(vals["py"] - vals["cpp"]) > 1e-6:
            mismatches[f] += 1
    rows.append(row)
    if (i + 1) % 1000 == 0:
        print(f"{i + 1}/{len(pairs)} done", file=sys.stderr, flush=True)

n = len(pairs)
print(f"pairs: {n}")
print(f"preprocessing: MolFromSmiles {t_mol:.2f}s, construct_graph (py only, on top) {t_graph:.2f}s")
print(f"{'filter':8} {'python [s]':>11} {'C++ [s]':>9} {'py/pair [us]':>13} {'C++/pair [us]':>14} {'speedup':>8} {'mismatches':>11}")
for f in FILTERS:
    tp, tc = totals[("py", f)], totals[("cpp", f)]
    print(f"{f:8} {tp:11.2f} {tc:9.2f} {tp / n * 1e6:13.1f} {tc / n * 1e6:14.1f} {tp / tc:7.1f}x {mismatches[f]:11d}")

if args.pairs_out:
    with open(args.pairs_out, "w") as f:
        f.write("id," + ",".join(f"{flt}_{impl}_{k}" for flt in FILTERS for impl in impls for k in ("value", "time")) + "\n")
        for row in rows:
            f.write(",".join(map(str, row)) + "\n")

git = lambda *a: subprocess.run(["git", "-C", str(REPO), *a], capture_output=True, text=True).stdout.strip()
commit = git("rev-parse", "--short", "HEAD") + ("-dirty" if git("status", "--porcelain", "--untracked-files=no", "cpp/src", "src") else "")
summary = Path(args.summary)
fields = {
    "date": datetime.date.today().isoformat(),
    "commit": commit,
    "dataset": Path(args.pairs).name,
    "n_pairs": n,
}
for f in FILTERS:
    tp, tc = totals[("py", f)], totals[("cpp", f)]
    fields[f"{f}_py_us"] = f"{tp / n * 1e6:.1f}"
    fields[f"{f}_cpp_us"] = f"{tc / n * 1e6:.1f}"
    fields[f"{f}_speedup"] = f"{tp / tc:.1f}"
    fields[f"{f}_mismatches"] = mismatches[f]
write_header = not summary.exists()
with open(summary, "a") as f:
    if write_header:
        f.write(",".join(fields) + "\n")
    f.write(",".join(map(str, fields.values())) + "\n")
print(f"summary appended to {summary}")
