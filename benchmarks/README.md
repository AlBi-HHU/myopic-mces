# Benchmarks

Compares the runtime of filter0/1/2 between the C++ (`myopic_mces_cpp`) and pure Python implementations
and checks that both return the same values.

```sh
curl -LO https://raw.githubusercontent.com/boecker-lab/myopic-mces-data/main/subsampled_instances_20k.csv
python benchmarks/bench_filters.py subsampled_instances_20k.csv
```

Each run appends one line to `results.csv` (times are per pair, in µs). Use `-n` for a quick run on fewer pairs.

## Log

C++ runtime per pair on `subsampled_instances_20k.csv` (called from Python), speedup over Python in parentheses.
C++ and Python values agree on all pairs in every entry.

| commit | filter0 | filter1 | filter2 | change |
|---|---|---|---|---|
| [aab50cc](https://github.com/AlBi-HHU/myopic-mces/commit/aab50cc) | 11.5 µs (3.3×) | 19.5 µs (8.3×) | 1016 µs (10.4×) | baseline |
| [9f45d4b](https://github.com/AlBi-HHU/myopic-mces/commit/9f45d4b) | 11.7 µs (2.8×) | 19.8 µs (8.0×) | 96.9 µs (108×) | filter2: per-atom neighborhoods precomputed, O(n³) Hungarian |
| [2f95545](https://github.com/AlBi-HHU/myopic-mces/commit/2f95545) | 5.9 µs (5.6×) | 7.1 µs (22.9×) | 91.6 µs (116×) | integer atom types, filter1 as sorted merge |
| [71c696d](https://github.com/AlBi-HHU/myopic-mces/commit/71c696d) | 5.7 µs (5.3×) | 6.5 µs (24.1×) | 38.1 µs (266×) | filter2: LAPJV instead of Hungarian |
