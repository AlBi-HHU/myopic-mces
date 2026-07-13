# Nextflow workflows for myopic MCES

Two workflows built on shared modules:

- `one_dataset.nf` — all pairwise MCES distances within a single SMILES list.
- `two_datasets.nf` — all cross MCES distances between two SMILES lists.

Both run the same four steps: prepare input into batched HDF5, compute each
batch in parallel, combine batches into one HDF5, and (optionally) sanity-check
a random sample of distances against a fresh re-computation.

## Usage

```bash
# one-dataset: all pairwise distances within one list
nextflow run nextflow/one_dataset.nf --smiles smiles.txt

# two-dataset: cross distances between two lists
nextflow run nextflow/two_datasets.nf --smiles_a dataset_a.txt --smiles_b dataset_b.txt

# with db lookup (needs PGHOST/PGDATABASE/... env vars)
nextflow run nextflow/two_datasets.nf --smiles_a dataset_a.txt --smiles_b dataset_b.txt --use_db_lookup

# with CPLEX
nextflow run nextflow/one_dataset.nf --smiles smiles.txt --cplex_home /opt/ibm/ILOG/CPLEX_Studio
```

## Parameters

Defaults are defined in `nextflow.config`; override any on the command line.

| Parameter | Default | Description |
|---|---|---|
| `smiles` | — | One-dataset mode: path to a `.txt` file (one SMILES per line) |
| `smiles_a` | — | Two-dataset mode: main dataset (passed as `input_file` to `prepare_input`) |
| `smiles_b` | — | Two-dataset mode: extra dataset (passed as `--hdf5_extra_input_file`) |
| `threshold` | `10` | Distance threshold for the computation |
| `batch_size` | `224960000` | Number of instances per HDF5 batch |
| `use_db_lookup` | `false` | Reuse precomputed distances from mcesdb (needs PG env vars) |
| `cplex_home` | `null` | CPLEX install dir; when unset, falls back to `COIN_CMD` |
| `choose_bound_dynamically` | `false` | Use a faster, potentially weaker lower bound when already above threshold |
| `cpus` | `8` | Cores per `COMPUTE_BATCH` task. Set to your node's total core count to run one batch at a time, each using all cores. Lower values = more batches in parallel, fewer cores each. |
| `sanity_check` | `true` | Run a random-sample sanity check after combining |
| `sanity_samples` | `10000` | Number of pairs to recompute in the sanity check |
| `out` | `results` | Output directory for `combined.hdf5` and `sanity_report.txt` |

## Output

- `combined.hdf5` — combined MCES distances (condensed `mces` for one-dataset
  mode, `(n_b, n_a)` matrix for two-dataset mode).
- `sanity_report.txt` — `True`/`False` from the sanity check.
