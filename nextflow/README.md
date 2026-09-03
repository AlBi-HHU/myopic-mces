# Nextflow workflows for myopic MCES

Two workflows built on shared modules:

- `one_dataset.nf` — all pairwise MCES distances within a single SMILES list.
- `two_datasets.nf` — all cross MCES distances between two SMILES lists.

Both run the same steps: prepare input into batched HDF5, compute each
batch in parallel, combine batches into one HDF5, (optionally) sanity-check
a random sample of distances against a fresh re-computation, and (optionally)
write the computed distances back into the mces database.

## Usage

`main.nf` at the repo root is the main entry point: it dispatches to
`one_dataset` or `two_datasets` based on which `--smiles*` params are set, so
`nextflow run .` works for either mode. It also lets you run the pipeline
straight from GitHub, on any branch, without cloning:

```bash
nextflow run AlBi-HHU/myopic-mces -r <branch> --smiles smiles.txt
```

The individual workflow files below remain runnable directly, too.

```bash
# one-dataset: all pairwise distances within one list
nextflow run nextflow/one_dataset.nf --smiles smiles.txt

# two-dataset: cross distances between two lists
nextflow run nextflow/two_datasets.nf --smiles_a dataset_a.txt --smiles_b dataset_b.txt

# with db lookup (needs PGHOST/PGDATABASE/... env vars)
nextflow run nextflow/two_datasets.nf --smiles_a dataset_a.txt --smiles_b dataset_b.txt --use_db_lookup

# write computed distances back to the mces database (needs PGHOST/PGDATABASE/... env vars)
nextflow run nextflow/one_dataset.nf --smiles smiles.txt --write_to_db

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
| `write_to_db` | `false` | Write computed distances back to mcesdb via `import-mces` (needs PG env vars) |
| `cplex_home` | `null` | CPLEX install dir; when unset, falls back to `COIN_CMD` |
| `solver_time_limit_seconds` | `null` | Per-instance ILP solver time limit; unset means no limit |
| `choose_bound_dynamically` | `false` | Use a faster, potentially weaker lower bound when already above threshold |
| `cpus` | `8` | Cores per `COMPUTE_BATCH` task. Set to your node's total core count to run one batch at a time, each using all cores. Lower values = more batches in parallel, fewer cores each. |
| `sanity_check` | `true` | Run a random-sample sanity check after combining |
| `sanity_samples` | `10000` | Number of pairs to recompute in the sanity check |
| `out` | `results` | Output directory for `combined.hdf5` and `sanity_report.txt` |

## Output

- `combined.hdf5` — combined MCES distances (condensed `mces` for one-dataset
  mode, `(n_b, n_a)` matrix for two-dataset mode).
- `sanity_report.txt` — `True`/`False` from the sanity check.
- `db_import.log` — output of `import-mces` when `--write_to_db` is set.

## Database integration

`--use_db_lookup` and `--write_to_db` both talk to a Postgres database managed
by the `mces-database` (`mcesdb`) package. Install it alongside this
package's `hdf5` extra, e.g.:

```bash
pip install -e /path/to/mces-database[importer]
pip install myopic_mces[hdf5]
```

and set the standard PostgreSQL environment variables (`PGHOST`, `PGDATABASE`,
`PGUSER`, `PGPASSWORD`, ...) before running Nextflow. `--write_to_db` shells
out to the `import-mces` CLI, which is only added to `PATH` once `mcesdb` is
installed.
