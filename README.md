# Computation of myopic MCES distances

Implementation of myopic MCES distance computation, see the preprint at [doi:10.1101/2023.03.27.534311](https://doi.org/10.1101/2023.03.27.534311) for details.

## Usage

Input and Output file are in csv format. Every line in the input-file is one comparison:

input-file: `index,SMILES1,SMILES2`

output-file: `index,time taken,myopic MCES distance,status (1 if exact distance, 2/4 if lower bound)`

Download via pip and execute:
```bash
pip install myopic-mces
myopic_mces input-file output-file
```

Alternatively, to install directly from this repository:
```bash
pip install -e .
```

For usage in Python:
```
from myopic_mces import MCES
MCES('CC(=O)OC1=CC=CC=C1C(=O)O', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C')
```

See [the PuLP documentation](https://coin-or.github.io/pulp/guides/how_to_configure_solvers.html) on how to configure ILP solvers. By default, the PuLP-provided COIN-OR solver will be used.

## Optional Arguments

General options
```
--threshold  int         Threshold for the comparison.
                         Exact distance is only calculated if the distance is lower than the threshold.
                         If set to -1 the exact disatnce is always calculated.

--solver string          Solver used for solving the ILP. Examples:'CPLEX_CMD', 'GUROBI_CMD', 'GLPK_CMD'

--num_jobs int           Number of jobs; instances to run in parallel.
                         By default this is set to the number of (logical) CPU cores.

--hdf5_mode              Input will be read from `input-file` in HDF5 format; output will be appended to this file.
                         See "Prepare HDF5 input" below
```

Options for the ILP solver
```
--solver_onethreaded    Limit ILP solver to one thread, likely resulting in faster
                        performance with parallel computations (not available for all solvers).

--solver_no_msg         Prevent solver from logging (not available for all solvers)

```

Experimental options for myopic MCES distance computation
```
--no_ilp_threshold           If set, do not add threshold as constraint to ILP,
                             resulting in longer runtimes and potential violations of the triangle equation.

--choose_bound_dynamically   If set, a potentially weaker but faster lower bound will be computed and used
                             when this bound is already greater than the threshold. By default (without
                             this option), always the strongest lower bound will be computed and used.
```

## Recommended settings

To speed up computations and save space, use the CPLEX solver, HDF5-mode (see [below](#prepare-hdf5-input)) and enable `--choose_bound_dynamically`:
```bash
PATH=$CPLEX_HOME/bin/x86-64_linux/:$PATH python -m myopic_mces.myopic_mces --threshold 10 --solver CPLEX_CMD --solver_onethreaded --solver_no_msg --hdf5_mode input-file.hdf5 tmpout
```

The `PATH`-variable has to be adapted to contain the directory of the CPLEX executable (see [the PuLP documentation](https://coin-or.github.io/pulp/guides/how_to_configure_solvers.html#cplex)).

## Dependencies and installation

Python packages required are:
```
rdkit(==2022.09.5)
networkx(==3.0)
pulp(==2.7.0)
scipy(==1.10.1)
joblib(==1.2.0)
```
Version numbers in braces correspond to an exemplary tested configuration (under Python version 3.11.0).
The program can be run on any standard operating system, tested on Windows 10 64 bit and Arch-Linux@linux-6.2.7 64 bit.

The recommended method of installation is via `pip`.
```bash
pip install myopic_mces
```

Dependencies can also be installed via [conda](https://docs.conda.io/en/latest/miniconda.html) or [mamba](https://github.com/mamba-org/mamba):
Download this repository, navigate to the download location and execute the following commands (replacing `conda` with `mamba` when using mamba):
```bash
conda env create -f conda_env.yml
# to activate the created enironment:
conda activate myopic_mces
```

A typical installation time should not exceed 5 minutes, mostly depending on the internet connection speed to download all required packages.

## Example data

The example provided in [example/example_data.csv](example/example_data.csv) can be run with:

```bash
pip install myopic-mces
myopic_mces example/example_data.csv example/example_data_out.csv
```

Typical runtime is about 10s on Windows 10 with all default options, running on 4 cores with 8GB RAM. Exemplary output is provided in [example/example_data_out.csv](example/example_data_out.csv).

## Utilities

### [`prepare_input.py`](src/myopic_mces/prepare_input.py)

For big datasets it is recommended to divide the input into batches, which can be done with this script:

```bash
python -m myopic_mces.prepare_input input-file.csv output-folder/ --batch_size 50_000_000
```

`input-file.csv` has to be formatted as shown above. This creates `output_folder` with the subdirectory `data` puts all batches (named `batch$i.csv`) inside.

#### Prepare HDF5 input

To conserve space, input for myopic MCES computation can now be provided as a HDF5-file containing the following "Datasets":

- `smiles`: list of all unique SMILES
- `computation_indices`: matrix with the shape (n, 3), each row representing one instance to compute. The first column contains the index of the pair, the second and third the indices of the two SMILES, respectively

Example:
```
smiles = [smiles_a, smiles_b, smiles_c]
computation_indices = [[0, 0, 1], # computes MCES for a vs. b
                       [1, 0, 2], # computes MCES for a vs. c
                       [2, 1, 2], # computes MCES for b vs. c
]
```

Instead of preparing the HDF5-file manually, with the additional option to create batches, [`prepare_input.py`](src/myopic_mces/prepare_input.py) can be used in HDF5-mode:

```bash
python -m myopic_mces.prepare_input --batch_size 50_000_000 --hdf5_mode input-smiles.txt output-folder/
```

Created batches (`batch$i.hdf5`) are written directly to the `output-folder`.

### [`combine_hdf5_batches.py`](src/myopic_mces/combine_hdf5_batches.py)

Batches results can be combined into a single HDF5 file with this script:
```bash
python -m myopic_mces.combine_hdf5_batches --out output-file batches/batch*.hdf5
```

To get a square matrix of MCES distances for convenience within Python from this file, simply do:
```python
import h5py
from scipy.spatial.distance import squareform

with h5py.File(hdf5outputfile, 'r') as f:
    mces_square = squareform(f['mces'])
    smiles = f['mces_smiles_order'] # column/row labels
```

### [`filter_dataset.py`](src/myopic_mces/filter_dataset.py)

If you want to filter datasets for similar structures in a database, you can use this script. This can speed up computations considerably, as for each query structure computations are stopped, when a "match" is found.

```bash
python filter_dataset.py --input input-file.json --out output-file --threshold 10
```

Input is provided in the following format (JSON file):
```
{query_smiles1: [db_smiles1, db_smiles2, ...],
 query_smiles2: [db_smiles1, db_smiles2, ...]}
```

The database can be different for each query SMILES, allow pre-filtering depending on the query.
