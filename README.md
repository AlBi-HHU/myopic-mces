# Computation of myopic MCES distances

Implementation of myopic MCES distance computation, see the preprint at [doi:10.1101/2023.03.27.534311](https://doi.org/10.1101/2023.03.27.534311) for details.

## Usage

Input and Output file are in csv format. Every line in the input-file is one comparison:

input-file: `index,SMILES1,SMILES2`

output-file: `index,time taken,myopic MCES distance,status (1 if exact distance, 2/4 if lower bound)`

Run from the command line:
```bash
python myopic_mces.py input-file output-file
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

The recommended method of installation is via [conda](https://docs.conda.io/en/latest/miniconda.html) or [mamba](https://github.com/mamba-org/mamba):
Download this repository, navigate to the download location and execute the following commands (replacing `conda` with `mamba` when using mamba):
```bash
conda env create -f conda_env.yml
# to activate the created enironment:
conda activate myopic_mces
```

A PyPI-package is also available, install via:
```bash
pip install myopic_mces
```

A typical installation time should not exceed 5 minutes, mostly depending on the internet connection speed to download all required packages.

## Example data

The example provided in [example/example_data.csv](example/example_data.csv) can be run with:

```bash
python myopic_mces.py example/example_data.csv example/example_data_out.csv
```

Alternatively, if the package was installed via pip:
```bash
myopic_mces example/example_data.csv example/example_data_out.csv
```

Typical runtime is about 10s on Windows 10 with all default options, running on 4 cores with 8GB RAM. Exemplary output is provided in [example/example_data_out.csv](example/example_data_out.csv).
