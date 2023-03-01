# Computation of myopic MCES distances

## Usage

Input and Output file are in csv format. Every line in the input-file is one comparison:

input-file: `index,SMILES1,SMILES2`

output-file: `index,time taken,myopic MCES distance,status (1 if exact distance, 2/4 if lower bound)`

Run from the command line:
```bash
python molecule_comparison.py input-file output-file
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


## Dependencies

Python packages required are:
```
rdkit
networkx
pulp
scipy
joblib
```

A conda (or [mamba](https://github.com/mamba-org/mamba)) environment with all necessary packages installed can be created with
```bash
conda env create -f conda_env.yml
# to activate:
conda activate myopic_mces
```
