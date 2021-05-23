# molecule-comparison

# Dependencies

rdkit

networkx

pulp

# Usage

Input and Output file are in csv format. Every line in the input-file is one comparison:

Input-file: Index,Smile1,Smile2

Output-file: Index, Time(s), Difference, status (1 if exact distance, 2 if lower bound)

Run from commandline:
```sh
python molecule_comparison.py input-file output-file
```

# Optinal Arguments
```
--threshold  int    Threshold for the comparison. Exact distance is only calculated if the distance is lower than the threshold.
         	        If set to -1 the exact disatnce is always calculated.

--solver string     Solver used for solving the ILP. Examples:'CPLEX_CMD','GLPK_CMD','CPLEX_CMD','GUROBI_CMD'
```
