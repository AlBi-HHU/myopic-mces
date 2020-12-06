# molecule-comparison

# Dependencies

rdkit

networkx

pulp

# Run

python molecule_comparison.py input-file output-file

Input-file: Index,Smile1,Smile2

Output-file: Index, Time(s), Difference, status (1 if exact distance, 2 if lower bound)
