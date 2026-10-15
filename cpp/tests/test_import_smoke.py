from rdkit import Chem
import rdkit
import mmces_filters

mol1 = Chem.MolFromSmiles("CCO")
mol2 = Chem.MolFromSmiles("CCC")

print("RDKit version:", rdkit.__version__)
print("RDKit paths:", rdkit.__path__)
print(type(mol1))
print(mmces_filters.filter2.__doc__)

res0 = mmces_filters.filter0(mol1, mol2)
res1 = mmces_filters.filter1(mol1, mol2)
res2 = mmces_filters.filter2(mol1, mol2)

print(f"Result Filter 0: {res0}")
print(f"Result Filter 1: {res1}")
print(f"Result Filter 2: {res2}")
