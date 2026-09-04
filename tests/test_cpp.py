from rdkit import Chem

import mmces_filters

def test_cpp():
    mol1 = Chem.MolFromSmiles("CCO")
    mol2 = Chem.MolFromSmiles("CCC")

    res0 = mmces_filters.filter0(mol1, mol2)
    res1 = mmces_filters.filter1(mol1, mol2)
    res2 = mmces_filters.filter2(mol1, mol2)

    print(f"Ergebnis Filter 0: {res0}")
    print(f"Ergebnis Filter 1: {res1}")
    print(f"Ergebnis Filter 2: {res2}")
