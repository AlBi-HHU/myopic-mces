from rdkit import DataStructs, Chem
from rdkit.Chem import MACCSkeys, Draw
import molecule_comparison

inp="./training_smiles.txt"
structs=[]
smiles=[]

with open(inp,"r+") as myfile: # read input
    data= myfile.read().split("\n")
    for line in data:
        structs.append(Chem.MolFromSmiles(line))
        smiles.append(line)

fps=[MACCSkeys.GenMACCSKeys(x) for x in structs]  # Generate MACCS fpts

for i in range(0,len(fps)-1):
    tanimoto = DataStructs.FingerprintSimilarity(fps[i],fps[i+1]) # Calc tanimoto
    mces = molecule_comparison.MCES(1,smiles[i], smiles[i+1], 5, "GLPK_CMD") # Calc MCES


    if mces[1]<=3 and tanimoto<0.7:    # Print out interesting pairs, for example "good" MCES and "bad" tanimoto
        print(str(mces[1])+" "+str(tanimoto)+" "+smiles[i]+" "+smiles[i+1])
       # im = Chem.Draw.MolsToGridImage( [structs[i],structs[i+1]]) # This draws the molecule pair, without this you can parallize
        #im.show()