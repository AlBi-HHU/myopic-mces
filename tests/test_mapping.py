from rdkit import Chem
from rdkit.Chem.Draw import MolsMatrixToGridImage, rdMolDraw2D
import pandas as pd
from myopic_mces.graph import construct_graph
from myopic_mces.MCES_ILP import construct_ILP, add_MCES_to_molgraphs
from myopic_mces import MCES
import ast
import pytest
import json
from joblib import Parallel, delayed

def calc_mapping(smiles1, smiles2):
    """
    Calculates the MCEs structure mapping between two molecules.
    """
    G1 = construct_graph(smiles1, save_mol=True, num_smiles=True)
    G2 = construct_graph(smiles2, save_mol = True, num_smiles=True)
    ILP = construct_ILP(G1, G2, threshold=-1)
    ILP.solve()
    mapping = add_MCES_to_molgraphs(ILP, G1, G2)
    return G1.num_smiles, G2.num_smiles, mapping


@pytest.mark.parametrize('smiles1, smiles2, mapping, name',
                            [("C-[C:1](=[O:2])-[c:3]1:[cH:4]:[cH:5]:[c:6](-[I:8]):[s:7]:1", \
                            "C-[C:1]1(-[CH3:17])-[CH:2]2-[CH2:3]-[CH:4]=[C:5]3-[CH2:6]-[O:7]-[CH:8]4-[CH:9]-3-[C:10]-2(-[CH2:11]-[CH2:12]-[C:13]-1=[O:14])-[CH2:15]-[O:16]-4", \
                            "{(0,1):(12,13),(1,2):(13,14),(1,3):(1,13),(3,4):(1,2),(4,5):(2,10),(5,6):(10,15)}", \
                            "id84972339"),

                            ("c1:[cH:1]:[cH:2]:[c:3](-[CH2:6]-[CH2:7]-[C:8](=[N:9]-[NH:10]-[c:11]2:[c:12](-[N+:20](=[O:21])-[O-:22]):[cH:13]:[c:14](-[N+:17](=[O:18])-[O-:19]):[cH:15]:[cH:16]:2)-[CH:23]=[CH:24]-[CH:25]=[CH:26]-[c:27]2:[cH:28]:[cH:29]:[cH:30]:[cH:31]:[cH:32]:2):[cH:4]:[cH:5]:1", \
                            "c1:[cH:1]:[cH:2]:[c:3](-[C:6](=[O:7])-[NH:8]-[c:9]2:[cH:10]:[cH:11]:[cH:12]:[c:13]3:[c:14]:2-[C:15](=[O:16])-[c:17]2:[c:18](:[c:21](-[NH:25]-[c:26]4:[n:27]:[c:28](-[NH:58]-[c:59]5:[cH:60]:[cH:61]:[cH:62]:[c:63]6:[c:64]:5-[C:65](=[O:66])-[c:67]5:[c:68](:[c:71](-[NH:75]-[C:76](=[O:77])-[c:78]7:[cH:79]:[cH:80]:[cH:81]:[cH:82]:[cH:83]:7):[cH:72]:[cH:73]:[cH:74]:5)-[C:69]-6=[O:70]):[n:29]:[c:30](-[NH:32]-[c:33]5:[cH:34]:[cH:35]:[cH:36]:[c:37]6:[c:38]:5-[C:39](=[O:40])-[c:41]5:[c:42](:[c:45](-[NH:49]-[C:50](=[O:51])-[c:52]7:[cH:53]:[cH:54]:[cH:55]:[cH:56]:[cH:57]:7):[cH:46]:[cH:47]:[cH:48]:5)-[C:43]-6=[O:44]):[n:31]:4):[cH:22]:[cH:23]:[cH:24]:2)-[C:19]-3=[O:20]):[cH:4]:[cH:5]:1", \
                            "{(0,1):(9,10),(0,5):(10,11),(1,2):(9,14),(11,12):(59,64),(11,16):(59,60),(12,13):(63,64),(13,14):(62,63),(14,15):(61,62),(15,16):(60,61),(2,3):(13,14),(23,24):(22,23),(24,25):(23,24),(25,26):(17,24),(27,28):(4,5),(27,32):(0,5),(28,29):(3,4),(29,30):(2,3),(3,4):(12,13),(3,6):(13,19),(31,32):(0,1),(4,5):(11,12),(6,7):(18,19),(7,8):(18,21),(8,23):(21,22),(8,9):(21,25)}", \
                            "id29332577")])
def test_draw_struct(smiles1, smiles2, mapping, name):
    """
    Draws MCES structure for two numbered SMILES strings and their mapping.
    Renumbers atoms based on the atommap encoded in the SMILES metadata.
    """
    raw_mol1 = Chem.MolFromSmiles(smiles1)
    raw_mol2 = Chem.MolFromSmiles(smiles2)

    dopts = rdMolDraw2D.MolDrawOptions()
    dopts.addAtomIndices = True
    dopts.setHighlightColour((.53, .87, .53, 1))
    dopts.noAtomLabels = True

    svg = MolsMatrixToGridImage(
                        [[raw_mol1, raw_mol2]],
                        drawOptions=dopts,
                        subImgSize=(400,400),
                        useSVG=True,
                        )

    svg_all_white = svg.replace('<svg ', '<svg style="background-color: white;" ')
    with open(f'testdata/raw_struct_{name}.svg', 'w') as out:
        out.write(svg_all_white)

    # renumber them with their atommap
    order1 = [atom.GetIdx() for atom in sorted(raw_mol1.GetAtoms(), key=lambda x: x.GetAtomMapNum())]
    mol1 = Chem.RenumberAtoms(raw_mol1, order1)
    
    order2 = [atom.GetIdx() for atom in sorted(raw_mol2.GetAtoms(), key=lambda x: x.GetAtomMapNum())]
    mol2 = Chem.RenumberAtoms(raw_mol2, order2)

    dopts = rdMolDraw2D.MolDrawOptions()
    dopts.addAtomIndices = True
    dopts.setHighlightColour((.53, .87, .53, 1))
    dopts.noAtomLabels = True

    # in highlight bonds, highlight only the ones that are in the MCES struct derived from add_MCES_to_molgraphs
    svg = MolsMatrixToGridImage(
                        [[mol1, mol2]],
                        drawOptions=dopts,
                        subImgSize=(400,400),
                        useSVG=True,
                        )

    svg_all_white = svg.replace('<svg ', '<svg style="background-color: white;" ')
    with open(f'testdata/reordered_struct_{name}.svg', 'w') as out:
        out.write(svg_all_white)

    G1_edges = []
    G2_edges = []
    mapping = ast.literal_eval(str(mapping))

    for atom_pair_G1 in mapping.keys():
        G1_edges.append(mol1.GetBondBetweenAtoms(*atom_pair_G1).GetIdx())

    for atom_pair_G2 in mapping.values():
        G2_edges.append(mol2.GetBondBetweenAtoms(*atom_pair_G2).GetIdx())
        
    dopts = rdMolDraw2D.MolDrawOptions()
    dopts.addAtomIndices = True
    dopts.setHighlightColour((.53, .87, .53, 1))
    dopts.noAtomLabels = True

    # in highlight bonds, highlight only the ones that are in the MCES struct derived from add_MCES_to_molgraphs
    svg = MolsMatrixToGridImage(
                        [[mol1, mol2]],
                        drawOptions=dopts,
                        subImgSize=(400,400),
                        highlightBondListsMatrix=[[G1_edges, G2_edges]],
                        useSVG=True,
                        )

    svg_all_white = svg.replace('<svg ', '<svg style="background-color: white;" ')
    with open(f'testdata/mces_struct_{name}.svg', 'w') as out:
        out.write(svg_all_white)

def test_mapping_reader():
    # aufruf auf MCES
    data = pd.read_csv("testdata/test_smiles.csv", header=None, names=['index', 'smiles1', 'smiles2'], nrows=10)
    reference = pd.read_csv("testdata/test_struct.csv", header=None, names=['i', 'dist', 't', 'mode', 'mapping', 'num_smiles1', 'num_smiles2'], keep_default_na=False, nrows=10)

    solver_options = {"timeLimit": 10, "msg": False}

    results = Parallel(n_jobs=-2)(
        delayed(MCES)(
            s1,
            s2,
            solver="COIN_CMD",
            threshold=-1,
            structure=True,
            num_smiles=True,
            solver_options=solver_options,
        )
        for s1, s2 in zip(data.smiles1, data.smiles2)
    )

    for result, ref in zip(results, reference.itertuples(index=False)):
        (
            i,
            dist,
            t,
            mode,
            mapping,
            num_smiles1,
            num_smiles2,
        ) = result

    def parse_mapping(x):
        return None if x == "" else ast.literal_eval(x)

    reference["mapping"] = reference["mapping"].map(parse_mapping)

    for result, ref in zip(results, reference.itertuples(index=False)):
        _, _, _, _, mapping, n1, n2 = result

        assert mapping == ref.mapping
        assert n1 == ref.num_smiles1
        assert n2 == ref.num_smiles2