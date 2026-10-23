# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:18:56 2020

@author: seipp
"""
from rdkit import Chem
import networkx as nx

CHIRAL_TAGS = {Chem.ChiralType.CHI_TETRAHEDRAL_CW: 'CW', Chem.ChiralType.CHI_TETRAHEDRAL_CCW: 'CCW'}
TRANS = (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOTRANS)
CIS = (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOCIS)

def construct_graph(s, save_mol=False, num_smiles=False, stereo=False):
    """
    Converts a SMILE into a graph

    Parameters
    ----------
    s : SMILES of molecule or rdkit molecule directly
    save_mol : whether to save the rdkit molecule as an attribute in the graph
    num_smiles : whether to save the numbered smiles as an attribute in the graph
    stereo : whether to save stereo information as attributes (None if unknown/unspecified):
             node attributes "chiral_tag" (CW/CCW) and "chiral_nbrs" (neighbour order the tag refers to),
             edge attribute "ez" (begin atom, end atom, stereo atoms, whether the stereo atoms are trans)

    Returns:
    -------
    networkx.classes.graph.Graph
        Graph that represents the molecule.
        The bond types are represented as edge weights.
        The atom types are represented as atom attributes of the nodes.
        The molecules are represented as rdkit molecules.
        The numbered_smiles are represented as SMILES with atommap metadata.

    """
    if (isinstance(s, Chem.rdchem.Mol)):
        m = s
    else:
        m = Chem.MolFromSmiles(s)
    # convert the molecule into a graph
    # The bond and atom types are converted to node/edge attributes
    G=nx.Graph()
    for atom in m.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
        G.add_node(atom.GetIdx(),atom=atom.GetSymbol())    
        if stereo:
            # chiral tag relative to the neighbour order of the atom,
            # implicit H / lone pair as placeholder (its position does not matter for comparisons)
            tag = CHIRAL_TAGS.get(atom.GetChiralTag())
            G.nodes[atom.GetIdx()]['chiral_tag'] = tag
            if tag is not None:
                nbrs = [b.GetOtherAtomIdx(atom.GetIdx()) for b in atom.GetBonds()]
                G.nodes[atom.GetIdx()]['chiral_nbrs'] = tuple(nbrs + ['H'] * (4 - len(nbrs)))
    for bond in m.GetBonds():
        G.add_edge(bond.GetBeginAtom().GetIdx(),bond.GetEndAtom().GetIdx(),weight=bond.GetBondTypeAsDouble(),idx=bond.GetIdx())
        if stereo:
            # (begin atom, end atom, stereo atom begin, stereo atom end, stereo atoms trans?)
            ez = None
            if bond.GetStereo() in TRANS + CIS:
                sa, sb = bond.GetStereoAtoms()
                ez = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), sa, sb, bond.GetStereo() in TRANS)
            G.edges[bond.GetBeginAtom().GetIdx(),bond.GetEndAtom().GetIdx()]['ez'] = ez
    if save_mol:
        G.mol = m
    if num_smiles:
        # canonical false to keep atom map numbering
        G.num_smiles = Chem.MolToSmiles(m, canonical=False, allBondsExplicit=True)
    return G
