# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:18:56 2020

@author: seipp
"""
from rdkit import Chem
import networkx as nx

def construct_graph(s, save_mol=False, num_smiles=False):
    """
    Converts a SMILE into a graph

    Parameters
    ----------
    s : SMILES of molecule or rdkit molecule directly
    save_mol : whether to save the rdkit molecule as an attribute in the graph
    num_smiles : whether to save the numbered smiles as an attribute in the graph

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
    for bond in m.GetBonds():
        G.add_edge(bond.GetBeginAtom().GetIdx(),bond.GetEndAtom().GetIdx(),weight=bond.GetBondTypeAsDouble(),idx=bond.GetIdx())
    if save_mol:
        G.mol = m
    if num_smiles:
        # canonical false to keep atom map numbering
        G.num_smiles = Chem.MolToSmiles(m, canonical=False, allBondsExplicit=True)
    return G
