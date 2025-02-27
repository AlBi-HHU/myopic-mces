# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:18:56 2020

@author: seipp
"""
from rdkit import Chem
import networkx as nx

def construct_graph(s):
    """
    Converts a SMILE into a graph

    Parameters
    ----------
    s : SMILES of molecule or rdkit molecule directly

    Returns:
    -------
    networkx.classes.graph.Graph
        Graph that represents the molecule.
        The bond types are represented as edge weights.
        The atom types are represented as atom attributes of the nodes.
    """
    if (isinstance(s, str)):
        #read the smile
        m = Chem.MolFromSmiles(s)
    else:
        m = s
    # convert the molecule into a graph
    # The bond and atom types are converted to node/edge attributes
    G=nx.Graph()
    for atom in m.GetAtoms():
        G.add_node(atom.GetIdx(),atom=atom.GetSymbol())
    for bond in m.GetBonds():
        G.add_edge(bond.GetBeginAtom().GetIdx(),bond.GetEndAtom().GetIdx(),weight=bond.GetBondTypeAsDouble())
    return G
