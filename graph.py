# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:18:56 2020

@author: seipp
"""
from rdkit import Chem
import networkx as nx

def construct_graph(s):
    m = Chem.MolFromSmiles(s)
    Chem.Kekulize(m)
    G=nx.Graph()
    l={}
    e={}
    for atom in m.GetAtoms():
        G.add_node(atom.GetIdx())
        l[atom.GetIdx()]=atom.GetSymbol()
    for bond in m.GetBonds():
        G.add_edge(bond.GetBeginAtom().GetIdx(),bond.GetEndAtom().GetIdx())
        e[tuple([bond.GetBeginAtom().GetIdx(),bond.GetEndAtom().GetIdx()])]=bond.GetBondTypeAsDouble()
        e[tuple([bond.GetEndAtom().GetIdx(),bond.GetBeginAtom().GetIdx()])]=bond.GetBondTypeAsDouble()
    return G,l,e