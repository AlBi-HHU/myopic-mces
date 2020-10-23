# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 17:59:05 2020

@author: seipp
"""

import networkx as nx

def filter1(G1,G2,l1,l2):
    
    atom_types1=[]
    for i in G1.nodes:
        if l1[i] not in atom_types1:
            atom_types1.append(l1[i])
    type_map1={}
    for i in atom_types1:
        type_map1[i]=list(filter(lambda x: i==l1[x],G1.nodes))    
        
    atom_types2=[]
    for i in G2.nodes:
        if l2[i] not in atom_types2:
            atom_types2.append(l2[i])
    type_map2={}
    for i in atom_types2:
        type_map2[i]=list(filter(lambda x: i==l2[x],G2.nodes))

    difference=0
    for i in atom_types1:
        if i in atom_types2:
            n=min(len(type_map1[i]),len(type_map2[i]))
            degreelist1=sorted(type_map1[i],key=lambda x:G1.degree[x])
            degreelist2=sorted(type_map2[i],key=lambda x:G2.degree[x])
            for j in range(n):
                difference+=G1.degree[degreelist1[j]]+G2.degree[degreelist2[j]]-2*min(G1.degree[degreelist1[j]],G2.degree[degreelist2[j]])
            if len(degreelist1)>n:
                for j in range(n,len(degreelist1)):
                    difference+=G1.degree[degreelist1[j]]
            if len(degreelist2)>n:
                for j in range(n,len(degreelist2)):
                    difference+=G2.degree[degreelist2[j]]
        else:
            for j in type_map1[i]:
                difference+=G1.degree[j]
    for i in atom_types2:
        if i not in atom_types1:
            for j in type_map2[i]:
                difference+=G2.degree[j]
    return difference

def apply_filter(G1,G2,l1,l2,e1,e2,threshold):
    d=filter1(G1,G2,l1,l2)/2
    if d<threshold:
        return True
    else:
        return False
    
            
                                              