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
    return difference/2

def get_cost(G1,G2,l1,l2,e1,e2,i,j):
    atom_types1=[]
    for k in G1.neighbors(i):
        if l1[k] not in atom_types1:
            atom_types1.append(l1[k])
    type_map1={}
    for k in atom_types1:
        type_map1[k]=list(filter(lambda x: k==l1[x],G1.neighbors(i))) 
   # print(atom_types1)    
        
    atom_types2=[]
    for k in G2.neighbors(j):
        if l2[k] not in atom_types2:
            atom_types2.append(l2[k])
    type_map2={}
    for k in atom_types2:
        type_map2[k]=list(filter(lambda x: k==l2[x],G2.neighbors(j)))
    #print(type_map1)
    #print(type_map2)
    difference=0.
    for k in atom_types1:
        if k in atom_types2:
            n=min(len(type_map1[k]),len(type_map2[k]))
            edgelist1=sorted(type_map1[k],key=lambda x:e1[tuple([x,i])],reverse=True)
            edgelist2=sorted(type_map2[k],key=lambda x:e2[tuple([x,j])],reverse=True)
            for l in range(n):
                difference+=(max(e1[tuple([edgelist1[l],i])],e2[tuple([edgelist2[l],j])])-min(e1[tuple([edgelist1[l],i])],e2[tuple([edgelist2[l],j])]))/2
            if len(edgelist1)>n:
                for l in range(n,len(edgelist1)):
                    difference+=e1[tuple([edgelist1[l],i])]/2
            if len(edgelist2)>n:
                for l in range(n,len(edgelist2)):
                    difference+=e2[tuple([edgelist2[l],j])]/2
        else:
            for l in type_map1[k]:
                difference+=e1[tuple([i,l])]/2
    for k in atom_types2:
        if k not in atom_types1:
            for l in type_map2[k]:
                difference+=e2[tuple([j,l])]/2
    return difference

def filter2(G1,G2,l1,l2,e1,e2):
    atom_types1=[]
    for i in G1.nodes:
        if l1[i] not in atom_types1:
            atom_types1.append(l1[i])
       
    atom_types2=[]
    for i in G2.nodes:
        if l2[i] not in atom_types2:
            atom_types2.append(l2[i])
    
    atom_types=atom_types1

    for i in atom_types2:
        if i not in atom_types:
            atom_types.append(i)
   
    res=0
    for i in atom_types:
        nodes1=list(filter(lambda x: i==l1[x],G1.nodes))
        nodes2=list(filter(lambda x: i==l2[x],G2.nodes))
        G=nx.Graph()
        for j in nodes1:
            G.add_node(tuple([1,j]))
        for j in nodes2:
            G.add_node(tuple([2,j]))
        for j in nodes1:
            for k in nodes2:
                if l1[j]==l2[k]:
                    G.add_edge(tuple([1,j]),tuple([2,k]),weight=get_cost(G1,G2,l1,l2,e1,e2,j,k))
        if len(nodes1)<len(nodes2):
            diff=len(nodes2)-len(nodes1)
            for j in range(1,diff+1):
                G.add_node(tuple([1,-j]))
                for k in nodes2:
                    G.add_edge(tuple([1,-j]),tuple([2,k]),weight=sum([e2[tuple([l,k])] for l in G2.neighbors(k)])/2)
        if len(nodes2)<len(nodes1):
            diff=len(nodes1)-len(nodes2)
            for j in range(1,diff+1):
                G.add_node(tuple([2,-j]))
                for k in nodes1:
                    G.add_edge(tuple([1,k]),tuple([2,-j]),weight=sum([e1[tuple([l,k])] for l in G1.neighbors(k)])/2)
        h=nx.bipartite.minimum_weight_full_matching(G)
        for k in h:
            if k[0]==1:
                res=res+G[k][h[k]]["weight"]
            
    return res

def apply_filter(G1,G2,l1,l2,e1,e2,threshold,ind):
    d=filter1(G1,G2,l1,l2)
    if d<=threshold:
        #print(ind)
        d=filter2(G1,G2,l1,l2,e1,e2)
        if d<=threshold:
            return d
    
    return d
    
            
                                              