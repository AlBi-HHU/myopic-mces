# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:17:41 2020

@author: seipp
"""
import pulp
import networkx as nx

def MCES_ILP(G1,l1,e1,G2,l2,e2):
    ILP=pulp.LpProblem("MCES", pulp.LpMaximize)
    nodepairs=[]
    for i in G1.nodes:
        for j in G2.nodes:
            nodepairs.append(tuple([i,j]))
    y=pulp.LpVariable.dicts('nodepairs', nodepairs, 
                            lowBound = 0,
                            upBound = 1,
                            cat = pulp.LpInteger)
    edgepairs=[]
    for i in G1.edges:
        for j in G2.edges:
            edgepairs.append(tuple([i,j]))
            c=pulp.LpVariable.dicts('edgepairs', edgepairs, 
                            lowBound = 0,
                            upBound = 1,
                            cat = pulp.LpInteger)
            
    h1=[] 
    h2=[]
    for i in edgepairs:
        if e1[i[0]]==2.0 and e2[i[1]]==2.0:
            h1.append(i)
        else:
            h2.append(i)
    ILP += pulp.lpSum([ 2*c[i] for i in h1])+pulp.lpSum([c[i] for i in h2])

    
    
    
    
    for i in G1.nodes:
        h=[]
        for j in G2.nodes:
            h.append(tuple([i,j]))
        ILP+=pulp.lpSum([y[k] for k in h])<=1
        
    for i in G2.nodes:
        h=[]
        for j in G1.nodes:
            h.append(tuple([j,i]))
        ILP+=pulp.lpSum([y[k] for k in h])<=1
        
    for i in G1.edges:
        ls=[]
        rs=[]
        for j in G2.edges:
            ls.append(tuple([i,j]))
        for j in G2.nodes:
            rs.append(tuple([i[0],j]))
        ILP+=pulp.lpSum([c[k] for k in ls])<=pulp.lpSum([y[k] for k in rs])
        
    for i in G2.edges:
        ls=[]
        rs=[]
        for j in G1.edges:
            ls.append(tuple([j,i]))
        for j in G1.nodes:
            rs.append(tuple([j,i[1]]))
        ILP+=pulp.lpSum([c[k] for k in ls])<=pulp.lpSum([y[k] for k in rs])
        
    for i in G1.nodes:
        for j in G2.edges:
            ls=[]
            for k in G1.neighbors(i):
                if tuple([tuple([i,k]),j]) in c:
                    ls.append(tuple([tuple([i,k]),j]))
                else:
                    ls.append(tuple([tuple([k,i]),j]))
            ILP+=pulp.lpSum([c[k] for k in ls])<=y[tuple([i,j[0]])]+y[tuple([i,j[1]])]
            
    for i in G2.nodes:
        for j in G1.edges:
            ls=[]
            for k in G2.neighbors(i):
                if tuple([j,tuple([i,k])]) in c:
                    ls.append(tuple([j,tuple([i,k])]))
                else:
                    ls.append(tuple([j,tuple([k,i])]))
            ILP+=pulp.lpSum([c[k] for k in ls])<=y[tuple([j[0],i])]+y[tuple([j[1],i])]
            
    for i in G1.nodes:
        for j in G2.nodes:
            if l1[i]!=l2[j]:
                ILP+=y[i,j]==0
    
    #for i in G1.edges:
    #    for j in G2.edges:
    #        if e1[i]!=e2[j]:
    #            ILP+=c[i,j]==0
    #n1=len(G1.edges)
    #n2=len(G2.edges)
    #ILP += pulp.lpSum([ c[i] for i in edgepairs])>=(n1+n2)/2-5
    
    solver=pulp.PULP_CBC_CMD(msg=False)
    ILP.solve(solver)
    if ILP.status==1:
        return int(ILP.objective.value())
    else:
        return -1