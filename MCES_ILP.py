# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:17:41 2020

@author: seipp
"""
import pulp
import networkx as nx

def MCES_ILP(G1,G2,threshold,solver):
    """ 
     Calculates the exact distance between two molecules using an ILP
     
     Parameters
     ----------
     G1 : networkx.classes.graph.Graph 
         Graph representing the first molecule.
     G2 : networkx.classes.graph.Graph 
         Graph representing the second molecule.
     threshold : int
         Threshold for the comparison. Exact distance is only calculated if the distance is lower than the threshold.
     solver: string
         ILP-solver used for solving MCES. Example:GUROBI_CMD
         
     Returns:
     -------
     float
         Distance between the molecules
     int
         Type of Distance:
             1 : Exact Distance
             2 : Lower bound (If the exact distance is above the threshold)
     
    """
    
    ILP=pulp.LpProblem("MCES", pulp.LpMinimize)
    
    #Variables for nodepairs
    nodepairs=[]
    for i in G1.nodes:
        for j in G2.nodes:
            nodepairs.append(tuple([i,j]))
    y=pulp.LpVariable.dicts('nodepairs', nodepairs, 
                            lowBound = 0,
                            upBound = 1,
                            cat = pulp.LpInteger)
    #variables for edgepairs and weight
    edgepairs=[]
    w={}
    for i in G1.edges:
        for j in G2.edges:
            edgepairs.append(tuple([i,j]))
            w[tuple([i,j])]=max(G1[i[0]][i[1]]["weight"],G2[j[0]][j[1]]["weight"])-min(G1[i[0]][i[1]]["weight"],G2[j[0]][j[1]]["weight"])
    
    #variables for not mapping an edge
    for i in G1.edges:
        edgepairs.append(tuple([i,-1]))
        w[tuple([i,-1])]=G1[i[0]][i[1]]["weight"]
    for j in G2.edges:
        edgepairs.append(tuple([-1,j]))
        w[tuple([-1,j])]=G2[j[0]][j[1]]["weight"]
    c=pulp.LpVariable.dicts('edgepairs', edgepairs, 
                            lowBound = 0,
                            upBound = 1,
                            cat = pulp.LpInteger)
        
            
    #objective function
    ILP += pulp.lpSum([ w[i]*c[i] for i in edgepairs])

    
    
    
    #Every node in G1 can only be mapped to at most one in G2
    for i in G1.nodes:
        h=[]
        for j in G2.nodes:
            h.append(tuple([i,j]))
        ILP+=pulp.lpSum([y[k] for k in h])<=1
      
    #Every node in G1 can only be mapped to at most one in G1    
    for i in G2.nodes:
        h=[]
        for j in G1.nodes:
            h.append(tuple([j,i]))
        ILP+=pulp.lpSum([y[k] for k in h])<=1
    
    #Every edge in G1 has to be mapped to an edge in G2 or the variable for not mapping has to be 1
    for i in G1.edges:
        ls=[]
        rs=[]
        for j in G2.edges:
            ls.append(tuple([i,j]))
        for j in G2.nodes:
            rs.append(tuple([i[0],j]))
        ILP+=pulp.lpSum([c[k] for k in ls])+c[tuple([i,-1])]==1
     
    #Every edge in G2 has to be mapped to an edge in G1 or the variable for not mapping has to be 1    
    for i in G2.edges:
        ls=[]
        rs=[]
        for j in G1.edges:
            ls.append(tuple([j,i]))
        for j in G1.nodes:
            rs.append(tuple([j,i[1]]))
        ILP+=pulp.lpSum([c[k] for k in ls])+c[tuple([-1,i])]==1
    
    #The mapping of the edges has to match the mapping of the nodes    
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
    
    #nodes can only be mapped if they represent the same atom
    for i in G1.nodes:
        for j in G2.nodes:
            if G1.nodes[i]["atom"]!=G2.nodes[j]["atom"]:
                ILP+=y[i,j]==0
    
    #constraint for the threshold
    ILP +=pulp.lpSum([ w[i]*c[i] for i in edgepairs])<=threshold
    
    #solve the ILP
    if solver=="default":
        ILP.solve()
    else:
        sol=pulp.getSolver(solver)
        ILP.solve(sol)
    if ILP.status==1:
        return float(ILP.objective.value()),1
    else:
        return 10,2