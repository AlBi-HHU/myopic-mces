# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 00:53:19 2020

@author: seipp
"""
from pycryptosat import Solver
import networkx as nx
def MCES_SAT(G1,l1,e1,G2,l2,e2):
    res=-1
    diff=10
    for v in range(20,-1,-1):
        s = Solver()
        diff=v/2
        
        nodepairs=[]
        nodepair_dict={}
        index=1
        for i in G1.nodes:
            for j in G2.nodes:
                nodepairs.append(tuple([i,j]))
                nodepair_dict[tuple([i,j])]=index
                index+=1
                
                
        edgepairs=[]
        w={}
        edgedict={}        
        Q3_var={}
        Q4_var={}
        for i in G1.edges:
            for j in G2.edges:
                edgepairs.append(tuple([i,j]))
                w[tuple([i,j])]=max(e1[i],e2[j])-min(e1[i],e2[j])
                edgedict[tuple([i,j])]=index
                index+=1
                Q3_var[tuple([i,j])]=index
                index+=1
                Q4_var[tuple([i,j])]=index
                index+=1
                
        for i in G1.edges:
            edgepairs.append(tuple([i,-1]))
            w[tuple([i,-1])]=e1[i]
            edgedict[tuple([i,-1])]=index
            index+=1
        for j in G2.edges:
            edgepairs.append(tuple([-1,j]))
            w[tuple([-1,j])]=e2[j]
            edgedict[tuple([-1,j])]=index
            index+=1
        
        T_var={}
        Q1_var={}
        Q2_var={}

        for i in range(v+1):
            T_var[tuple([0,i/2])]=index
            s.add_clause([index])
            index+=1

        for i in range(v+1):
            for j in range(1,len(edgepairs)+1):
                T_var[tuple([j,i/2])]=index
                index+=1
                Q1_var[tuple([j,i/2])]=index
                index+=1
                Q2_var[tuple([j,i/2])]=index
                index+=1

            
        for i in G1.nodes:
            for j in range(len(G2.nodes)):
                for k  in range(j+1,len(G2.nodes)):
                    s.add_clause([-nodepair_dict[tuple([i,list(G2.nodes)[j]])],-nodepair_dict[tuple([i,list(G2.nodes)[k]])]])
    
        for i in G2.nodes:
            for j in range(len(G1.nodes)):
                for k  in range(j+1,len(G1.nodes)):
                    s.add_clause([-nodepair_dict[tuple([list(G1.nodes)[j],i])],-nodepair_dict[tuple([list(G1.nodes)[k],i])]])
        
        for i in G1.edges:
            h=[]
            for j in G2.edges:
                h.append(edgedict[tuple([i,j])])
            h.append(edgedict[tuple([i,-1])])
            s.add_clause(h)
            
        for i in G2.edges:
            h=[]
            for j in G1.edges:
                h.append(edgedict[tuple([j,i])])
            h.append(edgedict[tuple([-1,i])])
            s.add_clause(h)
            
        for i in G1.edges:
            for j in range(len(G2.edges)):
                for k  in range(j+1,len(G2.edges)):
                    s.add_clause([-edgedict[tuple([i,list(G2.edges)[j]])],-edgedict[tuple([i,list(G2.edges)[k]])]])
            
        for i in G2.edges:
            for j in range(len(G1.edges)):
                for k  in range(j+1,len(G1.edges)):
                    s.add_clause([-edgedict[tuple([list(G1.edges)[j],i])],-edgedict[tuple([list(G1.edges)[k],i])]])
        
        for i in G1.edges:
            for j in G2.edges:
                s.add_clause([-edgedict[i,j],nodepair_dict[tuple([i[0],j[0]])],nodepair_dict[tuple([i[0],j[1]])]])
                s.add_clause([-edgedict[i,j],nodepair_dict[tuple([i[1],j[0]])],nodepair_dict[tuple([i[1],j[1]])]])
                s.add_clause([edgedict[i,j],Q3_var[i,j],Q4_var[i,j]])
                
                s.add_clause([Q3_var[i,j],nodepair_dict[tuple([i[0],j[0]])],nodepair_dict[tuple([i[0],j[1]])]])
                s.add_clause([-Q3_var[i,j],-nodepair_dict[tuple([i[0],j[0]])]])
                s.add_clause([-nodepair_dict[tuple([i[0],j[1]])],-Q3_var[i,j]])
                
                s.add_clause([Q4_var[i,j],nodepair_dict[tuple([i[1],j[0]])],nodepair_dict[tuple([i[1],j[1]])]])
                s.add_clause([-Q4_var[i,j],-nodepair_dict[tuple([i[1],j[0]])]])
                s.add_clause([-nodepair_dict[tuple([i[1],j[1]])],-Q4_var[i,j]])                              
                                                                                      

                

                
        for i in range(v+1):
            s.add_clause([T_var[tuple([0,i/2])]])
        
        for i in range(1,len(edgepairs)+1):
            for j in range(v+1):
                s.add_clause([-T_var[tuple([i-1,j/2])],edgedict[edgepairs[i-1]],Q1_var[tuple([i,j/2])]])
                s.add_clause([T_var[tuple([i-1,j/2])],-Q1_var[tuple([i,j/2])]])
                s.add_clause([-edgedict[edgepairs[i-1]],-Q1_var[tuple([i,j/2])]]) 
                h=[]
                h.append(-T_var[tuple([i,j/2])])
                h.append(Q1_var[tuple([i,j/2])])
                s.add_clause([-Q1_var[tuple([i,j/2])],T_var[tuple([i,j/2])]])
                if j>0:
                    h.append(T_var[tuple([i,(j-1)/2])])
                    s.add_clause([T_var[tuple([i,j/2])],-T_var[tuple([i,(j-1)/2])]])
                if j/2>=w[edgepairs[i-1]]:
                    s.add_clause([-T_var[tuple([i-1,j/2-w[edgepairs[i-1]]])],-edgedict[edgepairs[i-1]],Q2_var[tuple([i,j/2])]])
                    s.add_clause([T_var[tuple([i-1,j/2-w[edgepairs[i-1]]])],-Q2_var[tuple([i,j/2])]])
                    s.add_clause([edgedict[edgepairs[i-1]],-Q2_var[tuple([i,j/2])]])
                    h.append(Q2_var[tuple([i,j/2])])
                    s.add_clause([-Q2_var[tuple([i,j/2])],T_var[tuple([i,j/2])]])
                s.add_clause(h) 
                
        for i in G1.nodes:
            for j in G2.nodes:
                if l1[i]!=l2[j]:
                    s.add_clause([-nodepair_dict[i,j]])
        h=[]
        for i in range(int(2*diff)+1):
            h.append(T_var[tuple([len(edgepairs),i/2])])
        s.add_clause(h)
        
        sat, solution = s.solve()
        
        if not sat:
            return res
        else:
            res=diff
    return 0