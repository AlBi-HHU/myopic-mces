# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:16:05 2020

@author: seipp
"""
import sys
import time
from graph import construct_graph
from MCES_ILP import MCES_ILP
import pulp
import networkx as nx
from filter_MCES import apply_filter

from joblib import Parallel, delayed
import multiprocessing

def MCES(ind,s1,s2):
     start=time.time()
     G1,l1,e1=construct_graph(s1)
     G2,l2,e2=construct_graph(s2)       
     
     n=0
     for i in G1.edges:
         n=n+e1[i]
     for i in G2.edges:
         n=n+e2[i]
     if not apply_filter(G1,G2,l1,l2,e1,e2,10):
        end=time.time()
        total_time=str(end-start)
        return ind,-1,total_time,n    
         
     res=MCES_ILP(G1,l1,e1,G2,l2,e2,n)
     end=time.time()
     total_time=str(end-start)
     #print(ind)
     return ind,res,total_time,n
 
    
if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Invalid Numbers of Arguments. Script will be terminated.')
    else:
        num_cores = multiprocessing.cpu_count()
        F=sys.argv[1]
        F2=sys.argv[2]
        f=open(F,"r")
        
        inputs=[]
        for line in f:
            
            args=line.split(",")
            inputs.append(tuple([args[0],args[1],args[2]]))
        f.close()
        results = Parallel(n_jobs=num_cores)(delayed(MCES)(i[0],i[1],i[2]) for i in inputs)
        
        out=open(F2,"w")
        for i in results:
            if i[1]==-1:
                out.write(i[0]+","+i[2]+","+str(-1)+"\n")
            else:
                out.write(i[0]+","+i[2]+","+str(i[3]-2*i[1])+"\n")
        out.close()