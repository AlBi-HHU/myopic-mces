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

import argparse
def MCES(ind,s1,s2,threshold,solver):
     start=time.time()
     G1,l1,e1=construct_graph(s1)
     G2,l2,e2=construct_graph(s2)       
     d=apply_filter(G1,G2,l1,l2,e1,e2,threshold,ind)
     if d>threshold:
        end=time.time()
        total_time=str(end-start)
        return ind,d,total_time,2   
         
     res=MCES_ILP(G1,l1,e1,G2,l2,e2,threshold,solver)
     end=time.time()
     total_time=str(end-start)
     return ind,res[0],total_time,res[1]
 
    
if __name__ == '__main__':
    
    parser=argparse.ArgumentParser()
    parser.add_argument("input",help="input file in the format: id,smile1,smile2")
    parser.add_argument("output",help="output file")
    parser.add_argument("--threshold",type=int,default=10,action="store",help="threshold for the distance")
    parser.add_argument("--solver", type=str,default="default",action="store",help="Solver for the ILP. example:GUROBI_CMD")
    args = parser.parse_args()

    threshold=args.threshold
    
    num_cores = multiprocessing.cpu_count()
    F=args.input
    F2=args.output
    f=open(F,"r")
    solver=args.solver
       
    inputs=[]
    for line in f:
            
        args=line.split(",")
        inputs.append(tuple([args[0],args[1],args[2]]))
    f.close()
    results = Parallel(n_jobs=num_cores)(delayed(MCES)(i[0],i[1],i[2],threshold,solver) for i in inputs)
        
    out=open(F2,"w")
    for i in results:
        out.write(i[0]+","+i[2]+","+str(i[1])+","+str(i[3])+"\n")
    out.close()
