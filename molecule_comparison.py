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
def MCES(ind,s1,s2,threshold):
     start=time.time()
     G1,l1,e1=construct_graph(s1)
     G2,l2,e2=construct_graph(s2)       
     
     if not apply_filter(G1,G2,l1,l2,e1,e2,threshold,ind):
        end=time.time()
        total_time=str(end-start)
        return ind,-1,total_time   
         
     res=MCES_ILP(G1,l1,e1,G2,l2,e2,threshold)
     end=time.time()
     total_time=str(end-start)
     return ind,res,total_time
 
    
if __name__ == '__main__':
    
    parser=argparse.ArgumentParser()
    parser.add_argument("input",help="input file in the format: id,smile1,smile2")
    parser.add_argument("output",help="output file")
    parser.add_argument("--threshold",type=int,default=10,action="store",help="threshold for the distance")
    args = parser.parse_args()

    threshold=args.threshold
    
    num_cores = multiprocessing.cpu_count()
    F=args.input
    F2=args.output
    f=open(F,"r")
       
    inputs=[]
    for line in f:
            
        args=line.split(",")
        inputs.append(tuple([args[0],args[1],args[2]]))
    f.close()
    results = Parallel(n_jobs=num_cores)(delayed(MCES)(i[0],i[1],i[2],threshold) for i in inputs)
        
    out=open(F2,"w")
    for i in results:
        out.write(i[0]+","+i[2]+","+str(i[1])+"\n")
    out.close()
