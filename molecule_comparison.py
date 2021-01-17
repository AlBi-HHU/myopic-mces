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
     """ 
     Calculates the distance between two molecules
     
     Parameters
     ----------
     s1 : str 
         SMILE of the first molecule
     s2 : str 
         SMILE of the second molecule
     threshold : int
         Threshold for the comparison. Exact distance is only calculated if the distance is lower than the threshold
     solver: string
         ILP-solver used for solving MCES. Example:GUROBI_CMD
         
     Returns:
     -------
     float
         Distance between the molecules
     float
         Time taken for the calculation
     int
         Type of Distance:
             1 : Exact Distance
             2 : Lower bound (If the exact distance is above the threshold)
     
     """
     start=time.time()
     #construct graph for both smiles. 
     G1=construct_graph(s1)
     G2=construct_graph(s2)
     # filter out if distance is above the threshold 
     d=apply_filter(G1,G2,threshold)
     if d>threshold:
        end=time.time()
        total_time=str(end-start)
        return ind,d,total_time,2   
     # calculate MCES
     res=MCES_ILP(G1,G2,threshold,solver)
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