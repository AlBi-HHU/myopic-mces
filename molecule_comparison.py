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

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Invalid Numbers of Arguments. Script will be terminated.')
    else:
        F=sys.argv[1]
        F2=sys.argv[2]
        f=open(F,"r")
        out=open(F2,"w")
        for line in f:
            start=time.time()
            args=line.split(",")
            #print(args)
            s1=args[1]
            s2=args[2]
            G1,l1,e1=construct_graph(s1)
            G2,l2,e2=construct_graph(s2)
            res=MCES_ILP(G1,l1,e1,G2,l2,e2)
            
            n=0
            for i in G1.edges:
                if e1[i]==2.0:
                    n=n+2
                else:
                    n=n+1
            for i in G2.edges:
                if e2[i]==2.0:
                    n=n+2
                else:
                    n=n+1
            
            
            end=time.time()
            total_time=str(end-start)
            if res==-1:
                out.write(args[0]+","+total_time+","+str(-1)+"\n")
            else:
                out.write(args[0]+","+total_time+","+str(n-2*res)+"\n")
            print(args[0])
        f.close()
        out.close()