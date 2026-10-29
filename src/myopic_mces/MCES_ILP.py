# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:17:41 2020

@author: seipp
"""
import pulp
import networkx as nx
from myopic_mces.filter_MCES import ComputationMode
import re

def construct_ILP(G1, G2, threshold, no_ilp_threshold=False):
    """
     Constructs the ILP for calculating the exact distance between two molecules

     Parameters
     ----------
     G1 : networkx.classes.graph.Graph
         Graph representing the first molecule.
     G2 : networkx.classes.graph.Graph
         Graph representing the second molecule.
     threshold : float
         Threshold for the comparison. Exact distance is only calculated if the distance is lower than the threshold.
     no_ilp_threshold: bool
         if true, always return exact distance even if it is below the threshold (slower)

     Returns:
     -------
     pulp.LpProblem
         ILP
    """

    ILP=pulp.LpProblem("MCES", pulp.LpMinimize)

    #Variables for nodepairs
    nodepairs=[]
    for i in G1.nodes:
        for j in G2.nodes:
            if G1.nodes[i]["atom"]==G2.nodes[j]["atom"]:
                nodepairs.append(tuple([i,j]))
    y=ILP.add_variable_dicts('nodepairs', nodepairs,
                               lowBound = 0,
                               upBound = 1,
                               cat = pulp.LpInteger)
    #variables for edgepairs and weight
    edgepairs=[]
    w={}
    for i in G1.edges:
        for j in G2.edges:
            if (G1.nodes[i[0]]["atom"]==G2.nodes[j[0]]["atom"] and G1.nodes[i[1]]["atom"]==G2.nodes[j[1]]["atom"]) or (G1.nodes[i[1]]["atom"]==G2.nodes[j[0]]["atom"] and G1.nodes[i[0]]["atom"]==G2.nodes[j[1]]["atom"]):
                edgepairs.append(tuple([i,j]))
                w[tuple([i,j])]=max(G1[i[0]][i[1]]["weight"],G2[j[0]][j[1]]["weight"])-min(G1[i[0]][i[1]]["weight"],G2[j[0]][j[1]]["weight"])

    #variables for not mapping an edge
    for i in G1.edges:
        edgepairs.append(tuple([i,-1]))
        w[tuple([i,-1])]=G1[i[0]][i[1]]["weight"]
    for j in G2.edges:
        edgepairs.append(tuple([-1,j]))
        w[tuple([-1,j])]=G2[j[0]][j[1]]["weight"]
    c=ILP.add_variable_dicts('edgepairs', edgepairs,
                               lowBound = 0,
                               upBound = 1,
                               cat = pulp.LpInteger)


    #objective function
    ILP += pulp.lpSum([ w[i]*c[i] for i in edgepairs])




    #Every node in G1 can only be mapped to at most one in G2
    for i in G1.nodes:
        h=[]
        for j in G2.nodes:
            if G1.nodes[i]["atom"]==G2.nodes[j]["atom"]:
                h.append(tuple([i,j]))
        ILP+=pulp.lpSum([y[k] for k in h])<=1

    #Every node in G1 can only be mapped to at most one in G1
    for i in G2.nodes:
        h=[]
        for j in G1.nodes:
            if G1.nodes[j]["atom"]==G2.nodes[i]["atom"]:
                h.append(tuple([j,i]))
        ILP+=pulp.lpSum([y[k] for k in h])<=1

    #Every edge in G1 has to be mapped to an edge in G2 or the variable for not mapping has to be 1
    for i in G1.edges:
        ls=[]
        rs=[]
        for j in G2.edges:
            if (G1.nodes[i[0]]["atom"]==G2.nodes[j[0]]["atom"] and G1.nodes[i[1]]["atom"]==G2.nodes[j[1]]["atom"]) or (G1.nodes[i[1]]["atom"]==G2.nodes[j[0]]["atom"] and G1.nodes[i[0]]["atom"]==G2.nodes[j[1]]["atom"]):
                ls.append(tuple([i,j]))
        ILP+=pulp.lpSum([c[k] for k in ls])+c[tuple([i,-1])]==1

    #Every edge in G2 has to be mapped to an edge in G1 or the variable for not mapping has to be 1
    for i in G2.edges:
        ls=[]
        rs=[]
        for j in G1.edges:
            if (G1.nodes[j[0]]["atom"]==G2.nodes[i[0]]["atom"] and G1.nodes[j[1]]["atom"]==G2.nodes[i[1]]["atom"]) or (G1.nodes[j[1]]["atom"]==G2.nodes[i[0]]["atom"] and G1.nodes[j[0]]["atom"]==G2.nodes[i[1]]["atom"]):
                ls.append(tuple([j,i]))
        ILP+=pulp.lpSum([c[k] for k in ls])+c[tuple([-1,i])]==1

    #The mapping of the edges has to match the mapping of the nodes
    for i in G1.nodes:
        for j in G2.edges:
            ls=[]
            for k in G1.neighbors(i):
                if tuple([tuple([i,k]),j]) in c:
                    ls.append(tuple([tuple([i,k]),j]))
                else:
                    if  tuple([tuple([k,i]),j]) in c:
                        ls.append(tuple([tuple([k,i]),j]))
            rs=[]
            if G1.nodes[i]["atom"]==G2.nodes[j[0]]["atom"]:
                rs.append(tuple([i,j[0]]))
            if G1.nodes[i]["atom"]==G2.nodes[j[1]]["atom"]:
                rs.append(tuple([i,j[1]]))
            ILP+=pulp.lpSum([c[k] for k in ls])<=pulp.lpSum([y[k] for k in rs])


    for i in G2.nodes:
        for j in G1.edges:
            ls=[]
            for k in G2.neighbors(i):
                if tuple([j,tuple([i,k])]) in c:
                    ls.append(tuple([j,tuple([i,k])]))
                else:
                    if tuple([j,tuple([k,i])]) in c:
                        ls.append(tuple([j,tuple([k,i])]))
            rs=[]
            if G2.nodes[i]["atom"]==G1.nodes[j[0]]["atom"]:
                rs.append(tuple([j[0],i]))
            if G2.nodes[i]["atom"]==G1.nodes[j[1]]["atom"]:
                rs.append(tuple([j[1],i]))
            ILP+=pulp.lpSum([c[k] for k in ls])<=pulp.lpSum(y[k] for k in rs)

    #constraint for the threshold
    if threshold!=-1 and not no_ilp_threshold:
        ILP +=pulp.lpSum([ w[i]*c[i] for i in edgepairs])<=threshold

    return ILP


# only for these solvers, a time limit without solution can be told apart from an infeasible ILP
# (e.g. CPLEX_CMD reports both as infeasible, which would wrongly be returned as above threshold)
TIME_LIMIT_SOLVERS = ('COIN_CMD', 'CPLEX_PY')
# distances are multiples of 0.5 (bond orders), so a solution within these gaps is still exact
# (for the relative gap: as long as the distance is below 5000)
MAX_GAP_ABS = 0.5
MAX_GAP_REL = 1e-4

def check_solver_options(solver, solver_options):
    """
    Raises a ValueError for solver options with which the returned computation modes could be wrong.
    """
    if solver_options.get('timeLimit') is not None and solver not in TIME_LIMIT_SOLVERS:
        raise ValueError(f'A time limit is only supported for the solvers {", ".join(TIME_LIMIT_SOLVERS)} (got {solver}). '
                         'With other solvers, a time limit without solution can be mistaken for a distance above the threshold. '
                         'Switch to a supported solver or remove the time limit.')
    if (solver_options.get('gapAbs') or 0) >= MAX_GAP_ABS:
        raise ValueError(f'gapAbs has to be below {MAX_GAP_ABS}, otherwise solutions reported as exact may not be exact.')
    if (solver_options.get('gapRel') or 0) > MAX_GAP_REL:
        raise ValueError(f'gapRel must not be above {MAX_GAP_REL}, otherwise solutions reported as exact may not be exact.')


def MCES_ILP(G1, G2, threshold, solver='COIN_CMD', solver_options={}, no_ilp_threshold=False, structure=False):
    """
     Calculates the exact distance between two molecules using an ILP

     Parameters
     ----------
     G1 : networkx.classes.graph.Graph
         Graph representing the first molecule.
     G2 : networkx.classes.graph.Graph
         Graph representing the second molecule.
     threshold : float
         Threshold for the comparison. Exact distance is only calculated if the distance is lower than the threshold.
     solver: string
         ILP-solver used for solving MCES. Example: CPLEX_PY
     solver_options: dict
         additional options to pass to solvers. Example: threads=1, msg=False for better multi-threaded performance
     no_ilp_threshold: bool
         if true, always return exact distance even if it is below the threshold (slower)
     structure: bool
         if true, calculate pairwise MCES structures

     Returns:
     -------
     float
         Distance between the molecules
     int
         Type of Distance:
             1 : Exact Distance
             2 : Lower bound (If the exact distance is above the threshold)
     (dict 
         Mapping of MCES structure, in the format (atom1_mol1, atom2_mol1): (atom1_mol2, atom2_mol2) to describe that the bond
         between atom1 and atom2 in molecule1 was matched to the bond of atom1 and atom2 in molecule2.
         Only returned if structure is set to True.)

    """

    check_solver_options(solver, solver_options)
    ILP = construct_ILP(G1, G2, threshold=threshold, no_ilp_threshold=no_ilp_threshold)

    # if structure is set, only calculate it if the returned mode is exact
    def result(dist, mode):
        if structure:
            if mode == ComputationMode.EXACT.value:
                mapping = add_MCES_to_molgraphs(ILP, G1, G2)
                return dist, mode, mapping
            else:
                return dist, mode, None
        return dist, mode

    # solve the ILP and return the correct computation mode
    sol=pulp.getSolver(solver, **solver_options)
    stats = ILP.solve(sol)
    S = pulp.LpSolveStatus

    # "Optimal (within gap tolerance)" (GapLimit) was reported as optimal before PuLP 4
    if stats.status in (S.Optimal, S.GapLimit) and stats.has_solution:
        return result(float(ILP.objective.value()), ComputationMode.EXACT.value)
    elif stats.status == S.TimeLimit:
        if stats.has_solution:
            # hit time limit, but still found a solution
            return result(float(ILP.objective.value()), ComputationMode.TIMEOUT_ILP_UNPROVEN.value)
        # hit time limit, no solution
        return result(-1, ComputationMode.TIMEOUT_BOUND.value)       # TODO: now we should use a filter
    elif stats.status == S.Infeasible:
        return result(threshold, ComputationMode.ABOVE_THRESHOLD.value)
    else:
        raise Exception('unknown ILP status: ', stats.status_str)

def add_MCES_to_molgraphs(ILP, G1, G2):
    mapping_G1_to_G2 = {}
    for i, v in enumerate(ILP.variables()):
        if (v.name.startswith('nodepairs_') and v.value() == 1):
            i, j = map(int, re.findall(r'\d+', v.name.split('_', 1)[1]))
            G1.nodes[i]['in_mces'] = True
            G2.nodes[j]['in_mces'] = True
        elif (v.name.startswith('edgepairs_') and len(indices:=list(map(int, re.findall(r'\d+', v.name))))==4 and v.value() == 1):
            i1, j1, i2, j2 = indices
            G1.edges[(i1, j1)]['in_mces'] = True
            G2.edges[(i2, j2)]['in_mces'] = True
            # # in G1, save all of the corresponding edges in G2
            # G1.edges[(i1, j1)]['corresponding_edge'] = [(i2, j2)]

            # frozenset would be more correct, because edges are undirected, but its too long
            #mapping_G1_to_G2[frozenset((i1, j1))] = frozenset((i2, j2))
            mapping_G1_to_G2[(i1, j1)] = (i2, j2)
    return mapping_G1_to_G2

