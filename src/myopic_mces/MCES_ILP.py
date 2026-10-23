# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:17:41 2020

@author: seipp
"""
import pulp
import networkx as nx
from myopic_mces.filter_MCES import ComputationMode
import math
import re
from itertools import combinations, permutations

def permutation_parity(perm):
    """0 for even, 1 for odd permutations of range(len(perm))"""
    perm, swaps = list(perm), 0
    for i in range(len(perm)):
        while perm[i] != i:
            j = perm[i]
            perm[i], perm[j] = perm[j], perm[i]
            swaps += 1
    return swaps % 2


def same_chirality(n1, tag1, n2, tag2, mapping):
    """
     Whether two tetrahedral centres have the same configuration under a (partial) mapping of their neighbours

     n1, n2: neighbours (4, incl. placeholders for implicit H) in the order the chiral tags refer to
     mapping: dict from 3 neighbours in n1 to neighbours in n2, the remaining neighbours are paired up
    """
    rest = [x for x in n2 if x not in mapping.values()][0]
    perm = [n2.index(mapping.get(x, rest)) for x in n1]
    return (tag1 == tag2) == (permutation_parity(perm) == 0)


def ez_trans(ez, i, a, b):
    """
     Whether substituent a (bonded to double bond atom i) and substituent b (bonded to the other double bond atom)
     are trans, for a double bond with stereo info `ez` (see construct_graph)
    """
    begin, end, sa, sb, trans = ez
    if i != begin:
        a, b = b, a
    return trans ^ (a != sa) ^ (b != sb)


def stereo_terms(ILP, G1, G2, y, stereo_cost, stereo_unknown_cost):
    """
     Adds stereo penalty variables and constraints to the ILP, returns the objective terms.

     A tetrahedral centre pair i->k is penalized if at least 3 substituents are mapped onto substituents and
     the permutation parity of the mapping does not match the chiral tags; a double bond pair if one substituent
     on each side is mapped and the cis/trans relation differs. If stereo is specified only for one of the two,
     stereo_unknown_cost is charged once the geometry is defined (>= 3 resp. 1+1 mapped substituents).

     Note: the resulting distance satisfies the triangle inequality only if stereo_unknown_cost >= stereo_cost / 2
    """
    terms = []

    def add_penalty(name, cost, lower_bounds):
        # penalty variable q >= each lower bound, only added if needed
        if len(lower_bounds) > 0:
            q = ILP.add_variable(name, lowBound=0)
            for lb in lower_bounds:
                ILP.addConstraint(q >= lb)
            terms.append(cost * q)

    # tetrahedral centres
    for i, k in y:
        t1, t2 = G1.nodes[i].get("chiral_tag"), G2.nodes[k].get("chiral_tag")
        if t1 is None and t2 is None:
            continue
        cost = stereo_cost if (t1 is not None and t2 is not None) else stereo_unknown_cost
        n1, n2 = list(G1.neighbors(i)), list(G2.neighbors(k))
        if cost == 0 or len(n1) < 3 or len(n2) < 3:
            continue
        lower_bounds = []
        for T in combinations(n1, 3):
            if t1 is not None and t2 is not None:
                if len(G1.nodes[i]["chiral_nbrs"]) != 4 or len(G2.nodes[k]["chiral_nbrs"]) != 4:
                    break
                for img in permutations(n2, 3):
                    pairs = list(zip(T, img))
                    if all(p in y for p in pairs) and not same_chirality(
                            G1.nodes[i]["chiral_nbrs"], t1, G2.nodes[k]["chiral_nbrs"], t2, dict(pairs)):
                        lower_bounds.append(y[(i, k)] + pulp.lpSum(y[p] for p in pairs) - 3)
            else:
                # only one centre is specified: penalize once the geometry is defined
                mapped = [y[(a, b)] for a in T for b in n2 if (a, b) in y]
                if len(mapped) > 0:
                    lower_bounds.append(y[(i, k)] + pulp.lpSum(mapped) - 3)
        add_penalty(f"stereoatom_{i}_{k}", cost, lower_bounds)
    # double bonds
    for e1 in G1.edges:
        if G1.edges[e1]["weight"] != 2:
            continue
        for e2 in G2.edges:
            if G2.edges[e2]["weight"] != 2:
                continue
            z1, z2 = G1.edges[e1].get("ez"), G2.edges[e2].get("ez")
            if z1 is None and z2 is None:
                continue
            cost = stereo_cost if (z1 is not None and z2 is not None) else stereo_unknown_cost
            if cost == 0:
                continue
            lower_bounds = []
            k, l = e2
            for i, j in (e1, e1[::-1]):
                if (i, k) not in y or (j, l) not in y:
                    continue
                for a in G1.neighbors(i):
                    for b in G1.neighbors(j):
                        if a in (i, j) or b in (i, j) or a == b:
                            continue
                        for a2 in G2.neighbors(k):
                            for b2 in G2.neighbors(l):
                                if a2 in (k, l) or b2 in (k, l) or a2 == b2 or (a, a2) not in y or (b, b2) not in y:
                                    continue
                                if z1 is not None and z2 is not None and \
                                        ez_trans(z1, i, a, b) == ez_trans(z2, k, a2, b2):
                                    continue
                                lower_bounds.append(y[(i, k)] + y[(j, l)] + y[(a, a2)] + y[(b, b2)] - 3)
            add_penalty(f"stereobond_{e1[0]}_{e1[1]}_{e2[0]}_{e2[1]}", cost, lower_bounds)
    return terms


def construct_ILP(G1, G2, threshold, no_ilp_threshold=False, stereo=False, stereo_cost=0.5, stereo_unknown_cost=0.25):
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
     stereo: bool
         if true, add penalties for mapping tetrahedral centres and double bonds with different stereochemistry.
         The geometry is compared under the mapping (permutation parity of the mapped neighbours, cis/trans relation
         of mapped substituents), and only charged once enough substituents are mapped (3 for centres, 1 on each side
         for double bonds), see `stereo_terms`.
         Graphs have to be constructed with `construct_graph(..., stereo=True)`.
         As all penalties are non-negative, the non-stereo distance is a lower bound.
     stereo_cost: float
         penalty for mapping two elements with different known stereochemistry
     stereo_unknown_cost: float
         penalty for mapping two elements where stereochemistry is known for only one of them

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
    objective=[ w[i]*c[i] for i in edgepairs]
    #penalty for mapping stereo centres / double bonds with different configurations
    if stereo:
        objective += stereo_terms(ILP, G1, G2, y, stereo_cost, stereo_unknown_cost)
    ILP += pulp.lpSum(objective)




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
        ILP +=pulp.lpSum(objective)<=threshold

    return ILP


def MCES_ILP(G1, G2, threshold, solver='COIN_CMD', solver_options={}, no_ilp_threshold=False, structure=False,
             stereo=False, stereo_cost=0.5, stereo_unknown_cost=0.25):
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
     stereo, stereo_cost, stereo_unknown_cost:
         stereo-aware distance, see `construct_ILP`

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

    ILP = construct_ILP(G1, G2, threshold=threshold, no_ilp_threshold=no_ilp_threshold,
                        stereo=stereo, stereo_cost=stereo_cost, stereo_unknown_cost=stereo_unknown_cost)

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
    if solver=="COIN_CMD":
        sol=pulp.getSolver(solver="COIN_CMD", **solver_options)
        ILP.solve(sol)

        ilp_code = ILP.status
        ilp_sol_code = ILP.sol_status

        if ilp_code == pulp.constants.LpStatusOptimal:

            if ilp_sol_code == pulp.constants.LpSolutionOptimal:

                return result(float(ILP.objective.value()), ComputationMode.EXACT.value)
            
            if (solver_options.get('timeLimit') is not None and ilp_sol_code == pulp.constants.LpSolutionIntegerFeasible): 
                # timelimit hit, so mode 5
                return result(float(ILP.objective.value()), ComputationMode.TIMEOUT_ILP_UNPROVEN.value)
                # if this is zero, we use filter
        
        elif (ilp_code == pulp.constants.LpStatusInfeasible or ilp_sol_code == pulp.constants.LpSolutionInfeasible):
            
            if (solver_options.get('timeLimit') is not None\
                and math.isclose(ILP.solutionCpuTime, solver_options.get('timeLimit'), abs_tol=1) \
                and float(ILP.objective.value()) == 0): 
                return result(-1, ComputationMode.TIMEOUT_BOUND.value)
            
            return result(threshold, ComputationMode.ABOVE_THRESHOLD.value)
        
        elif (ilp_code == pulp.constants.LpStatusNotSolved):
            # must be timelimit
            return result(-1, ComputationMode.TIMEOUT_BOUND.value)
        else:
            raise Exception('unknown ILP status: ', ILP.status, pulp.constants.LpStatus[ILP.status])
    else:
        sol=pulp.getSolver(solver, **solver_options)
        ILP.solve(sol)

        ilp_code = ILP.status
        ilp_sol_code = ILP.sol_status

        try:
            ilp_code_detailed = ILP.solverModel.solution.get_status()
            ilp_code_time_limit_feasible = ILP.solverModel.solution.status.MIP_time_limit_feasible
            ilp_code_time_limit_infeasible = ILP.solverModel.solution.status.MIP_time_limit_infeasible
        except AttributeError:
            # only works for some solvers, namely CPLEX_PY
            ilp_code_detailed = ilp_code_time_limit_feasible = ilp_code_time_limit_infeasible = None

        if ilp_code == pulp.constants.LpStatusOptimal:
            if ilp_code_detailed is not None and ilp_code_detailed == ilp_code_time_limit_feasible:
                # hit time limit, but still found a solution
                return result(float(ILP.objective.value()), ComputationMode.TIMEOUT_ILP_UNPROVEN.value)
            return result(float(ILP.objective.value()), ComputationMode.EXACT.value)
        elif ilp_code == pulp.constants.LpStatusInfeasible:
            if ilp_code_detailed is not None and ilp_code_detailed == ilp_code_time_limit_infeasible:
                # hit time limit, no solution
                return result(-1, ComputationMode.TIMEOUT_BOUND.value)       # TODO: now we should use a filter
            return result(threshold, ComputationMode.ABOVE_THRESHOLD.value)
        # elif ilp_code == pulp.constants.LpStatusNotSolved:
        #     # must be time limit
        #     return float(ILP.objective.value()), ComputationMode.TIMEOUT_EXACT_SOLUTION.value
        else:
            raise Exception('unknown ILP status: ', ILP.status, pulp.constants.LpStatus[ILP.status])

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

