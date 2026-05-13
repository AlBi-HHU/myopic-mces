import pulp
import networkx as nx

import myopic_mces.MCES_ILP
from myopic_mces.MCES_ILP import construct_ILP
from myopic_mces.graph import construct_graph

def test_ILP_status():

    longs1 = 'CCCCCCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCC=CCCCCCCCC)OC(=O)CCCCCCCCCCCCCC'
    longs2 = 'CCCCCCCCCCCCCCCCCCCCCC(=O)OC(COC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCC=CCCCCCCCC'

    for s1, s2, thr, time_limit, status_exp in [('CCC', 'CCCCCCCC', 2., None, pulp.constants.LpStatusInfeasible),
                                                ('CCC', 'CCCCCCCC', 10., None, pulp.constants.LpStatusOptimal),
                                                (longs1, longs2, 10., 10, pulp.constants.LpStatusInfeasible),
                                                (longs1, longs2, 10., 150, pulp.constants.LpStatusOptimal),
                                    ]:
            sol=pulp.getSolver('CPLEX_PY', threads=1, timeLimit=time_limit)
            ILP = construct_ILP(construct_graph(s1), construct_graph(s2), thr)
            ILP.solve(sol)
            st = ILP.status
            st_detailed = ILP.solverModel.solution.get_status()
            print('status expected:', pulp.constants.LpStatus[status_exp], 'status actual:', pulp.constants.LpStatus[st],
                  'status detailed:', st_detailed, 'MIP_time_limit_feasible?', st_detailed == ILP.solverModel.solution.status.MIP_time_limit_feasible,
                  'MIP_time_limit_infeasible?', st_detailed == ILP.solverModel.solution.status.MIP_time_limit_infeasible)
            
