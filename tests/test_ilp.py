import pulp
import networkx as nx
import shutil

import myopic_mces.MCES_ILP
from myopic_mces.MCES_ILP import construct_ILP
from myopic_mces.graph import construct_graph
from myopic_mces import MCES

# CPLEX_PY:
# 1: exact → problem status optimal und solution status optimal
# 5: timelimit reached, solution → problem status optimal und solution status integer feasible
# 6: timelimit reached, no solution → problem status infeasible und solution status integer infeasible 
def test_ILP_CPLEX_PY():

    longs1 = 'CCCCCCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCC=CCCCCCCCC)OC(=O)CCCCCCCCCCCCCC'
    longs2 = 'CCCCCCCCCCCCCCCCCCCCCC(=O)OC(COC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCC=CCCCCCCCC'

    for s1, s2, thr, time_limit, status_exp, solution_exp in [  # infeasible
                                                                ('CCC', 'CCCCCCCC', 2., None, pulp.constants.LpStatusInfeasible, pulp.constants.LpSolutionInfeasible),
                                                                # infeasible in timelimit
                                                                ('CCC', 'CCCCCCCC', 2., 0.2, pulp.constants.LpStatusInfeasible, pulp.constants.LpSolutionInfeasible),
                                                                # optimal
                                                                ('CCC', 'CCCCCCCC', 10., None, pulp.constants.LpStatusOptimal, pulp.constants.LpSolutionOptimal),
                                                                # infeasible in timelimit
                                                                (longs1, longs2, 10., 10, pulp.constants.LpStatusInfeasible, pulp.constants.LpSolutionInfeasible),
                                                                # optimal in timelimit
                                                                (longs1, longs2, 10., 150, pulp.constants.LpStatusOptimal, pulp.constants.LpSolutionIntegerFeasible),
                                                            ]:
            # solver_options = dict(threads=1, timeLimit=time_limit)
            # MCES_ILP(s1, s2, threshold=thr, solver='CPLEX_PY', **solver_options)
            sol=pulp.CPLEX_PY(threads=1, timeLimit=time_limit, msg=True)
            ILP = construct_ILP(construct_graph(s1), construct_graph(s2), thr)
            ILP.solve(sol)
            st = ILP.status
            st_detailed = ILP.solverModel.solution.get_status()
            st_sol = ILP.sol_status
            print('status expected:', pulp.constants.LpStatus[status_exp], 'status actual:', pulp.constants.LpStatus[st],
                  'status detailed:', st_detailed, 'MIP_time_limit_feasible?', st_detailed == ILP.solverModel.solution.status.MIP_time_limit_feasible,
                  'MIP_time_limit_infeasible?', st_detailed == ILP.solverModel.solution.status.MIP_time_limit_infeasible)
            assert pulp.constants.LpStatus[status_exp] == pulp.constants.LpStatus[st]
            assert pulp.constants.LpSolution[solution_exp] == pulp.constants.LpSolution[st_sol]

# CBC via COIN_CMD:
# 1: exact → problem status optimal und solution status optimal [1,1]
# now pulp bugs:
# 5: timelimit reached, solution → problem status optimal und solution status optimal [1,1]
# 6: timelimit reached, no solution → problem status optimal und solution status integer feasible [1,2] (see lines 329ff. in pulp/apis/coin_api.py)
def test_ILP_COIN_CMD():
    longs1 = 'CCCCCCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCC=CCCCCCCCC)OC(=O)CCCCCCCCCCCCCC'
    longs2 = 'CCCCCCCCCCCCCCCCCCCCCC(=O)OC(COC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCC=CCCCCCCCC'
    # TODO: change logic if new pulp version fixes computation modes
    for s1, s2, thr, time_limit, status_exp, solution_exp in [  # infeasible
                                                                ('CCC', 'CCCCCCCC', 2., None, pulp.constants.LpStatusInfeasible, pulp.constants.LpSolutionInfeasible),
                                                                # optimal
                                                                ('CCC', 'CCCCCCCC', 10., None, pulp.constants.LpStatusOptimal, pulp.constants.LpSolutionOptimal),
                                                                # we cannot check for timelimit (in)feasible
                                                                # infeasible in timelimit
                                                                (longs1, longs2, 10., 10, pulp.constants.LpStatusOptimal, pulp.constants.LpSolutionIntegerFeasible),
                                                                # optimal in timelimit
                                                                (longs1, longs2, 10., 800, pulp.constants.LpStatusOptimal, pulp.constants.LpSolutionOptimal),
                                                            ]:

            cbc_path = ".venv/lib/python3.12/site-packages/pulp/solverdir/cbc/linux/i64/cbc"
            sol = pulp.COIN_CMD(path=cbc_path, threads=1, timeLimit=time_limit, timeMode='cpu', msg=True)
            ILP = construct_ILP(construct_graph(s1), construct_graph(s2), thr)
            ILP.solve(sol)
            st = ILP.status
            st_sol = ILP.sol_status
            print(f'time limit: {time_limit}')
            print(f'expected status {status_exp}, solution {solution_exp}')
            print(f'status {st}, solution {st_sol}')
            assert pulp.constants.LpStatus[status_exp] == pulp.constants.LpStatus[st]
            assert pulp.constants.LpSolution[solution_exp] == pulp.constants.LpSolution[st_sol]
