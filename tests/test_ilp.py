import pulp
import pytest

from myopic_mces.MCES_ILP import construct_ILP
from myopic_mces.graph import construct_graph

longs1 = 'CCCCCCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCC=CCCCCCCCC)OC(=O)CCCCCCCCCCCCCC'
longs2 = 'CCCCCCCCCCCCCCCCCCCCCC(=O)OC(COC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCC=CCCCCCCCC'


@pytest.mark.parametrize('s1, s2, thr, time_limit, status_exp, solution_exp',  
                        [   # infeasible
                            ('CCC', 'CCCCCCCC', 2., None, pulp.constants.LpStatusInfeasible, pulp.constants.LpSolutionInfeasible),
                            # optimal
                            ('CCC', 'CCCCCCCC', 10., None, pulp.constants.LpStatusOptimal, pulp.constants.LpSolutionOptimal),
                            # infeasible in timelimit
                            pytest.param(longs1, longs2, 10., 10, pulp.constants.LpStatusInfeasible, pulp.constants.LpSolutionInfeasible, marks = pytest.mark.timelimit()),
                            # optimal in timelimit
                            pytest.param(longs1, longs2, 10., 150, pulp.constants.LpStatusOptimal, pulp.constants.LpSolutionIntegerFeasible, marks=pytest.mark.timelimit()),
                        ])
def test_CPLEX_PY(s1, s2, thr, time_limit, status_exp, solution_exp):
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

@pytest.mark.parametrize('s1, s2, thr, time_limit, status_exp, solution_exp',   
                        [   # infeasible → 2: above threshold → problem status infeasible, solution status infeasible [-1,-1]
                            ('CCC', 'CCCCCCCC', 2., None, pulp.constants.LpStatusInfeasible, pulp.constants.LpSolutionInfeasible),
                            # optimal → 1: exact → problem status optimal und solution status optimal [1,1]
                            ('CCC', 'CCCCCCCC', 10., None, pulp.constants.LpStatusOptimal, pulp.constants.LpSolutionOptimal),
                            # infeasible in timelimit → 6: timelimit reached, no solution → problem status not solved, solution status no solution found [0,0]
                            pytest.param(longs1, longs2, 10., 10, pulp.constants.LpStatusNotSolved, pulp.constants.LpSolutionNoSolutionFound, marks=pytest.mark.timelimit()),
                            # optimal in timelimit → 5: timelimit reached, solution → problem status optimal und solution integer feasible [1,2]
                            pytest.param(longs1, longs2, 10., 150, pulp.constants.LpStatusOptimal, pulp.constants.LpSolutionIntegerFeasible, marks=pytest.mark.timelimit())
                        ])
def test_COIN_CMD(s1, s2, thr, time_limit, status_exp, solution_exp):
    sol = pulp.COIN_CMD(threads=1, timeLimit=time_limit, timeMode='cpu', msg=True)
    ILP = construct_ILP(construct_graph(s1), construct_graph(s2), thr)
    ILP.solve(sol)
    st = ILP.status
    st_sol = ILP.sol_status
    print(f'time limit: {time_limit}')
    print(f'expected status {status_exp}, solution {solution_exp}')
    print(f'status {st}, solution {st_sol}')
    assert pulp.constants.LpStatus[status_exp] == pulp.constants.LpStatus[st]
    assert pulp.constants.LpSolution[solution_exp] == pulp.constants.LpSolution[st_sol]