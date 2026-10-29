import pulp
import pytest

from myopic_mces.MCES_ILP import construct_ILP
from myopic_mces.graph import construct_graph

longs1 = 'CCCCCCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCC=CCCCCCCCC)OC(=O)CCCCCCCCCCCCCC'
longs2 = 'CCCCCCCCCCCCCCCCCCCCCC(=O)OC(COC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCC=CCCCCCCCC'


@pytest.mark.parametrize('s1, s2, thr, time_limit, status_exp, has_solution_exp',  
                        [   # infeasible
                            ('CCC', 'CCCCCCCC', 2., None, pulp.LpSolveStatus.Infeasible, False),
                            # optimal
                            ('CCC', 'CCCCCCCC', 10., None, pulp.LpSolveStatus.Optimal, True),
                            # infeasible in timelimit
                            pytest.param(longs1, longs2, 10., 10, pulp.LpSolveStatus.TimeLimit, False, marks = pytest.mark.timelimit()),
                            # optimal in timelimit
                            pytest.param(longs1, longs2, 10., 150, pulp.LpSolveStatus.TimeLimit, True, marks=pytest.mark.timelimit()),
                        ])
def test_CPLEX_PY(s1, s2, thr, time_limit, status_exp, has_solution_exp):
    sol=pulp.CPLEX_PY(threads=1, timeLimit=time_limit, msg=True)
    ILP = construct_ILP(construct_graph(s1), construct_graph(s2), thr)
    stats = ILP.solve(sol)
    print(f'expected status {status_exp.name}, solution {has_solution_exp}')
    print(f'status {stats.status_str}, solution {stats.has_solution}')
    assert stats.status == status_exp
    assert stats.has_solution == has_solution_exp

@pytest.mark.parametrize('s1, s2, thr, time_limit, status_exp, has_solution_exp',   
                        [   # infeasible → 2: above threshold → status infeasible, no solution
                            ('CCC', 'CCCCCCCC', 2., None, pulp.LpSolveStatus.Infeasible, False),
                            # optimal → 1: exact → status optimal, solution
                            ('CCC', 'CCCCCCCC', 10., None, pulp.LpSolveStatus.Optimal, True),
                            # infeasible in timelimit → 6: timelimit reached, no solution → status time limit, no solution
                            pytest.param(longs1, longs2, 10., 10, pulp.LpSolveStatus.TimeLimit, False, marks=pytest.mark.timelimit()),
                            # optimal in timelimit → 5: timelimit reached, solution → status time limit, solution
                            pytest.param(longs1, longs2, 10., 150, pulp.LpSolveStatus.TimeLimit, True, marks=pytest.mark.timelimit())
                        ])
def test_COIN_CMD(s1, s2, thr, time_limit, status_exp, has_solution_exp):
    sol = pulp.COIN_CMD(threads=1, timeLimit=time_limit, timeMode='cpu', msg=True)
    ILP = construct_ILP(construct_graph(s1), construct_graph(s2), thr)
    stats = ILP.solve(sol)
    print(f'time limit: {time_limit}')
    print(f'expected status {status_exp.name}, solution {has_solution_exp}')
    print(f'status {stats.status_str}, solution {stats.has_solution}')
    assert stats.status == status_exp
    assert stats.has_solution == has_solution_exp
