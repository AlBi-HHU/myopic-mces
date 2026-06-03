import pytest

from myopic_mces.MCES_ILP import ComputationMode, MCES_ILP
from myopic_mces.graph import construct_graph
from myopic_mces import MCES

longs1 = 'CCCCCCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCC=CCCCCCCCC)OC(=O)CCCCCCCCCCCCCC'
longs2 = 'CCCCCCCCCCCCCCCCCCCCCC(=O)OC(COC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCC=CCCCCCCCC'
longs3 = 'CC(C)OC(=O)C(CSCC(=O)NC1=C(C=C(C=C1Br)Br)CN(C)C2CCCCC2)NC(=O)C'
longs4 = 'CC(=O)N1CCC(CC1)N2CC3(C2)CN(CC4N3CCOC4)C(=O)C5=CC=C(C=C5)OC'
@pytest.mark.parametrize('s1, s2, thr, distance_exp, mode_exp, dynamic, bound_zero',
                        [   # optimal → 1: exact → problem status optimal und solution status optimal [1,1]
                            ('CCC', 'CCCCCCCC', -1, 5., ComputationMode.EXACT, True, True),
                            # infeasible → 2: above threshold → problem status infeasible, solution status infeasible [-1,-1]
                            ('CCCCCCCCCCCCCCCCCCCC', 'c1ccccc1N(=O)=O', 20., 20., ComputationMode.ABOVE_THRESHOLD, False, False),
                            # filter cases (20, 21, 22 depending on chosen filter)
                            ('CCC', longs1, 2, 63., ComputationMode.DYNAMIC_BOUND0, True, True),
                            ('CCC', longs1, 64, 67., ComputationMode.DYNAMIC_BOUND1, True, True),
                            (longs3, longs4, 10, 16., ComputationMode.DYNAMIC_BOUND2, True, True),
                            # always stronger bound case (4)
                            ('CCC', longs1, 10, 67., ComputationMode.STRONGEST_BOUND, False, False)
                        ])
def test_CPLEX_PY_mode(s1, s2, thr, distance_exp, mode_exp, dynamic, bound_zero):
    _, distance, _, compute_mode = MCES(s1, s2, thr, 0, 'CPLEX_PY', always_stronger_bound=not dynamic, use_bound_zero=bound_zero)
    assert distance == distance_exp
    assert ComputationMode(compute_mode) == mode_exp

@pytest.mark.timelimit('timelimit may differ for each machine')
# timeout is tested using MCES_ILP, otherwise the modes are overwritten by the filters
@pytest.mark.parametrize('s1, s2, thr, time_limit, distance_exp, mode_exp', 
                        [   # feasible under timelimit → 5: timeout exact
                            (longs1, longs2, 10, 150, 7., ComputationMode.TIMEOUT_ILP_UNPROVEN),
                            # fails on one thread under timelimit as infeasible → 6: timeout bound
                            (longs1, longs2, 10, 10, -1, ComputationMode.TIMEOUT_BOUND)
                        ])
def test_CPLEX_PY_timelimit(s1, s2, thr, time_limit, distance_exp, mode_exp):
    G1 = construct_graph(s1)
    G2 = construct_graph(s2)
    solver_options = dict(timeLimit=time_limit, threads=1, parallel=1)
    distance, compute_mode = MCES_ILP(G1, G2, thr, 'CPLEX_PY', solver_options)
    assert distance == distance_exp
    assert ComputationMode(compute_mode) == mode_exp

@pytest.mark.parametrize('s1, s2, thr, distance_exp, mode_exp, dynamic, bound_zero',
                        [   # optimal → 1: exact → problem status optimal und solution status optimal [1,1]
                            ('CCC', 'CCCCCCCC', -1, 5., ComputationMode.EXACT, True, True),
                            # infeasible → 2: above threshold → problem status infeasible, solution status infeasible [-1,-1]
                            ('CCCCCCCCCCCCCCCCCCCC', 'c1ccccc1N(=O)=O', 20., 20., ComputationMode.ABOVE_THRESHOLD, False, False),
                            # filter cases (20, 21, 22 depending on chosen filter)
                            ('CCC', longs1, 2, 63., ComputationMode.DYNAMIC_BOUND0, True, True),
                            ('CCC', longs1, 64, 67., ComputationMode.DYNAMIC_BOUND1, True, True),
                            (longs3, longs4, 10, 16., ComputationMode.DYNAMIC_BOUND2, True, True),
                            # always stronger bound case (4)
                            ('CCC', longs1, 10, 67., ComputationMode.STRONGEST_BOUND, False, False)
                        ])
def test_COIN_CMD_mode(s1, s2, thr, distance_exp, mode_exp, dynamic, bound_zero):
    _, distance, _, compute_mode = MCES(s1, s2, thr, 0, 'COIN_CMD',  always_stronger_bound=not dynamic, use_bound_zero=bound_zero)
    assert distance == distance_exp
    assert ComputationMode(compute_mode) == mode_exp

@pytest.mark.timelimit('timelimit may differ for each machine')
# timeout is tested using MCES_ILP, otherwise the modes are overwritten by the filters
@pytest.mark.parametrize('s1, s2, thr, time_limit, distance_exp, mode_exp', 
                        [   # feasible under timelimit → 5: timeout exact
                            (longs1, longs2, 10, 150, 7., ComputationMode.TIMEOUT_ILP_UNPROVEN),
                            # fails on one thread under timelimit as infeasible → 6: timeout bound
                            (longs1, longs2, 10, 10, -1, ComputationMode.TIMEOUT_BOUND)
                        ])
def test_COIN_CMD_timelimit(s1, s2, thr, time_limit, distance_exp, mode_exp):
    G1 = construct_graph(s1)
    G2 = construct_graph(s2)
    solver_options = dict(timeLimit=time_limit, threads=1, timeMode='cpu')
    distance, compute_mode = MCES_ILP(G1, G2, thr, 'COIN_CMD', solver_options)
    assert distance == distance_exp
    assert ComputationMode(compute_mode) == mode_exp