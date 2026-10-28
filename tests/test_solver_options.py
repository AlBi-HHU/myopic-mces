import pytest

from myopic_mces import MCES
from myopic_mces.MCES_ILP import check_solver_options
# solver options with which the computation modes could be wrong have to be rejected


@pytest.mark.parametrize('solver, solver_options', [
    ('CPLEX_CMD', dict(timeLimit=10)),
    ('GLPK_CMD', dict(timeLimit=10)),
    ('COIN_CMD', dict(gapAbs=0.5)),
    ('COIN_CMD', dict(gapRel=1e-3)),
])
def test_rejected(solver, solver_options):
    with pytest.raises(ValueError):
        check_solver_options(solver, solver_options)
    # also with catch_errors, which would otherwise turn the error into distance -1
    with pytest.raises(ValueError):
        MCES('CCC', 'CCCCCCCC', solver=solver, solver_options=solver_options, catch_errors=True)


@pytest.mark.parametrize('solver, solver_options', [
    ('COIN_CMD', dict(timeLimit=10)),
    ('CPLEX_PY', dict(timeLimit=10)),
    ('CPLEX_CMD', dict(threads=1)),
    ('COIN_CMD', dict(gapAbs=0.1, gapRel=1e-4)),
])
def test_accepted(solver, solver_options):
    check_solver_options(solver, solver_options)
