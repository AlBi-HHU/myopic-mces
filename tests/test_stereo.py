import itertools

import pulp
import pytest
from rdkit import Chem

from myopic_mces import MCES
from myopic_mces.filter_MCES import ComputationMode
from myopic_mces.graph import construct_graph

SOLVER = 'CPLEX_PY' if 'CPLEX_PY' in pulp.listSolvers(onlyAvailable=True) else 'COIN_CMD'
SOLVER_OPTIONS = dict(msg=False, threads=1)


def dist(s1, s2, **kwargs):
    return MCES(s1, s2, threshold=-1, solver=SOLVER, solver_options=SOLVER_OPTIONS, **kwargs)[1]


STEREO_PAIRS = [
    ('C/C=C/C', 'C/C=C\\C'),                                # E vs Z
    ('C/C=C/C', 'CC=CC'),                                   # E vs unknown
    ('C[C@H](N)O', 'C[C@@H](N)O'),                          # enantiomers
    ('C[C@H](N)O', 'CC(N)O'),                               # specified vs unknown
    ('C[C@@H](O)[C@H](C)O', 'C[C@@H](O)[C@@H](C)O'),        # meso vs (R,R)
    ('N[C@@H](C)C(=O)O', 'N[C@@H](CS)C(=O)O'),              # L-Ala vs L-Cys (CIP label flip)
    ('C[C@@H](O)CCCCCC', 'C[C@H](N)CCCCCC'),
    ('C/C=C/C[C@H](N)O', 'C/C=C\\C[C@@H](N)O'),
    ('CCO', 'CCCO'),                                        # no stereo at all
    ('C1=CC=C(C=C1)SC2=CC=CC=C2', 'C1=CC=C(C=C1)NC2=CC=CC=C2'),
]


def test_graph_attributes():
    G = construct_graph('C/C=C/C[C@H](N)O', stereo=True)
    assert G.edges[1, 2]['ez'] is not None
    assert G.edges[0, 1]['ez'] is None
    assert G.nodes[4]['chiral_tag'] in ('CW', 'CCW')
    assert len(G.nodes[4]['chiral_nbrs']) == 4
    assert G.nodes[0]['chiral_tag'] is None
    # without stereo, no attributes
    G = construct_graph('C/C=C/C[C@H](N)O')
    assert 'ez' not in G.edges[1, 2] and 'chiral_tag' not in G.nodes[4]


@pytest.mark.parametrize('s1, s2', STEREO_PAIRS)
def test_2d_is_lower_bound(s1, s2):
    assert dist(s1, s2) <= dist(s1, s2, stereo=True) + 1e-9


@pytest.mark.parametrize('s1, s2, kwargs, expected', [
    ('C/C=C/C', 'C/C=C\\C', {}, 0.5),
    ('C/C=C/C', 'C/C=C\\C', dict(stereo_cost=1.), 1.),
    ('C/C=C/C', 'C/C=C/C', {}, 0.),
    ('F/C=C/Cl', 'Cl/C=C/F', {}, 0.),
    ('C/C=N/O', 'C/C=N\\O', {}, 0.5),                            # oxime
    ('C/N=N/C', 'C/N=N\\C', {}, 0.5),                            # azo
    ('CC=CC', 'CC=CC', {}, 0.),
    ('C/C=C/C', 'CC=CC', {}, 0.25),                            # default: u = s/2
    ('C/C=C/C', 'CC=CC', dict(stereo_unknown_cost=0.5), 0.5),
    ('C[C@H](N)O', 'C[C@@H](N)O', {}, 0.5),
    ('C[C@H](N)O', 'C[C@H](N)O', {}, 0.),
    ('C[C@H](N)O', 'O[C@@H](N)C', {}, 0.),                        # same molecule, different atom order
    ('C[C@H](N)O', 'CC(N)O', {}, 0.25),
    ('C[C@H](N)O', 'CC(N)O', dict(stereo_unknown_cost=0.5), 0.5),
    ('C[C@@H](O)[C@H](C)O', 'C[C@@H](O)[C@@H](C)O', {}, 0.5),   # meso vs (R,R): one of two centres differs
    ('O[C@H]([C@H](O)C(=O)O)C(=O)O', 'O[C@@H]([C@@H](O)C(=O)O)C(=O)O', {}, 0.),  # meso vs its mirror image
    ('C/C=C/C[C@H](N)O', 'C/C=C\\C[C@@H](N)O', {}, 1.),
    ('C[C@H]1CCCC[C@@H]1C', 'C[C@H]1CCCC[C@H]1C', {}, 0.5),      # ring cis/trans, 1,2
    ('C[C@H]1CC[C@@H](C)CC1', 'C[C@H]1CC[C@H](C)CC1', {}, 0.5),  # ring cis/trans, 1,4
    # same geometry, different CIP labels: no penalty
    ('N[C@@H](C)C(=O)O', 'N[C@@H](CS)C(=O)O', {}, 1.),
    ('C/C=C(/CC)C', 'C/C=C(/CC)CO', {}, 1.),
    # different geometry, same CIP labels: penalty
    ('N[C@@H](C)C(=O)O', 'N[C@H](CS)C(=O)O', {}, 1.5),
    # stereocentre mapped, but not enough substituents: no penalty
    ('C[C@@H](O)CCCCCC', 'C[C@H](N)CCCCCC', {}, 2.),
])
def test_exact_values(s1, s2, kwargs, expected):
    assert dist(s1, s2, stereo=True, **kwargs) == pytest.approx(expected)
    assert dist(s2, s1, stereo=True, **kwargs) == pytest.approx(expected)


@pytest.mark.parametrize('s', ['C[C@@H]1CC[C@H](C(C)C)[C@@H](O)C1', 'OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O',
                               'C/C=C/C=C/C(=O)O', 'C/C(=C\\CO)C(C)C'])
def test_random_smiles(s):
    m = Chem.MolFromSmiles(s)
    for seed in range(5):
        s_rand = Chem.MolToSmiles(m, canonical=False, doRandom=True)
        assert dist(s, s_rand, stereo=True) == pytest.approx(0.)


@pytest.mark.parametrize('s1, s2', [p for p in STEREO_PAIRS if '@' not in ''.join(p) and '/' not in ''.join(p)])
def test_no_stereo_unchanged(s1, s2):
    assert dist(s1, s2, stereo=True) == pytest.approx(dist(s1, s2))


def test_threshold_includes_stereo():
    # 2D distance 0 passes the filters, stereo distance 0.5 is above the threshold
    _, d, _, mode = MCES('C/C=C/C', 'C/C=C\\C', threshold=0.25, solver=SOLVER, solver_options=SOLVER_OPTIONS,
                         stereo=True)
    assert mode == ComputationMode.ABOVE_THRESHOLD.value
    _, d, _, mode = MCES('C/C=C/C', 'C/C=C\\C', threshold=1, solver=SOLVER, solver_options=SOLVER_OPTIONS,
                         stereo=True)
    assert mode == ComputationMode.EXACT.value and d == pytest.approx(0.5)


def test_structure_mapping():
    # structure output still works with stereo penalty variables in the ILP
    _, d, _, mode, mapping = MCES('C[C@H](N)O', 'C[C@@H](N)O', threshold=-1, solver=SOLVER,
                                  solver_options=SOLVER_OPTIONS, stereo=True, structure=True)
    assert d == pytest.approx(0.5) and len(mapping) == 3


@pytest.mark.parametrize('triple', [
    ('C[C@H](N)O', 'CC(N)O', 'C[C@@H](N)O'),
    ('C/C=C/C', 'CC=CC', 'C/C=C\\C'),
    ('C[C@@H](O)[C@H](C)O', 'CC(O)C(C)O', 'C[C@@H](O)[C@@H](C)O'),
])
def test_triangle_inequality(triple):
    for a, b, c in itertools.permutations(triple):
        assert dist(a, c, stereo=True) <= dist(a, b, stereo=True) + dist(b, c, stereo=True) + 1e-9
