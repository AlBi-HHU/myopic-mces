import pytest
import pandas as pd

from myopic_mces.graph import construct_graph
from myopic_mces._cpp import mmces_filters, HAS_CPP, CPP_IMPORT_ERROR
# this is testing the c++ filters against the same reference values as the python filters

pytestmark = pytest.mark.cpp


def test_cpp_available():
    # fail instead of silently falling back to the python filters
    assert HAS_CPP, f'C++ filters not available: {CPP_IMPORT_ERROR}'


@pytest.mark.parametrize('filter_nr', [0, 1, 2])
def test_cpp_filter(filter_nr):
    assert HAS_CPP, f'C++ filters not available: {CPP_IMPORT_ERROR}'
    filter_fn = getattr(mmces_filters, f'filter{filter_nr}')
    ref = pd.read_csv(f'testdata/f{filter_nr}_res.csv', header=None, names=['i', 'mces', 't', 'mode'])
    smiles = pd.read_csv('testdata/test_smiles.csv', header=None, names=['i', 'smiles1', 'smiles2'])
    df = smiles.merge(ref[['i', 'mces']], on='i')

    df['filter_distance'] = df.apply(
        lambda row: filter_fn(construct_graph(row.smiles1, save_mol=True).mol,
                              construct_graph(row.smiles2, save_mol=True).mol),
        axis=1)
    assert (df['filter_distance'] == df['mces']).all()
