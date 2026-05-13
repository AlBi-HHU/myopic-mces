import time
import pytest
import pandas as pd

from myopic_mces.graph import construct_graph
from myopic_mces.filter_MCES import filter0, filter1, filter2
from myopic_mces import MCES

def test_filter0_tl():
    start = time.time()
    f0_data = pd.read_csv('tests/testdata/f0_res.csv', header=None, names=['i', 'mces', 't', 'mode'])
    smiles = pd.read_csv('tests/testdata/test_smiles.csv', header=None, names=['i', 'smiles1', 'smiles2'])
    df = smiles.merge(f0_data[['i', 'mces']], on='i')
    
    df['filter_distance'] = df.apply(
        lambda row: filter0(construct_graph(row.smiles1),
                            construct_graph(row.smiles2)),
        axis=1
    )
    print(time.time() - start)
    assert (df['filter_distance'] == df['mces']).all()

def test_filter1_tl():
    start = time.time()
    f1_data = pd.read_csv('tests/testdata/f1_res.csv', header=None, names=['i', 'mces', 't', 'mode'])
    smiles = pd.read_csv('tests/testdata/test_smiles.csv', header=None, names=['i', 'smiles1', 'smiles2'])
    df = smiles.merge(f1_data[['i', 'mces']], on='i')
    
    df['filter_distance'] = df.apply(
        lambda row: filter1(construct_graph(row.smiles1),
                            construct_graph(row.smiles2)),
        axis=1
    )
    print(time.time() - start)
    assert (df['filter_distance'] == df['mces']).all()


def test_filter2_tl():
    start = time.time()
    f2_data = pd.read_csv('tests/testdata/f2_res.csv', header=None, names=['i', 'mces', 't', 'mode'])
    smiles = pd.read_csv('tests/testdata/test_smiles.csv', header=None, names=['i', 'smiles1', 'smiles2'])
    df = smiles.merge(f2_data[['i', 'mces']], on='i')
    
    df['filter_distance'] = df.apply(
        lambda row: filter2(construct_graph(row.smiles1),
                            construct_graph(row.smiles2)),
        axis=1)
    print(time.time() - start)
    assert (df['filter_distance'] == df['mces']).all()

def test_all_filters(): 
    # assert that the three filters do not return the same result
    start = time.time()
    smiles = pd.read_csv('tests/testdata/test_smiles.csv', header=None, names=['i', 'smiles1', 'smiles2'])
    
    smiles['filter0_distance'] = smiles.apply(
        lambda row: filter0(construct_graph(row.smiles1),
                            construct_graph(row.smiles2)),
        axis=1
    )
    smiles['filter1_distance'] = smiles.apply(
        lambda row: filter1(construct_graph(row.smiles1),
                            construct_graph(row.smiles2)),
        axis=1
    )
    smiles['filter2_distance'] = smiles.apply(
        lambda row: filter2(construct_graph(row.smiles1),
                            construct_graph(row.smiles2)),
        axis=1
    )
    
    smiles.set_index('i')
    print(time.time() - start)
    

    # filter0 is not always smaller than filter1, so only assert that filter1 is smaller than filter2 in all cases
    assert (smiles['filter1_distance'] <= smiles['filter2_distance']).all()

    # specific cases in which filter0 is not smaller than filter1
    special_idx = pd.read_csv('tests/testdata/special_idx.csv', header=None, names=['i'])
    idx = special_idx['i'].tolist()
    filter0_vals_special = smiles['filter0_distance'].loc[smiles['i'].isin(idx)]
    filter1_vals_special = smiles['filter1_distance'].loc[smiles['i'].isin(idx)]
    assert (filter0_vals_special >= filter1_vals_special).all()
    
    # and assert all other cases 
    filter0_vals = smiles['filter0_distance'].loc[~smiles['i'].isin(idx)]
    filter1_vals = smiles['filter1_distance'].loc[~smiles['i'].isin(idx)]
    assert (filter0_vals <= filter1_vals).all()

    # violation_01 = smiles['filter0_distance'] > smiles['filter1_distance']
    # smiles.loc[violation_01].set_index('i').to_csv('tests/testdata/special.csv', header=None)
    

def test_same_as_main():
    start = time.time()
    main_df = pd.read_csv('tests/testdata/bio_dists_main_mces.csv', header=None, names=['i', 'mces'])
    smiles = pd.read_csv('tests/testdata/smiles_subset.csv', header=None, names=['i', 'smiles1', 'smiles2'])
    smiles['i'], smiles['mces'], smiles['t'], smiles['mode'] = zip(*smiles.apply(
        lambda row: MCES(row.smiles1, row.smiles2),
        axis=1
    ))
    print(time.time() - start)
    assert (main_df['mces'] == smiles['mces']).all()

