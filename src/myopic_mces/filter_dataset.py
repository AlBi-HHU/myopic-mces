"""For each structure of the query dataset, searchs for a match in the 'database' structures.
Match means: structure with MCES distance smaller than the provided threshold.
Output: List of structures with a match.
"""

from myopic_mces.graph import construct_graph
from myopic_mces.MCES_ILP import MCES_ILP
from myopic_mces.filter_MCES import apply_filter
from joblib import Parallel, delayed
import multiprocessing
import json
import argparse
import time

def MCES(G1, G2, threshold, solver, solver_options={}, always_stronger_bound=False):
    d, filter_id = apply_filter(G1, G2, threshold, always_stronger_bound=always_stronger_bound)
    if d > threshold:
        return d, filter_id
    # calculate MCES
    d_exact, filter_id = MCES_ILP(G1, G2, threshold, solver, solver_options=solver_options)
    return round(d_exact, 1), filter_id # ILP output can have varying precision

def MCES_query(i, smiles_query, db_list, threshold, solver, solver_options={}, always_stronger_bound=False):
    t0 = time.time()
    G1 = construct_graph(smiles_query)
    for j, smiles_db in enumerate(db_list):
        G2 = construct_graph(smiles_db)
        d, filter_id = MCES(G1, G2, threshold, solver, solver_options=solver_options, always_stronger_bound=always_stronger_bound)
        if (d <= threshold) and (filter_id == 1): # has to be solved exactly
            return i, smiles_query, smiles_db, time.time() - t0, d, filter_id
    return None

def MCES_query_shared_graphs(graphs, i, smiles_query, db_list, threshold, solver, solver_options={}, always_stronger_bound=False):
    t0 = time.time()
    G1 = construct_graph(smiles_query)
    for j, smiles_db in enumerate(db_list):
        G2 = graphs[smiles_db]
        d, filter_id = MCES(G1, G2, threshold, solver, solver_options=solver_options, always_stronger_bound=always_stronger_bound)
        if (d <= threshold) and (filter_id == 1): # has to be solved exactly
            return i, smiles_query, smiles_db, time.time() - t0, d, filter_id
    return None

# started converting this funciton for mpire usage
def get_mces_fun(threshold, solver, solver_options={}, always_stronger_bound=False):
    return lambda args: MCES_query(
        args[0], args[1][0], args[1][1],
        threshold=threshold, solver=solver, solver_options=solver_options, always_stronger_bound=always_stronger_bound)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input', help='json dictionary with format {query: [db_list]} (SMILES)')
    parser.add_argument('--out', help='structures of input with match')
    parser.add_argument('--threshold', type=float, help='threshold for the distance')
    parser.add_argument('--solver', type=str, default='default',
                        action='store', help='Solver for the ILP. example:CPLEX_CMD')
    parser.add_argument('--solver_onethreaded', action='store_true',
                        help='limit ILP solver to one thread, resulting in faster '
                        'performance with parallel computations (not available for all solvers)')
    parser.add_argument('--solver_no_msg', action='store_true',
                        help='prevent solver from logging (not available for all solvers)')
    parser.add_argument('--num_jobs', type=int, help='Number of jobs; instances to run in parallel. '
                        'By default this is set to the number of (logical) CPU cores.')
    parser.add_argument('--mpire', action='store_true')

    args = parser.parse_args()
    num_jobs = multiprocessing.cpu_count() if args.num_jobs is None else args.num_jobs

    solver_options = {}
    if args.solver_onethreaded:
        solver_options['threads'] = 1
    if args.solver_no_msg:
        solver_options['msg'] = False

    # preprocessing DB: SMILES -> Graphs
    # t0 = time.time()
    # db = {l.strip(): construct_graph(l.strip()) for l in open(args.db).readlines()}
    # print(f'contructing graphs done in {time.time() - t0:.1}s')

    print('reading input...')
    t0 = time.time()
    with open(args.input) as f:
        input_data = json.load(f)
    print(f'read input in{time.time() - t0:.1f}s, {len(input_data)} query compounds')

    if (args.mpire):
        fun = get_mces_fun(threshold=args.threshold, solver=args.solver, solver_options=solver_options, always_stronger_bound=False)
        from mpire import WorkerPool
        from mpire.utils import make_single_arguments
        with WorkerPool(n_jobs=num_jobs) as pool:
            results = pool.map_unordered(fun, make_single_arguments(enumerate(input_data.items())), progress_bar=True, iterable_len=len(input_data))
    else:
        results = Parallel(n_jobs=num_jobs, verbose=5)(
            delayed(MCES_query)(i, query, db_list, args.threshold, args.solver, solver_options=solver_options, always_stronger_bound=False)
            for i, (query, db_list) in enumerate(input_data.items()))
    with open(args.out, 'w') as out:
        for res in results:
            if res is None:
                continue
            out.write(','.join([str(x) for x in res]) + '\n')
