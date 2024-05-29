# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:16:05 2020
@author: seipp
"""
import time
from joblib import Parallel, delayed
import multiprocessing
import numpy as np
import argparse
from myopic_mces.graph import construct_graph
from myopic_mces.MCES_ILP import MCES_ILP
from myopic_mces.filter_MCES import apply_filter

def MCES(smiles1, smiles2, threshold=10, i=0, solver='default', solver_options={}, no_ilp_threshold=False, always_stronger_bound=True, catch_errors=False):
    """
    Calculates the distance between two molecules

    Parameters
    ----------
    smiles1 : str
        SMILES of the first molecule
    smiles2 : str
        SMILES of the second molecule
    threshold : float
        Threshold for the comparison. Exact distance is only calculated if the distance is lower than the threshold.
        If set to -1 the exact distance is always calculated.
    i : int
        index, mainly for parallelization
    solver: string
        ILP-solver used for solving MCES. Example:CPLEX_CMD
    solver_options: dict
        additional options to pass to solvers. Example: threads=1 for better multi-threaded performance
    no_ilp_threshold: bool
        if true, always return exact distance even if it is below the threshold (slower)
    always_stronger_bound: bool
        if true, always compute and use the second stronger bound

    Returns:
    -------
    int
        index
    float
        Distance between the molecules
    float
        Time taken for the calculation
    int
        Type of Distance:
            1 : Exact Distance
            2 : Lower bound (if the exact distance is above the threshold; bound chosen dynamically)
            4 : Lower bound (second lower bound was used)

    """
    start = time.time()
    # construct graph for both smiles.
    G1 = construct_graph(smiles1)
    G2 = construct_graph(smiles2)
    if threshold != -1:         # with `-1` always compute exact distance
        # filter out if distance is above the threshold
        try:
            distance, compute_mode = apply_filter(G1, G2, threshold, always_stronger_bound=always_stronger_bound)
            if distance > threshold:
                return i, distance, time.time() - start, compute_mode
        except Exception as e:
            print(e)
            if (catch_errors):
                distance = -1
                compute_mode = 2
            else:
                raise e
    # calculate MCES
    try:
        distance, compute_mode = MCES_ILP(G1, G2, threshold, solver, solver_options=solver_options,
                                          no_ilp_threshold=no_ilp_threshold)
    except Exception as e:
        print(e)
        if (catch_errors):
            distance = -1
            compute_mode = 1
        else:
            raise e
    return i, distance, time.time() - start, compute_mode

def hdf5_input(file_path):
    import h5py
    print('reading hdf5 input')
    t0 = time.time()
    with h5py.File(file_path, 'r') as f:
        smiles = np.asarray(f['smiles'])
        indices = f['computation_indices']
        inputs_raw = [indices[:, 0], smiles[indices[:, 1]], smiles[indices[:, 2]]]
    print(f'done, took {(time.time() - t0) / 60:.1f}min')
    return zip(*inputs_raw)

def hdf5_output(results, file_path, write_times=True, write_modes=True, args={}):
    import h5py
    print('writing hdf5 output')
    t0 = time.time()
    with h5py.File(file_path, 'a') as f:
        indices = f['computation_indices']
        assert list(indices[:, 0]) == [row[0] for row in results], 'something went wrong with the index order' # TODO: this could also be fixed, but just shouldn't happen
        f.create_dataset('mces', data=[row[1] for row in results], compression='gzip')
        if (write_times):
            f.create_dataset('computation_times', data=[row[2] for row in results], compression='gzip')
        if (write_modes):
            f.create_dataset('computation_modes', data=[row[3] for row in results], dtype='uint8', compression='gzip')
        comp_args = f.create_group('computation_args')
        for k, v in args:
            if v is not None:
                comp_args[k] = v
    print(f'done, took {(time.time() - t0) / 60:.1f}min')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input", help="input file in the format: id,smiles1,smiles2 OR hdf5 file when using hdf5 mode")
    parser.add_argument("output", help="output file")
    parser.add_argument("--threshold", type=float, default=10.,
                        action="store", help="threshold for the distance")
    parser.add_argument("--no_ilp_threshold", action="store_true",
                        help="(experimental) if set, do not add threshold as constraint to ILP, "
                        "resulting in longer runtimes and potential violations of the triangle equation")
    parser.add_argument("--choose_bound_dynamically", action="store_true",
                        help="if this is set, compute and use potentially weaker but faster lower bound if "
                        "already greater than the threshold. Otherwise (default), the strongest lower bound "
                        "is always computed and used. Enabling this can lead to massive speedups.")
    parser.add_argument("--solver", type=str, default="default",
                        action="store", help="Solver for the ILP. example:CPLEX_CMD")
    parser.add_argument("--solver_onethreaded", action="store_true",
                        help="limit ILP solver to one thread, resulting in faster "
                        "performance with parallel computations (not available for all solvers)")
    parser.add_argument("--solver_no_msg", action="store_true",
                        help="prevent solver from logging (not available for all solvers)")
    parser.add_argument("--num_jobs", type=int, help="Number of jobs; instances to run in parallel. "
                        "By default this is set to the number of (logical) CPU cores.")
    parser.add_argument('--hdf5_mode', action='store_true',
                        help='more time and space efficient mode of input/output using the `input` hdf5-file. '
                        'Has to contain one nx3 array with indices (`computation_indices`): (id_, smiles1_i, smiles2_i) '
                        'and another array with the corresponding SMILES (`smiles`). Output will be written to the same file.')
    parser.add_argument('--hide_rdkit_warnings', action='store_true')
    parser.add_argument('--catch_computation_errors', action='store_true')
    parser.add_argument('--jobs_batch_size', type=int, default=1000)
    parser.add_argument('--jobs_dispatch', default='10*n_jobs')
    args = parser.parse_args()

    if (args.hide_rdkit_warnings):
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
        # TODO: does not work
    num_jobs = multiprocessing.cpu_count() if args.num_jobs is None else args.num_jobs
    additional_mces_options = dict(no_ilp_threshold=args.no_ilp_threshold, solver_options=dict(),
                                   always_stronger_bound=not args.choose_bound_dynamically,
                                   catch_errors=args.catch_computation_errors)
    if (args.solver_onethreaded):
        additional_mces_options['solver_options']['threads'] = 1
    if (args.solver_no_msg):
        additional_mces_options['solver_options']['msg'] = False

    if (args.hdf5_mode):
        inputs = hdf5_input(args.input)
    else:
        with open(args.input) as in_handle:
            inputs = [line.strip().split(',')[:3] for line in in_handle]

    if (num_jobs > 1):
        results = Parallel(n_jobs=num_jobs, verbose=5, batch_size=args.jobs_batch_size, pre_dispatch=args.jobs_dispatch)(
            delayed(MCES)(smiles1, smiles2, args.threshold, i, args.solver, **additional_mces_options) for i, smiles1, smiles2 in inputs)
    else:
        results = [MCES(smiles1, smiles2, args.threshold, i, args.solver, **additional_mces_options) for i, smiles1, smiles2 in inputs]

    if (args.hdf5_mode):
        hdf5_output(results, args.input, write_times=True, write_modes=True, args=args._get_kwargs())
    else:
        with open(args.output, 'w') as out_handle:
            for ind, distance, duration, compute_mode in results:
                out_handle.write(f'{ind},{distance},{duration},{compute_mode}\n')

if __name__ == '__main__':
    main()
