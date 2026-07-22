# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:16:05 2020
@author: seipp
"""
import time
from joblib import Parallel, delayed
import numpy as np
import argparse
import sys
from myopic_mces.graph import construct_graph
from myopic_mces.MCES_ILP import MCES_ILP
from myopic_mces.filter_MCES import apply_filter, ComputationMode

def MCES(smiles1, smiles2, threshold=10, i=0, solver='COIN_CMD', solver_options={},
         no_ilp_threshold=False, always_stronger_bound=True, use_bound_zero=False,
         catch_errors=False, structure=False, num_smiles=False):
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
        ILP-solver used for solving MCES. Example: CPLEX_PY
    solver_options: dict
        additional options to pass to solvers. Example: threads=1 for better multi-threaded performance
    no_ilp_threshold: bool
        if true, always return exact distance even if it is below the threshold (slower)
    always_stronger_bound: bool
        if true, always compute and use the second stronger bound
    use_bound_zero : bool
        if true, also use an additional weak (molecular formula-based) filter
    catch_errors : bool
         if true, return distance -1 when errors are encountered

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
            6 : Timelimit reached, no solution found

    """
    start = time.time()
    # construct graph for both smiles. if we want to calculate the MCES structure mapping, we need save_mol
    G1 = construct_graph(smiles1, save_mol=structure, num_smiles=num_smiles)
    G2 = construct_graph(smiles2, save_mol=structure, num_smiles=num_smiles)
    if num_smiles:
        num_smiles1 = G1.num_smiles
        num_smiles2 = G2.num_smiles
    distance = np.inf
    compute_mode = ComputationMode.UNKNOWN
    if threshold != -1:         # with `-1` always compute exact distance
        # filter out if distance is above the threshold
        try:
            distance, compute_mode = apply_filter(G1, G2, threshold, always_stronger_bound=always_stronger_bound,
                                                  use_bound_zero=use_bound_zero)
            if distance > threshold:
                if structure:
                    if num_smiles:
                        return i, distance, time.time() - start, compute_mode, None, num_smiles1, num_smiles2
                    else:
                        return i, distance, time.time() - start, compute_mode, None
                return i, distance, time.time() - start, compute_mode
        except Exception as e:
            print('ERROR:', smiles1, smiles2, 'filter', e, file=sys.stderr)
            if (catch_errors):
                distance = -1
                compute_mode = ComputationMode.ABOVE_THRESHOLD.value
            else:
                raise e
    # store the filter results
    distance_filter = distance
    compute_mode_filter = compute_mode

    # calculate MCES
    try:
        if not structure:
            distance, compute_mode = MCES_ILP(G1, G2, threshold, solver, solver_options=solver_options,
                                            no_ilp_threshold=no_ilp_threshold, structure=structure)
        elif structure:
            distance, compute_mode, mapping = MCES_ILP(G1, G2, threshold, solver, solver_options=solver_options,
                                            no_ilp_threshold=no_ilp_threshold, structure=structure)
            
    except Exception as e:
        print('ERROR:', smiles1, smiles2, 'exact', e, file=sys.stderr)
        if (catch_errors):
            distance = -1
            compute_mode = ComputationMode.EXACT.value
        else:
            raise e
    # if ILP computation does not have a result because the time limit was reached, use the filter
    # supported for CPLEX and CBC
    if (threshold != -1):
        if (distance == -1 and compute_mode == ComputationMode.TIMEOUT_BOUND.value and distance_filter > 0):
            distance = distance_filter
            compute_mode = compute_mode_filter
    if num_smiles:
        if structure:
            return i, distance, time.time() - start, compute_mode, mapping, num_smiles1, num_smiles2
        else:
            return i, distance, time.time() - start, compute_mode, num_smiles1, num_smiles2
    elif structure:
        return i, distance, time.time() - start, compute_mode, mapping
    else: 
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
            f.create_dataset('computation_times', data=[row[2] if len(row) > 2 else -1 for row in results], compression='gzip')
        if (write_modes):
            f.create_dataset('computation_modes', data=[row[3] if len(row) > 3 else -1 for row in results], compression='gzip')
        comp_args = f.create_group('computation_args')
        for k, v in args:
            if v is not None:
                comp_args[k] = v
    print(f'done, took {(time.time() - t0) / 60:.1f}min')

def filter_inputs(inputs,dmatrix_file,threshold=None):
    from scipy.spatial.distance import squareform
    import h5py
    print('filtering for precomputed mces')
    t0 = time.time()
    filtered_inputs = []
    precomputed_mces = []
    with h5py.File(dmatrix_file,'r') as hf:
            mces = hf["mces"][:]
            all_smiles =hf['mces_smiles_order'][:]  #Decoding will result in not finding any value since hdf input process wont decode either! [s.decode() for s in hf['mces_smiles_order'][:]]
    mces = squareform(mces)
    smiles_index = {}
    for i, smiles in enumerate(all_smiles):
        smiles_index[smiles] = i
    for i, s1, s2 in inputs:
        idx1 = smiles_index.get(s1)
        idx2 = smiles_index.get(s2)
        
        if idx1 is not None and idx2 is not None:
            val = mces[idx1][idx2]
            if val != -1:
                if threshold is None or val < threshold:
                    precomputed_mces.append((i, val,-1,-1))
                    continue
        
        filtered_inputs.append((i, s1, s2))
    print(f"done, took {(time.time() - t0) / 60:.1f}min and found {len(precomputed_mces)} in libary, {len(filtered_inputs)} to go")
    return filtered_inputs, precomputed_mces

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input', help='input file in the format: index,smiles1,smiles2 OR hdf5 file when using hdf5 mode')
    parser.add_argument('output', help='output file in the format: index,myopic MCES distance,computation time in seconds,computation mode')

    # general options
    parser.add_argument('--threshold', type=float, default=10.,
                        action='store', help='threshold for the distance. exact distance is only calculated if the distance '\
                        'is lower than the threshold. if set to -1, exact distance is always calculated. ')
    parser.add_argument('--solver', type=str, default='COIN_CMD',
                        action='store', help='Solver for the ILP. example: CPLEX_PY')
    parser.add_argument('--num_jobs', type=int, help='Number of jobs; instances to run in parallel. '
                        'By default this is set to the number of (logical) CPU cores.', default=-1)
    parser.add_argument('--hdf5_mode', action='store_true',
                        help='more time and space efficient mode of input/output using the `input` hdf5-file. '
                        'Has to contain one nx3 array with indices (`computation_indices`): (id_, smiles1_i, smiles2_i) '
                        'and another array with the corresponding SMILES (`smiles`). Output will be written to the same file.')
    parser.add_argument('--hide_rdkit_warnings', action='store_true', help='attempts to suppress RDKit warning')
    
    # ilp solver options
    parser.add_argument('--solver_onethreaded', action='store_true',
                        help='limit ILP solver to one thread, resulting in faster '
                        'performance with parallel computations (not available for all solvers)')
    parser.add_argument('--solver_no_msg', action='store_true',
                        help='prevent solver from logging (not available for all solvers)')
    parser.add_argument('--solver_time_limit_seconds', type=float, default=None,
                        help='EXPERIMENTAL: set a time limit for the ILP solver. Computations below the threshold are not guaranteed to be exact anymore!'
                         ' Supported solver is CPLEX_PY, for others, correctness cannot be guaranteed.')


    # experimental options
    parser.add_argument('--no_ilp_threshold', action='store_true',
                        help='(experimental) if set, do not add threshold as constraint to ILP, '
                        'resulting in longer runtimes and potential violations of the triangle equation')
    parser.add_argument('--choose_bound_dynamically', action='store_true',
                        help='if this is set, compute and use potentially weaker but faster lower bound if '
                        'already greater than the threshold. Otherwise (default), the strongest lower bound '
                        'is always computed and used. Enabling this can lead to massive speedups.')
    parser.add_argument('--use_bound_zero', action='store_true',
                        help='if this is set, compute and use an additional weak molecular formula-based '
                        'lower bound. Use in conjunction with `choose_bound_dynamically`.')
    parser.add_argument('--catch_computation_errors', action='store_true', help='(experimental) instead of aborting '
                        'the computations, instances that failed to compute receive distance "-1"')
    parser.add_argument('--jobs_batch_size', type=int, default=32, help='(experimental) batch size for parallelization')
    parser.add_argument('--jobs_dispatch', default='10*n_jobs', help='(experimental) pre-dispatch of jobs for parallelization')
    parser.add_argument('--use_matrix_lookup', help='(experimental) Use with the '
                        'path to a HDF5 file with precomputed MCES distances. Computation for these instances will be '
                        'skipped, using the provided values. HDF5 has to contain distances (key `mces`) and SMILES '
                        '(`mces_smiles_order`), like the HDF5 files produced by this script. '
                        'NOTE: When used in combination with `prepare_input`, only use with `--no_shuffle`', action='store_true')
    parser.add_argument('--lookup_threshold', help='(experimental) Use with `--use_matrix_lookup`: '
                        'Precomputed values equal or greater than the threshold will be ignored; these '
                        'instances will be recomputed', default=None, type=float)
    parser.add_argument('--num_smiles', action='store_true', help='(experimental) if set, calculates and saves SMILES with numbered atoms ' \
                        'to the output for each input SMILES')
    parser.add_argument('--structure', action='store_true', help='(experimental) if set, saves a mapping of the pairwise MCES structure ' \
                        'to the output. Only works for molecules for which the MCES distance has been proven exact by the ILP.')
    args = parser.parse_args()

    if (args.hide_rdkit_warnings):
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
        # TODO: does not work with forks (preferred parallelization)

    additional_mces_options = dict(no_ilp_threshold=args.no_ilp_threshold, solver_options=dict(),
                                   always_stronger_bound=not args.choose_bound_dynamically,
                                   use_bound_zero=args.use_bound_zero,
                                   catch_errors=args.catch_computation_errors,
                                   num_smiles=args.num_smiles,
                                   structure=args.structure)
    
    # TODO: disclaimer that mapping does not work with hdf5 format yet
    if (args.solver != 'CPLEX_PY' and args.solver_time_limit_seconds is not None):
        print('Unsupported solver for timelimit flag, correctness not guaranteed!', file=sys.stderr)
        if (args.solver == 'COIN_CMD'):
            additional_mces_options['solver_options']['timeMode'] = 'cpu'

    if (args.solver_onethreaded):
        additional_mces_options['solver_options']['threads'] = 1
    if (args.solver_no_msg):
        additional_mces_options['solver_options']['msg'] = False
    if (args.solver_time_limit_seconds is not None):
        additional_mces_options['solver_options']['timeLimit'] = args.solver_time_limit_seconds

    if (args.hdf5_mode):
        inputs = hdf5_input(args.input)
    else:
        with open(args.input) as in_handle:
            inputs = [line.strip().split(',')[:3] for line in in_handle] # ignores extra input columns
    
    if args.use_matrix_lookup:
        inputs_to_process, results = filter_inputs(inputs=inputs, dmatrix_file=args.use_matrix_lookup,threshold=args.lookup_threshold)
    else:
        inputs_to_process, results = inputs, []

    results += Parallel(n_jobs=args.num_jobs, verbose=5, batch_size=args.jobs_batch_size, pre_dispatch=args.jobs_dispatch)(
            delayed(MCES)(smiles1, smiles2, args.threshold, i, args.solver, **additional_mces_options) for i, smiles1, smiles2 in inputs_to_process)

    if args.use_matrix_lookup:
        results.sort(key=lambda x: x[0])

    if (args.hdf5_mode):
        hdf5_output(results, args.input, write_times=True, write_modes=True, args=args._get_kwargs())
    else:
        with open(args.output, 'w') as out_handle:
            if args.structure:

                if mapping is not None:
                        mapping_str = repr(mapping).replace(', ', ',').replace(': ', ':')
                else:
                    mapping_str = ''

                if args.num_smiles:
                    for ind, distance, duration, compute_mode, mapping, num_smiles1, num_smiles2 in results:
                        out_handle.write(
                            f'{ind},{distance},{duration},{compute_mode},{mapping_str},{num_smiles1},{num_smiles2}\n'
                            )
                else:
                    for ind, distance, duration, compute_mode, mapping in results:
                        out_handle.write(
                            f'{ind},{distance},{duration},{compute_mode},{mapping_str}\n'
                            )
            elif args.num_smiles:                        
                for ind, distance, duration, compute_mode, num_smiles1, num_smiles2 in results:
                        out_handle.write(
                            f'{ind},{distance},{duration},{compute_mode},{num_smiles1},{num_smiles2}\n'
                            )
            else:
                for ind, distance, duration, compute_mode in results:
                        out_handle.write(
                            f'{ind},{distance},{duration},{compute_mode}\n'
                            )


if __name__ == '__main__':
    main()
