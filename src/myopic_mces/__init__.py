import os
import sys

import rdkit # has to be loaded before cpp filters are
from rdkit import Chem, Geometry

from ._cpp import mmces_filters, HAS_CPP, CPP_IMPORT_ERROR
# only report once: the env var is inherited by worker processes (joblib/multiprocessing)
if not HAS_CPP and not os.environ.get('MYOPIC_MCES_NO_CPP_HINT'):
    os.environ['MYOPIC_MCES_NO_CPP_HINT'] = '1'
    if isinstance(CPP_IMPORT_ERROR, ModuleNotFoundError): # extra simply not installed
        print('myopic_mces: much faster computations possible with C++ bounds '
              '(>250x speedup for bound2): pip install "myopic_mces[cpp]"', file=sys.stderr)
    else: # installed, but broken or version mismatch
        import warnings
        warnings.warn(f'C++ filters not available: {CPP_IMPORT_ERROR}')

from .myopic_mces import MCES
from .MCES_ILP import MCES_ILP
from .filter_MCES import apply_filter, filter1, filter2
from .filter_MCES import apply_CPP_filter
from .graph import construct_graph

__all__ = ['MCES', 'MCES_ILP', 'apply_filter', 'filter1', 'filter2', 'apply_CPP_filter', 'construct_graph']
