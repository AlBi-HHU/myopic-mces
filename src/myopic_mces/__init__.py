import rdkit # has to be loaded before cpp filters are
from rdkit import Chem, Geometry

from ._cpp import mmces_filters, HAS_CPP, CPP_IMPORT_ERROR
if not HAS_CPP:
    import warnings
    warnings.warn(f'C++ filters not available: {CPP_IMPORT_ERROR}')

from .myopic_mces import MCES
from .MCES_ILP import MCES_ILP
from .filter_MCES import apply_filter, filter1, filter2
from .filter_MCES import apply_CPP_filter
from .graph import construct_graph

__all__ = ['MCES', 'MCES_ILP', 'apply_filter', 'filter1', 'filter2', 'apply_CPP_filter', 'construct_graph']
