import rdkit # has to be loaded before cpp filters are

try: 
    import mmces_filters
    HAS_CPP = True
except ImportError as e:
    mmces_filters = None
    HAS_CPP = False
    import warnings
    warnings.warn(f'C++ filters not available: {e}')

from .myopic_mces import MCES
from .MCES_ILP import MCES_ILP
from .filter_MCES import apply_filter, filter1, filter2
from .filter_MCES import apply_CPP_filter
from .graph import construct_graph

__all__ = ['MCES', 'MCES_ILP', 'apply_filter', 'filter1', 'filter2', 'apply_CPP_filter', 'construct_graph']
