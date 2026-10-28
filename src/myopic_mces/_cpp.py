"""
Import of the optional C++ filters (`myopic_mces_cpp`, module `mmces_filters`).

`pip install -U myopic_mces` (without the `[cpp]` extra) does not update the
C++ extension, so the installed version is checked against the one pinned in
the `cpp` extra. On a mismatch the C++ filters are disabled (python fallback).
"""

import re
from importlib.metadata import version, requires, PackageNotFoundError

import rdkit # has to be loaded before cpp filters are


def _required_cpp_version():
    """version of myopic_mces_cpp pinned in the `cpp` extra, None if not determinable"""
    try:
        reqs = requires('myopic_mces') or []
    except PackageNotFoundError: # e.g. run from source without installation
        return None
    for req in reqs:
        m = re.match(r'myopic[-_.]mces[-_.]cpp\s*==\s*([^\s;,]+)', req, re.IGNORECASE)
        if m:
            return m.group(1)
    return None


def _installed_cpp_version():
    try:
        return version('myopic_mces_cpp')
    except PackageNotFoundError: # e.g. module built manually and not installed as package
        return None


def load_cpp_filters():
    """returns (module, None) if the C++ filters are usable, otherwise (None, ImportError)"""
    try:
        import mmces_filters
    except ImportError as e:
        return None, e

    required = _required_cpp_version()
    installed = _installed_cpp_version()
    if required is not None and installed is not None and installed != required:
        return None, ImportError(f'myopic_mces_cpp {installed} is installed, but this version of myopic_mces '
                                 f'requires myopic_mces_cpp {required}. '
                                 'Update with: pip install -U "myopic_mces[cpp]"')
    return mmces_filters, None


mmces_filters, CPP_IMPORT_ERROR = load_cpp_filters()
HAS_CPP = mmces_filters is not None
