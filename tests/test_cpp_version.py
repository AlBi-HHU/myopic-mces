import sys
import types

import pytest

from myopic_mces import _cpp
# this is testing the version check of the C++ extension against the pinned `cpp` extra


@pytest.fixture
def fake_cpp_module(monkeypatch):
    # works without the C++ extension being installed
    monkeypatch.setitem(sys.modules, 'mmces_filters', types.ModuleType('mmces_filters'))


def test_required_cpp_version_parsing(monkeypatch):
    monkeypatch.setattr(_cpp, 'requires', lambda _: ["rdkit>=2022.09.5", "myopic-mces-cpp==1.3.2; extra == 'cpp'"])
    assert _cpp._required_cpp_version() == '1.3.2'


def test_matching_version(monkeypatch, fake_cpp_module):
    monkeypatch.setattr(_cpp, '_required_cpp_version', lambda: '1.3.2')
    monkeypatch.setattr(_cpp, '_installed_cpp_version', lambda: '1.3.2')
    module, error = _cpp.load_cpp_filters()
    assert module is not None and error is None


def test_mismatching_version(monkeypatch, fake_cpp_module):
    monkeypatch.setattr(_cpp, '_required_cpp_version', lambda: '1.3.3')
    monkeypatch.setattr(_cpp, '_installed_cpp_version', lambda: '1.3.2')
    module, error = _cpp.load_cpp_filters()
    assert module is None
    assert isinstance(error, ImportError)
    assert '1.3.2' in str(error) and '1.3.3' in str(error)


def test_unknown_version_is_accepted(monkeypatch, fake_cpp_module):
    # e.g. running from source or a manually built module: no metadata -> no check
    monkeypatch.setattr(_cpp, '_required_cpp_version', lambda: '1.3.2')
    monkeypatch.setattr(_cpp, '_installed_cpp_version', lambda: None)
    module, error = _cpp.load_cpp_filters()
    assert module is not None and error is None
