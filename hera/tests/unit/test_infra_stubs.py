"""The stub layer must let hera's externally-dependent modules import."""
import sys

import pytest


@pytest.mark.unit
def test_stub_modules_are_registered():
    """install() is called by conftest at collection time."""
    for name in ("PyFoam", "paraview", "evtk", "PyFoam.Basics.DataStructures"):
        assert name in sys.modules, f"{name} was not stubbed"


@pytest.mark.unit
def test_namespace_stub_supports_submodule_import():
    """A namespace stub must carry __path__ or `import a.b` fails."""
    assert sys.modules["PyFoam"].__path__ == []


@pytest.mark.unit
def test_install_is_idempotent():
    from hera.tests.unit import _stubs

    before = sys.modules["PyFoam"]
    _stubs.install()
    assert sys.modules["PyFoam"] is before, "install() replaced an existing stub"


@pytest.mark.unit
def test_openfoam_toolkit_imports_under_stubs():
    """The whole point: a module that needs PyFoam+evtk must import."""
    from hera.simulations.openFoam.toolkit import OFToolkit

    assert OFToolkit is not None
