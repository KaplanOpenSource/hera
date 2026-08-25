"""The heavy simulation modules, reached through the stub layer.

These are the ones the design deferred to last because they import OpenFOAM,
paraview, hermes and evtk at module level.  The Phase 0 stub layer exists so
they can be imported anyway; this file is the check that it actually delivers
that, since every future test of those modules depends on it.
"""
import importlib
import pathlib

import pytest

STUBBED_DEPENDENTS = [
    "hera.simulations.openFoam.toolkit",
    "hera.simulations.openFoam.OFWorkflow",
    "hera.simulations.LSM.toolkit",
    "hera.simulations.LSM.template",
    "hera.simulations.LSM.hermesWorkflowToolkit",
    "hera.simulations.hermesWorkflowToolkit",
    "hera.simulations.WRF",
]


@pytest.mark.unit
@pytest.mark.parametrize("moduleName", STUBBED_DEPENDENTS)
def test_the_stub_layer_lets_it_import(moduleName):
    assert importlib.import_module(moduleName) is not None


@pytest.mark.unit
class TestTheStubsAreDoingTheWork:
    def test_hermes_is_present_as_a_stub_or_a_real_module(self):
        import sys

        assert "hermes" in sys.modules

    def test_pyfoam_is_stubbed_as_a_namespace_package(self):
        import sys

        assert sys.modules["PyFoam"].__path__ == []

    def test_the_openfoam_toolkit_class_is_reachable(self):
        """Importing is not enough; the class has to be usable as a symbol."""
        from hera.simulations.openFoam.toolkit import OFToolkit

        assert OFToolkit is not None

    def test_the_lsm_template_class_is_reachable(self):
        from hera.simulations.LSM.template import LSMTemplate

        assert callable(LSMTemplate.prepareParams)


@pytest.mark.unit
class TestTorchIsDeliberatelyNotStubbed:
    """Phase 0 recorded this exclusion; the reason is checked here."""

    def test_the_stub_module_says_so(self):
        source = pathlib.Path("hera/tests/unit/_stubs.py").read_text(encoding="utf-8")
        assert "Deliberately NOT stubbed: ``torch``" in source

    def test_a_leaf_magicmock_cannot_carry_submodules(self):
        """The real reason, checked rather than asserted.

        A MagicMock in sys.modules has no __path__, so `import torch.utils`
        raises "'torch' is not a package".  Stubbing torch would need the
        namespace-package form PyFoam gets, not a leaf mock.
        """
        import sys
        from unittest.mock import MagicMock

        sys.modules.setdefault("heraProbePackage", MagicMock())
        with pytest.raises(ModuleNotFoundError, match="is not a package"):
            __import__("heraProbePackage.submodule")

    def test_the_namespace_form_does_carry_submodules(self):
        """Which is exactly what PyFoam gets, and why that one works."""
        import sys

        assert hasattr(sys.modules["PyFoam"], "__path__")
        assert "PyFoam.Basics.DataStructures" in sys.modules

    def test_the_module_is_skipped_rather_than_failed_when_torch_is_absent(self):
        pytest.importorskip("torch")
        assert importlib.import_module(
            "hera.simulations.machineLearningDeepLearning.toolkit"
        ) is not None


@pytest.mark.unit
class TestWRFImportMessage:
    def test_it_imports_without_python_wrf(self):
        """It prints a notice rather than raising, so callers can still import."""
        assert importlib.import_module("hera.simulations.WRF") is not None

    def test_the_notice_is_a_bare_print_not_a_warning(self):
        """Characterisation: a warnings.warn would be filterable and routable.

        Two bare prints in wrfDatalayer.py put text on stdout at import time,
        which is the same pattern as hill2stl (B42) in a milder form.
        """
        source = pathlib.Path("hera/simulations/WRF/wrfDatalayer.py").read_text(
            encoding="utf-8", errors="replace"
        )
        assert 'print("You must install python-wrf' in source
        assert "warnings.warn" not in source
