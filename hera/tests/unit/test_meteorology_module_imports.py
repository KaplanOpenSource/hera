"""Which meteorology modules can be imported."""
import importlib
import pathlib

import pytest

IMPORTABLE = [
    "hera.measurements.meteorology.lowfreqdata.analysis",
    "hera.measurements.meteorology.lowfreqdata.toolkit",
    "hera.measurements.meteorology.highfreqdata.toolkit",
    "hera.measurements.meteorology.highfreqdata.analysis.turbulencestatistics",
    "hera.measurements.meteorology.highfreqdata.analysis.meandatacalculator",
    "hera.measurements.meteorology.highfreqdata.analysis.abstractcalculator",
    "hera.measurements.meteorology.highfreqdata.parsers.TOA5",
    "hera.measurements.meteorology.highfreqdata.parsers.CampbellBinary",
]


@pytest.mark.unit
@pytest.mark.parametrize("moduleName", IMPORTABLE)
def test_the_module_imports(moduleName):
    assert importlib.import_module(moduleName) is not None


@pytest.mark.unit
class TestGFS:
    """Not a defect: scikit-learn is declared, just absent from this env.

    GFS.py also imports `from osgeo import gdal` at module level -- a
    second, independent heavy dependency the module needs regardless of
    scikit-learn's presence.
    """

    def test_scikit_learn_is_declared_in_requirements(self):
        requirements = pathlib.Path("requirements.txt").read_text(encoding="utf-8")
        assert "scikit-learn" in requirements

    def test_the_module_imports_when_the_dependency_is_present(self):
        pytest.importorskip("sklearn")
        pytest.importorskip("osgeo")
        assert importlib.import_module(
            "hera.measurements.meteorology.GFS"
        ) is not None
