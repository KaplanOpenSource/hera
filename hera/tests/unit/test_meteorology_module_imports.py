"""Which meteorology modules can be imported.

Three of them used to fail for reasons that had nothing to do with a missing
optional dependency: one named a datalayer class that no longer exists (B40),
one read a hard-coded absolute path belonging to a particular developer's home
directory at import time (B41), and one ran module-level scratch code (B333).
"""
import importlib
import pathlib

import pytest

IMPORTABLE = [
    "hera.measurements.meteorology.analysis",
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
class TestRadiosonde:
    """B40, fixed: it used to subclass a datalayer class that was removed."""

    def test_the_named_base_class_is_gone(self):
        import hera.datalayer as datalayer

        assert not hasattr(datalayer, "ProjectMultiDBPublic")

    def test_the_source_no_longer_refers_to_it(self):
        source = pathlib.Path(
            "hera/measurements/meteorology/radiosonde.py"
        ).read_text(encoding="utf-8")
        assert "ProjectMultiDBPublic" not in source

    def test_the_module_imports(self):
        assert importlib.import_module(
            "hera.measurements.meteorology.radiosonde"
        ) is not None

    def test_the_datalayer_is_a_project(self):
        import hera.datalayer as datalayer
        from hera.measurements.meteorology import radiosonde

        assert issubclass(radiosonde.DataLayer, datalayer.Project)


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
