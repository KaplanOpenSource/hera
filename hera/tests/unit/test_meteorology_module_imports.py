"""Which meteorology modules can be imported.

Two of them cannot, for reasons that have nothing to do with a missing
optional dependency: one names a datalayer class that no longer exists, and
one reads a hard-coded absolute path belonging to a particular developer's
home directory at import time.
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
    """B40: it subclasses a datalayer class that was removed."""

    def test_the_named_base_class_is_gone(self):
        import hera.datalayer as datalayer

        assert not hasattr(datalayer, "ProjectMultiDBPublic")

    def test_the_source_still_refers_to_it(self):
        source = pathlib.Path(
            "hera/measurements/meteorology/radiosonde.py"
        ).read_text(encoding="utf-8")
        assert "ProjectMultiDBPublic" in source

    @pytest.mark.xfail(
        strict=True,
        reason="B40: radiosonde.py declares "
               "`class DataLayer(datalayer.ProjectMultiDBPublic)`, but that class "
               "no longer exists in hera.datalayer, so the module raises "
               "AttributeError on import. It is dead in any installation. "
               "See the consolidated findings issue.",
    )
    def test_the_module_imports(self):
        assert importlib.import_module(
            "hera.measurements.meteorology.radiosonde"
        ) is not None


@pytest.mark.unit
class TestHighFreqMain:
    """B41, fixed: the personal absolute path is gone and the module imports.

    It read /home/ilay/hera_unittest_data/... at module level, so importing it
    raised FileNotFoundError anywhere but one machine.  The path is now an
    argument and the read happens under __main__.
    """

    OLD_HARDCODED = "/home/ilay/hera_unittest_data"

    def test_no_absolute_path_remains(self):
        source = pathlib.Path(
            "hera/measurements/meteorology/highfreqdata/__main__.py"
        ).read_text(encoding="utf-8")
        assert self.OLD_HARDCODED not in source

    def test_the_work_is_behind_a_main_guard(self):
        source = pathlib.Path(
            "hera/measurements/meteorology/highfreqdata/__main__.py"
        ).read_text(encoding="utf-8")
        assert 'if __name__ == "__main__":' in source

    def test_it_reports_usage_when_given_no_path(self):
        from hera.measurements.meteorology.highfreqdata.__main__ import main

        with pytest.raises(SystemExit, match="usage:"):
            main([] and None or __import__("sys").argv.clear() or None)

    def test_the_module_imports(self):
        assert importlib.import_module(
            "hera.measurements.meteorology.highfreqdata.__main__"
        ) is not None


@pytest.mark.unit
class TestGFS:
    """Not a defect: scikit-learn is declared, just absent from this env."""

    def test_scikit_learn_is_declared_in_requirements(self):
        requirements = pathlib.Path("requirements.txt").read_text(encoding="utf-8")
        assert "scikit-learn" in requirements

    def test_the_module_imports_when_the_dependency_is_present(self):
        pytest.importorskip("sklearn")
        assert importlib.import_module(
            "hera.measurements.meteorology.GFS"
        ) is not None
