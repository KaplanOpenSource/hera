"""The experiment parser module, and the pydoc factory convention it would use.

CLAUDE.md requires that Calculator, InjuryLevel, Injury and Action subclasses
carry pydoc.locate-compatible names, because the factories build a class path
from a string.  That convention is checked here, together with the state of
hera/measurements/experiment/parsers.py -- which turns out to be orphaned.
"""
import pathlib
import pydoc

import pytest

from hera.measurements.experiment.parsers import (
    Parser_CampbellBinary,
    Parser_OldStyleMetaDataParquet,
    Parser_TOA5,
)

CALCULATOR_PATH = "hera.riskassessment.agents.effects.Calculator.Calculator%s"


@pytest.mark.unit
class TestPydocFactoryConvention:
    """The naming rule CLAUDE.md enforces, exercised the way the factory does."""

    @pytest.mark.parametrize("suffix", ["Haber", "TenBerge", "MaxConcentration"])
    def test_each_calculator_resolves_by_name(self, suffix):
        assert pydoc.locate(CALCULATOR_PATH % suffix) is not None

    def test_an_unknown_name_resolves_to_none_rather_than_raising(self):
        """Injury.py checks the result, so None is the contract, not an error."""
        assert pydoc.locate(CALCULATOR_PATH % "NoSuchThing") is None

    def test_the_resolved_class_is_the_real_one(self):
        from hera.riskassessment.agents.effects.Calculator import CalculatorHaber

        assert pydoc.locate(CALCULATOR_PATH % "Haber") is CalculatorHaber

    def test_the_resolved_class_is_constructible(self):
        """Resolution is only useful if the class can then be used."""
        calculator = pydoc.locate(CALCULATOR_PATH % "Haber")()
        assert calculator.injuryBreathingRate is not None


@pytest.mark.unit
class TestParserModuleIsOrphaned:
    """Nothing imports this module, and its dispatch code is commented out."""

    SOURCE_ROOT = pathlib.Path("hera")

    def test_no_module_imports_it(self):
        importers = []
        for path in self.SOURCE_ROOT.rglob("*.py"):
            if "experiment/parsers.py" in str(path) or "tests" in str(path):
                continue
            text = path.read_text(encoding="utf-8", errors="replace")
            if "experiment.parsers" in text or "experiment import parsers" in text:
                importers.append(str(path))
        assert importers == []

    def test_the_dispatch_that_would_use_it_is_commented_out(self):
        """lowfreqdata/toolkit.py once built 'Parser_{name}' paths; now inert."""
        source = pathlib.Path(
            "hera/measurements/meteorology/lowfreqdata/toolkit.py"
        ).read_text(encoding="utf-8")
        assert "#         parserPath = f\"{className}.parsers.Parser_{parser}\"" in source

    def test_the_live_parsers_live_elsewhere(self):
        """highfreqdata/toolkit.py imports its parsers from its own package."""
        source = pathlib.Path(
            "hera/measurements/meteorology/highfreqdata/toolkit.py"
        ).read_text(encoding="utf-8")
        assert "from .parsers.CampbellBinary import Parser" in source
        assert "from .parsers.TOA5 import ASCIIParser" in source

    def test_the_classes_still_import_cleanly(self):
        """Orphaned, but not broken -- worth knowing before deleting."""
        assert Parser_OldStyleMetaDataParquet is not None
        assert Parser_CampbellBinary is not None
        assert Parser_TOA5 is not None


@pytest.mark.unit
class TestTOA5Stub:
    """B44: Parser_TOA5.parse is an empty body that silently returns None."""

    def test_the_class_can_be_constructed(self):
        assert Parser_TOA5() is not None

    def test_a_working_toa5_parser_exists_elsewhere(self):
        """So the stub is a duplicate, not the only implementation."""
        from hera.measurements.meteorology.highfreqdata.parsers.TOA5 import ASCIIParser

        assert hasattr(ASCIIParser, "getPandasFromFile")
        assert hasattr(ASCIIParser, "getPandasFromDir")

    @pytest.mark.xfail(
        strict=True,
        reason="B44: Parser_TOA5.parse has an empty body, so it returns None for "
               "any input instead of parsing or refusing. Its docstring documents "
               "the `file` parameter and promises nothing about the return, and a "
               "real TOA5 parser exists at "
               "meteorology/highfreqdata/parsers/TOA5.ASCIIParser. "
               "See the consolidated findings issue.",
    )
    def test_parsing_either_returns_data_or_refuses(self, tmp_path):
        """A stub should raise NotImplementedError, not hand back None."""
        target = tmp_path / "sample.dat"
        target.write_text("TOA5,station\n", encoding="utf-8")

        result = Parser_TOA5().parse(str(target))
        assert result is not None

    def test_it_returns_none_even_for_a_path_that_does_not_exist(self):
        """Confirms the body is unreachable code rather than a silent failure path."""
        assert Parser_TOA5().parse("/no/such/file.dat") is None
