"""Which simulation modules can be imported at all.

Half of this batch's nominal scope -- evaporation and deposition -- turned out
to be unreachable, so importability itself is what these tests assert.  A
module that cannot be imported has no coverage to gain and no behaviour to
test; the gap is a defect, not a testing shortfall.
"""
import importlib

import pytest

IMPORTABLE = [
    "hera.simulations.hydrodynamics.nearWallFlow",
    "hera.simulations.windProfile.toolkit",
    "hera.simulations.utils.canopyWindProfile",
    "hera.simulations.utils.interpolations",
    "hera.simulations.utils.coordinateHandler",
    "hera.simulations.utils.inputForModelsCreation",
    "hera.simulations.analysis.errorCalculation",
]

# B31, fixed: these did `from ..utils import tonumber, tounit`, which resolved
# to hera.simulations.utils -- an empty package.  The names live in hera.utils,
# and the imports now say so.
FIXED_BY_THE_IMPORT_CORRECTION = [
    "hera.simulations.evaporation",
    "hera.simulations.evaporation.models",
    "hera.simulations.deposition",
    "hera.simulations.deposition.models",
]

# monaghan.py needs pyriskassessment, a private package declared in
# requirements.txt:372 and installed from an internal index.  Absent here, so
# skipped rather than failed -- the same treatment GFS.py gets for sklearn.
NEEDS_A_PRIVATE_PACKAGE = ["hera.simulations.evaporation.monaghan"]


@pytest.mark.unit
@pytest.mark.parametrize("moduleName", IMPORTABLE)
def test_the_module_imports(moduleName):
    assert importlib.import_module(moduleName) is not None


@pytest.mark.unit
class TestTheRelativeImportFix:
    """B31: the mechanism, kept so a regression is recognisable."""

    def test_the_names_live_in_the_top_level_utils(self):
        """Where the three broken modules should have imported from."""
        from hera.utils import tonumber, tounit

        assert callable(tonumber)
        assert callable(tounit)

    def test_the_simulations_utils_package_exports_nothing(self):
        """hera/simulations/utils/__init__.py is empty, so `..utils` has no names."""
        import hera.simulations.utils as simulationsUtils

        assert not hasattr(simulationsUtils, "tonumber")
        assert not hasattr(simulationsUtils, "tounit")

    def test_no_module_still_uses_the_wrong_level(self):
        """Two dots would resolve to the empty hera.simulations.utils again."""
        import pathlib

        offenders = [
            str(path)
            for path in pathlib.Path("hera/simulations").rglob("*.py")
            if "from ..utils import" in path.read_text(errors="replace")
        ]
        assert offenders == []

    @pytest.mark.parametrize("moduleName", FIXED_BY_THE_IMPORT_CORRECTION)
    def test_the_physics_packages_import(self, moduleName):
        assert importlib.import_module(moduleName) is not None

    @pytest.mark.parametrize("moduleName", NEEDS_A_PRIVATE_PACKAGE)
    def test_the_rest_imports_once_its_private_dependency_is_present(self, moduleName):
        pytest.importorskip("pyriskassessment")
        assert importlib.import_module(moduleName) is not None
