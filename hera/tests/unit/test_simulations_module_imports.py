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

# Blocked by B31: `from ..utils import tonumber, tounit` resolves to
# hera.simulations.utils, whose __init__ is empty.  The names live in
# hera.utils, three levels up.
BROKEN_BY_RELATIVE_IMPORT = [
    "hera.simulations.evaporation",
    "hera.simulations.evaporation.models",
    "hera.simulations.evaporation.monaghan",
    "hera.simulations.deposition",
    "hera.simulations.deposition.models",
]


@pytest.mark.unit
@pytest.mark.parametrize("moduleName", IMPORTABLE)
def test_the_module_imports(moduleName):
    assert importlib.import_module(moduleName) is not None


@pytest.mark.unit
class TestTheBrokenRelativeImport:
    """B31, isolated to its mechanism so the fix is unambiguous."""

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

    def test_a_sibling_module_shows_the_correct_form(self):
        """DropletCloud.py uses three dots and works; these three use two."""
        source = __import__("pathlib").Path(
            "hera/simulations/gaussian/DropletCloud.py"
        ).read_text(encoding="utf-8")
        assert "from ...utils import tounit,tonumber" in source

    @pytest.mark.xfail(
        strict=True,
        reason="B31: evaporation/models.py:3, evaporation/monaghan.py:5 and "
               "deposition/models.py:2 do `from ..utils import tonumber, tounit`, "
               "which resolves to hera.simulations.utils -- an empty package. The "
               "names are in hera.utils, so the import needs three dots, as "
               "DropletCloud.py already uses. Both packages are therefore "
               "unimportable in any installation. "
               "See the consolidated findings issue.",
    )
    @pytest.mark.parametrize("moduleName", BROKEN_BY_RELATIVE_IMPORT)
    def test_the_physics_packages_import(self, moduleName):
        assert importlib.import_module(moduleName) is not None
