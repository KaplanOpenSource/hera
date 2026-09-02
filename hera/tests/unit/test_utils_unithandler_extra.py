"""utils/unitHandler.py: the unum<->pint bridge helpers not covered
elsewhere. All deprecated but still live code paths."""
import warnings

import pytest
import unum.units as u

from hera.utils import ureg
from hera.utils import unitHandler as uh


@pytest.fixture(autouse=True)
def _ignore_deprecation():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


@pytest.mark.unit
class TestUnumToStr:
    def test_a_unum_object_is_rendered_as_value_star_unit(self):
        assert uh.unumToStr(5 * u.m) == "5*m"

    def test_a_plain_value_is_stringified_directly(self):
        assert uh.unumToStr(42) == "42"


@pytest.mark.unit
class TestStrToUnum:
    def test_a_unum_object_passes_through_unchanged(self):
        val = 5 * u.m
        assert uh.strToUnum(val) is val

    def test_a_non_unum_value_passes_through_unchanged(self):
        assert uh.strToUnum("5 m") == "5 m"


@pytest.mark.unit
class TestExtractUnumUnitsFromPint:
    def test_it_returns_the_unum_equivalent_unit(self):
        result = uh.extractUnumUnitsFromPint(1 * ureg.m)
        assert str(result) == "1.0 [m]"

    def test_an_unmapped_unit_raises(self):
        with pytest.raises(ValueError, match="not mapped to Unum"):
            uh.extractUnumUnitsFromPint(1 * ureg.knot)


@pytest.mark.unit
class TestUnumToBaseUnits:
    def test_it_converts_to_mks_base_units(self):
        result = uh.unumToBaseUnits(5 * u.km)
        assert str(result) == "5000.0 [m]"
