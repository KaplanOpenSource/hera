"""Unit handling: pint registry, custom units, and the unum bridge.

Every conversion below is checked against its definition -- a dunam is
1000 m2 by definition, mmH2O is defined in the module as atm/10197.162129779,
0 degC is 273.15 K -- rather than against what the module returns.
"""
import pytest

from hera.utils import unitHandler as uh
from hera.utils.unitHandler import Quantity, tonumber, tounit, ureg


@pytest.mark.unit
class TestCustomUnits:
    """The two units hera defines on top of pint's registry."""

    def test_dunam_is_a_thousand_square_metres(self):
        assert (1 * ureg.dunam).to(ureg.m**2).magnitude == pytest.approx(1000.0)

    def test_dunam_has_area_dimensions(self):
        assert (1 * ureg.dunam).check("[length]**2")

    def test_one_atmosphere_in_mmh2o(self):
        """Defined as atm / 10197.162129779, so the inverse must hold exactly."""
        assert (1 * ureg.atm).to(ureg.mmH2O).magnitude == pytest.approx(10197.162129779)

    def test_mmh2o_has_pressure_dimensions(self):
        assert (1 * ureg.mmH2O).check("[pressure]")

    def test_ten_dunam_is_a_hectare(self):
        assert (10 * ureg.dunam).to(ureg.hectare).magnitude == pytest.approx(1.0)


@pytest.mark.unit
class TestRegistryBasics:
    def test_length_conversion(self):
        assert (1 * ureg.km).to(ureg.m).magnitude == pytest.approx(1000.0)

    def test_time_conversion(self):
        assert (1 * ureg.hour).to(ureg.s).magnitude == pytest.approx(3600.0)

    def test_incompatible_conversion_raises(self):
        from pint.errors import DimensionalityError

        with pytest.raises(DimensionalityError):
            (1 * ureg.m).to(ureg.s)

    def test_unknown_unit_raises(self):
        from pint.errors import UndefinedUnitError

        with pytest.raises(UndefinedUnitError):
            ureg.parse_expression("bananas").to(ureg.m)


@pytest.mark.unit
class TestOffsetUnits:
    """Temperature needs Quantity(), not multiplication.

    Pinned because this is the trap for callers of this module: `5 * ureg.degC`
    looks like every other unit expression and raises instead.
    """

    def test_multiplication_by_an_offset_unit_is_ambiguous(self):
        from pint.errors import OffsetUnitCalculusError

        with pytest.raises(OffsetUnitCalculusError):
            (0 * ureg.degC).to(ureg.kelvin)

    def test_quantity_constructor_converts_correctly(self):
        assert Quantity(0, "degC").to("kelvin").magnitude == pytest.approx(273.15)

    def test_absolute_zero(self):
        assert Quantity(-273.15, "degC").to("kelvin").magnitude == pytest.approx(0.0)

    def test_kelvin_multiplication_is_fine(self):
        """Kelvin has no offset, so the usual form works."""
        assert (300 * ureg.kelvin).to(ureg.kelvin).magnitude == pytest.approx(300.0)


@pytest.mark.unit
class TestToNumber:
    """Strips units, returning a plain number in the requested unit."""

    def test_pint_quantity_is_converted_then_stripped(self):
        assert tonumber(1 * ureg.km, ureg.m) == pytest.approx(1000.0)

    def test_a_plain_number_passes_through_untouched(self):
        """No units to strip, so the value is returned as-is."""
        assert tonumber(5, ureg.m) == 5

    def test_result_is_not_a_quantity(self):
        assert not hasattr(tonumber(1 * ureg.km, ureg.m), "magnitude")

    @pytest.mark.skipif(not uh.unumSupport, reason="unum is not installed")
    def test_unum_quantity_is_converted(self):
        import unum.units as unum_units

        assert tonumber(1000 * unum_units.m, unum_units.km) == pytest.approx(1.0)


@pytest.mark.unit
class TestToUnit:
    """Attaches or converts units, keeping the value a Quantity."""

    def test_plain_number_and_string_unit_become_a_quantity(self):
        result = tounit(5, "m")
        assert isinstance(result, Quantity)
        assert result.magnitude == pytest.approx(5.0)

    @pytest.mark.xfail(
        strict=True,
        reason="B12: tounit builds Quantity(x, unit) from the bare pint import "
               "rather than ureg.Quantity, so the result belongs to pint's default "
               "application registry. See the consolidated findings issue.",
    )
    def test_result_belongs_to_heras_registry(self):
        """A value produced by tounit has to be usable with hera's units.

        Two consequences of the current behaviour, both verified:
          * tounit(5, "m") + 5 * ureg.m raises
            "Cannot operate with Quantity and Quantity of different registries"
          * tounit(1, "dunam") raises UndefinedUnitError -- the custom unit
            hera's registry exists to provide is exactly the one it cannot make
        """
        result = tounit(5, "m")
        assert result._REGISTRY is ureg
        assert (result + 5 * ureg.m).to(ureg.m).magnitude == pytest.approx(10.0)

    @pytest.mark.xfail(
        strict=True,
        reason="B12: tounit uses the default pint registry, which has no 'dunam'. "
               "See the consolidated findings issue.",
    )
    def test_custom_units_are_reachable_through_tounit(self):
        assert tounit(1, "dunam").to(ureg.m**2).magnitude == pytest.approx(1000.0)

    def test_quantity_is_converted(self):
        result = tounit(1 * ureg.km, ureg.m)
        assert result.magnitude == pytest.approx(1000.0)

    def test_tounum_is_an_alias_of_tounit(self):
        assert uh.tounum is tounit

    def test_incompatible_conversion_raises(self):
        from pint.errors import DimensionalityError

        with pytest.raises(DimensionalityError):
            tounit(1 * ureg.m, ureg.s)


@pytest.mark.unit
@pytest.mark.skipif(not uh.unumSupport, reason="unum is not installed")
class TestUnumBridge:
    def test_unum_to_pint_preserves_the_value(self):
        import unum.units as unum_units

        converted = uh.unumToPint(1000 * unum_units.m)
        assert isinstance(converted, Quantity)
        assert converted.to(ureg.km).magnitude == pytest.approx(1.0)

    def test_a_mapped_unit_round_trips(self):
        import unum.units as unum_units

        original = 5 * unum_units.m
        assert uh.unumToPint(original).to(ureg.m).magnitude == pytest.approx(5.0)

    def test_an_unmapped_pint_unit_raises_rather_than_guessing(self):
        """PINT_TO_UNUM_MAP is explicit; a miss must be loud, not approximate."""
        with pytest.raises(ValueError, match="not mapped to Unum"):
            uh.pintToUnum(1 * ureg.km)

    def test_a_mapped_pint_unit_converts(self):
        assert uh.pintToUnum(5 * ureg.m) is not None


@pytest.mark.unit
class TestTemperatureConstants:
    """celsius and K are module-level constants used by older call sites."""

    def test_both_constants_exist(self):
        assert uh.celsius is not None
        assert uh.K is not None

    @pytest.mark.xfail(
        strict=True,
        reason="B8: the comment above the override states 'celsius and K are NOT "
               "overridden -- pint versions are kept', and the code immediately "
               "below overrides both with unum objects. The comment and the code "
               "disagree. See the consolidated findings issue.",
    )
    @pytest.mark.skipif(not uh.unumSupport, reason="the conflict needs unum present")
    def test_constants_are_the_pint_units_the_comment_promises(self):
        assert type(uh.celsius).__module__.startswith("pint")
        assert type(uh.K).__module__.startswith("pint")
