r"""Atmospheric profiles for the Gaussian dispersion toolkit.

Each quantity is checked against the physics it claims to implement, using
values that can be looked up independently:

* pressure   -- barometric, 760 mmHg at sea level, e-folding height 1/1.186e-4
* temperature-- the standard environmental lapse rate, 6.5 K per km
* density    -- dry air at 20 C is 1.204 kg/m3
* viscosity  -- dry air at 20 C is about 1.81e-5 Pa s
* wind       -- a log profile must reproduce u_ref exactly at the reference
                height, and vanish at the roughness length
"""
import warnings

import numpy as np
import pytest

from hera.simulations.gaussian.Meteorology import (
    MeteorologyFactory,
    StandardMeteorolgyConstant_log,
    StandardMeteorolgyConstant_powerLaw,
    StandardMeteorolgyConstant_uniformWind,
)
from hera.utils.unitHandler import ureg

STABILITIES = ["A", "B", "C", "D", "E", "F"]


@pytest.fixture()
def factory():
    return MeteorologyFactory()


@pytest.fixture()
def meteorology(factory):
    """Log profile, 5 m/s at 10 m, 20 C at the ground, z0 = 10 cm."""
    return factory.getMeteorologyFromU10(
        u10=5 * ureg.m / ureg.s,
        inversion=1000 * ureg.m,
        verticalProfileType="log",
    )


@pytest.mark.unit
class TestAirPressure:
    def test_sea_level_pressure_is_one_atmosphere(self, meteorology):
        pressure = meteorology.getAirPressure(0 * ureg.m)
        assert pressure.to(ureg.mmHg).magnitude == pytest.approx(760.0)
        assert pressure.to(ureg.atm).magnitude == pytest.approx(1.0, rel=1e-6)

    def test_the_e_folding_height_matches_the_coefficient(self, meteorology):
        """At z = 1/1.186e-4 = 8432 m the pressure must fall to 1/e."""
        scale_height = 1.0 / 1.186e-4
        ratio = (
            meteorology.getAirPressure(scale_height * ureg.m)
            / meteorology.getAirPressure(0 * ureg.m)
        ).magnitude
        assert ratio == pytest.approx(1 / np.e, rel=1e-4)

    def test_pressure_falls_with_height(self, meteorology):
        values = [
            meteorology.getAirPressure(z * ureg.m).to(ureg.mmHg).magnitude
            for z in (0, 100, 1000, 5000)
        ]
        assert values == sorted(values, reverse=True)

    def test_a_height_in_kilometres_is_converted(self, meteorology):
        from_km = meteorology.getAirPressure(1 * ureg.km).to(ureg.mmHg).magnitude
        from_m = meteorology.getAirPressure(1000 * ureg.m).to(ureg.mmHg).magnitude
        assert from_km == pytest.approx(from_m)


@pytest.mark.unit
class TestAirTemperature:
    def test_ground_temperature_is_what_was_configured(self, meteorology):
        assert meteorology.getAirTemperature(0 * ureg.m).to(
            ureg.kelvin
        ).magnitude == pytest.approx(293.15)

    def test_the_lapse_rate_is_six_and_a_half_kelvin_per_kilometre(self, meteorology):
        ground = meteorology.getAirTemperature(0 * ureg.m).to(ureg.kelvin).magnitude
        aloft = meteorology.getAirTemperature(1000 * ureg.m).to(ureg.kelvin).magnitude
        assert ground - aloft == pytest.approx(6.5)

    def test_the_lapse_rate_is_linear(self, meteorology):
        first = meteorology.getAirTemperature(0 * ureg.m).to(ureg.kelvin).magnitude
        second = meteorology.getAirTemperature(500 * ureg.m).to(ureg.kelvin).magnitude
        third = meteorology.getAirTemperature(1000 * ureg.m).to(ureg.kelvin).magnitude
        assert (first - second) == pytest.approx(second - third)

    def test_temperature_falls_with_height(self, meteorology):
        values = [
            meteorology.getAirTemperature(z * ureg.m).to(ureg.kelvin).magnitude
            for z in (0, 500, 1000, 2000)
        ]
        assert values == sorted(values, reverse=True)


@pytest.mark.unit
class TestAirDensity:
    def test_matches_dry_air_at_twenty_celsius(self, meteorology):
        """Textbook value: 1.204 kg/m3 at 20 C and one atmosphere."""
        density = meteorology.getAirDensity(0 * ureg.m)
        assert density.to(ureg.kg / ureg.m**3).magnitude == pytest.approx(1.204, rel=1e-3)

    def test_the_result_carries_mks_units(self, meteorology):
        density = meteorology.getAirDensity(0 * ureg.m)
        assert density.units == ureg.kg / ureg.m**3

    def test_density_falls_with_height(self, meteorology):
        values = [
            meteorology.getAirDensity(z * ureg.m).to(ureg.kg / ureg.m**3).magnitude
            for z in (0, 500, 1000, 3000)
        ]
        assert values == sorted(values, reverse=True)

    def test_colder_ground_air_is_denser(self, factory):
        cold = factory.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s,
            inversion=1000 * ureg.m,
            temperature=ureg.Quantity(0, ureg.degC),
        )
        warm = factory.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s,
            inversion=1000 * ureg.m,
            temperature=ureg.Quantity(30, ureg.degC),
        )
        assert cold.getAirDensity(0 * ureg.m) > warm.getAirDensity(0 * ureg.m)


@pytest.mark.unit
class TestAirViscosity:
    def test_is_within_a_few_percent_of_the_textbook_value(self, meteorology):
        """Dry air at 20 C is about 1.81e-5 Pa s.

        The correlation used here gives 1.854e-5, about 2.4% high, so the
        tolerance is loose on purpose -- the test asserts the right order of
        magnitude and sign of the temperature dependence, not the correlation.
        """
        viscosity = meteorology.getAirDynamicViscosity(0 * ureg.m)
        assert viscosity.to(ureg.Pa * ureg.s).magnitude == pytest.approx(
            1.81e-5, rel=0.05
        )

    def test_viscosity_rises_with_temperature(self, factory):
        """Unlike liquids, gas viscosity increases with temperature."""
        cold = factory.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s,
            inversion=1000 * ureg.m,
            temperature=ureg.Quantity(0, ureg.degC),
        ).getAirDynamicViscosity(0 * ureg.m)
        warm = factory.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s,
            inversion=1000 * ureg.m,
            temperature=ureg.Quantity(30, ureg.degC),
        ).getAirDynamicViscosity(0 * ureg.m)
        assert warm > cold

    def test_viscosity_falls_with_height(self, meteorology):
        """Because it follows temperature, which falls with height."""
        ground = meteorology.getAirDynamicViscosity(0 * ureg.m)
        aloft = meteorology.getAirDynamicViscosity(2000 * ureg.m)
        assert aloft < ground


@pytest.mark.unit
class TestTurbulentKineticEnergy:
    def test_equals_the_documented_expression(self, meteorology):
        """TKE = (3.25 u*)^2, with u* = 0.3 m/s by default."""
        expected = (3.25 * 0.3) ** 2
        assert meteorology.getTKE(0 * ureg.m).magnitude == pytest.approx(expected)

    def test_the_height_argument_is_ignored(self, meteorology):
        """Documented as "currently implemented only neutral conditions".

        Pinned so that adding a height dependence is a deliberate change
        rather than an accident -- callers today may rely on a constant.
        """
        at_ground = meteorology.getTKE(0 * ureg.m)
        aloft = meteorology.getTKE(250 * ureg.m)
        assert at_ground == aloft

    def test_scales_with_the_square_of_ustar(self, factory):
        one = factory.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s,
            inversion=1000 * ureg.m,
            ustar=0.2 * ureg.m / ureg.s,
        ).getTKE(0 * ureg.m)
        two = factory.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s,
            inversion=1000 * ureg.m,
            ustar=0.4 * ureg.m / ureg.s,
        ).getTKE(0 * ureg.m)
        assert (two / one).magnitude == pytest.approx(4.0)


@pytest.mark.unit
class TestLogWindProfile:
    def test_reproduces_the_reference_wind_at_the_reference_height(self, meteorology):
        """The defining property of the calibration u*/kappa = u_ref/ln(z_ref/z0)."""
        at_reference = meteorology.getWindVelocity(meteorology.refHeight * ureg.m)
        assert at_reference.to(ureg.m / ureg.s).magnitude == pytest.approx(5.0)

    def test_the_wind_vanishes_at_the_roughness_length(self, meteorology):
        assert meteorology.getWindVelocity(meteorology.z0).to(
            ureg.m / ureg.s
        ).magnitude == pytest.approx(0.0)

    def test_the_wind_is_zero_below_the_roughness_length(self, meteorology):
        assert meteorology.getWindVelocity(0.01 * ureg.m).to(
            ureg.m / ureg.s
        ).magnitude == pytest.approx(0.0)

    def test_the_wind_increases_with_height(self, meteorology):
        values = [
            meteorology.getWindVelocity(z * ureg.m).to(ureg.m / ureg.s).magnitude
            for z in (1, 5, 10, 50, 100)
        ]
        assert values == sorted(values)

    def test_a_rougher_surface_gives_more_shear_above_the_reference(self, factory):
        r"""With u10 pinned, a larger z0 raises the wind aloft.

        u(z) = u10 * ln(z/z0) / ln(10/z0).  For z > 10 m that ratio increases
        with z0, because the normalising ln(10/z0) shrinks faster than
        ln(z/z0) does.  So rough terrain, holding the 10 m wind fixed, implies
        a steeper profile and a stronger wind above.  Checked with numbers:
        z0 = 1 m gives 10.0 m/s at 100 m, z0 = 1 cm gives 6.67 m/s.
        """
        smooth = factory.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s, inversion=1000 * ureg.m, z0=0.01 * ureg.m
        ).getWindVelocity(100 * ureg.m)
        rough = factory.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s, inversion=1000 * ureg.m, z0=1.0 * ureg.m
        ).getWindVelocity(100 * ureg.m)

        assert rough > smooth
        assert rough.to(ureg.m / ureg.s).magnitude == pytest.approx(10.0)
        assert smooth.to(ureg.m / ureg.s).magnitude == pytest.approx(6.6667, rel=1e-4)

    def test_the_profile_is_clamped_at_three_hundred_metres(self, meteorology):
        """Undocumented in the docstring, but real: height is clipped to [0, 300].

        Pinned because a caller asking for 1 km would otherwise silently get the
        300 m value.
        """
        at_300 = meteorology.getWindVelocity(300 * ureg.m)
        assert meteorology.getWindVelocity(400 * ureg.m) == at_300
        assert meteorology.getWindVelocity(1000 * ureg.m) == at_300


@pytest.mark.unit
class TestOtherWindProfiles:
    def test_uniform_wind_ignores_height(self, factory):
        uniform = factory.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s,
            inversion=1000 * ureg.m,
            verticalProfileType="uniformWind",
        )
        assert uniform.getWindVelocity(1 * ureg.m) == uniform.getWindVelocity(200 * ureg.m)

    def test_uniform_wind_returns_the_reference_value(self, factory):
        uniform = factory.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s,
            inversion=1000 * ureg.m,
            verticalProfileType="uniformWind",
        )
        assert uniform.getWindVelocity(50 * ureg.m).to(
            ureg.m / ureg.s
        ).magnitude == pytest.approx(5.0)

    def test_power_law_reproduces_the_reference_wind(self, factory):
        power = factory.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s,
            inversion=1000 * ureg.m,
            verticalProfileType="powerLaw",
        )
        assert power.getWindVelocity(10 * ureg.m).to(
            ureg.m / ureg.s
        ).magnitude == pytest.approx(5.0)

    def test_power_law_increases_with_height(self, factory):
        power = factory.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s,
            inversion=1000 * ureg.m,
            verticalProfileType="powerLaw",
        )
        values = [
            power.getWindVelocity(z * ureg.m).to(ureg.m / ureg.s).magnitude
            for z in (1, 10, 50, 100)
        ]
        assert values == sorted(values)

    @pytest.mark.parametrize(
        "stability, exponent",
        [("A", 0.07), ("B", 0.07), ("C", 0.1), ("D", 0.15), ("E", 0.35), ("F", 0.55)],
    )
    def test_hot_spot_uses_the_published_exponents(self, factory, stability, exponent):
        """u(z) = u10 (z/10)^p, with p rising from unstable to stable air."""
        meteorology = factory.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s, inversion=1000 * ureg.m, stability=stability
        )
        computed = meteorology.getWindVelocity_hotSpot(100 * ureg.m)
        expected = 5.0 * (100.0 / 10.0) ** exponent
        assert computed.to(ureg.m / ureg.s).magnitude == pytest.approx(expected)

    def test_hot_spot_returns_u10_at_ten_metres(self, meteorology):
        assert meteorology.getWindVelocity_hotSpot(10 * ureg.m).to(
            ureg.m / ureg.s
        ).magnitude == pytest.approx(5.0)


@pytest.mark.unit
class TestFactory:
    @pytest.mark.parametrize(
        "profileType, expected",
        [
            ("powerLaw", StandardMeteorolgyConstant_powerLaw),
            ("log", StandardMeteorolgyConstant_log),
            ("uniformWind", StandardMeteorolgyConstant_uniformWind),
        ],
    )
    def test_each_profile_name_maps_to_its_class(self, factory, profileType, expected):
        built = factory.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s,
            inversion=1000 * ureg.m,
            verticalProfileType=profileType,
        )
        assert type(built) is expected

    def test_log_is_the_default_profile(self, factory):
        built = factory.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s, inversion=1000 * ureg.m
        )
        assert isinstance(built, StandardMeteorolgyConstant_log)

    def test_an_unknown_profile_name_raises(self, factory):
        with pytest.raises(KeyError):
            factory.getMeteorologyFromU10(
                u10=5 * ureg.m / ureg.s,
                inversion=1000 * ureg.m,
                verticalProfileType="parabolic",
            )

    def test_building_from_an_arbitrary_reference_height(self, factory):
        built = factory.getMeteorologyFromURefHeight(
            u=8 * ureg.m / ureg.s, inversion=1000 * ureg.m, refHeight=50 * ureg.m
        )
        assert built.getWindVelocity(50 * ureg.m).to(
            ureg.m / ureg.s
        ).magnitude == pytest.approx(8.0)

    @pytest.mark.parametrize("stability", STABILITIES)
    def test_every_stability_class_is_accepted(self, factory, stability):
        built = factory.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s, inversion=1000 * ureg.m, stability=stability
        )
        assert built.stability == stability

    @pytest.mark.parametrize("stability", STABILITIES)
    def test_the_power_law_exponent_is_resolved_for_every_class(self, factory, stability):
        """From Irwin (1984), interpolated over z0; it must be a real number."""
        built = factory.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s,
            inversion=1000 * ureg.m,
            stability=stability,
            verticalProfileType="powerLaw",
        )
        assert 0 < float(built.wind_p) < 1

    def test_more_stable_air_has_a_steeper_profile(self, factory):
        """The Irwin exponent grows monotonically from A to F."""
        exponents = []
        for stability in STABILITIES:
            built = factory.getMeteorologyFromU10(
                u10=5 * ureg.m / ureg.s,
                inversion=1000 * ureg.m,
                stability=stability,
                verticalProfileType="powerLaw",
            )
            exponents.append(float(built.wind_p))
        assert exponents == sorted(exponents)


@pytest.mark.unit
class TestRawStringPrefixOnTheWrongQuotes:
    r"""Five docstrings in gaussian/ carry the r prefix on the CLOSING quotes.

        def getAirDensity(self, height):
            \"\"\"
                ... \rho_{air} ...
            r\"\"\"

    The r therefore becomes the last character of a non-raw docstring, so the
    LaTeX backslashes are invalid escape sequences.  Two consequences follow,
    and the second is the interesting one.
    """

    def test_the_module_is_importable_despite_the_escapes(self):
        """The skip message in test_no_invalid_escapes.py claims otherwise.

        It reports "pre-existing SyntaxError, file is not importable", which is
        false: the module imports fine.  What it actually hit was a SyntaxWarning
        that its own filter had escalated into a SyntaxError.
        """
        import hera.simulations.gaussian.Meteorology as module

        assert module.MeteorologyFactory is not None

    def test_an_escalated_escape_warning_surfaces_as_a_syntax_error(self):
        """The mechanism behind B26, shown on a one-line probe.

        With the escape warning escalated to an error, compile() raises
        SyntaxError. test_no_invalid_escapes.py catches SyntaxError and
        calls pytest.skip, so the failure it exists to detect can only ever
        become a skip.

        Python <=3.11 emits DeprecationWarning for an invalid escape
        sequence; >=3.12 emits SyntaxWarning -- both must be escalated for
        this probe to behave the same way on every supported interpreter.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("error", SyntaxWarning)
            warnings.simplefilter("error", DeprecationWarning)
            with pytest.raises(SyntaxError):
                compile('x = "\\c"\n', "probe.py", "exec")

    def test_compiling_the_module_emits_no_warning(self):
        import pathlib

        source = pathlib.Path(
            "hera/simulations/gaussian/Meteorology.py"
        ).read_text(encoding="utf-8")

        with warnings.catch_warnings():
            warnings.simplefilter("error", SyntaxWarning)
            compile(source, "Meteorology.py", "exec")
