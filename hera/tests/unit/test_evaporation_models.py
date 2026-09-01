"""Evaporation models.

This module was unimportable until the relative-import fix (B31), so nothing
here had ever been exercised.  The correlations it implements have published
forms, and the assertions follow them:

    FSG diffusion    D = 1e-7 T^1.75 sqrt(1/Ma + 1/Mb) / (Va^1/3 + Vb^1/3)^2
    EPA diffusion    D = 4.09e-9 T^1.9 sqrt(1/Ma + 1/Mb) / Mb^1/3
    air viscosity    mu = 1.8205e-2 sqrt(T/293)

The named-model dispatch (`molecularDiffusion` -> `molecularDiffusion_<name>`)
is the same pydoc-style convention used elsewhere in hera.
"""
import numpy as np
import pytest

from hera.simulations.evaporation.models import evaporationModels
from hera.utils.unitHandler import ureg

AGENT = {
    "name": "test_agent",
    "effectParameters": {"tenbergeCoefficient": 1},
    # Note the capital P: PhysicalPropeties reads configJSON["physicalProperties"],
    # and molecularWeight goes through ureg(), so it must be a unit-bearing string.
    "physicalProperties": {
        "molecularWeight": "100 g/mol",
        "molecularVolume": 80.0,
    },
    "effects": {
        "inhalation": {
            "type": "Lognormal10",
            "calculator": {"TenBerge": {}},
            "units": "mg/m**3",
            "parameters": {
                "type": "Lognormal10DoseResponse",
                "levels": ["Severe"],
                "parameters": {"Severe": {"TL_50": 400, "sigma": 0.5}},
            },
        }
    },
}


@pytest.fixture()
def models():
    """Built around the broken constructor (B50).

    evaporationModels.__init__ calls RiskToolkit.getAgent(agent) on the CLASS,
    but getAgent is an instance method (self, nameOrDesc, version=None), so the
    agent lands in `self` and nameOrDesc is missing.  The class cannot be
    constructed at all.  Every correlation in it is nonetheless a pure function
    of the fields, so the object is assembled directly here and the constructor
    is pinned separately in TestConstruction.
    """
    from hera.riskassessment.agents.Agents import Agent

    built = evaporationModels.__new__(evaporationModels)
    built._agent = Agent(AGENT)
    built._Mair = 28.967
    built._Vair = 20.1
    built._Vagent = 80.0
    built._molecularDiffusionModel = "FSG"
    built._dinamicViscocityModel = "powerLaw"
    built._evaporationModel = "US"
    return built


@pytest.mark.unit
class TestConstruction:
    def test_the_documented_default_models_are_selected(self, models):
        assert models.evaporationModel == "US"
        assert models.dinamicViscocityModel == "powerLaw"
        assert models.molecularDiffusionModel == "FSG"

    def test_air_properties_are_the_standard_ones(self, models):
        """Air: 28.967 g/mol, FSG molar volume 20.1 cm3/mol."""
        assert models.Mair == pytest.approx(28.967)
        assert models.Vair == pytest.approx(20.1)

    def test_the_model_names_are_settable(self, models):
        models.molecularDiffusionModel = "EPA"
        assert models.molecularDiffusionModel == "EPA"

    def test_getagent_is_an_instance_method(self):
        """Characterisation of B50's mechanism."""
        import inspect

        from hera.riskassessment import RiskToolkit

        parameters = list(inspect.signature(RiskToolkit.getAgent).parameters)
        assert parameters[0] == "self"
        assert "nameOrDesc" in parameters

    @pytest.mark.xfail(
        strict=True,
        reason="B50: __init__ calls RiskToolkit.getAgent(agent) on the class, but "
               "getAgent is an instance method (self, nameOrDesc, version=None). "
               "The agent is bound to self and nameOrDesc is missing, so every "
               "construction raises TypeError -- the class cannot be instantiated "
               "at all. Undetected until now because the module was unimportable "
               "(B31). See the consolidated findings issue.",
    )
    def test_the_class_can_be_constructed(self):
        assert evaporationModels(AGENT) is not None


@pytest.mark.unit
class TestMolecularDiffusion:
    def test_fsg_matches_its_published_form(self, models):
        """D = 1e-7 T^1.75 sqrt(1/Ma + 1/Mb) / (Va^1/3 + Vb^1/3)^2."""
        temperature = 300.0
        expected = (
            1e-7
            * temperature**1.75
            * np.sqrt(1 / models.Mair + 1 / models.Magent)
            / (np.cbrt(models.Vair) + np.cbrt(models.Vagent)) ** 2
        )
        assert models.molecularDiffusion_FSG(temperature) == pytest.approx(expected)

    def test_epa_matches_its_published_form(self, models):
        """D = 4.09e-9 T^1.9 sqrt(1/Ma + 1/Mb) / Mb^1/3."""
        temperature = 300.0
        expected = (
            4.09e-9
            * temperature**1.9
            * np.sqrt(1 / models.Mair + 1 / models.Magent)
            / np.cbrt(models.Magent)
        )
        assert models.molecularDiffusion_EPA(temperature) == pytest.approx(expected)

    def test_the_named_model_is_dispatched(self, models):
        """molecularDiffusion defers to molecularDiffusion_<name>."""
        assert models.molecularDiffusion(300.0 * ureg.K) == pytest.approx(
            models.molecularDiffusion_FSG(300.0)
        )

    def test_switching_the_model_changes_the_dispatch(self, models):
        models.molecularDiffusionModel = "EPA"
        assert models.molecularDiffusion(300.0 * ureg.K) == pytest.approx(
            models.molecularDiffusion_EPA(300.0)
        )

    def test_an_unknown_model_name_raises(self, models):
        models.molecularDiffusionModel = "NoSuchMethod"
        with pytest.raises(AttributeError):
            models.molecularDiffusion(300.0 * ureg.K)

    def test_a_temperature_in_celsius_is_converted(self, models):
        """CLAUDE.md requires pint; the conversion must reach the correlation."""
        fromCelsius = models.molecularDiffusion(ureg.Quantity(26.85, ureg.degC))
        fromKelvin = models.molecularDiffusion(300.0 * ureg.K)
        assert fromCelsius == pytest.approx(fromKelvin, rel=1e-3)

    @pytest.mark.parametrize("temperature", [250.0, 300.0, 350.0, 400.0])
    def test_diffusion_grows_with_temperature(self, models, temperature):
        """Both correlations are positive powers of T."""
        assert models.molecularDiffusion_FSG(temperature) > 0
        assert models.molecularDiffusion_FSG(temperature + 50) > (
            models.molecularDiffusion_FSG(temperature)
        )

    def test_epa_has_the_steeper_temperature_dependence(self, models):
        """T^1.9 against T^1.75, so the EPA form must rise faster."""
        ratioEPA = models.molecularDiffusion_EPA(400.0) / models.molecularDiffusion_EPA(
            200.0
        )
        ratioFSG = models.molecularDiffusion_FSG(400.0) / models.molecularDiffusion_FSG(
            200.0
        )
        assert ratioEPA > ratioFSG

    def test_the_exponent_is_recoverable_from_two_points(self, models):
        """log-log slope must come out at 1.75 for FSG."""
        low, high = 200.0, 400.0
        slope = np.log(
            models.molecularDiffusion_FSG(high) / models.molecularDiffusion_FSG(low)
        ) / np.log(high / low)
        assert slope == pytest.approx(1.75, rel=1e-9)


@pytest.mark.unit
class TestAirViscosity:
    def test_it_matches_its_published_form(self, models):
        """mu = 1.8205e-2 sqrt(T/293)."""
        assert models.dynamicViscocityAir_powerLaw(300.0) == pytest.approx(
            1.8205e-2 * np.sqrt(300.0 / 293.0)
        )

    def test_the_reference_temperature_gives_the_reference_value(self, models):
        """At 293 K the square root is 1, so mu is exactly the coefficient."""
        assert models.dynamicViscocityAir_powerLaw(293.0) == pytest.approx(1.8205e-2)

    def test_viscosity_grows_with_temperature(self, models):
        """As a gas must, unlike a liquid."""
        assert models.dynamicViscocityAir_powerLaw(350.0) > (
            models.dynamicViscocityAir_powerLaw(250.0)
        )

    def test_the_named_model_is_dispatched(self, models):
        assert models.dynamicViscocityAir(300.0 * ureg.K) == pytest.approx(
            models.dynamicViscocityAir_powerLaw(300.0)
        )

    def test_the_square_root_law_is_recoverable(self, models):
        """Doubling T must multiply mu by sqrt(2)."""
        ratio = models.dynamicViscocityAir_powerLaw(
            600.0
        ) / models.dynamicViscocityAir_powerLaw(300.0)
        assert ratio == pytest.approx(np.sqrt(2.0))


@pytest.mark.unit
class TestDimensionlessGroups:
    def test_reynolds_grows_with_velocity(self, models):
        slow = models.Reynolds(diameter=1.0, velocity=1.0, temperature=300.0)
        fast = models.Reynolds(diameter=1.0, velocity=10.0, temperature=300.0)
        assert fast == pytest.approx(10 * slow)

    def test_reynolds_grows_with_the_length_scale(self, models):
        small = models.Reynolds(diameter=1.0, velocity=1.0, temperature=300.0)
        large = models.Reynolds(diameter=5.0, velocity=1.0, temperature=300.0)
        assert large == pytest.approx(5 * small)

    def test_reynolds_is_positive(self, models):
        assert models.Reynolds(diameter=1.0, velocity=1.0, temperature=300.0) > 0

    def test_schmidt_is_positive(self, models):
        assert models.Schmidt(300.0) > 0

    def test_schmidt_is_the_ratio_of_the_two_diffusivities(self, models):
        """Sc = nu / D, with nu = mu / rho."""
        temperature = 300.0
        density = models.Mair / (temperature * 8.205e-5)
        expected = models.dynamicViscocityAir(temperature * ureg.K) / (
            density * models.molecularDiffusion(temperature * ureg.K)
        )
        assert models.Schmidt(temperature) == pytest.approx(expected)


@pytest.mark.unit
class TestEvaporativeFlux:
    """B51: flux cannot run at all.

    flux_US strips the units off its temperature with
    `temperature = tonumber(temperature, ureg.K)` and then hands the bare float
    to agent.physicalproperties.vaporPressure, which starts with
    `unumToPint(temperature).m_as(ureg.kelvin)`.  A plain number is
    dimensionless, so that raises DimensionalityError for every input.
    Undetected until now because the module was unimportable (B31) and its
    constructor is broken (B50).
    """

    ARGS = dict(
        diameter=1 * ureg.m, velocity=2 * ureg.m / ureg.s, temperature=300 * ureg.K
    )

    def test_the_temperature_is_stripped_before_being_passed_on(self):
        """Characterisation of the mechanism, read off the source."""
        import inspect

        source = inspect.getsource(evaporationModels.flux_US)
        assert "temperature = tonumber(temperature, ureg.K)" in source
        assert "vaporPressure(temperature)" in source

    def test_vapour_pressure_requires_units(self):
        from hera.riskassessment.agents.Agents import PhysicalPropeties
        import inspect

        source = inspect.getsource(PhysicalPropeties.vaporPressure)
        assert "unumToPint(temperature).m_as(ureg.kelvin)" in source

    @pytest.mark.xfail(
        strict=True,
        reason="B51: flux_US passes a unit-stripped temperature to vaporPressure, "
               "which requires units, so every call raises DimensionalityError. "
               "See the consolidated findings issue.",
    )
    def test_it_is_positive(self, models):
        assert models.flux(**self.ARGS).magnitude > 0

    @pytest.mark.xfail(
        strict=True,
        reason="B51: same cause. See the consolidated findings issue.",
    )
    def test_it_carries_flux_units(self, models):
        assert models.flux(**self.ARGS).check("[mass]/([length]**2*[time])")

    @pytest.mark.xfail(
        strict=True,
        reason="B51: same cause. See the consolidated findings issue.",
    )
    def test_a_faster_wind_evaporates_more(self, models):
        slow = models.flux(diameter=1 * ureg.m, velocity=1 * ureg.m / ureg.s,
                           temperature=300 * ureg.K)
        fast = models.flux(diameter=1 * ureg.m, velocity=5 * ureg.m / ureg.s,
                           temperature=300 * ureg.K)
        assert fast > slow

    @pytest.mark.xfail(
        strict=True,
        reason="B51: same cause. See the consolidated findings issue.",
    )
    def test_the_output_units_are_configurable(self, models):
        assert models.flux(units=ureg.mg / (ureg.m**2 * ureg.s), **self.ARGS) is not None

    def test_the_named_model_is_dispatched(self, models):
        """The dispatch itself is fine -- an unknown name raises AttributeError."""
        models.evaporationModel = "NoSuchModel"
        with pytest.raises(AttributeError):
            models.flux(**self.ARGS)
