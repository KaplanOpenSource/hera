"""Agent and PhysicalPropeties: the hazardous-agent descriptor and its
physical-property correlations.

Constructing a full ``Agent`` also builds its ``PhysicalPropeties`` and its
``Injury`` effects (via ``InjuryFactory``), so the fixture below wires a
minimal but complete descriptor -- one effect, one level, one calculator --
using the exact shape ``RiskToolkit.analysis.getRiskAreas`` builds inline.
"""
import pytest

from hera.riskassessment.agents.Agents import Agent, PhysicalPropeties
from hera.utils.unitHandler import ureg


def _agent_descriptor(**physical_properties):
    desc = {
        "name": "TestAgent",
        "effectParameters": {"tenbergeCoefficient": 2},
        "effects": {
            "RegularPopulation": {
                "type": "Threshold",
                "calculator": {"TenBerge": {}},
                "parameters": {
                    "type": "Threshold",
                    "levels": ["Severe"],
                    "parameters": {"Severe": {"threshold": "10 mg/m**3"}},
                },
            }
        },
    }
    if physical_properties:
        desc["physicalProperties"] = physical_properties
    return desc


@pytest.mark.unit
class TestAgentConstruction:
    def test_effect_names_lists_the_configured_effects(self):
        agent = Agent(_agent_descriptor())
        assert agent.effectNames == ["RegularPopulation"]

    def test_an_effect_is_reachable_by_dict_style_access(self):
        from hera.riskassessment.agents.effects.Injury import Injury

        agent = Agent(_agent_descriptor())
        assert isinstance(agent["RegularPopulation"], Injury)

    def test_the_name_comes_from_the_descriptor(self):
        assert Agent(_agent_descriptor()).name == "TestAgent"

    def test_full_description_is_the_original_descriptor(self):
        desc = _agent_descriptor()
        assert Agent(desc).fullDescription is desc

    def test_effect_properties_is_the_effect_parameters_dict(self):
        desc = _agent_descriptor()
        desc["effectParameters"] = {"tenbergeCoefficient": 3}
        assert Agent(desc).effectproperties == {"tenbergeCoefficient": 3}

    def test_physical_properties_default_to_an_empty_shell_without_the_key(self):
        agent = Agent(_agent_descriptor())
        assert agent.physicalproperties.molecularWeight is None
        assert agent.physicalproperties.spreadFactor is None


@pytest.mark.unit
class TestAgentToJSON:
    def test_toJSON_includes_the_name_and_physical_properties(self):
        agent = Agent(_agent_descriptor(molecularWeight="108*g/mol"))
        j = agent.toJSON()
        assert j["name"] == "TestAgent"
        assert j["physicalProperties"] == {"molecularWeight": "108*g/mol"}

    def test_toJSON_includes_one_entry_per_effect(self):
        agent = Agent(_agent_descriptor())
        j = agent.toJSON()
        assert set(j["effect"]) == {"RegularPopulation"}
        assert j["effect"]["RegularPopulation"]["calculator"]["type"] == "tenBerge"

    def test_str_is_the_pretty_printed_json(self):
        import json

        agent = Agent(_agent_descriptor())
        assert str(agent) == json.dumps(agent.toJSON(), indent=4)


@pytest.mark.unit
class TestPhysicalPropetiesConstruction:
    def test_without_a_physical_properties_key_nothing_is_set(self):
        props = PhysicalPropeties({"name": "X"})
        assert props.molecularWeight is None
        assert props.sorptionCoefficient is None
        assert props.spreadFactor is None

    def test_molecular_weight_is_converted_to_grams_per_mole(self):
        props = PhysicalPropeties({"physicalProperties": {"molecularWeight": "0.108 kg/mol"}})
        assert props.molecularWeight.m_as(ureg.g / ureg.mol) == pytest.approx(108.0)

    def test_molecular_weight_defaults_to_one_gram_per_mole(self):
        props = PhysicalPropeties({"physicalProperties": {}})
        assert props.molecularWeight.m_as(ureg.g / ureg.mol) == pytest.approx(1.0)

    def test_sorption_coefficient_is_converted_to_centimetres_per_second(self):
        props = PhysicalPropeties({"physicalProperties": {"sorptionCoefficient": "0.01 m/s"}})
        assert props.sorptionCoefficient.m_as(ureg.cm / ureg.s) == pytest.approx(1.0)

    def test_spread_factor_defaults_to_one(self):
        props = PhysicalPropeties({"physicalProperties": {}})
        assert props.spreadFactor == pytest.approx(1.0)

    def test_spread_factor_is_coerced_to_float(self):
        props = PhysicalPropeties({"physicalProperties": {"spreadFactor": "2.5"}})
        assert props.spreadFactor == pytest.approx(2.5)


@pytest.mark.unit
class TestPhysicalPropetiesGetters:
    @pytest.fixture()
    def props(self):
        return PhysicalPropeties({
            "physicalProperties": {
                "molecularWeight": "108*g/mol",
                "sorptionCoefficient": "1*cm/s",
                "spreadFactor": 2,
            }
        })

    def test_get_molecular_weight_matches_the_property(self, props):
        assert props.getMolecularWeight() == props.molecularWeight

    def test_get_spread_factor_matches_the_property(self, props):
        assert props.getSpreadFactor() == props.spreadFactor

    def test_get_sorption_coefficient_matches_the_property(self, props):
        assert props.getSorptionCoefficient() == props.sorptionCoefficient

    def test_molecular_volume_reads_the_raw_params_dict(self):
        props = PhysicalPropeties({"physicalProperties": {"molecularVolume": 42}})
        assert props.molecularVolume == 42

    def test_molecular_volume_without_the_key_raises(self, props):
        with pytest.raises(KeyError, match="molecularVolume"):
            props.molecularVolume

    def test_toJSON_returns_the_original_params_dict(self):
        original = {"molecularWeight": "108*g/mol"}
        props = PhysicalPropeties({"physicalProperties": original})
        assert props.toJSON() == original


@pytest.mark.unit
class TestPhysicalPropetiesCorrelations:
    """Regression-style checks against the documented formulas, not against
    current output -- see the module docstring's a/b/c/d constants."""

    @pytest.fixture()
    def props(self):
        return PhysicalPropeties({
            "physicalProperties": {
                "molecularWeight": "100*g/mol",
                "volatilityConstants": [8.0, 2000.0, 200.0, 0],
                "densityConstants": [1.5, 0.001, 20.0],
                "vaporPressure": {"A": 5.0, "B": 1000.0, "C": 50.0, "units": "bar"},
            }
        })

    def test_density_decreases_linearly_above_the_reference_temperature(self, props):
        a, b, c = 1.5, 0.001, 20.0
        assert props.getDensity(ureg.Quantity(20, ureg.degC)).m_as(ureg.g / ureg.cm**3) == pytest.approx(a)
        assert props.getDensity(ureg.Quantity(120, ureg.degC)).m_as(ureg.g / ureg.cm**3) == pytest.approx(
            a - b * 100
        )

    @pytest.mark.xfail(
        strict=True,
        reason="B66: getVolatility multiplies MW (a Quantity in g/mol) "
               "straight into the return expression and then multiplies by "
               "g/cm**3 again, without ever consuming the mol**-1 -- the "
               "result comes out in g**2/(mol*cm**3), not the documented "
               "g/cm**3 vapor saturation concentration. "
               "See the consolidated findings issue.",
    )
    def test_volatility_is_a_mass_concentration(self, props):
        assert props.getVolatility(ureg.Quantity(25, ureg.degC)).m_as(ureg.g / ureg.cm**3) > 0

    def test_volatility_currently_carries_an_extra_inverse_mole(self, props):
        """Characterisation of B66."""
        volatility = props.getVolatility(ureg.Quantity(25, ureg.degC))
        assert volatility.units == ureg.g**2 / ureg.mol / ureg.cm**3

    def test_vapor_pressure_uses_the_antoine_style_constants(self, props):
        A, B, C = 5.0, 1000.0, 50.0
        temperature_k = 300.0
        expected = 10 ** (A - B / (temperature_k - C))
        assert props.vaporPressure(temperature_k * ureg.kelvin) == pytest.approx(expected)

    def test_vapor_pressure_defaults_d_e_f_to_zero_when_absent(self, props):
        """The D/E/F terms are optional; omitting them must not raise."""
        props.vaporPressure(300.0 * ureg.kelvin)
