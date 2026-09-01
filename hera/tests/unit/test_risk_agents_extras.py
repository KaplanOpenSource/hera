"""riskassessment: tenbergeCoefficient property and CalculatorMaxConcentration."""
import pandas
import pytest

from hera.riskassessment.agents.Agents import Agent
from hera.riskassessment.agents.effects.Calculator import CalculatorMaxConcentration
from hera.utils import ureg


@pytest.mark.unit
class TestTenbergeCoefficient:
    @pytest.fixture()
    def agent(self):
        a = Agent.__new__(Agent)
        a._effectParameters = {}
        a._agentconfig = {"effects": {}}
        a._effects = {}
        return a

    def test_defaults_to_one(self, agent):
        assert agent.tenbergeCoefficient == 1

    def test_setter_stores_as_float_and_rebuilds_effects(self, agent):
        agent.tenbergeCoefficient = "2.5"
        assert agent.tenbergeCoefficient == 2.5
        assert agent._effectParameters["tenbergeCoefficient"] == 2.5


@pytest.mark.unit
class TestCalculatorMaxConcentration:
    def test_it_computes_the_max_rolling_average(self):
        calc = CalculatorMaxConcentration(sampling=2, breathingRate=10 * ureg.L / ureg.min)
        df = pandas.DataFrame({"C": [1.0, 2.0, 3.0, 4.0]})
        result = calc.calculate(
            df, field="C", breathingRate=10 * ureg.L / ureg.min, inUnits=1 * ureg.mg / ureg.m**3
        )
        assert result["C"] == pytest.approx(3.5)

    def test_to_json_reports_its_type(self):
        calc = CalculatorMaxConcentration(sampling=2)
        assert calc.toJSON()["type"] == "maxConcentrations"
