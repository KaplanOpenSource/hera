"""The two InjuryLevel subclasses batch 5 didn't reach: Threshold's
serialisation and the exponential dose-response level end to end."""
import numpy
import pytest

from hera.riskassessment.agents.effects.InjuryLevel import InjuryLevelExponential, InjuryLevelThreshold
from hera.utils.unitHandler import ureg


@pytest.mark.unit
class TestThresholdToJSON:
    def test_to_json_names_its_type_and_the_threshold(self):
        level = InjuryLevelThreshold("Severe", threshold="10 mg/m**3")
        j = level.toJSON()
        assert j["type"] == "threshold"
        assert j["name"] == "Severe"
        assert "10" in j["threshold"]


@pytest.mark.unit
class TestExponentialConstruction:
    def test_k_is_coerced_to_a_float(self):
        level = InjuryLevelExponential("Severe", k="0.5")
        assert level.k == pytest.approx(0.5)
        assert isinstance(level.k, float)

    def test_missing_k_raises(self):
        with pytest.raises(KeyError):
            InjuryLevelExponential("Severe")


@pytest.mark.unit
class TestExponentialGetPercent:
    def test_zero_toxic_load_gives_zero_percent(self):
        level = InjuryLevelExponential("Severe", k=0.1)
        assert level.getPercent(0.0) == pytest.approx(0.0)

    def test_the_formula_is_one_minus_exp_minus_k_times_load(self):
        level = InjuryLevelExponential("Severe", k=0.1)
        toxic_load = 5.0
        assert level.getPercent(toxic_load) == pytest.approx(1 - numpy.exp(-0.1 * toxic_load))

    def test_the_response_is_monotonically_increasing_with_load(self):
        level = InjuryLevelExponential("Severe", k=0.2)
        loads = numpy.array([0.0, 1.0, 5.0, 20.0])
        percents = level.getPercent(loads)
        assert numpy.all(numpy.diff(percents) > 0)

    def test_a_large_load_approaches_full_injury(self):
        level = InjuryLevelExponential("Severe", k=1.0)
        assert level.getPercent(50.0) == pytest.approx(1.0, abs=1e-10)

    def test_to_json_names_its_type_and_k(self):
        level = InjuryLevelExponential("Severe", k=0.5)
        j = level.toJSON()
        assert j["type"] == "exponential"
        assert j["k"] == "0.5"
