"""Dose-response injury levels.

The log-normal level has a defining identity: at the toxic load TL_50 exactly
half the population is affected, and getToxicLoad is its inverse.  Those two
facts pin the whole probit machinery, so most of the assertions here are about
them and about the shape properties any CDF must have.
"""
import numpy as np
import pytest

from hera.riskassessment.agents.effects.InjuryLevel import (
    InjuryLevel,
    InjuryLevelLognormal10DoseResponse,
    InjuryLevelThreshold,
)
from hera.utils.unitHandler import ureg


@pytest.fixture()
def severe():
    """Half the population affected at a toxic load of 100, probit slope 0.5."""
    return InjuryLevelLognormal10DoseResponse("Severe", TL_50=100.0, sigma=0.5)


@pytest.mark.unit
class TestLognormalDoseResponse:
    def test_half_the_population_is_affected_at_tl50(self, severe):
        """The definition of TL_50, and the anchor for everything else."""
        assert severe.getPercent(100.0) == pytest.approx(0.5)

    def test_the_inverse_returns_tl50_at_one_half(self, severe):
        assert severe.getToxicLoad(0.5) == pytest.approx(100.0)

    @pytest.mark.parametrize("fraction", [0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99])
    def test_the_two_directions_are_inverse(self, severe, fraction):
        assert severe.getPercent(severe.getToxicLoad(fraction)) == pytest.approx(
            fraction
        )

    def test_the_response_is_monotonic(self, severe):
        values = [float(severe.getPercent(load)) for load in (1, 10, 50, 100, 200, 1000)]
        assert values == sorted(values)

    def test_the_response_is_a_probability(self, severe):
        for load in (1e-9, 1.0, 100.0, 1e9):
            assert 0.0 <= float(severe.getPercent(load)) <= 1.0

    def test_a_decade_either_side_is_two_sigma(self, severe):
        """With sigma = 0.5 in log10, TL/10 and 10 TL are +-2 sigma.

        The normal CDF at -2 is 0.02275, so those two points pin both the base
        of the logarithm and the meaning of sigma.
        """
        assert float(severe.getPercent(10.0)) == pytest.approx(0.02275, abs=1e-4)
        assert float(severe.getPercent(1000.0)) == pytest.approx(0.97725, abs=1e-4)

    def test_a_smaller_sigma_makes_a_sharper_transition(self):
        sharp = InjuryLevelLognormal10DoseResponse("A", TL_50=100.0, sigma=0.1)
        broad = InjuryLevelLognormal10DoseResponse("B", TL_50=100.0, sigma=1.0)

        assert float(sharp.getPercent(200.0)) > float(broad.getPercent(200.0))
        assert float(sharp.getPercent(50.0)) < float(broad.getPercent(50.0))

    def test_both_still_cross_at_tl50(self):
        """Sigma changes the slope, never the midpoint."""
        for sigma in (0.1, 0.5, 1.0, 2.0):
            level = InjuryLevelLognormal10DoseResponse("L", TL_50=100.0, sigma=sigma)
            assert level.getPercent(100.0) == pytest.approx(0.5)

    def test_a_larger_tl50_shifts_the_curve_right(self):
        weak = InjuryLevelLognormal10DoseResponse("A", TL_50=100.0, sigma=0.5)
        strong = InjuryLevelLognormal10DoseResponse("B", TL_50=1000.0, sigma=0.5)
        assert float(strong.getPercent(100.0)) < float(weak.getPercent(100.0))

    def test_an_array_of_loads_is_accepted(self, severe):
        result = severe.getPercent(np.array([10.0, 100.0, 1000.0]))
        assert len(result) == 3
        assert result[1] == pytest.approx(0.5)


@pytest.mark.unit
class TestLognormalConstruction:
    def test_tl50_is_required(self):
        with pytest.raises(ValueError, match="Must supply TL_50"):
            InjuryLevelLognormal10DoseResponse("x", sigma=0.5)

    def test_sigma_is_required(self):
        with pytest.raises(ValueError, match="Must supply sigma"):
            InjuryLevelLognormal10DoseResponse("x", TL_50=100.0)

    def test_a_plain_number_is_given_the_default_dosage_units(self):
        level = InjuryLevelLognormal10DoseResponse("x", TL_50=100.0, sigma=0.5)
        assert level.TL_50.check("[mass]/[length]**3")

    def test_the_name_is_kept(self, severe):
        assert severe.name == "Severe"

    def test_it_is_an_injury_level(self, severe):
        assert isinstance(severe, InjuryLevel)

    def test_it_serialises_its_parameters(self, severe):
        payload = severe.toJSON()
        assert "TL_50" in str(payload) or "sigma" in str(payload)


@pytest.mark.unit
class TestThresholdLevel:
    def test_a_unit_bearing_string_is_accepted(self):
        level = InjuryLevelThreshold("T", threshold="50 mg/m**3")
        assert level.threshold.to(ureg.mg / ureg.m**3).magnitude == pytest.approx(50.0)

    def test_below_the_threshold_nobody_is_affected(self):
        """Passed as a hera Quantity, which is the only form that works."""
        level = InjuryLevelThreshold("T", threshold="50 mg/m**3")
        assert float(level.getPercent(10 * ureg.mg / ureg.m**3)) == 0

    def test_above_the_threshold_everybody_is(self):
        level = InjuryLevelThreshold("T", threshold="50 mg/m**3")
        assert float(level.getPercent(100 * ureg.mg / ureg.m**3)) == 1

    def test_the_step_is_at_the_threshold(self):
        level = InjuryLevelThreshold("T", threshold="50 mg/m**3")
        assert float(level.getPercent(49.9 * ureg.mg / ureg.m**3)) == 0
        assert float(level.getPercent(50.1 * ureg.mg / ureg.m**3)) == 1

    @pytest.mark.xfail(
        strict=True,
        reason="B39 (a downstream victim of B12): getPercent evaluates "
               "tounit(x, self.units) > threshold. For a plain number tounit "
               "builds a Quantity in pint's DEFAULT registry while the threshold "
               "lives in hera's, so the comparison raises "
               "'Cannot operate with Quantity and Quantity of different "
               "registries'. A numeric toxic load -- the ordinary calling form -- "
               "cannot be evaluated at all. See the consolidated findings issue.",
    )
    @pytest.mark.parametrize("load, expected", [(10.0, 0), (100.0, 1)])
    def test_a_plain_numeric_toxic_load_is_evaluated(self, load, expected):
        level = InjuryLevelThreshold("T", threshold="50 mg/m**3")
        assert float(np.atleast_1d(level.getPercent(load))[0]) == expected

    def test_the_threshold_is_required(self):
        with pytest.raises(ValueError, match="Cannot find the threshold"):
            InjuryLevelThreshold("T")

    def test_a_units_argument_changes_the_stored_units(self):
        level = InjuryLevelThreshold("T", units=ureg.g / ureg.m**3, threshold="50 mg/m**3")
        assert level.threshold.to(ureg.g / ureg.m**3).magnitude == pytest.approx(0.05)

    def test_a_bare_numeric_string_is_rejected_for_being_dimensionless(self):
        """At least this failure names the real problem."""
        from pint.errors import DimensionalityError

        with pytest.raises(DimensionalityError):
            InjuryLevelThreshold("T", threshold="50")

    @pytest.mark.xfail(
        strict=True,
        reason="B38: the threshold is passed straight to ureg(), which parses "
               "strings only, so a plain number fails with pint's internal "
               "AttributeError: 'float' object has no attribute 'replace'. Its "
               "sibling InjuryLevelLognormal10DoseResponse accepts TL_50=100.0 "
               "through tounit(), and the docstring here says only 'Must include "
               "threshold'. See the consolidated findings issue.",
    )
    def test_a_plain_number_is_accepted_like_its_sibling_takes_tl50(self):
        level = InjuryLevelThreshold("T", threshold=50.0)
        assert level.threshold.to(ureg.mg / ureg.m**3).magnitude == pytest.approx(50.0)

    def test_the_sibling_does_accept_a_plain_number(self):
        """Establishes the inconsistency B38 is about."""
        sibling = InjuryLevelLognormal10DoseResponse("x", TL_50=50.0, sigma=0.5)
        assert sibling.TL_50.to(ureg.mg / ureg.m**3).magnitude == pytest.approx(50.0)
