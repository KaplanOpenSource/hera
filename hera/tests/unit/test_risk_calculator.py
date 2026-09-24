r"""Toxic-load calculators.

Two dose-response laws, with definitions that make the expected numbers
computable by hand:

    Haber      D(T) = \int_0^T C dt
    ten Berge  D(T) = \int_0^T C^n dt      (Haber is the n = 1 case)

Both are scaled by the ratio of the exposed population's breathing rate to the
rate at which the injury level was established.  A constant concentration
therefore gives D = C x t exactly, which is what most of these assert.
"""
import numpy as np
import pandas as pd
import pytest
import xarray as xr

from hera.riskassessment.agents.effects.Calculator import (
    AbstractCalculator,
    CalculatorHaber,
    CalculatorMaxConcentration,
    CalculatorTenBerge,
)
from hera.utils.unitHandler import ureg

MG_PER_M3 = ureg.mg / ureg.m**3


def constant_field(concentration=10.0, steps=5, minutes=1.0, units=MG_PER_M3):
    """A uniform concentration over `steps` equal time steps."""
    dataset = xr.Dataset(
        {"C": ("datetime", np.full(steps, float(concentration)))},
        coords={"datetime": np.arange(steps)},
    )
    dataset.attrs["C"] = 1 * units
    dataset.attrs["dt"] = minutes * ureg.min
    return dataset


@pytest.mark.unit
class TestHaberLaw:
    def test_a_constant_exposure_accumulates_linearly(self):
        """C = 10 mg/m3 for 1-minute steps gives 10, 20, 30, ... mg.min/m3."""
        dose = CalculatorHaber().calculate(
            constant_field(), field="C", time="datetime"
        )
        assert np.asarray(dose).tolist() == pytest.approx([10.0, 20.0, 30.0, 40.0, 50.0])

    def test_the_dose_scales_with_concentration(self):
        single = CalculatorHaber().calculate(
            constant_field(concentration=1.0), field="C", time="datetime"
        )
        double = CalculatorHaber().calculate(
            constant_field(concentration=2.0), field="C", time="datetime"
        )
        assert np.asarray(double)[-1] == pytest.approx(2 * np.asarray(single)[-1])

    def test_the_dose_scales_with_the_time_step(self):
        short = CalculatorHaber().calculate(
            constant_field(minutes=1.0), field="C", time="datetime"
        )
        long = CalculatorHaber().calculate(
            constant_field(minutes=2.0), field="C", time="datetime"
        )
        assert np.asarray(long)[-1] == pytest.approx(2 * np.asarray(short)[-1])

    def test_a_zero_concentration_gives_no_dose(self):
        dose = CalculatorHaber().calculate(
            constant_field(concentration=0.0), field="C", time="datetime"
        )
        assert np.asarray(dose).tolist() == pytest.approx([0.0] * 5)

    def test_the_dose_never_decreases(self):
        """A cumulative integral of a non-negative field is monotonic."""
        values = np.asarray(
            CalculatorHaber().calculate(constant_field(), field="C", time="datetime")
        ).tolist()
        assert values == sorted(values)

    def test_a_breathing_rate_above_the_reference_raises_the_dose(self):
        """The ratio is exposed rate / reference rate, so double means double."""
        reference = CalculatorHaber().calculate(
            constant_field(), field="C", time="datetime"
        )
        doubled = CalculatorHaber().calculate(
            constant_field(),
            field="C",
            time="datetime",
            breathingRate=20 * ureg.L / ureg.min,
        )
        assert np.asarray(doubled)[-1] == pytest.approx(2 * np.asarray(reference)[-1])

    def test_the_concentration_units_are_converted(self):
        """1 g/m3 is 1000 mg/m3, so the dose must be a thousand times larger."""
        milligrams = CalculatorHaber().calculate(
            constant_field(units=MG_PER_M3), field="C", time="datetime"
        )
        grams = CalculatorHaber().calculate(
            constant_field(units=ureg.g / ureg.m**3), field="C", time="datetime"
        )
        assert np.asarray(grams)[-1] == pytest.approx(
            1000 * np.asarray(milligrams)[-1]
        )


@pytest.mark.unit
class TestTenBergeLaw:
    def test_the_n_equals_one_case_reduces_to_haber(self):
        """The defining relationship between the two laws."""
        field = constant_field()
        haber = CalculatorHaber().calculate(field, field="C", time="datetime")
        tenberge = CalculatorTenBerge(tenbergeCoefficient=1.0).calculate(
            field, field="C", time="datetime"
        )
        assert np.asarray(tenberge).tolist() == pytest.approx(
            np.asarray(haber).tolist()
        )

    def test_the_exponent_applies_to_the_concentration(self):
        """C = 10, n = 2 gives 100 per minute, so 100, 200, 300, ..."""
        dose = CalculatorTenBerge(tenbergeCoefficient=2.0).calculate(
            constant_field(), field="C", time="datetime"
        )
        assert np.asarray(dose).tolist() == pytest.approx(
            [100.0, 200.0, 300.0, 400.0, 500.0]
        )

    @pytest.mark.parametrize("exponent", [0.5, 1.0, 1.5, 2.0, 3.0])
    def test_the_first_step_is_the_concentration_raised_to_n(self, exponent):
        dose = CalculatorTenBerge(tenbergeCoefficient=exponent).calculate(
            constant_field(concentration=4.0), field="C", time="datetime"
        )
        assert np.asarray(dose)[0] == pytest.approx(4.0**exponent)

    def test_a_larger_exponent_punishes_high_concentrations_harder(self):
        """The reason ten Berge exists: brief high peaks matter more than n=1."""
        low = CalculatorTenBerge(tenbergeCoefficient=1.0).calculate(
            constant_field(concentration=10.0), field="C", time="datetime"
        )
        high = CalculatorTenBerge(tenbergeCoefficient=2.0).calculate(
            constant_field(concentration=10.0), field="C", time="datetime"
        )
        assert np.asarray(high)[-1] > np.asarray(low)[-1]

    def test_the_coefficient_is_recorded_on_the_instance(self):
        assert CalculatorTenBerge(tenbergeCoefficient=2.5).n == 2.5


@pytest.mark.unit
class TestSerialisation:
    def test_the_reference_breathing_rate_defaults_to_ten_litres_a_minute(self):
        assert CalculatorHaber().injuryBreathingRate.to(
            ureg.L / ureg.min
        ).magnitude == pytest.approx(10.0)

    def test_haber_names_its_type(self):
        assert CalculatorHaber().toJSON()["type"] == "haber"

    def test_the_breathing_rate_survives_serialisation(self):
        calculator = CalculatorHaber(breathingRate=15 * ureg.L / ureg.min)
        assert "15" in calculator.toJSON()["breathingRate"]

    def test_str_is_valid_json(self):
        import json

        assert isinstance(json.loads(str(CalculatorHaber())), dict)


@pytest.mark.unit
class TestInputValidation:
    def test_an_unsupported_container_is_rejected(self):
        with pytest.raises(ValueError, match="not a pandas.DataFrame or xarray.Dataset"):
            CalculatorHaber().calculate(
                [1, 2, 3], field="C", time="datetime", inUnits=1 * MG_PER_M3
            )

    def test_the_abstract_calculator_has_no_calculate(self):
        assert not hasattr(AbstractCalculator, "calculate")

    def test_max_concentration_requires_a_sampling_argument(self):
        """Unlike its siblings, its constructor takes a mandatory `sampling`."""
        with pytest.raises(TypeError):
            CalculatorMaxConcentration()


@pytest.mark.unit
class TestPandasInput:
    """The pandas branch is documented in detail and cannot produce a number.

    Its docstring even promises "For pandas, we do not assume that the time is
    equispaced", so it is meant to be a supported path.
    """

    @staticmethod
    def _frame(concentration=10.0, steps=5):
        index = pd.date_range("2020-01-01", periods=steps, freq="1min")
        frame = pd.DataFrame({"C": [concentration] * steps}, index=index)
        frame.index.name = "datetime"
        return frame

    @pytest.mark.xfail(
        strict=True,
        reason="B36: the default-units branch tests hasattr(field, 'attrs'), which "
               "modern pandas DataFrames satisfy, so it takes the xarray path and "
               "evaluates df.attrs[None] -> KeyError: None. The documented pandas "
               "default of mg/m3 is unreachable. "
               "See the consolidated findings issue.",
    )
    def test_the_documented_default_units_apply(self):
        dose = CalculatorHaber().calculate(self._frame(), time="datetime")
        assert not np.isnan(np.asarray(dose["C"])).any()

    @pytest.mark.xfail(
        strict=True,
        reason="B37: even with inUnits supplied, the result is all NaN. "
               "concentrationField[:-1] keeps the DatetimeIndex while dt_min[1:] "
               "carries the RangeIndex produced by reset_index(), so the "
               "multiplication aligns on disjoint indices. "
               "See the consolidated findings issue.",
    )
    def test_a_constant_exposure_accumulates_linearly(self):
        dose = CalculatorHaber().calculate(
            self._frame(), time="datetime", inUnits=1 * MG_PER_M3
        )
        assert np.asarray(dose["C"]).tolist() == pytest.approx(
            [10.0, 20.0, 30.0, 40.0, 50.0]
        )

    def test_the_two_indices_really_do_not_overlap(self):
        """Characterisation of B37's mechanism, so the cause is unambiguous."""
        frame = self._frame()
        elapsed = frame.reset_index()["datetime"].diff().apply(lambda x: x.seconds) / 60.0

        assert isinstance(frame[:-1].index, pd.DatetimeIndex)
        assert not isinstance(elapsed[1:].index, pd.DatetimeIndex)
        assert np.isnan((frame[:-1]["C"] * elapsed[1:]).to_numpy()).all()
