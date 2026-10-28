"""Injury and the calculators, exercised on pandas.DataFrame input.

Batch 5 exercised the calculators against xarray only; ``Injury`` also
supports a plain pandas.DataFrame ("each point is a column, time is the
index"), which turns out to have three independent problems:

* B62: ``Injury.calculateToxicLoads`` calls
  ``self.calculator.calculate(concentrationField, field, breathingRate=...,
  time=...)`` -- passing ``field`` positionally. That matches
  ``CalculatorTenBerge``/``CalculatorMaxConcentration``, whose second
  positional parameter is ``field``, but not ``CalculatorHaber``, whose
  signature is ``(self, concentrationField, breathingRate=..., time=...,
  field=None, ...)``. Any ``Injury`` configured with the Haber calculator
  raises ``TypeError`` the moment it is used.
* B63: all three calculators decide whether to read units from
  ``concentrationField.attrs[field]`` by checking
  ``hasattr(concentrationField, "attrs")`` -- meant to detect an xarray
  object. ``pandas.DataFrame`` has carried an (empty) ``.attrs`` dict since
  pandas 1.0, so the same branch runs for pandas input too, and
  ``df.attrs[field]`` raises ``KeyError``. Calling the calculator directly
  works around it with an explicit ``inUnits=``, but
  ``Injury.calculateToxicLoads``/``calculatePointWiseFractionInjured`` have
  no ``inUnits`` parameter to forward one through -- so through the
  ``Injury`` API, every pandas call hits this unconditionally, regardless
  of which calculator is configured.
* B64: past both of those, ``CalculatorHaber.calculate``'s pandas branch
  computes ``concentrationField[:-1].fillna(0) * dt_min[1:]`` -- a
  DataFrame multiplied by a Series without ``axis=0``, so pandas aligns the
  Series' index against the DataFrame's *columns* rather than its rows.
  With no matching labels the result is an all-NaN frame with the Series'
  index values as bogus extra columns. ``CalculatorTenBerge``'s equivalent
  line avoids exactly this by converting the same ``dt_min`` to
  ``.values.reshape(...)`` before multiplying.
"""
import pandas
import pytest

from hera.riskassessment.agents.Agents import Agent
from hera.riskassessment.agents.effects.Calculator import CalculatorHaber, CalculatorTenBerge
from hera.utils.unitHandler import ureg


def _series_dataframe(values, freq="1min"):
    times = pandas.date_range("2020-01-01", periods=len(values), freq=freq)
    return pandas.DataFrame({"P1": values}, index=pandas.Index(times, name="datetime"))


def _agent_with_calculator(calculator_name, **calculator_params):
    desc = {
        "name": "TestAgent",
        "effectParameters": {},
        "effects": {
            "RegularPopulation": {
                "type": "Threshold",
                "calculator": {calculator_name: calculator_params},
                "parameters": {
                    "type": "Threshold",
                    "levels": ["Severe"],
                    "parameters": {"Severe": {"threshold": "10 mg/m**3"}},
                },
            }
        },
    }
    return Agent(desc)["RegularPopulation"]


@pytest.mark.unit
class TestCalculatorHaberOnPandas:
    @pytest.mark.xfail(
        strict=True,
        reason="B64: the pandas branch multiplies a DataFrame by a Series "
               "without axis=0, aligning on columns instead of rows, so the "
               "cumulative dose comes out all NaN with bogus extra columns. "
               "See the consolidated findings issue.",
    )
    def test_a_constant_exposure_accumulates_linearly(self):
        df = _series_dataframe([10.0, 10.0, 10.0, 10.0, 10.0])
        dose = CalculatorHaber().calculate(
            df, time="datetime", field=None, inUnits=1 * ureg.mg / ureg.m**3
        )
        assert dose["P1"].iloc[-1] == pytest.approx(30.0)

    def test_the_result_is_currently_all_nan_with_bogus_extra_columns(self):
        """Characterisation of B64."""
        df = _series_dataframe([10.0, 10.0, 10.0, 10.0, 10.0])
        dose = CalculatorHaber().calculate(
            df, time="datetime", field=None, inUnits=1 * ureg.mg / ureg.m**3
        )
        assert set(dose.columns) == {"P1", 1, 2, 3, 4}
        assert dose.isna().all().all()

    def test_default_inunits_raises_because_pandas_now_has_attrs_too(self):
        """B63: hasattr(df, 'attrs') is true for pandas.DataFrame since
        pandas 1.0, so the xarray-only branch runs and df.attrs[field]
        raises."""
        df = _series_dataframe([10.0, 10.0])
        with pytest.raises(KeyError):
            CalculatorHaber().calculate(df, time="datetime", field=None)


@pytest.mark.unit
class TestCalculatorTenBergeOnPandas:
    """TenBerge sidesteps B64 (it reshapes dt_min to a bare array first),
    so it is the calculator to reach for on pandas input today."""

    def test_a_constant_exposure_with_n_one_accumulates_linearly(self):
        df = _series_dataframe([10.0, 10.0, 10.0, 10.0, 10.0])
        dose = CalculatorTenBerge(tenbergeCoefficient=1).calculate(
            df, field=None, time="datetime", inUnits=1 * ureg.mg / ureg.m**3
        )
        assert dose["P1"].iloc[-1] == pytest.approx(40.0)

    def test_default_inunits_raises_because_pandas_now_has_attrs_too(self):
        """B63 again, on TenBerge."""
        df = _series_dataframe([10.0, 10.0])
        with pytest.raises(KeyError):
            CalculatorTenBerge(tenbergeCoefficient=1).calculate(df, field=None, time="datetime")

    def test_to_json_names_its_type_and_records_n(self):
        j = CalculatorTenBerge(tenbergeCoefficient=2.5).toJSON()
        assert j["type"] == "tenBerge"
        assert j["n"] == "2.5"


@pytest.mark.unit
class TestInjuryCalculateToxicLoadsOnPandas:
    """Neither calculator choice survives the Injury-level API on pandas
    input: TenBerge hits B63 (no inUnits passthrough), Haber hits B62
    first (wrong positional slot) and would still hit B63 after that."""

    def test_toxic_loads_via_haber_raises_on_the_positional_mismatch(self):
        """Characterisation of B62."""
        injury = _agent_with_calculator("Haber")
        df = _series_dataframe([0.0, 5.0, 20.0])
        with pytest.raises(TypeError, match="breathingRate"):
            injury.calculateToxicLoads(df, time="datetime")

    def test_toxic_loads_via_tenberge_raises_on_the_missing_attrs_key(self):
        """Characterisation of B63: calculateToxicLoads has no inUnits
        parameter, so there is no way around the attrs lookup at all."""
        injury = _agent_with_calculator("TenBerge", tenbergeCoefficient=1)
        df = _series_dataframe([0.0, 5.0, 20.0])
        with pytest.raises(KeyError):
            injury.calculateToxicLoads(df, time="datetime")

    @pytest.mark.xfail(
        strict=True,
        reason="B63: Injury.calculateToxicLoads has no inUnits parameter to "
               "forward to the calculator, so concentrationField.attrs[field] "
               "is always attempted for pandas input and always raises "
               "KeyError, regardless of which calculator is configured. "
               "See the consolidated findings issue.",
    )
    def test_toxic_loads_via_tenberge_should_work_on_pandas_input(self):
        injury = _agent_with_calculator("TenBerge", tenbergeCoefficient=1)
        df = _series_dataframe([0.0, 5.0, 20.0])
        toxic = injury.calculateToxicLoads(df, time="datetime")
        assert toxic["P1"].iloc[-1] > 0

    def test_calculate_point_wise_fraction_injured_rejects_non_pandas_input(self):
        injury = _agent_with_calculator("TenBerge", tenbergeCoefficient=1)
        with pytest.raises(ValueError, match="not implemented"):
            injury.calculatePointWiseFractionInjured({"not": "a dataframe"}, time="datetime")

    @pytest.mark.xfail(
        strict=True,
        reason="B63, same root cause, on the documented pandas entry point.",
    )
    def test_calculate_point_wise_fraction_injured_should_work_on_pandas_input(self):
        injury = _agent_with_calculator("TenBerge", tenbergeCoefficient=1)
        df = _series_dataframe([0.0, 5.0, 20.0, 20.0, 20.0])
        frac = injury.calculatePointWiseFractionInjured(df, time="datetime")
        assert set(frac["level"]) == {"Severe"}

    def test_calculate_point_wise_fraction_injured_currently_raises(self):
        """Characterisation of B63 on calculatePointWiseFractionInjured."""
        injury = _agent_with_calculator("TenBerge", tenbergeCoefficient=1)
        df = _series_dataframe([0.0, 5.0, 20.0, 20.0, 20.0])
        with pytest.raises(KeyError):
            injury.calculatePointWiseFractionInjured(df, time="datetime")

    def test_the_deprecated_calculate_warns_before_hitting_the_same_wall(self):
        injury = _agent_with_calculator("TenBerge", tenbergeCoefficient=1)
        df = _series_dataframe([0.0, 5.0, 20.0])
        with pytest.warns(UserWarning, match="obselete"):
            with pytest.raises(KeyError):
                injury.calculate(df, field=None, time="datetime")

    def test_to_json_nests_levels_and_the_calculator(self):
        injury = _agent_with_calculator("TenBerge", tenbergeCoefficient=1)
        j = injury.toJSON()
        assert set(j["levels"]) == {"Severe"}
        assert j["calculator"]["type"] == "tenBerge"


@pytest.mark.unit
class TestInjuryFactory:
    def test_the_factory_name_is_never_assigned_anywhere(self):
        """Characterisation: _name has no setter and __init__ never runs,
        so the singleton's `name` is always None."""
        from hera.riskassessment.agents.effects import injuryfactory

        assert injuryfactory.name is None
