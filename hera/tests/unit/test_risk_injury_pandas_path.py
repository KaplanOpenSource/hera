"""Injury and the calculators, exercised on pandas.DataFrame input.

Batch 5 exercised the calculators against xarray only; ``Injury`` also
supports a plain pandas.DataFrame ("each point is a column, time is the
index"), which turns out to have three independent problems:

* B62, fixed: ``Injury.calculateToxicLoads`` called
  ``self.calculator.calculate(concentrationField, field, breathingRate=...,
  time=...)`` -- passing ``field`` positionally. That matched
  ``CalculatorTenBerge``/``CalculatorMaxConcentration``, whose second
  positional parameter is ``field``, but not ``CalculatorHaber``, whose
  signature is ``(self, concentrationField, breathingRate=..., time=...,
  field=None, ...)``, so any ``Injury`` configured with the Haber calculator
  raised ``TypeError`` the moment it was used. It is passed by keyword now.
* B36/B63, fixed: all three calculators decided whether to read units from
  ``concentrationField.attrs[field]`` by checking
  ``hasattr(concentrationField, "attrs")`` -- meant to detect an xarray
  object. ``pandas.DataFrame`` has carried an (empty) ``.attrs`` dict since
  pandas 1.0, so the same branch ran for pandas input too and
  ``df.attrs[field]`` raised ``KeyError``. Calling the calculator directly
  worked around it with an explicit ``inUnits=``, but
  ``Injury.calculateToxicLoads``/``calculatePointWiseFractionInjured`` had
  no ``inUnits`` parameter to forward one through -- so through the
  ``Injury`` API every pandas call hit it unconditionally. The branch now
  tests ``isinstance(concentrationField, xarray.Dataset)``, like the
  computation below it.
* B64, fixed: ``CalculatorHaber.calculate``'s pandas branch computed
  ``concentrationField[:-1].fillna(0) * dt_min[1:]`` -- a DataFrame
  multiplied by a Series without ``axis=0``, so pandas aligned the Series'
  index against the DataFrame's *columns* rather than its rows. With no
  matching labels the result was an all-NaN frame with the Series' index
  values as bogus extra columns. It now converts ``dt_min`` to
  ``.values.reshape(...)`` before multiplying, exactly as
  ``CalculatorTenBerge``'s equivalent line always did.
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
    def test_a_constant_exposure_accumulates_linearly(self):
        """Five samples a minute apart span four known intervals, so the
        dose is 40 -- the same number ten Berge with n = 1 gives below."""
        df = _series_dataframe([10.0, 10.0, 10.0, 10.0, 10.0])
        dose = CalculatorHaber().calculate(
            df, time="datetime", field=None, inUnits=1 * ureg.mg / ureg.m**3
        )
        assert dose["P1"].iloc[-1] == pytest.approx(40.0)

    def test_the_input_columns_are_the_only_columns(self):
        """B64 surfaced as bogus integer columns and an all-NaN frame, from
        the Series index being aligned against the DataFrame's columns."""
        df = _series_dataframe([10.0, 10.0, 10.0, 10.0, 10.0])
        dose = CalculatorHaber().calculate(
            df, time="datetime", field=None, inUnits=1 * ureg.mg / ureg.m**3
        )
        assert set(dose.columns) == {"P1"}
        assert not dose.isna().any().any()

    def test_the_documented_default_units_apply(self):
        """B36, fixed: hasattr(df, 'attrs') is true for pandas.DataFrame
        since pandas 1.0, so the xarray-only branch used to run for pandas
        too and df.attrs[field] raised KeyError."""
        df = _series_dataframe([10.0, 10.0])
        dose = CalculatorHaber().calculate(df, time="datetime", field=None)
        assert dose["P1"].iloc[-1] == pytest.approx(10.0)


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

    def test_the_documented_default_units_apply(self):
        """B36 again, on TenBerge."""
        df = _series_dataframe([10.0, 10.0])
        dose = CalculatorTenBerge(tenbergeCoefficient=1).calculate(df, field=None, time="datetime")
        assert dose["P1"].iloc[-1] == pytest.approx(10.0)

    def test_to_json_names_its_type_and_records_n(self):
        j = CalculatorTenBerge(tenbergeCoefficient=2.5).toJSON()
        assert j["type"] == "tenBerge"
        assert j["n"] == "2.5"


@pytest.mark.unit
class TestInjuryCalculateToxicLoadsOnPandas:
    """The Injury-level API on pandas input.  Haber used to hit B62 (the
    wrong positional slot) and both calculators used to hit B36 (the attrs
    lookup) before that; the pointwise entry point still runs into the
    recorded B39/B12 cross-registry defect in InjuryLevelThreshold."""

    def test_toxic_loads_via_haber_work_on_pandas_input(self):
        """B62, fixed: `field` used to be passed positionally, which landed
        in CalculatorHaber's `breathingRate` slot."""
        injury = _agent_with_calculator("Haber")
        df = _series_dataframe([0.0, 5.0, 20.0])
        toxic = injury.calculateToxicLoads(df, time="datetime")
        assert toxic["P1"].tolist() == pytest.approx([0.0, 5.0])

    def test_toxic_loads_via_tenberge_work_on_pandas_input(self):
        """B36, fixed: the attrs lookup used to run for pandas input too,
        and calculateToxicLoads had no inUnits to route around it."""
        injury = _agent_with_calculator("TenBerge", tenbergeCoefficient=1)
        df = _series_dataframe([0.0, 5.0, 20.0])
        toxic = injury.calculateToxicLoads(df, time="datetime")
        assert toxic["P1"].iloc[-1] > 0

    def test_an_explicit_inunits_is_forwarded_to_the_calculator(self):
        """B63, fixed: the Injury API had no inUnits parameter at all, so a
        pandas field in any other unit could not be rescued."""
        injury = _agent_with_calculator("TenBerge", tenbergeCoefficient=1)
        df = _series_dataframe([10.0, 10.0, 10.0])

        inGrams = injury.calculateToxicLoads(
            df, time="datetime", inUnits=1 * ureg.g / ureg.m**3)
        default = injury.calculateToxicLoads(df, time="datetime")

        assert inGrams["P1"].iloc[-1] == pytest.approx(1000 * default["P1"].iloc[-1])

    def test_calculate_point_wise_fraction_injured_rejects_non_pandas_input(self):
        injury = _agent_with_calculator("TenBerge", tenbergeCoefficient=1)
        with pytest.raises(ValueError, match="not implemented"):
            injury.calculatePointWiseFractionInjured({"not": "a dataframe"}, time="datetime")

    @pytest.mark.xfail(
        strict=True,
        reason="B39/B12 (recorded, and pinned in test_utils_units.py and "
               "test_risk_injurylevel.py): the toxic loads now reach "
               "InjuryLevelThreshold.getPercent as plain numbers, and tounit "
               "builds those in pint's default registry while the threshold "
               "lives in hera's, so the comparison raises. "
               "See the consolidated findings issue.",
    )
    def test_calculate_point_wise_fraction_injured_should_work_on_pandas_input(self):
        injury = _agent_with_calculator("TenBerge", tenbergeCoefficient=1)
        df = _series_dataframe([0.0, 5.0, 20.0, 20.0, 20.0])
        frac = injury.calculatePointWiseFractionInjured(df, time="datetime")
        assert set(frac["level"]) == {"Severe"}

    def test_calculate_point_wise_fraction_injured_currently_raises(self):
        """Characterisation of B39/B12 on calculatePointWiseFractionInjured:
        the dose is computed, and only the threshold comparison fails."""
        injury = _agent_with_calculator("TenBerge", tenbergeCoefficient=1)
        df = _series_dataframe([0.0, 5.0, 20.0, 20.0, 20.0])
        with pytest.raises(ValueError, match="different registries"):
            injury.calculatePointWiseFractionInjured(df, time="datetime")

    def test_the_deprecated_calculate_warns_and_is_xarray_only(self):
        """calculateRegionOfInjured selects with .sel()/.isel(), so the
        deprecated entry point is reachable but xarray-only."""
        injury = _agent_with_calculator("TenBerge", tenbergeCoefficient=1)
        df = _series_dataframe([0.0, 5.0, 20.0])
        with pytest.warns(UserWarning, match="obselete"):
            with pytest.raises(AttributeError, match="no attribute 'sel'"):
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
