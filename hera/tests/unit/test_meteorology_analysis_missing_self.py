"""measurements/meteorology/analysis.py: addDatesColumns and calcHourlyDist
are defined inside the analysis class with no `self` parameter -- same
bug class as B93/B100. Calling them the normal way (on an instance) binds
the instance to the first parameter and breaks; they only work called
directly on the class.
"""
import pandas
import pytest

from hera.measurements.meteorology.analysis import analysis


@pytest.fixture()
def df():
    return pandas.DataFrame({"v": [1.0, 2.0, 3.0]}, index=pandas.date_range("2020-01-01", periods=3, freq="4h"))


@pytest.mark.unit
class TestAnalysisConstruction:
    def test_datalayer_is_stored(self):
        a = analysis("myDataLayer")
        assert a.datalayer == "myDataLayer"


@pytest.mark.unit
class TestAddDatesColumnsIsBroken:
    @pytest.mark.xfail(
        strict=True,
        reason="B104: addDatesColumns has no `self` parameter (first "
               "positional is `data`), so calling it on an instance binds "
               "the instance to `data` instead. See the consolidated "
               "findings issue.",
    )
    def test_calling_it_on_an_instance_should_add_the_columns(self, df):
        a = analysis(None)
        a.addDatesColumns(df)

    def test_calling_it_on_an_instance_currently_raises(self, df):
        """Characterisation of B104."""
        a = analysis(None)
        with pytest.raises(AttributeError):
            a.addDatesColumns(df)

    def test_it_only_works_called_directly_on_the_class(self, df):
        """Characterisation of B104."""
        result = analysis.addDatesColumns(df)
        assert "yearonly" in result.columns
        assert "season" in result.columns


@pytest.mark.unit
class TestCalcHourlyDistIsBroken:
    @pytest.mark.xfail(
        strict=True,
        reason="B104: same missing-`self` bug as addDatesColumns. "
               "See the consolidated findings issue.",
    )
    def test_calling_it_on_an_instance_should_compute_a_distribution(self, df):
        a = analysis(None)
        a.calcHourlyDist(df, "v")

    def test_calling_it_on_an_instance_currently_raises(self, df):
        """Characterisation of B104."""
        a = analysis(None)
        with pytest.raises(AttributeError):
            a.calcHourlyDist(df, "v")
