"""measurements/meteorology/analysis.py: addDatesColumns and calcHourlyDist
used to be defined inside the analysis class with no `self` parameter -- same
bug class as B93/B100 -- so calling them the normal way (on an instance) bound
the instance to the first parameter and broke.  B104 is fixed: both now take
`self` like their lowfreqdata counterparts.
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
class TestAddDatesColumns:
    def test_calling_it_on_an_instance_adds_the_columns(self, df):
        a = analysis(None)
        result = a.addDatesColumns(df)
        assert "yearonly" in result.columns
        assert "season" in result.columns

    def test_the_date_column_can_be_named_explicitly(self, df):
        a = analysis(None)
        withDate = df.assign(when=df.index)
        result = a.addDatesColumns(withDate, datecolumn="when")
        assert list(result["yearonly"]) == [2020, 2020, 2020]


@pytest.mark.unit
class TestCalcHourlyDist:
    def test_calling_it_on_an_instance_computes_a_distribution(self, df):
        a = analysis(None)
        x_mid, y_mid, M = a.calcHourlyDist(df, "v", bins=3)
        assert len(x_mid) == 3
        assert len(y_mid) == 3
        assert M.shape == (3, 3)
