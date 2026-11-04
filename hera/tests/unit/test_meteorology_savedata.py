"""The data-persistence dispatcher shared by the high-frequency calculators.

B69: ``SaveDataHandler.getSaveData_HDF`` calls ``data.to_HDF(path, key)`` --
capitalised -- but pandas only has ``to_hdf`` (lowercase). Saving to HDF5
always raises ``AttributeError``.
"""
import pandas
import pytest

from hera.measurements.meteorology.highfreqdata.analysis.abstractcalculator import (
    AbstractCalculator,
    SaveDataHandler,
    getSaveData,
)


@pytest.mark.unit
class TestAbstractCalculatorJoinMethod:
    def test_the_join_method_defaults_to_left(self):
        calc = AbstractCalculator(pandas.DataFrame({"u": [1.0]}), {"samplingWindow": "30min"})
        assert calc.JoinMethod == "left"


@pytest.mark.unit
class TestGetSaveDataDispatch:
    def test_it_dispatches_to_the_matching_handler_by_name(self):
        df = pandas.DataFrame({"a": [1, 2, 3]})
        assert getSaveData(dataFormat="JSON_pandas", data=df) == df.to_json()

    def test_an_unknown_format_raises_attribute_error(self):
        with pytest.raises(AttributeError):
            getSaveData(dataFormat="NoSuchFormat", data=pandas.DataFrame())


@pytest.mark.unit
class TestSaveDataHandlerJSON:
    def test_without_a_path_it_returns_the_json_string_directly(self):
        df = pandas.DataFrame({"a": [1, 2]})
        assert SaveDataHandler.getSaveData_JSON_pandas(df) == df.to_json()

    def test_with_a_path_it_writes_the_file_and_returns_the_path(self, tmp_path):
        df = pandas.DataFrame({"a": [1, 2]})
        target = str(tmp_path / "out.json")
        result = SaveDataHandler.getSaveData_JSON_pandas(df, path=target)
        assert result == target
        assert (tmp_path / "out.json").read_text() == df.to_json()


@pytest.mark.unit
class TestSaveDataHandlerParquet:
    def test_a_dask_dataframe_is_written_to_the_target_directory(self, tmp_path):
        dask_dataframe = pytest.importorskip("dask.dataframe")
        ddf = dask_dataframe.from_pandas(pandas.DataFrame({"a": [1, 2, 3]}), npartitions=1)
        target = str(tmp_path / "parquet_out")
        result = SaveDataHandler.getSaveData_parquet(ddf, target)
        assert result == target
        assert (tmp_path / "parquet_out").is_dir()
        assert list((tmp_path / "parquet_out").glob("*.parquet"))


@pytest.mark.unit
class TestSaveDataHandlerHDFIsBroken:
    @pytest.mark.xfail(
        strict=True,
        reason="B69: getSaveData_HDF calls data.to_HDF(path, key) -- "
               "capitalised -- but pandas.DataFrame only defines the "
               "lowercase to_hdf. Saving to HDF5 always raises "
               "AttributeError. See the consolidated findings issue.",
    )
    def test_hdf_saving_should_write_a_file(self, tmp_path):
        df = pandas.DataFrame({"a": [1, 2]})
        SaveDataHandler.getSaveData_HDF(df, str(tmp_path / "out.h5"), "k")

    def test_hdf_saving_currently_raises(self, tmp_path):
        """Characterisation of B69."""
        df = pandas.DataFrame({"a": [1, 2]})
        with pytest.raises(AttributeError, match="to_HDF"):
            SaveDataHandler.getSaveData_HDF(df, str(tmp_path / "out.h5"), "k")
