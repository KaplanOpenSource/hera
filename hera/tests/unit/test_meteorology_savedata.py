"""The data-persistence dispatcher shared by the high-frequency calculators.

B69, fixed: ``SaveDataHandler.getSaveData_HDF`` used to call
``data.to_HDF(path, key)`` -- capitalised -- but pandas only has ``to_hdf``
(lowercase, and ``key`` keyword-only since 2.0), so saving to HDF5 always
raised ``AttributeError``.
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
class TestSaveDataHandlerHDF:
    def test_hdf_saving_writes_a_file_and_returns_path_and_key(self, tmp_path):
        pytest.importorskip("tables")
        df = pandas.DataFrame({"a": [1, 2]})
        target = str(tmp_path / "out.h5")
        result = SaveDataHandler.getSaveData_HDF(df, target, "k")
        assert result == dict(path=target, key="k")
        assert (tmp_path / "out.h5").is_file()
        pandas.testing.assert_frame_equal(pandas.read_hdf(target, key="k"), df)

    def test_it_no_longer_calls_the_capitalised_name(self):
        import inspect

        source = inspect.getsource(SaveDataHandler.getSaveData_HDF)
        assert "to_HDF" not in source
        assert "to_hdf" in source
