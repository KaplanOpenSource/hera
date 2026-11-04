"""Normalisation helpers shared by the high-frequency parsers."""
import pandas
import pytest

from hera.measurements.meteorology.highfreqdata.parsers import (
    detect_device_type,
    normalise_sonic_df,
    normalise_trh_df,
)


@pytest.mark.unit
class TestNormaliseSonicDf:
    def test_uppercase_uvw_columns_are_lowercased(self):
        df = pandas.DataFrame({
            "TIMESTAMP": pandas.date_range("2020-01-01", periods=2, freq="1s"),
            "U": [1.0, 2.0], "V": [0.5, 0.6], "W": [0.1, 0.2], "T": [20.0, 20.1],
        })
        normalised, _ = normalise_sonic_df(df)
        assert list(normalised.columns) == ["u", "v", "w", "T"]

    def test_the_timestamp_column_becomes_the_index(self):
        df = pandas.DataFrame({
            "TIMESTAMP": pandas.date_range("2020-01-01", periods=2, freq="1s"),
            "U": [1.0, 2.0], "V": [0.5, 0.6], "W": [0.1, 0.2],
        })
        normalised, _ = normalise_sonic_df(df)
        assert isinstance(normalised.index, pandas.DatetimeIndex)
        assert normalised.index.name == "Time"

    def test_height_station_and_instrument_are_pulled_into_metadata(self):
        df = pandas.DataFrame({
            "TIMESTAMP": pandas.date_range("2020-01-01", periods=2, freq="1s"),
            "U": [1.0, 2.0], "V": [0.5, 0.6], "W": [0.1, 0.2],
            "height": [10, 10], "station": ["A", "A"],
        })
        normalised, metadata = normalise_sonic_df(df)
        assert metadata["height"] == 10
        assert metadata["station"] == "A"
        assert "height" not in normalised.columns
        assert "station" not in normalised.columns

    def test_multiple_distinct_metadata_values_are_kept_as_a_list(self):
        df = pandas.DataFrame({
            "TIMESTAMP": pandas.date_range("2020-01-01", periods=2, freq="1s"),
            "U": [1.0, 2.0], "V": [0.5, 0.6], "W": [0.1, 0.2],
            "height": [10, 20],
        })
        _, metadata = normalise_sonic_df(df)
        assert sorted(metadata["height"]) == [10, 20]

    def test_record_column_is_dropped(self):
        df = pandas.DataFrame({
            "TIMESTAMP": pandas.date_range("2020-01-01", periods=2, freq="1s"),
            "U": [1.0, 2.0], "V": [0.5, 0.6], "W": [0.1, 0.2], "RECORD": [1, 2],
        })
        normalised, _ = normalise_sonic_df(df)
        assert "RECORD" not in normalised.columns

    def test_non_numeric_values_in_a_data_column_become_nan(self):
        df = pandas.DataFrame({
            "TIMESTAMP": pandas.date_range("2020-01-01", periods=2, freq="1s"),
            "U": ["1.0", "not_a_number"], "V": [0.5, 0.6], "W": [0.1, 0.2],
        })
        normalised, _ = normalise_sonic_df(df)
        assert normalised["u"].iloc[0] == pytest.approx(1.0)
        assert pandas.isna(normalised["u"].iloc[1])

    def test_the_device_type_is_recorded_as_sonic(self):
        df = pandas.DataFrame({
            "TIMESTAMP": pandas.date_range("2020-01-01", periods=1, freq="1s"),
            "U": [1.0], "V": [0.5], "W": [0.1],
        })
        _, metadata = normalise_sonic_df(df)
        assert metadata["deviceType"] == "sonic"

    def test_an_already_indexed_datetime_frame_keeps_its_index(self):
        times = pandas.date_range("2020-01-01", periods=2, freq="1s")
        df = pandas.DataFrame({"U": [1.0, 2.0], "V": [0.5, 0.6], "W": [0.1, 0.2]}, index=times)
        normalised, _ = normalise_sonic_df(df)
        assert list(normalised.index) == list(times)
        assert normalised.index.name == "Time"

    def test_an_unparseable_index_is_left_as_is_without_raising(self):
        df = pandas.DataFrame({"U": [1.0], "V": [0.5], "W": [0.1]}, index=["not-a-date"])
        normalised, _ = normalise_sonic_df(df)
        assert list(normalised.index) == ["not-a-date"]


@pytest.mark.unit
class TestNormaliseTrhDf:
    def test_tct_is_renamed_to_tc_t(self):
        df = pandas.DataFrame({
            "TIMESTAMP": pandas.date_range("2020-01-01", periods=2, freq="1s"),
            "TcT": [20.0, 20.1], "TRH": [50.0, 51.0], "RH": [45.0, 46.0],
        })
        normalised, _ = normalise_trh_df(df)
        assert "TC_T" in normalised.columns
        assert "TcT" not in normalised.columns

    def test_the_device_type_is_recorded_as_trh(self):
        df = pandas.DataFrame({
            "TIMESTAMP": pandas.date_range("2020-01-01", periods=1, freq="1s"),
            "TcT": [20.0],
        })
        _, metadata = normalise_trh_df(df)
        assert metadata["deviceType"] == "TRH"

    def test_only_known_trh_columns_are_kept(self):
        df = pandas.DataFrame({
            "TIMESTAMP": pandas.date_range("2020-01-01", periods=1, freq="1s"),
            "TcT": [20.0], "extraneous": ["x"],
        })
        normalised, _ = normalise_trh_df(df)
        assert "extraneous" not in normalised.columns


@pytest.mark.unit
class TestDetectDeviceType:
    def test_uvw_columns_are_detected_as_sonic(self):
        assert detect_device_type(pandas.DataFrame(columns=["U", "V", "W"])) == "sonic"

    def test_tc_t_is_detected_as_trh(self):
        assert detect_device_type(pandas.DataFrame(columns=["TC_T", "RH"])) == "TRH"

    def test_tct_lowercased_is_also_detected_as_trh(self):
        assert detect_device_type(pandas.DataFrame(columns=["TcT"])) == "TRH"

    def test_unrelated_columns_are_unknown(self):
        assert detect_device_type(pandas.DataFrame(columns=["foo", "bar"])) == "unknown"

    def test_detection_is_case_insensitive(self):
        assert detect_device_type(pandas.DataFrame(columns=["u", "v", "w"])) == "sonic"
