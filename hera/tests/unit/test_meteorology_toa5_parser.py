"""highfreqdata/parsers/TOA5.py: ASCIIParser, against synthetic .dat files
in the real Campbell Scientific TOA5 layout (4 header rows -- station
info, real column names, units, aggregation type -- then data rows)."""
import csv

import pandas
import pytest

from hera.measurements.meteorology.highfreqdata.parsers.TOA5 import ASCIIParser


def _write_raw_sonic_file(path, start_record=1, start_second=0, n=3):
    rows = [
        ["TOA5", "Station", "CR1000", "SN", "OS", "Prog", "Sig", "Raw_Sonic"],
        ["TIMESTAMP", "RECORD", "U_1", "V_1", "W_1", "T_1"],
        ["TS", "RN", "m/s", "m/s", "m/s", "C"],
        ["", "", "Smp", "Smp", "Smp", "Smp"],
    ]
    for i in range(n):
        rows.append([
            f"2020-01-01 00:00:{start_second + i:02d}",
            str(start_record + i),
            str(1.0 + i), str(2.0 + i), str(3.0 + i), str(20.0 + i),
        ])
    with open(path, "w", newline="") as f:
        csv.writer(f).writerows(rows)


@pytest.fixture()
def parser():
    return ASCIIParser()


@pytest.mark.unit
class TestGetColumns:
    def test_metadata_identifies_raw_sonic_and_counts_devices(self, parser, tmp_path):
        path = tmp_path / "a.dat"
        _write_raw_sonic_file(str(path))
        cols, n_devices, has_metadata = parser.get_columns(str(path))
        assert cols == ["U", "V", "W", "T"]
        assert n_devices == 1
        assert has_metadata is True

    def test_metadata_identifies_tct_trh(self, parser, tmp_path):
        path = tmp_path / "b.dat"
        with open(path, "w", newline="") as f:
            csv.writer(f).writerows([
                ["TOA5", "St", "CR", "SN", "OS", "Prog", "Sig", "Tct_Trh"],
                ["TIMESTAMP", "RECORD", "TC_T1", "TRH1", "RH1"],
            ])
        cols, n_devices, has_metadata = parser.get_columns(str(path))
        assert cols == ["TC_T", "TRH", "RH"]
        assert has_metadata is True

    def test_no_metadata_falls_back_to_the_filename(self, parser, tmp_path):
        path = tmp_path / "log_Raw_Sonic.dat"
        with open(path, "w", newline="") as f:
            csv.writer(f).writerows([["a", "b", "U_1", "V_1", "W_1", "T_1"], ["x"]])
        cols, n_devices, has_metadata = parser.get_columns(str(path))
        assert cols == ["U", "V", "W", "T"]
        assert has_metadata is False

    def test_no_metadata_and_no_filename_hint_raises(self, parser, tmp_path):
        path = tmp_path / "unlabeled.dat"
        with open(path, "w", newline="") as f:
            csv.writer(f).writerows([["a", "b", "c"], ["x"]])
        with pytest.raises(ValueError, match="No Device could be detected"):
            parser.get_columns(str(path))


@pytest.mark.unit
class TestGetPandasFromFile:
    def test_it_renames_columns_and_drops_the_header_rows(self, parser, tmp_path):
        path = tmp_path / "a.dat"
        _write_raw_sonic_file(str(path))
        dfs = parser.getPandasFromFile(str(path), None, None)
        assert list(dfs.keys()) == ["Raw_Sonic_1"]
        df = dfs["Raw_Sonic_1"]
        assert list(df.columns) == ["TIMESTAMP", "RECORD", "U", "V", "W", "T"]
        assert len(df) == 3
        assert df["U"].tolist() == ["1.0", "2.0", "3.0"]

    def test_parse_dispatches_a_dat_file_to_getpandasfromfile(self, parser, tmp_path):
        path = tmp_path / "a.dat"
        _write_raw_sonic_file(str(path))
        dfs = parser.parse(str(path))
        assert "Raw_Sonic_1" in dfs


@pytest.mark.unit
class TestGetPandasFromDir:
    def test_it_concatenates_and_time_sorts_all_files_in_a_directory(self, parser, tmp_path):
        _write_raw_sonic_file(str(tmp_path / "a.dat"), start_record=1, start_second=0)
        _write_raw_sonic_file(str(tmp_path / "b.dat"), start_record=4, start_second=3)
        dfs = parser.getPandasFromDir(str(tmp_path), None, None)
        assert len(dfs["Raw_Sonic_1"]) == 6
        assert dfs["Raw_Sonic_1"]["RECORD"].tolist() == ["1", "2", "3", "4", "5", "6"]

    def test_parse_dispatches_a_directory_to_getpandasfromdir(self, parser, tmp_path):
        _write_raw_sonic_file(str(tmp_path / "a.dat"))
        dfs = parser.parse(str(tmp_path))
        assert "Raw_Sonic_1" in dfs


@pytest.mark.unit
class TestFromTimeToTimeHandler:
    @pytest.fixture()
    def df(self):
        return pandas.DataFrame({
            "TIMESTAMP": pandas.date_range("2020-01-01", periods=5, freq="1s"),
            "V": range(5),
        })

    def test_both_bounds_narrow_to_the_interval(self, parser, df):
        result = parser.fromTime_toTime_handler(df, "2020-01-01 00:00:01", "2020-01-01 00:00:03")
        assert result["V"].tolist() == [1, 2, 3]

    def test_only_fromtime_extends_to_the_end(self, parser, df):
        result = parser.fromTime_toTime_handler(df, "2020-01-01 00:00:02", None)
        assert result["V"].tolist() == [2, 3, 4]

    def test_only_totime_extends_from_the_start(self, parser, df):
        result = parser.fromTime_toTime_handler(df, None, "2020-01-01 00:00:02")
        assert result["V"].tolist() == [0, 1, 2]
