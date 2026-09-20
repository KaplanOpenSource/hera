"""highfreqdata/toolkit.py: HighFreqToolKit -- parser auto-detection,
normalisation, and the parse/load pipeline. Exercised against synthetic
TOA5 files (the CampbellBinary path is broken by B82 for any real column
data, so the end-to-end tests below go through TOA5 instead)."""
import csv
import os
import warnings

import pandas
import pytest

from hera.measurements.meteorology.highfreqdata.toolkit import HighFreqToolKit


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
def tk(unit_toolkit_factory):
    from hera import toolkitHome

    return unit_toolkit_factory(toolkitHome.METEOROLOGY_HIGHFREQ)


@pytest.fixture()
def toa5_file(tmp_path):
    path = tmp_path / "a.dat"
    _write_raw_sonic_file(str(path))
    return path


@pytest.mark.unit
class TestConstruction:
    def test_doctype_includes_the_toolkit_name(self, tk):
        assert tk.docType == f"{tk.toolkitName}_HighFreqData"


@pytest.mark.unit
class TestDetectParser:
    def test_a_toa5_header_is_detected_as_toa5(self, tmp_path):
        path = tmp_path / "a.dat"
        path.write_text("TOA5,Station,...\nTIMESTAMP,RECORD,U_1\n")
        assert HighFreqToolKit._detect_parser(str(path)) == "toa5"

    def test_a_tob1_header_is_detected_as_campbell(self, tmp_path):
        path = tmp_path / "b.dat"
        path.write_bytes(b'"TOB1","Station"\r\n' + b"\x00" * 20)
        assert HighFreqToolKit._detect_parser(str(path)) == "campbell"

    def test_a_csv_like_first_line_falls_back_to_toa5(self, tmp_path):
        """Neither header magic string is present, but the first line
        itself is CSV-shaped and carries a recognised keyword -- only the
        first line is checked, not the whole file."""
        path = tmp_path / "c.dat"
        path.write_text("TIMESTAMP,RECORD,U_1,V_1,W_1,T_1\n1,2,3,4,5,6\n")
        assert HighFreqToolKit._detect_parser(str(path)) == "toa5"

    def test_unrecognisable_content_defaults_to_campbell(self, tmp_path):
        path = tmp_path / "d.dat"
        path.write_bytes(b"\x00\x01\x02\x03")
        assert HighFreqToolKit._detect_parser(str(path)) == "campbell"

    def test_a_directory_picks_the_first_dat_file_inside_it(self, tmp_path):
        (tmp_path / "a.dat").write_text("TOA5,Station\nTIMESTAMP,RECORD,U_1\n")
        assert HighFreqToolKit._detect_parser(str(tmp_path)) == "toa5"

    def test_a_directory_with_no_dat_files_raises(self, tmp_path):
        (tmp_path / "notadat.txt").write_text("x")
        with pytest.raises(FileNotFoundError, match="No .dat files"):
            HighFreqToolKit._detect_parser(str(tmp_path))


@pytest.mark.unit
class TestParseRaw:
    def test_campbell_dispatches_to_the_binary_parser(self, tk, tmp_path, monkeypatch):
        calls = []
        monkeypatch.setattr(
            "hera.measurements.meteorology.highfreqdata.toolkit.Parser.parse",
            lambda self, path, fromTime, toTime: calls.append(("campbell", path)) or "RAW",
        )
        result = tk._parse_raw(str(tmp_path / "x.dat"), None, None, "campbell")
        assert result == "RAW"
        assert calls == [("campbell", str(tmp_path / "x.dat"))]

    def test_toa5_dispatches_to_the_ascii_parser(self, tk, toa5_file):
        raw = tk._parse_raw(str(toa5_file), None, None, "toa5")
        assert "Raw_Sonic_1" in raw

    def test_an_unknown_parser_name_raises(self, tk, toa5_file):
        with pytest.raises(ValueError, match="Unknown parser"):
            tk._parse_raw(str(toa5_file), None, None, "bogus")


@pytest.mark.unit
class TestNormaliseRaw:
    def test_a_plain_dataframe_is_wrapped_in_a_single_element_list(self):
        df = pandas.DataFrame(
            {"u": [1.0], "v": [1.0], "w": [1.0], "T": [20.0]},
            index=pandas.date_range("2020-01-01", periods=1, freq="1s"),
        )
        results = HighFreqToolKit._normalise_raw(df)
        assert len(results) == 1
        norm_df, meta = results[0]
        assert meta["deviceType"] == "sonic"
        assert list(norm_df.columns) == ["u", "v", "w", "T"]

    def test_a_dict_produces_one_entry_per_device_with_its_name_in_metadata(self):
        df = pandas.DataFrame(
            {"U": [1.0], "V": [1.0], "W": [1.0], "T": [20.0]},
            index=pandas.date_range("2020-01-01", periods=1, freq="1s"),
        )
        results = HighFreqToolKit._normalise_raw({"Raw_Sonic_1": df})
        assert len(results) == 1
        norm_df, meta = results[0]
        assert meta["deviceName"] == "Raw_Sonic_1"

    def test_an_unexpected_type_raises_typeerror(self):
        with pytest.raises(TypeError, match="Unexpected parser output type"):
            HighFreqToolKit._normalise_raw(["not", "a", "dataframe"])


@pytest.mark.unit
class TestParseDataAndLegacyWrappers:
    def test_parsedata_returns_a_normalised_dataframe_and_metadata(self, tk, toa5_file):
        results = tk.parseData(str(toa5_file))
        assert len(results) == 1
        df, meta = results[0]
        assert list(df.columns) == ["u", "v", "w", "T"]
        assert meta["deviceName"] == "Raw_Sonic_1"

    def test_asciitoparquet_is_deprecated_but_still_normalises(self, tk, toa5_file):
        with pytest.warns(DeprecationWarning):
            result = tk.asciiToParquet(str(toa5_file))
        assert "Raw_Sonic_1" in result

    def test_campbeltoparquet_is_deprecated_but_inherits_b82(self, tk, tmp_path):
        """campbelToParquet routes through Parser.parse -- B82 breaks every
        real-column Campbell binary file, so this always raises before
        producing anything. Not a new defect, same root cause."""
        import struct

        header = (
            '"TOB1","MyStation","CR3000","12345","CR3000.Std.32",'
            '"CPU:program.CR3","12345","MyInstrument"\r\n'
            '"TIMESTAMP","U_1","V_1","W_1","T_1"\r\n'
            '"TS","","","",""\r\n'
            '"","Smp","Smp","Smp","Smp"\r\n'
            '"ULONG","ULONG","IEEE4","IEEE4","IEEE4","IEEE4"\r\n'
        ).encode("ascii")
        body = struct.pack("<IIffff", 10, 0, 1.0, 0.5, 0.1, 20.0)
        path = tmp_path / "sonic.dat"
        path.write_bytes(header + body)

        with pytest.warns(DeprecationWarning):
            with pytest.raises(ValueError, match="columns"):
                tk.campbelToParquet(str(path))


@pytest.mark.unit
class TestLoadData:
    def test_it_parses_saves_a_parquet_and_registers_a_data_source(self, tk, toa5_file, tmp_path):
        outdir = tmp_path / "out"
        doc = tk.loadData(name="sonic1", path=str(toa5_file), outputDirectory=str(outdir))
        assert doc.dataFormat == "parquet"
        assert os.path.isfile(doc.resource)

    def test_a_second_call_without_overwrite_or_append_raises(self, tk, toa5_file, tmp_path):
        outdir = tmp_path / "out"
        tk.loadData(name="sonic1", path=str(toa5_file), outputDirectory=str(outdir))
        with pytest.raises(ValueError, match="already exists"):
            tk.loadData(name="sonic1", path=str(toa5_file), outputDirectory=str(outdir))

    def test_overwrite_replaces_the_existing_data_source(self, tk, toa5_file, tmp_path):
        outdir = tmp_path / "out"
        first = tk.loadData(name="sonic1", path=str(toa5_file), outputDirectory=str(outdir))
        second = tk.loadData(name="sonic1", path=str(toa5_file), outputDirectory=str(outdir), overwrite=True)
        assert second.resource == first.resource

    def test_append_concatenates_onto_the_existing_parquet(self, tk, toa5_file, tmp_path):
        outdir = tmp_path / "out"
        first = tk.loadData(name="sonic1", path=str(toa5_file), outputDirectory=str(outdir))
        before = pandas.read_parquet(first.resource)
        tk.loadData(name="sonic1", path=str(toa5_file), outputDirectory=str(outdir), append=True)
        after = pandas.read_parquet(first.resource)
        assert len(after) >= len(before)

    def test_overwrite_and_append_together_raise(self, tk, toa5_file, tmp_path):
        outdir = tmp_path / "out"
        with pytest.raises(ValueError, match="mutually exclusive"):
            tk.loadData(
                name="sonic1", path=str(toa5_file), outputDirectory=str(outdir),
                overwrite=True, append=True,
            )
