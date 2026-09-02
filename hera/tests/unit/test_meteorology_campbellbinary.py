"""CampbellBinary: the TOB1 binary sonic/TRH parser.

Tested against hand-built synthetic TOB1 files (5 ASCII header lines
terminated by ``\\n``, then packed binary records) rather than a real
capture, since none was available. The low-level reader
(``CampbellBinaryInterface``) works correctly against these fixtures; the
higher-level ``Parser`` does not:

* B82: ``Parser.getData`` slices each record's data with
  ``line[cbi.columnsIndexes[i][0] : cbi.columnsIndexes[i][1]]``, but
  ``columnsIndexes`` (e.g. ``[1, 5]`` for a 4-column sonic record) is
  computed as if it indexed into the *raw* record -- which starts with two
  leading timestamp fields (SECONDS, NANOSECONDS). ``line`` itself is
  ``_getDataFromStream``'s ``retval[2:]``, which has *already* dropped
  those two fields. Slicing an off-by-two range out of the already-trimmed
  list drops the first data column and returns one column short of what
  the header declares -- pandas then refuses to build the DataFrame at
  all. This breaks every call to ``getData``/``getPandasFromFile``/
  ``getPandasFromDir``/``parse``, for both the sonic (u,v,w,T) and TRH
  (TcT,TRH,RH) column layouts.
"""
import struct

import pandas
import pytest

from hera.measurements.meteorology.highfreqdata.parsers.CampbellBinary import (
    CampbellBinaryInterface,
    Parser,
)


def _write_tob1(path, header1, header4, fmt, records):
    header0 = (
        '"TOB1","MyStation","CR3000","12345","CR3000.Std.32",'
        '"CPU:program.CR3","12345","MyInstrument"\r\n'
    )
    header2 = '"TS",' + ",".join(['""'] * (len(records[0]) - 2)) + "\r\n"
    header3 = '"",' + ",".join(['"Smp"'] * (len(records[0]) - 2)) + "\r\n"
    header = (header0 + header1 + header2 + header3 + header4).encode("ascii")
    body = b"".join(struct.pack(fmt, *record) for record in records)
    path.write_bytes(header + body)


@pytest.fixture()
def sonic_file(tmp_path):
    path = tmp_path / "sonic.dat"
    _write_tob1(
        path,
        header1='"TIMESTAMP","U_1","V_1","W_1","T_1"\r\n',
        header4='"ULONG","ULONG","IEEE4","IEEE4","IEEE4","IEEE4"\r\n',
        fmt="<IIffff",
        records=[
            (10, 0, 1.0, 0.5, 0.1, 20.0),
            (11, 500000000, 1.1, 0.6, 0.2, 20.1),
            (12, 0, 1.2, 0.7, 0.3, 20.2),
        ],
    )
    return path


@pytest.fixture()
def trh_file(tmp_path):
    path = tmp_path / "trh.dat"
    _write_tob1(
        path,
        header1='"TIMESTAMP","TC_T1","TRH_1","RH_1"\r\n',
        header4='"ULONG","ULONG","IEEE4","IEEE4","IEEE4"\r\n',
        fmt="<IIfff",
        records=[(10, 0, 20.0, 50.0, 45.0), (11, 0, 20.1, 50.1, 45.1)],
    )
    return path


@pytest.mark.unit
class TestCampbellBinaryInterfaceHeaders:
    def test_the_station_and_instrument_come_from_the_first_header_line(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        assert cbi.station == "MyStation"
        assert cbi.instrument == "MyInstrument"

    def test_headers_size_stops_right_after_the_fifth_newline(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        raw = sonic_file.read_bytes()
        assert raw[cbi.headersSize - 1] == ord("\n")
        assert raw[: cbi.headersSize].count(b"\n") == 5

    def test_the_raw_format_matches_the_fifth_header_line(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        assert cbi.rawFormat == ["ULONG", "ULONG", "IEEE4", "IEEE4", "IEEE4", "IEEE4"]

    def test_the_struct_format_string_matches_the_field_types(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        assert cbi.format == "<IIffff"

    def test_record_size_is_the_packed_size_of_the_format(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        assert cbi.recordSize == struct.calcsize("<IIffff")

    def test_records_num_counts_whole_records_after_the_header(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        assert cbi.recordsNum == 3


@pytest.mark.unit
class TestColumnDetection:
    def test_sonic_columns_are_named_uvw_t(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        assert cbi.columnsNames == [["u", "v", "w", "T"]]

    def test_a_single_sonic_group_reports_height_ten(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        assert cbi.heights == [10]

    def test_trh_columns_are_named_tct_trh_rh(self, trh_file):
        cbi = CampbellBinaryInterface(file=str(trh_file))
        assert cbi.columnsNames == [["TcT", "TRH", "RH"]]


@pytest.mark.unit
class TestRecordAccess:
    def test_get_record_by_index_decodes_the_timestamp_from_1990_epoch(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        time, values = cbi.getRecordByIndex(0)
        assert time == pandas.Timestamp(1990, 1, 1) + pandas.Timedelta(seconds=10)

    def test_get_record_by_index_decodes_the_data_values(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        _, values = cbi.getRecordByIndex(0)
        assert values == pytest.approx([1.0, 0.5, 0.1, 20.0], abs=1e-5)

    def test_first_time_matches_record_zero(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        assert cbi.firstTime == cbi.getRecordByIndex(0)[0]

    def test_last_time_matches_the_final_record(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        assert cbi.lastTime == cbi.getRecordByIndex(cbi.recordsNum - 1)[0]

    def test_get_record_index_by_time_finds_the_exact_match(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        target = cbi.getRecordByIndex(2)[0]
        assert cbi.getRecordIndexByTime(target) == 2

    def test_get_time_by_record_index_matches_get_record_by_index(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        assert cbi.getTimeByRecordIndex(1) == cbi.getRecordByIndex(1)[0]


@pytest.mark.unit
class TestParserGetDataIsBroken:
    @pytest.mark.xfail(
        strict=True,
        reason="B82: columnsIndexes assumes `line` still has the two "
               "leading timestamp fields, but _getDataFromStream already "
               "strips them (retval[2:]). The resulting slice is one "
               "column short of what the header declares, and pandas "
               "refuses to build the frame. See the consolidated findings "
               "issue.",
    )
    def test_get_pandas_from_file_should_return_all_sonic_columns(self, sonic_file):
        df = Parser().getPandasFromFile(str(sonic_file), fromTime=None, toTime=None)
        assert list(df.columns[:4]) == ["u", "v", "w", "T"]

    def test_get_pandas_from_file_currently_raises_for_sonic_data(self, sonic_file):
        """Characterisation of B82."""
        with pytest.raises(ValueError, match="columns"):
            Parser().getPandasFromFile(str(sonic_file), fromTime=None, toTime=None)

    def test_get_pandas_from_file_currently_raises_for_trh_data_too(self, trh_file):
        """Same root cause, the other column layout."""
        with pytest.raises(ValueError, match="columns"):
            Parser().getPandasFromFile(str(trh_file), fromTime=None, toTime=None)

    def test_get_data_drops_the_first_data_column(self, sonic_file):
        """Direct look at getData's output, one level below the DataFrame
        construction that raises: the first value (u) is silently missing
        from every row."""
        ts, cols, values = Parser().getData(str(sonic_file), fromTime=None, toTime=None)
        assert cols == [["u", "v", "w", "T"]]
        assert len(values[10][0]) == 3  # v, w, T -- u is missing

    def test_parse_dispatches_to_get_pandas_from_file_for_a_file_path(self, sonic_file, monkeypatch):
        calls = []
        monkeypatch.setattr(
            Parser, "getPandasFromFile",
            lambda self, path, fromTime, toTime: calls.append(path) or "RESULT",
        )
        result = Parser().parse(str(sonic_file))
        assert calls == [str(sonic_file)]
        assert result == "RESULT"

    def test_parse_dispatches_to_get_pandas_from_dir_for_a_directory(self, tmp_path, monkeypatch):
        calls = []
        monkeypatch.setattr(
            Parser, "getPandasFromDir",
            lambda self, path, fromTime, toTime: calls.append(path) or "RESULT",
        )
        result = Parser().parse(str(tmp_path))
        assert calls == [str(tmp_path)]
        assert result == "RESULT"


@pytest.mark.unit
class TestGetPandasFromDir:
    def test_it_concatenates_one_dataframe_per_dat_file_in_the_directory(self, tmp_path, monkeypatch):
        """getPandasFromDir itself just globs *.dat and concatenates results
        from getPandasFromFile -- isolate that from B82 (which breaks
        getPandasFromFile's own column slicing) by stubbing it out, the
        same technique the parse() dispatch tests above use."""
        (tmp_path / "a.dat").write_bytes(b"")
        (tmp_path / "b.dat").write_bytes(b"")
        (tmp_path / "ignored.txt").write_bytes(b"")

        calls = []

        def fake_get_pandas_from_file(self, path, fromTime, toTime):
            calls.append(path)
            return pandas.DataFrame({"v": [len(calls)]})

        monkeypatch.setattr(Parser, "getPandasFromFile", fake_get_pandas_from_file)
        result = Parser().getPandasFromDir(str(tmp_path), fromTime=None, toTime=None)

        assert len(calls) == 2
        assert sorted(calls) == sorted(str(tmp_path / n) for n in ("a.dat", "b.dat"))
        assert len(result) == 2


@pytest.mark.unit
class TestCampbellBinaryInterfaceLowLevelHelpers:
    def test_bindata_returns_the_full_raw_file_contents(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        raw = sonic_file.read_bytes()
        assert cbi.binData == raw
        assert len(cbi.binData) == cbi.headersSize + cbi.recordsNum * cbi.recordSize

    def test_get_record_by_time_finds_the_same_record_as_by_index(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        expected_time, expected_line = cbi.getRecordByIndex(1)
        time, line = cbi.getRecordByTime(expected_time)
        assert time == expected_time
        assert line == expected_line

    def test_bytetostr_decodes_and_strips_trailing_nulls(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        assert cbi._byteToStr(b"AB\x00\x00") == "AB"

    def test_floatconvert_applies_the_sign_bit(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        assert cbi._floatConvert(0x00, 100) == pytest.approx(100.0)
        assert cbi._floatConvert(0x80, 100) == pytest.approx(-100.0)

    def test_floatconvert_applies_the_scale_factor_bits(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        assert cbi._floatConvert(0x60, 100) == pytest.approx(0.1)   # /1000
        assert cbi._floatConvert(0x40, 100) == pytest.approx(1.0)   # /100
        assert cbi._floatConvert(0x20, 100) == pytest.approx(10.0)  # /10
        assert cbi._floatConvert(0x00, 100) == pytest.approx(100.0)  # /1

    def test_newfloatconvert_caches_by_key(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        first = cbi._newfloatConvert(100)
        assert cbi._lut[100] == first
        second = cbi._newfloatConvert(100)
        assert second == first


@pytest.mark.unit
class TestNewFloatConvertNanSentinelIsInconsistent:
    """B107: the NaN sentinel branch (key == 65183) sets self._lut[key] =
    float('nan') but falls off the end of the if-block without returning
    it -- the implicit `return None` fires instead. A second call for the
    same key returns nan correctly, because by then the lookup at the top
    of the function (`return self._lut[key]`) short-circuits the branch
    entirely. The two calls disagree about what this key means."""

    @pytest.mark.xfail(
        strict=True,
        reason="B107: the nan-sentinel branch of _newfloatConvert does not "
               "return its own cached value -- the first call for key "
               "65183 returns None while every subsequent (cached) call "
               "returns nan. See the consolidated findings issue.",
    )
    def test_the_first_call_should_also_return_nan(self, sonic_file):
        import math

        cbi = CampbellBinaryInterface(file=str(sonic_file))
        first = cbi._newfloatConvert(65183)
        assert first is not None and math.isnan(first)

    def test_the_first_call_currently_returns_none(self, sonic_file):
        """Characterisation of B107."""
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        assert cbi._newfloatConvert(65183) is None

    def test_a_second_call_for_the_same_key_returns_nan_from_the_cache(self, sonic_file):
        """Characterisation of B107: same key, different call, different
        answer -- because the second call hits the cache lookup instead of
        re-running the broken branch."""
        import math

        cbi = CampbellBinaryInterface(file=str(sonic_file))
        cbi._newfloatConvert(65183)
        second = cbi._newfloatConvert(65183)
        assert second is not None and math.isnan(second)
