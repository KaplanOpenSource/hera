"""hera/measurements/experiment/parsers.py: everything the two existing
files for this module leave uncovered.

``test_experiment_parsers.py`` covers the pydoc factory convention, the
TOA5 stub (B44) and the fact that the module is orphaned;
``test_experiment_parsers_logic.py`` covers the two broken entry points
(B89 ``getData``'s off-by-two column slice, B90 the wrong
``DataFrame.from_dict`` orientation).  This file covers the rest:

* ``Parser_CampbellBinary.parse`` -- the only method in this file that is
  *not* a byte-for-byte copy of
  ``meteorology/highfreqdata/parsers/CampbellBinary.py``.  That copy's
  ``parse`` simply returns the DataFrame; this one additionally builds a
  nested station/instrument/height metadata dict and wraps the frame in a
  single-partition Dask DataFrame.  Because ``getPandasFromFile`` itself
  is broken (B89), the dispatch and the metadata construction are
  exercised against a stubbed reader -- the same technique
  ``test_meteorology_campbellbinary.py`` uses.
* ``Parser_CampbellBinary.getPandasFromDir`` -- globbing and concatenation.
* ``CampbellBinaryInterface``: ``binData``, ``firstTime``, ``lastTime``,
  ``getRecordByTime``, ``getTimeByRecordIndex``, ``getRecordIndexByTime``,
  ``_byteToStr``, ``_floatConvert``, ``_newfloatConvert``, the untaken
  branches of ``_getFormat``, and the multi-height branches of
  ``_getColumnNames`` / ``_getColumnIndexes`` / ``heights``.
* ``Parser_OldStyleMetaDataParquet.getExperimentDict`` / ``_getLists`` --
  reached only through a raising ``parse`` in the existing file.  They are
  exercised here on the one station-locations orientation the code does
  accept (attribute-major: ``{"lat": {station: value}, ...}``), which is
  the flip side of B90.

Fixtures are hand-built TOB1 files (five ASCII header lines, then packed
binary records), following ``test_meteorology_campbellbinary.py``; no real
capture was available.  Nothing here needs a Project, a toolkit or a
database, so no datalayer fixtures are used.

Bugs pinned here (all in the binary reader, and all present in the
highfreqdata copy of the same code too, since the two files are
independent duplicates):

* B115: ``getRecordByIndex`` turns the record's integer SECONDS field
  into a timestamp with ``pandas.Timedelta(days=lastSec / 86400.0)``.
  The float division does not round-trip: a record stamped SECONDS=11
  decodes to ``1990-01-01 00:00:10.999999999``.  77 of the first 600
  integer second values are affected, so roughly one record in eight is
  off by a nanosecond.  ``getRecordIndexByTime`` compares timestamps for
  *exact* equality, so asking for the record at the second it is actually
  stamped with raises ``IndexError`` for those records.
  ``pandas.Timedelta(seconds=lastSec)`` is exact.
* B116: ``_getDataFromStream``'s conversion loop is
  ``for i in range(3, len(retval))``, but ``retval[0]`` and ``retval[1]``
  are the two timestamp fields, so the first *data* field is ``retval[2]``
  and is skipped.  The first data column of every record is therefore
  returned unconverted -- a raw FP2 ``uint16`` instead of a float, or raw
  ``bytes`` instead of a string.  Two identical raw FP2 values decode to
  two different numbers depending only on which column they sit in.
* B117: ``_newfloatConvert`` calls
  ``self._floatConvert(int(key % 256), key / 256)``.  ``_floatConvert``
  documents its second argument as the low byte, but ``/`` is true
  division, so the value passed is ``lowbyte + hbyte/256`` -- the high
  byte leaks back in as a fractional part and every FP2 value comes out
  with a spurious offset of up to ``1/factor``.  ``key // 256`` is the
  low byte.
* B118: the ``key == 65183`` NaN-sentinel branch of
  ``_newfloatConvert`` stores ``float('nan')`` in the cache but returns
  ``None`` instead of it, so the first call for that key disagrees with
  every later (cached) call.  This is the same defect already pinned as
  B107 against the highfreqdata copy; it is repeated here verbatim, and
  fixing one file will not fix the other.

Deliberately not covered here: ``Parser_TOA5`` (a stub, already pinned as
B44), ``getData``'s column slicing and ``parse``'s metadata orientation
(B89/B90, already pinned), and ``_chunkSize`` / ``_dataContent`` /
``_basenum`` / ``byteSize``, which are assigned in ``__init__`` and then
never read by any code path in the module.
"""
import math
import os
import struct

import pandas
import pytest

from hera.measurements.experiment.parsers import (
    CampbellBinaryInterface,
    Parser_CampbellBinary,
    Parser_OldStyleMetaDataParquet,
)


def _write_tob1(path, header1, header4, fmt, records):
    """Write a minimal TOB1 file: five header lines, then packed records."""
    header0 = (
        '"TOB1","MyStation","CR3000","12345","CR3000.Std.32",'
        '"CPU:program.CR3","12345","MyInstrument"\r\n'
    )
    dataColumns = len(records[0]) - 2
    header2 = '"TS",' + ",".join(['""'] * dataColumns) + "\r\n"
    header3 = '"",' + ",".join(['"Smp"'] * dataColumns) + "\r\n"
    header = (header0 + header1 + header2 + header3 + header4).encode("ascii")
    body = b"".join(struct.pack(fmt, *record) for record in records)
    path.write_bytes(header + body)


@pytest.fixture()
def sonic_file(tmp_path):
    """Four sonic records, stamped at seconds 10, 11, 12 and 13 since 1990."""
    path = tmp_path / "sonic.dat"
    _write_tob1(
        path,
        header1='"TIMESTAMP","U_1","V_1","W_1","T_1"\r\n',
        header4='"ULONG","ULONG","IEEE4","IEEE4","IEEE4","IEEE4"\r\n',
        fmt="<IIffff",
        records=[
            (10, 0, 1.0, 0.5, 0.1, 20.0),
            (11, 0, 1.1, 0.6, 0.2, 20.1),
            (12, 0, 1.2, 0.7, 0.3, 20.2),
            (13, 0, 1.3, 0.8, 0.4, 20.3),
        ],
    )
    return path


@pytest.fixture()
def fp2_file(tmp_path):
    """One record whose two FP2 data fields hold the identical raw value.

    ``0x001f`` little-endian is the byte pair ``1f 00``; FP2 reads the
    first byte as sign + decimal locator + high mantissa (0x1f: positive,
    no decimal shift, mantissa high 31) and the second as the low
    mantissa, so both fields mean 31 * 256 + 0 = 7936.
    """
    path = tmp_path / "fp2.dat"
    _write_tob1(
        path,
        header1='"TIMESTAMP","U_1","V_1"\r\n',
        header4='"ULONG","ULONG","FP2","FP2"\r\n',
        fmt="<IIHH",
        records=[(10, 0, 0x001F, 0x001F)],
    )
    return path


@pytest.fixture()
def ascii_file(tmp_path):
    """One record whose two data fields are both four-byte ASCII strings."""
    path = tmp_path / "ascii.dat"
    _write_tob1(
        path,
        header1='"TIMESTAMP","U_1","V_1"\r\n',
        header4='"ULONG","ULONG","ASCII(4)","ASCII(4)"\r\n',
        fmt="<II4s4s",
        records=[(10, 0, b"AB", b"AB")],
    )
    return path


def _stub_reader(monkeypatch, methodName, frame):
    """Replace one of the parser's file readers, recording its calls."""
    calls = []

    def _fake(self, path, fromTime, toTime):
        calls.append((path, fromTime, toTime))
        return frame

    monkeypatch.setattr(Parser_CampbellBinary, methodName, _fake)
    return calls


PARSED_FRAME = pandas.DataFrame(
    {
        "u": [1.0, 2.0, 3.0],
        "height": [10, 10, 20],
        "station": ["STN1", "STN1", "STN2"],
        "instrument": ["sonic", "sonic", "sonic"],
    }
)


@pytest.mark.unit
class TestParserCampbellBinaryParse:
    """``parse`` -- the one method unique to this copy of the parser."""

    def test_a_file_path_is_routed_to_the_single_file_reader(self, sonic_file, monkeypatch):
        calls = _stub_reader(monkeypatch, "getPandasFromFile", PARSED_FRAME)

        Parser_CampbellBinary().parse(str(sonic_file), fromTime="A", toTime="B")

        assert calls == [(str(sonic_file), "A", "B")]

    def test_a_directory_path_is_routed_to_the_directory_reader(self, tmp_path, monkeypatch):
        calls = _stub_reader(monkeypatch, "getPandasFromDir", PARSED_FRAME)

        Parser_CampbellBinary().parse(str(tmp_path))

        assert calls == [(str(tmp_path), None, None)]

    def test_the_data_is_returned_as_a_single_partition_dask_frame(self, sonic_file, monkeypatch):
        _stub_reader(monkeypatch, "getPandasFromFile", PARSED_FRAME)

        loaded, _ = Parser_CampbellBinary().parse(str(sonic_file))

        assert loaded.npartitions == 1
        computed = loaded.compute()
        assert list(computed.columns) == list(PARSED_FRAME.columns)
        assert list(computed["u"]) == [1.0, 2.0, 3.0]

    def test_the_metadata_is_nested_by_station_then_instrument_then_height(
        self, sonic_file, monkeypatch
    ):
        _stub_reader(monkeypatch, "getPandasFromFile", PARSED_FRAME)

        _, metadata = Parser_CampbellBinary().parse(str(sonic_file))

        assert sorted(metadata) == ["STN1", "STN2"]
        assert sorted(metadata["STN1"]["sonic"]) == [10]
        assert sorted(metadata["STN2"]["sonic"]) == [20]

    def test_each_leaf_describes_its_own_station_instrument_and_height(
        self, sonic_file, monkeypatch
    ):
        _stub_reader(monkeypatch, "getPandasFromFile", PARSED_FRAME)

        _, metadata = Parser_CampbellBinary().parse(str(sonic_file))

        assert metadata["STN2"]["sonic"][20] == dict(
            station="STN2", instrument="sonic", height=20
        )

    def test_caller_supplied_metadata_is_merged_into_every_leaf(self, sonic_file, monkeypatch):
        _stub_reader(monkeypatch, "getPandasFromFile", PARSED_FRAME)

        _, metadata = Parser_CampbellBinary().parse(str(sonic_file), campaign="C1")

        for station, instruments in metadata.items():
            for instrument, heights in instruments.items():
                for leaf in heights.values():
                    assert leaf["campaign"] == "C1"

    def test_heights_are_keyed_by_int_even_when_the_column_is_float(
        self, sonic_file, monkeypatch
    ):
        frame = PARSED_FRAME.assign(height=[10.0, 10.0, 20.0])
        _stub_reader(monkeypatch, "getPandasFromFile", frame)

        _, metadata = Parser_CampbellBinary().parse(str(sonic_file))

        assert list(metadata["STN1"]["sonic"]) == [10]
        assert metadata["STN1"]["sonic"][10]["height"] == 10


@pytest.mark.unit
class TestGetPandasFromDir:
    """Globbing and concatenation, isolated from B89 with a stubbed reader."""

    def test_every_dat_file_in_the_directory_is_read(self, tmp_path, monkeypatch):
        for name in ("a.dat", "b.dat"):
            (tmp_path / name).write_bytes(b"")
        calls = []

        def _fake(self, path, fromTime, toTime):
            calls.append(path)
            return pandas.DataFrame({"v": [len(calls)]})

        monkeypatch.setattr(Parser_CampbellBinary, "getPandasFromFile", _fake)

        result = Parser_CampbellBinary().getPandasFromDir(
            str(tmp_path), fromTime=None, toTime=None
        )

        assert sorted(calls) == [str(tmp_path / "a.dat"), str(tmp_path / "b.dat")]
        assert len(result) == 2

    def test_files_that_are_not_dat_are_skipped(self, tmp_path, monkeypatch):
        (tmp_path / "a.dat").write_bytes(b"")
        (tmp_path / "notes.txt").write_bytes(b"")
        (tmp_path / "b.dat.bak").write_bytes(b"")
        calls = _stub_reader(
            monkeypatch, "getPandasFromFile", pandas.DataFrame({"v": [1]})
        )

        Parser_CampbellBinary().getPandasFromDir(str(tmp_path), fromTime=None, toTime=None)

        assert [path for path, _, _ in calls] == [str(tmp_path / "a.dat")]

    def test_a_directory_with_no_dat_files_is_refused(self, tmp_path):
        """Nothing to concatenate, so pandas -- not the parser -- objects."""
        with pytest.raises(ValueError, match="No objects to concatenate"):
            Parser_CampbellBinary().getPandasFromDir(
                str(tmp_path), fromTime=None, toTime=None
            )


@pytest.mark.unit
class TestFormatDescriptorParsing:
    """The fifth header line, which declares one Campbell type per field."""

    def test_every_campbell_type_maps_to_its_struct_code(self, tmp_path):
        path = tmp_path / "types.dat"
        _write_tob1(
            path,
            header1='"TIMESTAMP","U_1"\r\n',
            header4='"ULONG","ULONG","FP2","IEEE8","USHORT","LONG","BOOL","ASCII(4)"\r\n',
            fmt="<IIHdHl?4s",
            records=[(10, 0, 0x001F, 1.5, 7, -3, True, b"ab")],
        )
        cbi = CampbellBinaryInterface(file=str(path))

        assert cbi.format == "<IIHdHl?4s"
        assert cbi.recordSize == struct.calcsize("<IIHdHl?4s")

    def test_the_raw_type_names_are_kept_alongside_the_struct_format(self, tmp_path):
        path = tmp_path / "types.dat"
        _write_tob1(
            path,
            header1='"TIMESTAMP","U_1"\r\n',
            header4='"ULONG","ULONG","IEEE8","BOOL"\r\n',
            fmt="<IId?",
            records=[(10, 0, 1.5, False)],
        )
        cbi = CampbellBinaryInterface(file=str(path))

        assert cbi.rawFormat == ["ULONG", "ULONG", "IEEE8", "BOOL"]

    def test_an_ascii_field_takes_the_width_declared_in_its_name(self, tmp_path):
        path = tmp_path / "ascii12.dat"
        _write_tob1(
            path,
            header1='"TIMESTAMP","U_1"\r\n',
            header4='"ULONG","ULONG","ASCII(12)"\r\n',
            fmt="<II12s",
            records=[(10, 0, b"hello")],
        )
        cbi = CampbellBinaryInterface(file=str(path))

        assert cbi.format == "<II12s"

    def test_an_unrecognised_type_name_is_rejected(self, tmp_path):
        path = tmp_path / "unknown.dat"
        _write_tob1(
            path,
            header1='"TIMESTAMP","U_1"\r\n',
            header4='"ULONG","ULONG","NOSUCHTYPE"\r\n',
            fmt="<III",
            records=[(10, 0, 1)],
        )
        cbi = CampbellBinaryInterface(file=str(path))

        with pytest.raises(Exception, match="Unknown NOSUCHTYPE Format"):
            cbi.format

    def test_a_format_line_holding_a_single_field_is_rejected(self, tmp_path):
        """The comma check exists to catch a header line that is not a
        format descriptor at all."""
        path = tmp_path / "nocomma.dat"
        _write_tob1(
            path,
            header1='"TIMESTAMP"\r\n',
            header4='"ULONG"\r\n',
            fmt="<II",
            records=[(10, 0)],
        )
        cbi = CampbellBinaryInterface(file=str(path))

        with pytest.raises(Exception, match="Missing Format Descriptor"):
            cbi.recordSize


@pytest.mark.unit
class TestMultipleHeightColumnLayouts:
    """A single mast file can carry three instruments, one per height."""

    def test_one_sonic_group_is_a_single_layer_at_ten_metres(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        assert cbi.columnsNames == [["u", "v", "w", "T"]]
        assert cbi.heights == [10]

    def test_three_sonic_groups_each_get_their_own_four_columns(self, tmp_path):
        path = tmp_path / "sonic3.dat"
        _write_tob1(
            path,
            header1=(
                '"TIMESTAMP","U_1","V_1","W_1","T_1","U_2","V_2","W_2","T_2",'
                '"U_3","V_3","W_3","T_3"\r\n'
            ),
            header4='"ULONG","ULONG"' + ',"IEEE4"' * 12 + "\r\n",
            fmt="<II" + "f" * 12,
            records=[tuple([10, 0] + [float(i) for i in range(12)])],
        )
        cbi = CampbellBinaryInterface(file=str(path))

        assert cbi.columnsNames == [["u", "v", "w", "T"]] * 3
        assert cbi.heights == [6, 11, 16]

    def test_three_sonic_groups_slice_four_fields_each(self, tmp_path):
        path = tmp_path / "sonic3.dat"
        _write_tob1(
            path,
            header1=(
                '"TIMESTAMP","U_1","V_1","W_1","T_1","U_2","V_2","W_2","T_2",'
                '"U_3","V_3","W_3","T_3"\r\n'
            ),
            header4='"ULONG","ULONG"' + ',"IEEE4"' * 12 + "\r\n",
            fmt="<II" + "f" * 12,
            records=[tuple([10, 0] + [float(i) for i in range(12)])],
        )
        cbi = CampbellBinaryInterface(file=str(path))

        assert cbi.columnsIndexes == [[1, 5], [5, 9], [9, 13]]

    def test_three_thermocouple_groups_put_the_humidity_pair_on_the_last_one(
        self, tmp_path
    ):
        """TcT is one column per height; the TRH/RH pair is declared once,
        and the parser attaches it to the topmost group."""
        path = tmp_path / "tct3.dat"
        _write_tob1(
            path,
            header1='"TIMESTAMP","TC_T(1)","TC_T(2)","TC_T(3)","TRH","RH"\r\n',
            header4='"ULONG","ULONG"' + ',"IEEE4"' * 5 + "\r\n",
            fmt="<II" + "f" * 5,
            records=[(10, 0, 20.0, 20.1, 20.2, 50.0, 45.0)],
        )
        cbi = CampbellBinaryInterface(file=str(path))

        assert cbi.columnsNames == [["TcT"], ["TcT"], ["TcT", "TRH", "RH"]]
        assert cbi.heights == [6, 11, 16]

    def test_three_thermocouple_groups_slice_one_field_each_plus_the_pair(
        self, tmp_path
    ):
        path = tmp_path / "tct3.dat"
        _write_tob1(
            path,
            header1='"TIMESTAMP","TC_T(1)","TC_T(2)","TC_T(3)","TRH","RH"\r\n',
            header4='"ULONG","ULONG"' + ',"IEEE4"' * 5 + "\r\n",
            fmt="<II" + "f" * 5,
            records=[(10, 0, 20.0, 20.1, 20.2, 50.0, 45.0)],
        )
        cbi = CampbellBinaryInterface(file=str(path))

        assert cbi.columnsIndexes == [[1, 2], [2, 3], [3, 6]]


@pytest.mark.unit
class TestRecordNavigation:
    def test_bin_data_is_the_whole_file(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        assert cbi.binData == sonic_file.read_bytes()
        assert len(cbi.binData) == cbi.headersSize + cbi.recordsNum * cbi.recordSize

    def test_first_time_and_last_time_bracket_the_records(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        assert cbi.firstTime == cbi.getRecordByIndex(0)[0]
        assert cbi.lastTime == cbi.getRecordByIndex(cbi.recordsNum - 1)[0]
        assert cbi.firstTime < cbi.lastTime

    def test_first_time_and_last_time_are_computed_once_and_cached(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        assert cbi.firstTime is cbi.firstTime
        assert cbi.lastTime is cbi.lastTime

    def test_get_time_by_record_index_is_the_timestamp_half_of_the_record(
        self, sonic_file
    ):
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        for index in range(cbi.recordsNum):
            assert cbi.getTimeByRecordIndex(index) == cbi.getRecordByIndex(index)[0]

    def test_get_record_by_time_returns_the_record_stamped_with_that_time(
        self, sonic_file
    ):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        expectedTime, expectedLine = cbi.getRecordByIndex(2)

        time, line = cbi.getRecordByTime(expectedTime)

        assert time == expectedTime
        assert line == expectedLine

    def test_get_record_index_by_time_locates_every_record(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        for index in range(cbi.recordsNum):
            time = cbi.getTimeByRecordIndex(index)
            assert cbi.getRecordIndexByTime(time) == index

    def test_a_time_with_no_record_is_refused(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        with pytest.raises(IndexError, match="There is no record at"):
            cbi.getRecordIndexByTime(pandas.Timestamp(2020, 1, 1))

    def test_the_time_window_selects_only_the_records_inside_it(self, sonic_file):
        """``getData``'s window is inclusive at both ends; only the
        timestamps are checked here, because the values it returns are
        broken by B89."""
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        fromTime = cbi.getTimeByRecordIndex(1)
        toTime = cbi.getTimeByRecordIndex(2)

        ts, _, _ = Parser_CampbellBinary().getData(
            str(sonic_file), fromTime=str(fromTime), toTime=str(toTime)
        )

        assert ts == [fromTime, toTime]


@pytest.mark.unit
class TestRecordTimestampsLoseTheirLastNanosecond:
    """B115: seconds are converted through a float number of days."""

    @pytest.mark.xfail(
        strict=True,
        reason="B115: getRecordByIndex builds the timestamp with "
               "pandas.Timedelta(days=lastSec / 86400.0); the float "
               "division does not round-trip, so a record stamped with "
               "SECONDS=11 decodes to 1990-01-01 00:00:10.999999999 "
               "instead of ...:11. pandas.Timedelta(seconds=lastSec) is "
               "exact. See the consolidated findings issue.",
    )
    def test_an_integer_seconds_field_decodes_to_that_exact_second(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        assert cbi.getRecordByIndex(1)[0] == pandas.Timestamp(1990, 1, 1) + \
            pandas.Timedelta(seconds=11)

    def test_the_eleventh_second_currently_decodes_a_nanosecond_early(self, sonic_file):
        """Characterisation of B115."""
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        assert cbi.getRecordByIndex(1)[0] == pandas.Timestamp(
            "1990-01-01 00:00:10.999999999"
        )

    def test_neighbouring_records_are_not_affected(self, sonic_file):
        """Characterisation of B115: the error depends on the value of
        the seconds field, so it hits some records and not others -- which
        is what makes it hard to notice."""
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        for index, second in ((0, 10), (2, 12), (3, 13)):
            assert cbi.getRecordByIndex(index)[0] == pandas.Timestamp(1990, 1, 1) + \
                pandas.Timedelta(seconds=second)

    @pytest.mark.xfail(
        strict=True,
        reason="B115: getRecordIndexByTime matches timestamps exactly, so "
               "a record whose seconds field decoded inexactly cannot be "
               "found by the time it is actually stamped with. "
               "See the consolidated findings issue.",
    )
    def test_a_record_can_be_looked_up_by_its_own_wall_clock_second(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        assert cbi.getRecordIndexByTime(pandas.Timestamp("1990-01-01 00:00:11")) == 1

    def test_looking_that_record_up_currently_raises(self, sonic_file):
        """Characterisation of B115."""
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        with pytest.raises(IndexError, match="There is no record at"):
            cbi.getRecordIndexByTime(pandas.Timestamp("1990-01-01 00:00:11"))


@pytest.mark.unit
class TestFirstDataColumnIsNeverConverted:
    """B116: the conversion loop starts one field too late."""

    @pytest.mark.xfail(
        strict=True,
        reason="B116: _getDataFromStream converts retval[3:] but the data "
               "starts at retval[2] (retval[0]/[1] are SECONDS and "
               "NANOSECONDS), so the first data field of every record is "
               "returned unconverted. Two identical raw FP2 values decode "
               "to two different numbers. See the consolidated findings "
               "issue.",
    )
    def test_two_identical_raw_fp2_fields_decode_to_the_same_value(self, fp2_file):
        cbi = CampbellBinaryInterface(file=str(fp2_file))

        _, line = cbi.getRecordByIndex(0)

        assert line[0] == line[1]

    def test_the_first_fp2_field_is_still_a_raw_integer(self, fp2_file):
        """Characterisation of B116: 0x001f comes back as the integer 31
        rather than the FP2 value 7936 it encodes."""
        cbi = CampbellBinaryInterface(file=str(fp2_file))

        _, line = cbi.getRecordByIndex(0)

        assert line[0] == 0x001F
        assert isinstance(line[0], int)
        assert line[1] != 0x001F

    def test_the_first_ascii_field_is_still_raw_bytes(self, ascii_file):
        """Characterisation of B116: the same off-by-one skips the
        bytes-to-string decode, not just the FP2 conversion."""
        cbi = CampbellBinaryInterface(file=str(ascii_file))

        _, line = cbi.getRecordByIndex(0)

        assert line[0] == b"AB\x00\x00"
        assert line[1] == "AB"


@pytest.mark.unit
class TestFp2ConversionHelpers:
    def test_the_sign_bit_flips_the_value(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        assert cbi._floatConvert(0x00, 100) == pytest.approx(100.0)
        assert cbi._floatConvert(0x80, 100) == pytest.approx(-100.0)

    def test_the_decimal_locator_bits_choose_the_divisor(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        assert cbi._floatConvert(0x00, 100) == pytest.approx(100.0)
        assert cbi._floatConvert(0x20, 100) == pytest.approx(10.0)
        assert cbi._floatConvert(0x40, 100) == pytest.approx(1.0)
        assert cbi._floatConvert(0x60, 100) == pytest.approx(0.1)

    def test_the_mantissa_spans_both_bytes(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        assert cbi._floatConvert(0x1F, 255) == pytest.approx(31 * 256 + 255)

    def test_a_converted_value_is_cached_under_its_raw_key(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        first = cbi._newfloatConvert(0x1234)

        assert cbi._lut[0x1234] == first
        assert cbi._newfloatConvert(0x1234) == first

    def test_a_byte_sequence_decodes_without_its_padding_nulls(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        assert cbi._byteToStr(b"AB\x00\x00") == "AB"


@pytest.mark.unit
class TestFp2LowByteIsDividedInsteadOfShifted:
    """B117: ``key / 256`` where ``key // 256`` is the low byte."""

    @pytest.mark.xfail(
        strict=True,
        reason="B117: _newfloatConvert passes key / 256 as _floatConvert's "
               "lowbyte, but true division leaves the high byte in as a "
               "fractional part, so every FP2 value is offset by "
               "(key mod 256) / 256 / factor. key // 256 is the low byte. "
               "See the consolidated findings issue.",
    )
    def test_a_raw_fp2_key_converts_to_the_value_its_two_bytes_encode(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        assert cbi._newfloatConvert(0x001F) == cbi._floatConvert(0x1F, 0x00)

    def test_the_high_byte_currently_leaks_into_the_low_byte(self, sonic_file):
        """Characterisation of B117: 0x001f encodes 7936, but the low
        byte is passed as 31/256 rather than 0."""
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        assert cbi._newfloatConvert(0x001F) == pytest.approx(7936 + 31 / 256)

    def test_the_offset_shrinks_with_the_decimal_locator(self, sonic_file):
        """Characterisation of B117: the spurious term is divided by the
        same factor as the value, so it is 1/1000 of a count for a
        three-decimal reading and a whole count for an integer one."""
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        assert cbi._newfloatConvert(0x0061) == pytest.approx(
            (1 * 256 + 0x61 / 256) / 1000.0
        )


@pytest.mark.unit
class TestNanSentinelDoesNotReturnItsOwnCachedValue:
    """B118: the same defect as B107, duplicated in this module."""

    @pytest.mark.xfail(
        strict=True,
        reason="B118: the key == 65183 branch of _newfloatConvert stores "
               "float('nan') in the cache but returns None instead of it, "
               "so the first call disagrees with every cached call that "
               "follows. Identical to B107 in "
               "meteorology/highfreqdata/parsers/CampbellBinary.py, which "
               "is an independent copy of this code. See the consolidated "
               "findings issue.",
    )
    def test_the_sentinel_converts_to_nan_on_the_first_call_too(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        value = cbi._newfloatConvert(65183)

        assert value is not None and math.isnan(value)

    def test_the_first_call_currently_returns_none(self, sonic_file):
        """Characterisation of B118."""
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        assert cbi._newfloatConvert(65183) is None

    def test_the_cache_holds_nan_all_the_same(self, sonic_file):
        """Characterisation of B118: the value is computed and stored,
        just not handed back, so the second call answers differently."""
        cbi = CampbellBinaryInterface(file=str(sonic_file))

        cbi._newfloatConvert(65183)

        assert math.isnan(cbi._lut[65183])
        assert math.isnan(cbi._newfloatConvert(65183))


CAMPAIGN_METADATA = {
    "stations": {
        "STN1": {
            "instruments": {"sonic": ["10m", "20m"], "trh": ["3m"]},
            "averagedHeight": 5.0,
            "buildingHeight": 12.0,
        },
        "STN2": {
            "instruments": {"sonic": ["10m"]},
            "averagedHeight": 7.0,
            "buildingHeight": 15.0,
        },
    }
}

# The only orientation _getLists can consume: attribute-major, so that
# DataFrame.from_dict's default orient="columns" puts the station names in
# the index and the attribute names in the columns.  See B90.
STATION_LOCATIONS = {
    "lat": {"STN1": 32.1, "STN2": 32.2},
    "lon": {"STN1": 34.8, "STN2": 34.9},
    "station_code": {"STN1": "S1", "STN2": "S2"},
}


def _description(pathToData):
    return {
        "campaignName": "MYCAMPAIGN",
        "stationLocations": STATION_LOCATIONS,
        "toolkitHandler": "someHandler",
        "pathToData": pathToData,
    }


@pytest.mark.unit
class TestOldStyleMetaDataParquetExperimentDict:
    """``getExperimentDict`` / ``_getLists`` on input they can consume."""

    def test_it_reports_trials_devices_and_assets(self):
        result = Parser_OldStyleMetaDataParquet().getExperimentDict(
            CAMPAIGN_METADATA, _description("/campaign")
        )

        assert sorted(result) == ["assets", "devices", "trials"]

    def test_one_device_is_emitted_per_instrument_height(self):
        result = Parser_OldStyleMetaDataParquet().getExperimentDict(
            CAMPAIGN_METADATA, _description("/campaign")
        )

        assert [device["deviceName"] for device in result["devices"]] == [
            "sonic_1",
            "sonic_2",
            "trh_1",
            "sonic_3",
        ]

    def test_the_device_type_is_the_instrument_name(self):
        result = Parser_OldStyleMetaDataParquet().getExperimentDict(
            CAMPAIGN_METADATA, _description("/campaign")
        )
        byName = {device["deviceName"]: device for device in result["devices"]}

        assert byName["trh_1"]["deviceType"] == "trh"
        assert byName["sonic_3"]["deviceType"] == "sonic"

    def test_the_device_data_path_is_station_instrument_height_under_the_root(self):
        result = Parser_OldStyleMetaDataParquet().getExperimentDict(
            CAMPAIGN_METADATA, _description("/campaign")
        )
        byName = {device["deviceName"]: device for device in result["devices"]}

        assert byName["sonic_2"]["deviceDataPath"] == os.path.join(
            "/campaign", "STN1", "sonic", "20m"
        )
        assert byName["sonic_3"]["deviceDataPath"] == os.path.join(
            "/campaign", "STN2", "sonic", "10m"
        )

    def test_each_device_carries_its_station_coordinates_and_height(self):
        result = Parser_OldStyleMetaDataParquet().getExperimentDict(
            CAMPAIGN_METADATA, _description("/campaign")
        )
        byName = {device["deviceName"]: device for device in result["devices"]}

        design = byName["sonic_3"]["trials"]["measurement1"]["design"]
        assert design == dict(
            name="sonic_3",
            deviceLocation="STN2",
            deviceCoords=[32.2, 34.9],
            deviceHeight="10m",
        )

    def test_the_design_and_deploy_entries_describe_the_same_device(self):
        result = Parser_OldStyleMetaDataParquet().getExperimentDict(
            CAMPAIGN_METADATA, _description("/campaign")
        )
        trials = result["devices"][0]["trials"]["measurement1"]

        assert trials["design"] == trials["deploy"]

    def test_each_station_becomes_an_asset_with_its_building_properties(self):
        result = Parser_OldStyleMetaDataParquet().getExperimentDict(
            CAMPAIGN_METADATA, _description("/campaign")
        )

        assert sorted(result["assets"]) == ["STN1", "STN2"]
        assert result["assets"]["STN1"] == dict(
            name="STN1",
            codeName="S1",
            properties=dict(
                coords=[32.1, 34.8],
                avgAreaBuildingsHgt=5.0,
                buildingHgt=12.0,
            ),
            trials={
                "measurement1": dict(
                    entities=["sonic_1", "sonic_2", "trh_1"],
                    types=["sonic", "sonic", "trh"],
                )
            },
        )

    def test_the_asset_of_the_second_station_lists_only_its_own_devices(self):
        result = Parser_OldStyleMetaDataParquet().getExperimentDict(
            CAMPAIGN_METADATA, _description("/campaign")
        )

        assert result["assets"]["STN2"]["trials"]["measurement1"] == dict(
            entities=["sonic_3"], types=["sonic"]
        )

    def test_device_numbering_continues_across_stations(self):
        """Characterisation, not a documented guarantee: the instrument
        counter is campaign-wide, so STN2's only sonic is sonic_3 rather
        than sonic_1.  Device names are therefore not reusable as a
        per-station identifier."""
        result = Parser_OldStyleMetaDataParquet().getExperimentDict(
            CAMPAIGN_METADATA, _description("/campaign")
        )
        stn2 = [
            device
            for device in result["devices"]
            if device["trials"]["measurement1"]["design"]["deviceLocation"] == "STN2"
        ]

        assert [device["deviceName"] for device in stn2] == ["sonic_3"]

    def test_a_single_measurement_trial_is_always_declared(self):
        result = Parser_OldStyleMetaDataParquet().getExperimentDict(
            CAMPAIGN_METADATA, _description("/campaign")
        )

        assert result["trials"] == [
            dict(
                trialSet="measurement",
                trialName="measurement1",
                properties=dict(trialStnProps={}),
            )
        ]


@pytest.mark.unit
class TestOldStyleMetaDataParquetParse:
    """``parse`` reads the two JSON files and wraps the result."""

    @staticmethod
    def _writeCampaign(tmp_path):
        import json

        (tmp_path / "metadata.json").write_text(json.dumps(CAMPAIGN_METADATA))
        (tmp_path / "campaignDescription.json").write_text(
            json.dumps(_description(str(tmp_path)))
        )

    def test_everything_is_nested_under_the_campaign_name(self, tmp_path):
        self._writeCampaign(tmp_path)

        result = Parser_OldStyleMetaDataParquet().parse(str(tmp_path))

        assert list(result) == ["MYCAMPAIGN"]
        assert sorted(result["MYCAMPAIGN"]) == [
            "Stations",
            "assets",
            "devices",
            "toolkitHandler",
            "trials",
        ]

    def test_the_station_locations_are_passed_through_unchanged(self, tmp_path):
        self._writeCampaign(tmp_path)

        result = Parser_OldStyleMetaDataParquet().parse(str(tmp_path))

        assert result["MYCAMPAIGN"]["Stations"] == STATION_LOCATIONS

    def test_the_toolkit_handler_comes_from_the_campaign_description(self, tmp_path):
        self._writeCampaign(tmp_path)

        result = Parser_OldStyleMetaDataParquet().parse(str(tmp_path))

        assert result["MYCAMPAIGN"]["toolkitHandler"] == "someHandler"

    def test_the_devices_match_what_the_experiment_dict_builds(self, tmp_path):
        self._writeCampaign(tmp_path)

        result = Parser_OldStyleMetaDataParquet().parse(str(tmp_path))

        assert [
            device["deviceName"] for device in result["MYCAMPAIGN"]["devices"]
        ] == ["sonic_1", "sonic_2", "trh_1", "sonic_3"]

    def test_both_metadata_files_are_read_from_the_given_directory(self, tmp_path):
        """Neither file is looked up relative to the process CWD."""
        self._writeCampaign(tmp_path)
        elsewhere = tmp_path / "elsewhere"
        elsewhere.mkdir()

        with pytest.raises(FileNotFoundError, match="metadata.json"):
            Parser_OldStyleMetaDataParquet().parse(str(elsewhere))
