"""hera/measurements/experiment/parsers.py: the actual parsing logic that
test_experiment_parsers.py doesn't reach (that file only covers the pydoc
factory convention and documents the module as orphaned/dead code -- these
tests exercise what the parsers actually do when called).

* Parser_CampbellBinary / CampbellBinaryInterface -- a byte-for-byte
  duplicate of hera/measurements/meteorology/highfreqdata/parsers/CampbellBinary.py
  (same class names, same bodies), carrying the identical bug already
  pinned there as B82:

  B89: ``Parser_CampbellBinary.getData`` slices each record's data with
  ``line[cbi.columnsIndexes[i][0] : cbi.columnsIndexes[i][1]]``, but
  ``columnsIndexes`` is computed as if it indexed into the *raw* record --
  which starts with two leading timestamp fields (SECONDS, NANOSECONDS).
  ``line`` is ``_getDataFromStream``'s ``retval[2:]``, which has *already*
  dropped those two fields. The off-by-two slice drops the first data
  column and returns one column short of what the header declares --
  pandas then refuses to build the DataFrame. Because this file is an
  independent copy (not an import) of the other parser, fixing one does
  not fix the other -- and since the whole module is orphaned (see
  test_experiment_parsers.py::TestParserModuleIsOrphaned), nothing is
  currently exercising this path in production either way.

* Parser_OldStyleMetaDataParquet:

  B90: ``_getLists`` builds ``stationsData`` with
  ``pandas.DataFrame.from_dict(descriptionData['stationLocations'])`` --
  the default ``orient="columns"`` treats each station name as a *column*
  and its attribute names (``lat``/``lon``/``station_code``) as the *row
  index*. The rest of the method (``.reset_index().rename(columns=
  {"index": "station_name"})``, then ``stndata['lat'].item()``) only makes
  sense if stations are rows and attributes are columns, i.e.
  ``orient="index"``. With the current orientation, ``reset_index()``
  produces a ``station_name`` column containing the strings
  ``"lat"``/``"lon"``/``"station_code"`` instead of station names, and
  ``stndata['lat']`` raises ``KeyError`` for any real multi-field station
  dict. This breaks ``_getLists``/``getExperimentDict``/``parse`` for any
  campaign with more than a single bare value per station.
"""
import json
import struct

import pandas
import pytest

from hera.measurements.experiment.parsers import (
    CampbellBinaryInterface,
    Parser_CampbellBinary,
    Parser_OldStyleMetaDataParquet,
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
        ],
    )
    return path


@pytest.mark.unit
class TestCampbellBinaryInterfaceStillWorks:
    def test_the_station_and_instrument_come_from_the_first_header_line(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        assert cbi.station == "MyStation"
        assert cbi.instrument == "MyInstrument"

    def test_records_num_counts_whole_records_after_the_header(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        assert cbi.recordsNum == 2

    def test_get_record_by_index_decodes_the_timestamp_from_1990_epoch(self, sonic_file):
        cbi = CampbellBinaryInterface(file=str(sonic_file))
        time, _ = cbi.getRecordByIndex(0)
        assert time == pandas.Timestamp(1990, 1, 1) + pandas.Timedelta(seconds=10)


@pytest.mark.unit
class TestParserCampbellBinaryGetDataIsBroken:
    @pytest.mark.xfail(
        strict=True,
        reason="B89: same off-by-two columnsIndexes bug as the meteorology "
               "highfreqdata CampbellBinary parser (B82), duplicated "
               "independently in this file. See the consolidated findings issue.",
    )
    def test_get_pandas_from_file_should_return_all_sonic_columns(self, sonic_file):
        df = Parser_CampbellBinary().getPandasFromFile(str(sonic_file), fromTime=None, toTime=None)
        assert list(df.columns[:4]) == ["u", "v", "w", "T"]

    def test_get_pandas_from_file_currently_raises(self, sonic_file):
        """Characterisation of B89."""
        with pytest.raises(ValueError, match="columns"):
            Parser_CampbellBinary().getPandasFromFile(str(sonic_file), fromTime=None, toTime=None)

    def test_get_data_drops_the_first_data_column(self, sonic_file):
        """Characterisation of B89, one level below the DataFrame
        construction that raises."""
        ts, cols, values = Parser_CampbellBinary().getData(str(sonic_file), fromTime=None, toTime=None)
        assert cols == [["u", "v", "w", "T"]]
        assert len(values[10][0]) == 3  # v, w, T -- u is missing


def _write_campaign(tmp_path, metadata, station_locations):
    description = {
        "campaignName": "MyCampaign",
        "stationLocations": station_locations,
        "toolkitHandler": "someHandler",
        "pathToData": str(tmp_path),
    }
    (tmp_path / "metadata.json").write_text(json.dumps(metadata))
    (tmp_path / "campaignDescription.json").write_text(json.dumps(description))


@pytest.mark.unit
class TestParserOldStyleMetaDataParquetGetListsIsBroken:
    """B90: stationsData is built with the wrong DataFrame.from_dict
    orientation -- see the module docstring."""

    METADATA = {
        "stations": {
            "STN1": {
                "instruments": {"sonic": ["10m"]},
                "averagedHeight": 5.0,
                "buildingHeight": 12.0,
            }
        }
    }
    STATION_LOCATIONS = {"STN1": {"lat": 32.1, "lon": 34.8, "station_code": "S1"}}

    @pytest.mark.xfail(
        strict=True,
        reason="B90: DataFrame.from_dict(stationLocations) uses the default "
               "orient='columns', putting attribute names in the row index "
               "instead of station names -- stndata['lat'] raises KeyError "
               "for any real (multi-field) station dict. "
               "See the consolidated findings issue.",
    )
    def test_parse_should_wrap_everything_under_the_campaign_name(self, tmp_path):
        _write_campaign(tmp_path, self.METADATA, self.STATION_LOCATIONS)
        Parser_OldStyleMetaDataParquet().parse(str(tmp_path))

    def test_parse_currently_raises_a_keyerror_on_lat(self, tmp_path):
        """Characterisation of B90."""
        _write_campaign(tmp_path, self.METADATA, self.STATION_LOCATIONS)
        with pytest.raises(KeyError, match="lat"):
            Parser_OldStyleMetaDataParquet().parse(str(tmp_path))

    def test_it_fails_the_same_way_for_any_station_shape(self, tmp_path):
        """Characterisation of B90: the wrong orientation means the
        'station_name' column never actually holds station names, so the
        query-then-lookup pattern fails regardless of how many attributes
        each station carries."""
        _write_campaign(tmp_path, self.METADATA, {"STN1": {"only_field": 1}})
        with pytest.raises(KeyError, match="lat"):
            Parser_OldStyleMetaDataParquet().parse(str(tmp_path))
