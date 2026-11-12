"""The DataHandler_* classes batch 2 didn't reach: time, csv, netcdf/xarray,
JSON-pandas, GeoJSON/GeoPackage geopandas, image, and plain-dict passthrough.
Each is a static save/load pair with real file I/O -- no DB involved.
"""
import geopandas
import numpy
import pandas
import pytest
import xarray
from shapely.geometry import Point

from hera.datalayer import datahandler


@pytest.mark.unit
class TestTimeHandler:
    def test_get_data_parses_the_resource_as_a_timestamp(self):
        result = datahandler.DataHandler_time.getData("2020-01-01T00:00:00")
        assert result == pandas.Timestamp("2020-01-01")

    def test_save_data_is_a_no_op_returning_an_empty_dict(self):
        assert datahandler.DataHandler_time.saveData(pandas.Timestamp.now(), "ignored") == {}


@pytest.mark.unit
class TestDictHandler:
    def test_get_data_returns_the_resource_unchanged(self):
        payload = {"a": 1, "b": [1, 2, 3]}
        assert datahandler.DataHandler_dict.getData(payload) == payload

    def test_save_data_is_a_no_op_returning_an_empty_dict(self):
        assert datahandler.DataHandler_dict.saveData({"a": 1}, "ignored") == {}


@pytest.mark.unit
class TestCsvPandasHandler:
    def test_round_trips_a_dataframe(self, tmp_path):
        target = tmp_path / "data.csv"
        frame = pandas.DataFrame({"a": [1, 2], "b": ["x", "y"]})
        datahandler.DataHandler_csv_pandas.saveData(frame, str(target), index=False)
        restored = datahandler.DataHandler_csv_pandas.getData(str(target))
        assert list(restored["a"]) == [1, 2]
        assert list(restored["b"]) == ["x", "y"]


@pytest.mark.unit
class TestJsonPandasHandler:
    def test_round_trips_a_dataframe(self, tmp_path):
        target = tmp_path / "data.json"
        frame = pandas.DataFrame({"a": [1, 2], "b": [3.0, 4.0]})
        meta = datahandler.DataHandler_JSON_pandas.saveData(frame, str(target))
        assert meta == {"usePandas": True}
        restored = datahandler.DataHandler_JSON_pandas.getData(str(target))
        assert list(restored["a"]) == [1, 2]

    def test_round_trips_a_series_with_the_series_flag(self, tmp_path):
        target = tmp_path / "series.json"
        series = pandas.Series([1, 2, 3], name="values")
        meta = datahandler.DataHandler_JSON_pandas.saveData(series, str(target))
        assert meta == {"pandasSeries": True, "usePandas": True}
        restored = datahandler.DataHandler_JSON_pandas.getData(str(target), pandasSeries=True)
        assert list(restored) == [1, 2, 3]


@pytest.mark.unit
class TestNetcdfXarrayHandler:
    def test_round_trips_a_dataset_and_its_attrs(self, tmp_path):
        """Note: attrs round-trip through JSONToConfiguration, which
        auto-hydrates any pint-parseable string value into a Quantity --
        so a plain non-unit string is used here to test the round trip
        itself, not that separate (and separately-tested) behaviour."""
        target = tmp_path / "data.nc"
        ds = xarray.Dataset({"temperature": ("x", [1.0, 2.0, 3.0])})
        ds.attrs["source"] = "unit test fixture"
        datahandler.DataHandler_netcdf_xarray.saveData(ds, str(target))
        restored = datahandler.DataHandler_netcdf_xarray.getData(str(target))
        assert list(restored["temperature"].values) == [1.0, 2.0, 3.0]
        assert restored.attrs["source"] == "unit test fixture"


@pytest.mark.unit
class TestJsonGeopandasHandlerIsBroken:
    """B87: saveData calls resource.to_json(fileName, **kwargs) -- but
    geopandas.GeoDataFrame.to_json's first positional parameter is `na`
    ('null'/'drop'/'keep'), and the method returns a JSON *string*, it
    never writes to a file at all. fileName lands in the na= slot and
    always fails validation."""

    @pytest.mark.xfail(
        strict=True,
        reason="B87: to_json's first positional argument is `na`, not a "
               "file path -- the method returns a string and never writes "
               "a file. fileName is passed as na= and rejected outright. "
               "See the consolidated findings issue.",
    )
    def test_save_data_should_write_a_geojson_file(self, tmp_path):
        target = tmp_path / "data.geojson"
        gdf = geopandas.GeoDataFrame({"name": ["A", "B"], "geometry": [Point(0, 0), Point(1, 1)]}, crs=4326)
        datahandler.DataHandler_JSON_geopandas.saveData(gdf, str(target))
        assert target.exists()

    def test_save_data_currently_raises_for_every_call(self, tmp_path):
        """Characterisation of B87."""
        target = tmp_path / "data.geojson"
        gdf = geopandas.GeoDataFrame({"name": ["A"], "geometry": [Point(0, 0)]}, crs=4326)
        with pytest.raises(ValueError, match="Unknown na method"):
            datahandler.DataHandler_JSON_geopandas.saveData(gdf, str(target))

    def test_get_data_itself_still_works_against_a_hand_written_file(self, tmp_path):
        """getData is independent of the broken saveData -- verified
        against a file written directly, the way it would exist if
        something else produced it correctly."""
        import json

        target = tmp_path / "data.geojson"
        target.write_text(json.dumps({
            "type": "FeatureCollection",
            "features": [
                {"type": "Feature", "properties": {"name": "A"},
                 "geometry": {"type": "Point", "coordinates": [0, 0]}},
            ],
        }))
        restored = datahandler.DataHandler_JSON_geopandas.getData(str(target))
        assert list(restored["name"]) == ["A"]

    def test_get_data_applies_the_crs_from_desc_if_given(self, tmp_path):
        import json

        target = tmp_path / "data.geojson"
        target.write_text(json.dumps({
            "type": "FeatureCollection",
            "features": [
                {"type": "Feature", "properties": {"name": "A"},
                 "geometry": {"type": "Point", "coordinates": [0, 0]}},
            ],
        }))
        restored = datahandler.DataHandler_JSON_geopandas.getData(str(target), desc={"crs": 4326})
        assert restored.crs.to_epsg() == 4326


@pytest.mark.unit
class TestGeopandasHandler:
    def test_round_trips_a_geodataframe_via_gpkg(self, tmp_path):
        target = tmp_path / "data.gpkg"
        gdf = geopandas.GeoDataFrame({"name": ["A", "B"], "geometry": [Point(0, 0), Point(1, 1)]}, crs=4326)
        meta = datahandler.DataHandler_geopandas.saveData(gdf, str(target))
        assert meta["crs"].to_epsg() == 4326
        restored = datahandler.DataHandler_geopandas.getData(str(target))
        assert list(restored["name"]) == ["A", "B"]


@pytest.mark.unit
class TestImageHandler:
    def test_round_trips_a_simple_rgb_array(self, tmp_path):
        target = tmp_path / "img.png"
        array = numpy.zeros((4, 4, 3), dtype=numpy.uint8)
        array[:, :, 0] = 255
        datahandler.DataHandler_image.saveData(array, str(target))
        restored = datahandler.DataHandler_image.getData(str(target))
        assert restored.shape[:2] == (4, 4)
        assert restored[0, 0, 0] > restored[0, 0, 1]
