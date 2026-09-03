"""hera/measurements/GIS/raster/topography.py -- the raster/elevation layer.

Covered here
------------
The module-level ``_open_raster`` and, on ``TopographyToolkit``:
``findElevationFile``, ``getPointElevation``, ``getPointListElevation``,
``convertPointsCRS``, ``create_xarray``, ``getElevation``,
``getElevationOfXarray``, ``createElevationSTL`` and ``getElevationSTL``.

How the inputs were produced
----------------------------
``_write_geotiff`` below synthesises every raster with ``rasterio``: a small
band, an explicit ``from_origin`` transform and an explicit EPSG code, written
into ``tmp_path``.  Two properties make the assertions exact rather than
golden:

* The tile is named as the code expects (``N32E035.hgt``) but its *transform*
  is chosen freely -- origin lon 35.0 / lat 33.0, pixel 0.0625 deg (= 1/16, an
  exact binary fraction).  Every latitude and longitude used below is
  therefore an exact binary fraction too, so the fractional pixel offsets the
  interpolator computes are exact and no test depends on floating-point
  rounding at a pixel boundary.
* The band is ``value(row, col) = 10 * row + 100 * col``, which is linear and
  hence reproduced *exactly* by the bilinear interpolation in
  ``getPointElevation``: the expected elevation is
  ``10 * rasterx + 100 * rastery``.  That gives an independent closed form to
  check against rather than a recorded number, and it also pins the
  orientation of every reshape (a transposed grid would not match).

``rasterio`` reads these ``.hgt``-named files with the GTiff driver (verified:
``src.driver == 'GTiff'``), so no real SRTM archive is needed.

Deliberately NOT tested, and why
--------------------------------
* ``topographyAnalysis.calculateStastics``: covered elsewhere and out of
  scope for this batch's target list.
* The ``osgeo``/gdal path: ``_gdal_available`` is False in the pinned
  environment (``No module named 'osgeo'``), and nothing in the file branches
  on it any more -- ``_open_raster`` is pure rasterio.
* FreeCAD: not reached.  ``getElevationSTL`` builds ASCII STL through
  ``stlFactory.rasterToSTL``, which is pure Python/numpy, so the STL is
  asserted on for real rather than against a MagicMock.
* ``stlFactory.rasterToSTL``'s own facet geometry: it has its own tests in
  ``test_gis_utils_stlfactory.py``; here only the toolkit's forwarding,
  orientation and shifting are checked.
"""
import math
import os

import numpy as np
import pytest

from hera.measurements.GIS.raster import topography as topography_module
from hera.measurements.GIS.raster.topography import _open_raster
from hera.measurements.GIS.utils import ITM, WSG84, convertCRS

# --------------------------------------------------------------------------
# raster synthesis
# --------------------------------------------------------------------------

_WEST = 35.0
_NORTH = 33.0
_RES = 0.0625  # 1/16 degree: exact in binary, so every offset below is exact.
_SIZE = 8


def _write_geotiff(path, array, west, north, resolution, crs="EPSG:4326"):
    """Write ``array`` as a one-band GeoTIFF at the given origin and return its path."""
    import rasterio
    from rasterio.transform import from_origin

    profile = dict(
        driver="GTiff",
        height=array.shape[0],
        width=array.shape[1],
        count=1,
        dtype=array.dtype,
        transform=from_origin(west, north, resolution, resolution),
    )
    if crs is not None:
        profile["crs"] = crs
    with rasterio.open(str(path), "w", **profile) as destination:
        destination.write(array, 1)
    return str(path)


def _ramp_band(size=_SIZE, dtype="float32"):
    """value(row, col) = 10 * row + 100 * col -- linear, so bilinear-exact."""
    rows, cols = np.mgrid[0:size, 0:size]
    return (10.0 * rows + 100.0 * cols).astype(dtype)


def _expected_elevation(lat, lon, west=_WEST, north=_NORTH, resolution=_RES):
    """The closed form the interpolator must reproduce on the ramp band."""
    rasterx = (lat - north) / -resolution
    rastery = (lon - west) / resolution
    return 10.0 * rasterx + 100.0 * rastery


@pytest.fixture()
def topo(unit_toolkit_factory):
    from hera import toolkitHome

    return unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)


@pytest.fixture()
def srtm_folder(tmp_path):
    """A folder holding the tile N32E035.hgt carrying the ramp band."""
    folder = tmp_path / "srtm"
    folder.mkdir()
    _write_geotiff(folder / "N32E035.hgt", _ramp_band(), _WEST, _NORTH, _RES)
    return str(folder)


@pytest.fixture()
def srtm(topo, srtm_folder):
    """The same folder, registered as a 'string' datasource, name returned."""
    topo.addDataSource("SRTM", srtm_folder, "string", version=(0, 0, 1))
    topo.setDataSourceDefaultVersion("SRTM", (0, 0, 1))
    return "SRTM"


# Points chosen so that (lat - north) / res and (lon - west) / res are exact.
# (lat, lon, rasterx, rastery)
_EXACT_POINTS = [
    (32.96875, 35.03125, 0.5, 0.5),
    (32.90625, 35.15625, 1.5, 2.5),
    (32.75000, 35.06250, 4.0, 1.0),
    (32.81250, 35.31250, 3.0, 5.0),
]

# A small box inside the tile, roughly 275 m by 235 m.  Note the argument
# names: for WSG84 the code reads minx/maxx as LATITUDE and miny/maxy as
# LONGITUDE -- see TestCreateXarrayArgumentOrder (B191).
_BOX = dict(minx=32.9000, miny=35.1500, maxx=32.9025, maxy=35.1525)


# --------------------------------------------------------------------------
# _open_raster
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestOpenRaster:
    """The one rasterio seam the whole module reads elevation through."""

    def test_it_returns_the_first_band_as_floats(self, srtm_folder):
        array, _ = _open_raster(os.path.join(srtm_folder, "N32E035.hgt"))
        assert array.dtype == np.float64
        np.testing.assert_allclose(array, _ramp_band().astype(float))

    def test_an_integer_band_is_promoted_to_float(self, tmp_path):
        path = _write_geotiff(tmp_path / "ints.tif",
                              np.arange(9, dtype="int16").reshape(3, 3), 35.0, 33.0, 0.5)
        array, _ = _open_raster(path)
        assert array.dtype == np.float64

    def test_the_geotransform_uses_the_gdal_convention(self, srtm_folder):
        """(xmin, xres, 0, ymax, 0, -yres) -- note the negative y resolution."""
        _, geotransform = _open_raster(os.path.join(srtm_folder, "N32E035.hgt"))
        assert len(geotransform) == 6
        assert geotransform[0] == pytest.approx(_WEST)
        assert geotransform[1] == pytest.approx(_RES)
        assert geotransform[2] == pytest.approx(0.0)
        assert geotransform[3] == pytest.approx(_NORTH)
        assert geotransform[4] == pytest.approx(0.0)
        assert geotransform[5] == pytest.approx(-_RES)

    def test_it_returns_a_tuple_not_an_affine(self, srtm_folder):
        _, geotransform = _open_raster(os.path.join(srtm_folder, "N32E035.hgt"))
        assert isinstance(geotransform, tuple)

    def test_a_missing_file_propagates_the_rasterio_error(self, tmp_path):
        with pytest.raises(Exception):
            _open_raster(str(tmp_path / "absent.tif"))


# --------------------------------------------------------------------------
# findElevationFile
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestFindElevationFile:
    def test_a_single_folder_string_is_treated_as_a_one_element_list(self, topo, srtm):
        found = topo.findElevationFile("N32E035.hgt", srtm)
        assert os.path.isfile(found)
        assert os.path.basename(found) == "N32E035.hgt"

    def test_it_joins_the_folder_and_the_filename_with_os_path_join(self, topo, srtm, srtm_folder):
        assert topo.findElevationFile("N32E035.hgt", srtm) == \
            os.path.join(srtm_folder, "N32E035.hgt")

    def test_it_searches_every_folder_of_a_list_datasource(self, topo, tmp_path, srtm_folder):
        empty = tmp_path / "empty"
        empty.mkdir()
        topo.addDataSource("MANY", [str(empty), srtm_folder], "JSON_dict", version=(0, 0, 1))
        assert topo.findElevationFile("N32E035.hgt", "MANY") == \
            os.path.join(srtm_folder, "N32E035.hgt")

    def test_it_returns_the_first_folder_that_has_the_tile(self, topo, tmp_path, srtm_folder):
        first = tmp_path / "first"
        first.mkdir()
        _write_geotiff(first / "N32E035.hgt", _ramp_band(), _WEST, _NORTH, _RES)
        topo.addDataSource("MANY", [str(first), srtm_folder], "JSON_dict", version=(0, 0, 1))
        assert topo.findElevationFile("N32E035.hgt", "MANY") == \
            os.path.join(str(first), "N32E035.hgt")

    def test_a_missing_tile_raises_file_not_found_naming_the_folders(self, topo, srtm, srtm_folder):
        with pytest.raises(FileNotFoundError) as raised:
            topo.findElevationFile("N99E099.hgt", srtm)
        assert "N99E099.hgt" in str(raised.value)
        assert srtm_folder in str(raised.value)

    def test_an_unknown_datasource_name_raises_type_error_on_the_none_it_gets(self, topo):
        """getDataSourceData returns None, which the folder loop cannot iterate."""
        with pytest.raises(TypeError):
            topo.findElevationFile("N32E035.hgt", "NO_SUCH_SOURCE")


# --------------------------------------------------------------------------
# getPointElevation
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestGetPointElevation:
    """Bilinear interpolation of a single point out of one SRTM tile."""

    @pytest.mark.parametrize("lat, lon, rasterx, rastery", _EXACT_POINTS)
    def test_it_bilinearly_interpolates_the_band(self, topo, srtm, lat, lon, rasterx, rastery):
        assert topo.getPointElevation(lat=lat, long=lon, dataSourceName=srtm) == \
            pytest.approx(10.0 * rasterx + 100.0 * rastery)

    def test_an_exact_pixel_centre_returns_that_pixels_value(self, topo, srtm):
        band = _ramp_band()
        # rasterx = 4, rastery = 1 -> the pixel itself, no interpolation weight.
        assert topo.getPointElevation(lat=32.75, long=35.0625, dataSourceName=srtm) == \
            pytest.approx(band[4, 1])

    def test_the_result_grows_eastwards_and_southwards(self, topo, srtm):
        """The ramp increases with both row and column, so both must show."""
        base = topo.getPointElevation(lat=32.90625, long=35.15625, dataSourceName=srtm)
        east = topo.getPointElevation(lat=32.90625, long=35.21875, dataSourceName=srtm)
        south = topo.getPointElevation(lat=32.84375, long=35.15625, dataSourceName=srtm)
        assert east > base
        assert south > base
        assert east - base == pytest.approx(100.0)
        assert south - base == pytest.approx(10.0)

    def test_at_the_southern_edge_it_degenerates_to_the_last_pixel(self, topo, srtm):
        """The +1 guard collapses all four corners onto height11."""
        band = _ramp_band()
        # rasterx = 7.5 -> int 7, and 7 + 1 == the row count, so the guard fires.
        assert topo.getPointElevation(lat=32.53125, long=35.15625, dataSourceName=srtm) == \
            pytest.approx(band[7, 2])

    def test_values_below_minus_one_thousand_are_treated_as_sea_level(self, topo, tmp_path):
        band = np.full((4, 4), 500.0, dtype="float32")
        band[0:2, 0:2] = -9999.0
        folder = tmp_path / "voids"
        folder.mkdir()
        _write_geotiff(folder / "N32E035.hgt", band, _WEST, _NORTH, _RES)
        topo.addDataSource("VOIDS", str(folder), "string", version=(0, 0, 1))
        assert topo.getPointElevation(lat=32.96875, long=35.03125,
                                      dataSourceName="VOIDS") == pytest.approx(0.0)

    def test_the_untouched_neighbourhood_keeps_its_value(self, topo, tmp_path):
        """The control for the void test: elsewhere the band is unchanged."""
        band = np.full((4, 4), 500.0, dtype="float32")
        band[0:2, 0:2] = -9999.0
        folder = tmp_path / "voids"
        folder.mkdir()
        _write_geotiff(folder / "N32E035.hgt", band, _WEST, _NORTH, _RES)
        topo.addDataSource("VOIDS", str(folder), "string", version=(0, 0, 1))
        assert topo.getPointElevation(lat=32.84375, long=35.15625,
                                      dataSourceName="VOIDS") == pytest.approx(500.0)

    def test_it_reads_the_configured_default_srtm_when_no_source_is_named(self, topo, srtm_folder):
        topo.addDataSource("SRTM", srtm_folder, "string", version=(0, 0, 1))
        topo.setConfig(defaultSRTM="SRTM")
        assert topo.getPointElevation(lat=32.90625, long=35.15625) == pytest.approx(265.0)

    def test_a_point_in_a_tile_that_is_not_there_raises_file_not_found(self, topo, srtm):
        with pytest.raises(FileNotFoundError, match=r"N40E040\.hgt"):
            topo.getPointElevation(lat=40.5, long=40.5, dataSourceName=srtm)


@pytest.mark.unit
class TestSRTMTileNaming:
    """B188: the tile name is built without padding or hemisphere letters."""

    @staticmethod
    def _requested_name(topo, srtm, lat, lon):
        """The filename the code asked for, read back out of its own error."""
        with pytest.raises(FileNotFoundError) as raised:
            topo.getPointElevation(lat=lat, long=lon, dataSourceName=srtm)
        return str(raised.value).split(" not found")[0]

    def test_the_longitude_is_padded_to_three_digits(self, topo, srtm):
        assert self._requested_name(topo, srtm, 40.5, 7.5).endswith("E007.hgt")

    def test_the_latitude_is_not_padded_at_all(self, topo, srtm):
        """Characterisation of B188: real SRTM tiles are N05..., not N5...."""
        assert self._requested_name(topo, srtm, 5.5, 7.5) == "N5E007.hgt"

    def test_a_southern_latitude_becomes_a_negative_n_tile(self, topo, srtm):
        """Characterisation of B188: no S hemisphere letter."""
        assert self._requested_name(topo, srtm, -1.5, 7.5) == "N-1E007.hgt"

    def test_a_western_longitude_becomes_a_negative_e_tile(self, topo, srtm):
        """Characterisation of B188: no W hemisphere letter."""
        assert self._requested_name(topo, srtm, 40.5, -2.5) == "N40E-02.hgt"

    def test_it_truncates_towards_zero_rather_than_flooring(self, topo, srtm):
        """Characterisation of B188: int(-1.5) is -1, but the tile is S02."""
        assert self._requested_name(topo, srtm, -1.5, 7.5) != "N-2E007.hgt"

    @pytest.mark.xfail(
        strict=True,
        reason="B188: getPointElevation and getPointListElevation build the "
               "tile name as 'N' + str(int(lat)) + 'E' + str(int(long)).zfill(3). "
               "The SRTM/HGT convention is two zero-padded latitude digits plus "
               "a hemisphere letter -- N05E007.hgt, S02W003.hgt -- so every "
               "latitude below 10 degrees, every southern latitude and every "
               "western longitude asks for a filename no SRTM archive contains "
               "('N5E007.hgt', 'N-1E007.hgt', 'N40E-02.hgt'). int() also "
               "truncates towards zero where the tile index needs a floor, so "
               "-1.5 must map to S02, not to 1. Israel sits at 29-33N / 34-36E, "
               "which is why this has never bitten. See the consolidated "
               "findings issue.",
    )
    @pytest.mark.parametrize(
        "lat, lon, expected",
        [(5.5, 7.5, "N05E007.hgt"), (-1.5, 7.5, "S02E007.hgt"), (40.5, -2.5, "N40W003.hgt")],
    )
    def test_the_tile_name_follows_the_srtm_convention(self, topo, srtm, lat, lon, expected):
        assert self._requested_name(topo, srtm, lat, lon) == expected


@pytest.mark.unit
class TestMissingDefaultSRTM:
    """B190: the 'datasource is not found' guard cannot fire on a missing key."""

    def test_the_guard_does_fire_when_the_config_stores_an_explicit_none(self, topo):
        topo.setConfig(defaultSRTM=None)
        with pytest.raises(ValueError, match="datasource is not found"):
            topo.getPointElevation(lat=32.90625, long=35.15625)

    def test_a_missing_config_key_raises_key_error_instead(self, topo):
        """Characterisation of B190."""
        with pytest.raises(KeyError, match="defaultSRTM"):
            topo.getPointElevation(lat=32.90625, long=35.15625)

    def test_the_point_list_method_has_the_same_shape_of_guard(self, topo):
        """Characterisation of B190 at the second call site."""
        import pandas as pd

        with pytest.raises(KeyError, match="defaultSRTM"):
            topo.getPointListElevation(pd.DataFrame({"lat": [32.90625], "lon": [35.15625]}))

    @pytest.mark.xfail(
        strict=True,
        reason="B190: getPointElevation and getPointListElevation both do "
               "`dataSourceName = self.getConfig()['defaultSRTM']` and only then "
               "check `if dataSourceName is None: raise ValueError('The "
               "datasource is not found!')`. A project that never configured a "
               "default has no such key, so the bare dict subscript raises "
               "KeyError('defaultSRTM') first and the deliberate, readable "
               "ValueError is unreachable except in the odd case of a config "
               "that stores None explicitly. getConfig().get('defaultSRTM') "
               "would make the guard do its job. See the consolidated findings "
               "issue.",
    )
    def test_a_missing_default_reports_the_intended_error(self, topo):
        with pytest.raises(ValueError, match="datasource is not found"):
            topo.getPointElevation(lat=32.90625, long=35.15625)


# --------------------------------------------------------------------------
# getPointListElevation
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestGetPointListElevationFromDataFrame:
    """The plain-DataFrame branch: lat/lon columns."""

    @pytest.fixture()
    def points(self):
        import pandas as pd

        return pd.DataFrame(
            {"lat": [lat for lat, _, _, _ in _EXACT_POINTS],
             "lon": [lon for _, lon, _, _ in _EXACT_POINTS]}
        )

    def test_it_adds_a_filename_and_an_elevation_column(self, topo, srtm, points):
        result = topo.getPointListElevation(points, srtm)
        assert {"filename", "elevation"} <= set(result.columns)
        assert list(result["filename"].unique()) == ["N32E035.hgt"]

    def test_every_elevation_matches_the_bilinear_closed_form(self, topo, srtm, points):
        result = topo.getPointListElevation(points, srtm)
        expected = [10.0 * rasterx + 100.0 * rastery for _, _, rasterx, rastery in _EXACT_POINTS]
        np.testing.assert_allclose(result["elevation"].values, expected)

    def test_it_agrees_with_the_single_point_method(self, topo, srtm, points):
        result = topo.getPointListElevation(points, srtm)
        for _, row in result.iterrows():
            assert row["elevation"] == pytest.approx(
                topo.getPointElevation(lat=row["lat"], long=row["lon"], dataSourceName=srtm)
            )

    def test_it_preserves_the_input_index_and_row_order(self, topo, srtm, points):
        shuffled = points.iloc[::-1]
        result = topo.getPointListElevation(shuffled, srtm)
        assert list(result.index) == list(shuffled.index)

    def test_it_returns_a_frame_of_the_same_length(self, topo, srtm, points):
        assert len(topo.getPointListElevation(points, srtm)) == len(points)

    def test_it_reads_the_configured_default_srtm(self, topo, srtm_folder, points):
        topo.addDataSource("SRTM", srtm_folder, "string", version=(0, 0, 1))
        topo.setConfig(defaultSRTM="SRTM")
        assert "elevation" in topo.getPointListElevation(points)

    def test_points_are_grouped_so_each_tile_is_opened_once(self, topo, tmp_path, monkeypatch):
        import pandas as pd

        folder = tmp_path / "two_tiles"
        folder.mkdir()
        _write_geotiff(folder / "N32E035.hgt", _ramp_band(), _WEST, _NORTH, _RES)
        _write_geotiff(folder / "N32E036.hgt", _ramp_band(), _WEST + 1.0, _NORTH, _RES)
        topo.addDataSource("TWO", str(folder), "string", version=(0, 0, 1))

        opened = []
        real_open = topography_module._open_raster

        def counting_open(fheight):
            opened.append(fheight)
            return real_open(fheight)

        monkeypatch.setattr(topography_module, "_open_raster", counting_open)

        points = pd.DataFrame(
            {"lat": [32.90625, 32.84375, 32.90625],
             "lon": [35.15625, 35.15625, 36.15625]}
        )
        result = topo.getPointListElevation(points, "TWO")
        assert sorted(os.path.basename(path) for path in opened) == \
            ["N32E035.hgt", "N32E036.hgt"]
        assert sorted(result["filename"].unique()) == ["N32E035.hgt", "N32E036.hgt"]

    def test_two_tiles_are_interpolated_independently(self, topo, tmp_path):
        import pandas as pd

        folder = tmp_path / "two_tiles"
        folder.mkdir()
        _write_geotiff(folder / "N32E035.hgt", _ramp_band(), _WEST, _NORTH, _RES)
        _write_geotiff(folder / "N32E036.hgt", _ramp_band(), _WEST + 1.0, _NORTH, _RES)
        topo.addDataSource("TWO", str(folder), "string", version=(0, 0, 1))
        points = pd.DataFrame({"lat": [32.90625, 32.90625], "lon": [35.15625, 36.15625]})
        result = topo.getPointListElevation(points, "TWO")
        # Both tiles carry the same band and both points sit at the same offset.
        np.testing.assert_allclose(result["elevation"].values, [265.0, 265.0])


@pytest.mark.unit
class TestGetPointListElevationFromGeometry:
    """The GeoSeries and GeoDataFrame branches."""

    @staticmethod
    def _geoseries():
        import geopandas as gpd
        from shapely.geometry import Point

        return gpd.GeoSeries(
            [Point(lon, lat) for lat, lon, _, _ in _EXACT_POINTS], crs=WSG84
        )

    def test_a_geoseries_gains_point_lon_lat_and_elevation(self, topo, srtm):
        result = topo.getPointListElevation(self._geoseries(), srtm)
        assert {"point", "filename", "lon", "lat", "elevation"} <= set(result.columns)

    def test_a_geoseries_reads_x_as_longitude_and_y_as_latitude(self, topo, srtm):
        result = topo.getPointListElevation(self._geoseries(), srtm)
        np.testing.assert_allclose(result["lat"].values,
                                   [lat for lat, _, _, _ in _EXACT_POINTS])
        np.testing.assert_allclose(result["lon"].values,
                                   [lon for _, lon, _, _ in _EXACT_POINTS])

    def test_a_geoseries_gets_the_same_elevations_as_a_dataframe(self, topo, srtm):
        result = topo.getPointListElevation(self._geoseries(), srtm)
        expected = [10.0 * rasterx + 100.0 * rastery for _, _, rasterx, rastery in _EXACT_POINTS]
        np.testing.assert_allclose(result["elevation"].values, expected)

    def test_a_geodataframe_with_a_point_column_is_accepted(self, topo, srtm):
        import geopandas as gpd

        frame = gpd.GeoDataFrame(
            {"id": list(range(len(_EXACT_POINTS))), "point": self._geoseries()},
            geometry="point", crs=WSG84,
        )
        result = topo.getPointListElevation(frame, srtm)
        expected = [10.0 * rasterx + 100.0 * rastery for _, _, rasterx, rastery in _EXACT_POINTS]
        np.testing.assert_allclose(result["elevation"].values, expected)
        assert list(result["id"]) == list(range(len(_EXACT_POINTS)))

    def test_a_geodataframe_without_a_point_column_raises_value_error(self, topo, srtm):
        import geopandas as gpd

        frame = gpd.GeoDataFrame({"id": [0]}, geometry=self._geoseries()[:1], crs=WSG84)
        with pytest.raises(ValueError, match="must contain a field 'point'"):
            topo.getPointListElevation(frame, srtm)


@pytest.mark.unit
class TestGetPointListElevationIgnoresInputCRS:
    """B189: the inputCRS parameter is accepted and never read."""

    @staticmethod
    def _wgs_frame():
        import pandas as pd

        return pd.DataFrame({"lat": [32.90625], "lon": [35.15625]})

    @staticmethod
    def _itm_frame():
        point = convertCRS([[35.15625, 32.90625]], inputCRS=WSG84, outputCRS=ITM)[0]
        import pandas as pd

        return pd.DataFrame({"lat": [point.y], "lon": [point.x]})

    def test_the_parameter_is_part_of_the_signature(self):
        import inspect

        signature = inspect.signature(
            topography_module.TopographyToolkit.getPointListElevation
        )
        assert signature.parameters["inputCRS"].default == WSG84

    def test_passing_it_changes_nothing_for_wgs84_numbers(self, topo, srtm):
        """Characterisation of B189."""
        default = topo.getPointListElevation(self._wgs_frame(), srtm)
        declared = topo.getPointListElevation(self._wgs_frame(), srtm, inputCRS=ITM)
        np.testing.assert_allclose(declared["elevation"].values, default["elevation"].values)

    def test_itm_coordinates_are_read_as_degrees_and_ask_for_a_nonsense_tile(self, topo, srtm):
        """Characterisation of B189: 200 000 m becomes latitude 200000."""
        with pytest.raises(FileNotFoundError, match=r"N\d{6}E\d{6}\.hgt"):
            topo.getPointListElevation(self._itm_frame(), srtm, inputCRS=ITM)

    @pytest.mark.xfail(
        strict=True,
        reason="B189: getPointListElevation declares inputCRS=WSG84 in its "
               "signature and then never uses the name -- no reprojection, no "
               "validation, not even a warning. A caller who hands it ITM "
               "eastings/northings and says inputCRS=ITM (the project's own "
               "constant) has the metres read as degrees: the tile name comes "
               "out as N753355E212021.hgt and the call dies in "
               "findElevationFile. Either honour the parameter by converting "
               "through convertPointsCRS, or drop it from the signature. See "
               "the consolidated findings issue.",
    )
    def test_the_same_physical_point_gives_the_same_elevation_in_either_crs(self, topo, srtm):
        expected = topo.getPointListElevation(self._wgs_frame(), srtm)["elevation"].iloc[0]
        got = topo.getPointListElevation(self._itm_frame(), srtm, inputCRS=ITM)
        assert got["elevation"].iloc[0] == pytest.approx(expected)


# --------------------------------------------------------------------------
# convertPointsCRS
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestConvertPointsCRS:
    """The toolkit's own CRS conversion, one branch per input type."""

    LON, LAT = 35.15625, 32.90625

    def _reference(self):
        return convertCRS([[self.LON, self.LAT]], inputCRS=WSG84, outputCRS=ITM)[0]

    def test_a_one_dimensional_array_is_a_single_point(self, topo):
        result = topo.convertPointsCRS(np.array([self.LON, self.LAT]), WSG84, ITM)
        assert len(result) == 1
        assert result.geometry.iloc[0].x == pytest.approx(self._reference().x)
        assert result.geometry.iloc[0].y == pytest.approx(self._reference().y)

    def test_a_two_dimensional_array_converts_every_row(self, topo):
        points = np.array([[self.LON, self.LAT], [self.LON + 0.01, self.LAT + 0.01]])
        result = topo.convertPointsCRS(points, WSG84, ITM)
        assert len(result) == 2
        assert result.geometry.iloc[1].x > result.geometry.iloc[0].x

    def test_a_dataframe_uses_the_x_and_y_columns_by_default(self, topo):
        import pandas as pd

        frame = pd.DataFrame({"x": [self.LON, self.LON + 0.01], "y": [self.LAT, self.LAT + 0.01]})
        result = topo.convertPointsCRS(frame, WSG84, ITM)
        assert result.geometry.iloc[0].x == pytest.approx(self._reference().x)

    def test_a_dataframe_honours_custom_column_names(self, topo):
        import pandas as pd

        frame = pd.DataFrame({"lon": [self.LON, self.LON + 0.01],
                              "lat": [self.LAT, self.LAT + 0.01]})
        result = topo.convertPointsCRS(frame, WSG84, ITM, x="lon", y="lat")
        assert result.geometry.iloc[0].x == pytest.approx(self._reference().x)

    def test_a_single_element_list_takes_the_scalar_transform_path(self, topo):
        result = topo.convertPointsCRS([[self.LON, self.LAT]], WSG84, ITM)
        assert len(result) == 1
        assert result.crs.to_epsg() == ITM
        assert result.geometry.iloc[0].x == pytest.approx(self._reference().x)

    def test_a_multi_element_list_converts_every_point(self, topo):
        result = topo.convertPointsCRS(
            [[self.LON, self.LAT], [self.LON + 0.01, self.LAT], [self.LON, self.LAT + 0.01]],
            WSG84, ITM,
        )
        assert len(result) == 3
        assert result.crs.to_epsg() == ITM

    def test_the_result_carries_the_output_crs(self, topo):
        for points in ([[self.LON, self.LAT]],
                       [[self.LON, self.LAT], [self.LON + 0.01, self.LAT]]):
            assert topo.convertPointsCRS(points, WSG84, ITM).crs.to_epsg() == ITM

    def test_it_agrees_with_the_shared_convertCRS_helper(self, topo):
        points = [[self.LON, self.LAT], [self.LON + 0.01, self.LAT + 0.01]]
        mine = topo.convertPointsCRS(points, WSG84, ITM)
        theirs = convertCRS(points, inputCRS=WSG84, outputCRS=ITM)
        for got, want in zip(mine.geometry, theirs):
            assert got.x == pytest.approx(want.x)
            assert got.y == pytest.approx(want.y)

    def test_a_round_trip_returns_the_original_degrees(self, topo):
        itm = topo.convertPointsCRS([[self.LON, self.LAT]], WSG84, ITM)
        back = topo.convertPointsCRS(
            [[itm.geometry.iloc[0].x, itm.geometry.iloc[0].y]], ITM, WSG84
        )
        assert back.geometry.iloc[0].x == pytest.approx(self.LON, abs=1e-7)
        assert back.geometry.iloc[0].y == pytest.approx(self.LAT, abs=1e-7)

    def test_an_unsupported_input_type_raises_value_error(self, topo):
        with pytest.raises(ValueError, match="Unsupported type"):
            topo.convertPointsCRS("35,32", WSG84, ITM)

    def test_a_tuple_is_not_supported(self, topo):
        with pytest.raises(ValueError, match="Unsupported type"):
            topo.convertPointsCRS(((self.LON, self.LAT),), WSG84, ITM)


# --------------------------------------------------------------------------
# create_xarray
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestCreateXarray:
    """The (i, j) grid builder, including its memory guard."""

    def test_it_returns_lat_and_lon_over_the_i_j_dimensions(self, topo):
        dataset = topo.create_xarray(dxdy=100, **_BOX)
        assert set(dataset.dims) == {"i", "j"}
        assert dataset["lat"].dims == ("i", "j")
        assert dataset["lon"].dims == ("i", "j")

    def test_lat_and_lon_are_data_variables_not_coordinates(self, topo):
        dataset = topo.create_xarray(dxdy=100, **_BOX)
        assert set(dataset.data_vars) == {"lat", "lon"}
        assert set(dataset.coords) == {"i", "j"}

    def test_the_index_coordinates_run_from_zero(self, topo):
        dataset = topo.create_xarray(dxdy=100, **_BOX)
        np.testing.assert_array_equal(dataset["i"].values, np.arange(dataset.sizes["i"]))
        np.testing.assert_array_equal(dataset["j"].values, np.arange(dataset.sizes["j"]))

    def test_a_projected_input_crs_uses_the_coordinates_as_given(self, topo):
        dataset = topo.create_xarray(minx=200000, miny=650000, maxx=200300, maxy=650300,
                                     dxdy=100, inputCRS=ITM)
        assert dataset["lat"].shape == (3, 3)

    def test_a_projected_grid_comes_back_in_degrees(self, topo):
        dataset = topo.create_xarray(minx=200000, miny=650000, maxx=200300, maxy=650300,
                                     dxdy=100, inputCRS=ITM)
        assert 31.0 < float(dataset["lat"].values.min()) < 33.0
        assert 34.0 < float(dataset["lon"].values.min()) < 36.0

    def test_the_first_row_is_the_northernmost(self, topo):
        """np.arange(...)[::-1] on the y axis, so i increases southwards."""
        dataset = topo.create_xarray(minx=200000, miny=650000, maxx=200500, maxy=650500,
                                     dxdy=100, inputCRS=ITM)
        latitudes = dataset["lat"].values
        assert latitudes[0, 0] > latitudes[-1, 0]

    def test_columns_run_west_to_east(self, topo):
        dataset = topo.create_xarray(minx=200000, miny=650000, maxx=200500, maxy=650500,
                                     dxdy=100, inputCRS=ITM)
        longitudes = dataset["lon"].values
        assert longitudes[0, 0] < longitudes[0, -1]

    def test_a_coarser_dxdy_gives_fewer_cells(self, topo):
        fine = topo.create_xarray(dxdy=50, **_BOX)
        coarse = topo.create_xarray(dxdy=100, **_BOX)
        assert fine["lat"].size > coarse["lat"].size

    def test_the_cell_spacing_is_dxdy_metres_on_the_ground(self, topo):
        dataset = topo.create_xarray(minx=200000, miny=650000, maxx=200500, maxy=650500,
                                     dxdy=100, inputCRS=ITM)
        first = convertCRS([[float(dataset["lon"].values[0, 0]),
                             float(dataset["lat"].values[0, 0])]],
                           inputCRS=WSG84, outputCRS=ITM)[0]
        second = convertCRS([[float(dataset["lon"].values[0, 1]),
                              float(dataset["lat"].values[0, 1])]],
                            inputCRS=WSG84, outputCRS=ITM)[0]
        assert second.x - first.x == pytest.approx(100.0, abs=0.1)

    def test_an_oversized_grid_is_refused_before_it_is_allocated(self, topo):
        with pytest.raises(MemoryError, match="Too many grid points"):
            topo.create_xarray(minx=0, miny=0, maxx=2000, maxy=2000, dxdy=1, inputCRS=ITM)

    def test_the_guard_reports_the_offending_shape(self, topo):
        with pytest.raises(MemoryError) as raised:
            topo.create_xarray(minx=0, miny=0, maxx=2000, maxy=2000, dxdy=1, inputCRS=ITM)
        assert "2000 x 2000" in str(raised.value)
        assert "Increase dxdy" in str(raised.value)

    def test_a_grid_just_under_the_guard_is_built(self, topo):
        dataset = topo.create_xarray(minx=0, miny=0, maxx=500, maxy=500, dxdy=1, inputCRS=ITM)
        assert dataset["lat"].shape == (500, 500)


@pytest.mark.unit
class TestCreateXarrayArgumentOrder:
    """B191: minx/miny mean (lat, lon) for WSG84 but (x, y) for anything else."""

    # An Israeli box, stated the way the parameter names read: x is longitude.
    LON_MIN, LON_MAX = 35.1500, 35.1525
    LAT_MIN, LAT_MAX = 32.9000, 32.9025

    def test_naming_them_lat_then_lon_puts_the_grid_where_asked(self, topo):
        """Characterisation of B191: minx is the latitude."""
        dataset = topo.create_xarray(minx=self.LAT_MIN, miny=self.LON_MIN,
                                     maxx=self.LAT_MAX, maxy=self.LON_MAX, dxdy=100)
        assert self.LON_MIN - 1e-3 <= float(dataset["lon"].values.min())
        assert float(dataset["lon"].values.max()) <= self.LON_MAX + 1e-3
        assert self.LAT_MIN - 1e-3 <= float(dataset["lat"].values.min())

    def test_naming_them_lon_then_lat_silently_relocates_the_grid(self, topo):
        """Characterisation of B191: the grid lands at lon 32, lat 35."""
        dataset = topo.create_xarray(minx=self.LON_MIN, miny=self.LAT_MIN,
                                     maxx=self.LON_MAX, maxy=self.LAT_MAX, dxdy=100)
        assert 32.0 < float(dataset["lon"].values.mean()) < 33.0
        assert 35.0 < float(dataset["lat"].values.mean()) < 36.0

    def test_the_projected_branch_reads_them_as_genuine_x_and_y(self, topo):
        """The inconsistency: with inputCRS=ITM, minx really is the easting."""
        dataset = topo.create_xarray(minx=200000, miny=650000, maxx=200300, maxy=650300,
                                     dxdy=100, inputCRS=ITM)
        corner = convertCRS([[float(dataset["lon"].values[-1, 0]),
                              float(dataset["lat"].values[-1, 0])]],
                            inputCRS=WSG84, outputCRS=ITM)[0]
        assert corner.x == pytest.approx(200000, abs=1.0)
        assert corner.y == pytest.approx(650000, abs=1.0)

    def test_the_class_uses_x_for_longitude_everywhere_else(self, topo):
        """convertPointsCRS in the same class takes (x, y) = (lon, lat)."""
        point = topo.convertPointsCRS([[self.LON_MIN, self.LAT_MIN]], WSG84, ITM)
        assert 100000 < point.geometry.iloc[0].x < 300000
        assert 400000 < point.geometry.iloc[0].y < 800000

    @pytest.mark.xfail(
        strict=True,
        reason="B191: TopographyToolkit.create_xarray's WSG84 branch calls "
               "self.convertPointsCRS(points=[[miny, minx]], ...), i.e. it "
               "feeds miny as the x (longitude) and minx as the y (latitude). "
               "The parameters named minx/maxx are therefore LATITUDES and "
               "miny/maxy LONGITUDES -- the opposite of their names, of "
               "convertPointsCRS's own (x, y) = (lon, lat) contract in the same "
               "class, and of the same method's `else` branch, which builds "
               "Point(minx, miny) and so does read minx as the easting. The "
               "docstring says only 'bounding box coordinates'. A caller "
               "following the names gets a grid silently placed at lon 32 / lat "
               "35 instead of lon 35 / lat 32, and getElevation and "
               "createElevationSTL inherit the trap. See the consolidated "
               "findings issue.",
    )
    def test_minx_is_the_longitude_its_name_implies(self, topo):
        dataset = topo.create_xarray(minx=self.LON_MIN, miny=self.LAT_MIN,
                                     maxx=self.LON_MAX, maxy=self.LAT_MAX, dxdy=100)
        assert self.LON_MIN - 1e-3 <= float(dataset["lon"].values.min())
        assert float(dataset["lon"].values.max()) <= self.LON_MAX + 1e-3


# --------------------------------------------------------------------------
# getElevation / getElevationOfXarray
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestGetElevation:
    """create_xarray, then one elevation per cell, reshaped back onto (i, j)."""

    def test_it_adds_an_elevation_coordinate_shaped_like_lat(self, topo, srtm):
        dataset = topo.getElevation(dxdy=100, dataSourceName=srtm, **_BOX)
        assert "elevation" in dataset.coords
        assert dataset["elevation"].dims == ("i", "j")
        assert dataset["elevation"].shape == dataset["lat"].shape

    def test_it_keeps_the_lat_and_lon_grid_untouched(self, topo, srtm):
        dataset = topo.getElevation(dxdy=100, dataSourceName=srtm, **_BOX)
        plain = topo.create_xarray(dxdy=100, **_BOX)
        np.testing.assert_allclose(dataset["lat"].values, plain["lat"].values)
        np.testing.assert_allclose(dataset["lon"].values, plain["lon"].values)

    def test_every_cell_matches_the_bilinear_closed_form_at_its_own_lat_lon(self, topo, srtm):
        """Pins the reshape orientation: a transposed grid would not match."""
        dataset = topo.getElevation(dxdy=100, dataSourceName=srtm, **_BOX)
        expected = _expected_elevation(dataset["lat"].values, dataset["lon"].values)
        np.testing.assert_allclose(dataset["elevation"].values, expected, rtol=1e-9)

    def test_it_equals_the_two_steps_applied_by_hand(self, topo, srtm):
        import pandas as pd

        dataset = topo.getElevation(dxdy=100, dataSourceName=srtm, **_BOX)
        grid = topo.create_xarray(dxdy=100, **_BOX)
        points = pd.DataFrame({"lat": grid["lat"].values.flatten(),
                               "lon": grid["lon"].values.flatten()})
        stepwise = topo.getPointListElevation(points, srtm)
        np.testing.assert_allclose(
            dataset["elevation"].values.flatten(), stepwise["elevation"].values
        )

    def test_it_reads_the_configured_default_srtm(self, topo, srtm_folder):
        topo.addDataSource("SRTM", srtm_folder, "string", version=(0, 0, 1))
        topo.setConfig(defaultSRTM="SRTM")
        assert "elevation" in topo.getElevation(dxdy=100, **_BOX).coords

    def test_the_memory_guard_still_applies(self, topo, srtm):
        with pytest.raises(MemoryError):
            topo.getElevation(minx=0, miny=0, maxx=2000, maxy=2000, dxdy=1,
                              inputCRS=ITM, dataSourceName=srtm)


@pytest.mark.unit
class TestGetElevationOfXarray:
    """The same, for a grid the caller already has."""

    def test_it_adds_elevation_to_an_existing_grid(self, topo, srtm):
        grid = topo.create_xarray(dxdy=100, **_BOX)
        result = topo.getElevationOfXarray(grid, srtm)
        assert "elevation" in result.coords
        assert result["elevation"].shape == grid["lat"].shape

    def test_it_agrees_with_getElevation_on_the_same_box(self, topo, srtm):
        grid = topo.create_xarray(dxdy=100, **_BOX)
        np.testing.assert_allclose(
            topo.getElevationOfXarray(grid, srtm)["elevation"].values,
            topo.getElevation(dxdy=100, dataSourceName=srtm, **_BOX)["elevation"].values,
        )

    def test_it_accepts_a_hand_built_dataset_of_any_shape(self, topo, srtm):
        import xarray as xr

        latitudes = np.array([[32.96875, 32.90625, 32.84375],
                              [32.75000, 32.68750, 32.62500]])
        longitudes = np.array([[35.03125, 35.15625, 35.06250],
                               [35.06250, 35.31250, 35.15625]])
        grid = xr.Dataset(
            coords={"i": np.arange(2), "j": np.arange(3)},
            data_vars={"lat": (["i", "j"], latitudes), "lon": (["i", "j"], longitudes)},
        )
        result = topo.getElevationOfXarray(grid, srtm)
        np.testing.assert_allclose(result["elevation"].values,
                                   _expected_elevation(latitudes, longitudes), rtol=1e-9)

    def test_it_leaves_the_original_lat_lon_in_place(self, topo, srtm):
        grid = topo.create_xarray(dxdy=100, **_BOX)
        result = topo.getElevationOfXarray(grid, srtm)
        np.testing.assert_allclose(result["lat"].values, grid["lat"].values)
        np.testing.assert_allclose(result["lon"].values, grid["lon"].values)


# --------------------------------------------------------------------------
# getElevationSTL / createElevationSTL
# --------------------------------------------------------------------------

def _stl_vertices(stl):
    """Every 'vertex x y z' triple in an ASCII STL string, in order."""
    return [
        tuple(float(token) for token in line.split()[1:])
        for line in stl.splitlines()
        if line.strip().startswith("vertex")
    ]


@pytest.mark.unit
class TestGetElevationSTL:
    """Elevation grid -> ASCII STL, in ITM metres."""

    @pytest.fixture()
    def elevation(self, topo, srtm):
        return topo.getElevation(dxdy=100, dataSourceName=srtm, **_BOX)

    def test_it_returns_an_ascii_stl_solid(self, topo, elevation):
        stl = topo.getElevationSTL(elevation)
        assert stl.startswith("solid Topography\n")
        assert stl.rstrip().endswith("endsolid Topography")

    def test_the_solid_name_is_honoured_at_both_ends(self, topo, elevation):
        stl = topo.getElevationSTL(elevation, solidName="Golan")
        assert stl.startswith("solid Golan\n")
        assert stl.rstrip().endswith("endsolid Golan")

    def test_every_facet_is_closed_and_has_three_vertices(self, topo, elevation):
        stl = topo.getElevationSTL(elevation)
        facets = stl.count("facet normal")
        assert facets > 0
        assert stl.count("endfacet") == facets
        assert stl.count("outer loop") == facets
        assert stl.count("endloop") == facets
        assert len(_stl_vertices(stl)) == 3 * facets

    def test_the_vertices_are_itm_metres_not_degrees(self, topo, elevation):
        vertices = _stl_vertices(topo.getElevationSTL(elevation))
        xs = [v[0] for v in vertices]
        ys = [v[1] for v in vertices]
        assert 100000 < min(xs) and max(xs) < 300000
        assert 400000 < min(ys) and max(ys) < 800000

    def test_the_heights_span_the_elevation_field(self, topo, elevation):
        vertices = _stl_vertices(topo.getElevationSTL(elevation))
        zs = [v[2] for v in vertices]
        assert max(zs) == pytest.approx(float(elevation["elevation"].values.max()))
        # rasterToSTL puts the base ten metres below the lowest point.
        assert min(zs) == pytest.approx(float(elevation["elevation"].values.min()) - 10.0)

    def test_shifting_translates_every_vertex_in_the_plane(self, topo, elevation):
        plain = _stl_vertices(topo.getElevationSTL(elevation))
        shifted = _stl_vertices(topo.getElevationSTL(elevation, shiftx=1000.0, shifty=2000.0))
        assert len(plain) == len(shifted)
        for (x0, y0, z0), (x1, y1, z1) in zip(plain, shifted):
            assert x1 == pytest.approx(x0 - 1000.0)
            assert y1 == pytest.approx(y0 - 2000.0)
            assert z1 == pytest.approx(z0)

    def test_a_zero_shift_is_the_default(self, topo, elevation):
        assert topo.getElevationSTL(elevation) == \
            topo.getElevationSTL(elevation, shiftx=0, shifty=0)

    def test_it_is_deterministic(self, topo, elevation):
        assert topo.getElevationSTL(elevation) == topo.getElevationSTL(elevation)


@pytest.mark.unit
class TestCreateElevationSTL:
    """The one-shot composition of getElevation and getElevationSTL."""

    def test_it_equals_the_two_steps_applied_by_hand(self, topo, srtm):
        combined = topo.createElevationSTL(dxdy=100, dataSourceName=srtm, **_BOX)
        stepwise = topo.getElevationSTL(
            topo.getElevation(dxdy=100, dataSourceName=srtm, **_BOX)
        )
        assert combined == stepwise

    def test_it_forwards_the_solid_name(self, topo, srtm):
        stl = topo.createElevationSTL(dxdy=100, dataSourceName=srtm, solidName="Carmel", **_BOX)
        assert stl.startswith("solid Carmel\n")

    def test_it_forwards_the_shifts(self, topo, srtm):
        plain = _stl_vertices(topo.createElevationSTL(dxdy=100, dataSourceName=srtm, **_BOX))
        shifted = _stl_vertices(
            topo.createElevationSTL(dxdy=100, dataSourceName=srtm, shiftx=500.0, shifty=0,
                                    **_BOX)
        )
        for (x0, _, _), (x1, _, _) in zip(plain, shifted):
            assert x1 == pytest.approx(x0 - 500.0)

    def test_the_default_solid_name_is_topography(self, topo, srtm):
        assert topo.createElevationSTL(dxdy=100, dataSourceName=srtm,
                                       **_BOX).startswith("solid Topography\n")

    def test_a_coarser_resolution_gives_a_smaller_mesh(self, topo, srtm):
        fine = topo.createElevationSTL(dxdy=50, dataSourceName=srtm, **_BOX)
        coarse = topo.createElevationSTL(dxdy=100, dataSourceName=srtm, **_BOX)
        assert fine.count("facet normal") > coarse.count("facet normal")

    def test_a_projected_input_crs_is_forwarded_to_the_grid_builder(self, topo, srtm):
        """Reaching the else-branch of create_xarray through the STL entry point."""
        origin = convertCRS([[35.15625, 32.90625]], inputCRS=WSG84, outputCRS=ITM)[0]
        stl = topo.createElevationSTL(
            minx=origin.x, miny=origin.y, maxx=origin.x + 250, maxy=origin.y + 250,
            dxdy=100, inputCRS=ITM, dataSourceName=srtm,
        )
        assert stl.startswith("solid Topography\n")
        assert stl.count("facet normal") > 0


# --------------------------------------------------------------------------
# module-level sanity
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestModuleLevelState:
    def test_the_module_defines_wsg84_as_the_shared_epsg_code(self):
        assert topography_module.WSG84 == WSG84

    def test_gdal_is_absent_in_this_environment(self):
        """Documents why the gdal branches are left untested here."""
        assert topography_module._gdal_available is False

    def test_the_toolkit_exposes_its_analysis_layer(self, topo):
        assert isinstance(topo.analysis, topography_module.topographyAnalysis)
        assert topo.analysis.datalayer is topo

    def test_the_toolkit_reports_its_own_name(self, topo):
        assert topo.toolkitName == "TopographyToolkit"

    def test_the_bilinear_closed_form_is_a_real_check_not_a_tautology(self):
        """Guard on the test helper itself: it must reproduce the band exactly."""
        band = _ramp_band()
        for row in range(_SIZE):
            for col in range(_SIZE):
                lat = _NORTH - row * _RES
                lon = _WEST + col * _RES
                assert _expected_elevation(lat, lon) == pytest.approx(band[row, col])
        assert math.isclose(_expected_elevation(_NORTH - 1.5 * _RES, _WEST + 2.5 * _RES), 265.0)
