"""hera/measurements/GIS/raster/landcover.py -- the retrieval, roughness and
presentation layers that ``test_gis_landcover_extras.py`` and
``test_gis_landcover_roughness.py`` leave uncovered.

Covered here
------------
``getLandCoverAtPoint``, ``getLandCover`` (and, through it, its nested
``vectorizeLandCoverCalc``), ``getRoughnessAtPoint``,
``getRoughnessFromLandcover``, ``getRoughness``,
``_getUrbanRoughnessFromLandCover``, ``_getRoughnessFromBuildingsDataFrame``,
and the whole ``presentation`` class: ``datalayer``, ``plotRoughness``,
``plotLandcover``, ``plotLambdas``, ``_plotWithRectangles``,
``_adddRectanglesToPlot``, ``_getLandcoverRectangles``,
``_getRoughnessRectangles`` and ``_getLambdasRectangles``.

How the inputs were produced
----------------------------
Every raster is synthesised with ``rasterio`` into ``tmp_path`` by
``_write_geotiff`` below: a handful of pixels, an explicit ``from_origin``
transform and an explicit EPSG code, so the pixel a given lon/lat lands in is
known by construction rather than by golden file.  The real
``rasterio.open`` -> ``read(1)`` -> ``transform`` path in the production code
runs unmodified.  Rasters are reached both ways the code supports: as a bare
filesystem path passed straight as ``dataSourceName``, and as a datasource
registered on the mongomock-backed toolkit with ``dataFormat="string"`` (that
handler returns the resource verbatim, which is exactly the "datasource stores
a file path" branch).

Land-cover xarrays are built with the module's own
``hera.measurements.GIS.utils.create_xarray`` so the dim names (``i``, ``j``)
and coord names (``lat``, ``lon``, ``dxdy``) are the real ones, then given a
``landcover`` / ``z0`` / lambda coord as the method under test expects.

Deliberately NOT tested, and why
--------------------------------
* The ``osgeo.gdal.Dataset`` branch of ``getLandCoverAtPoint`` and the
  ``hasattr(ds, "GetRasterBand")`` branch of ``getLandCover``.  ``osgeo`` is
  not importable in the pinned environment (``ModuleNotFoundError: No module
  named 'osgeo'``), so the first branch is skipped by its own ``except
  ImportError`` and the second can only be reached with a real
  ``gdal.Dataset``.  Faking them would test the fake, not the code.
* ``_handleType1`` itself and ``roughnesslength2sandgrainroughness``'s value:
  already covered entry-by-entry in ``test_gis_landcover_roughness.py`` and
  ``test_gis_landcover_extras.py``.  Only the *duplication* of the latter is
  noted here.
* The real ``GIS_Buildings`` toolkit.  ``_getUrbanRoughnessFromLandCover``
  reaches it through the module-level ``toolkitHome`` name, which is patched
  as a seam; ``LambdaFromBuildingData`` has its own tests.
"""
import inspect
import math

import numpy as np
import pytest

from hera.measurements.GIS.raster import landcover as landcover_module
from hera.measurements.GIS.utils import BETA, ITM, KARMAN, WSG84, convertCRS, create_xarray

# --------------------------------------------------------------------------
# raster synthesis
# --------------------------------------------------------------------------

def _write_geotiff(path, array, west, north, resolution, crs="EPSG:4326"):
    """Write ``array`` as a one-band GeoTIFF and return its path.

    The transform is ``from_origin(west, north, resolution, resolution)``, so
    pixel (row, col) spans lon ``west + col * resolution`` and lat
    ``north - row * resolution``.  ``crs=None`` writes a raster with no CRS,
    which is a branch the production code handles explicitly.
    """
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


# A 4x4 raster of the codes 0..15, one code per pixel, 0.001 deg pixels whose
# top-left corner is lon 35.000 / lat 32.004.
_WEST = 35.0
_NORTH = 32.004
_RES = 0.001


def _codes_raster(tmp_path, name="landcover.tif", crs="EPSG:4326"):
    array = np.arange(16, dtype="int16").reshape(4, 4)
    return array, _write_geotiff(tmp_path / name, array, _WEST, _NORTH, _RES, crs=crs)


def _pixel_centre(row, col):
    """The (lon, lat) at the centre of pixel (row, col) of the codes raster."""
    return _WEST + (col + 0.5) * _RES, _NORTH - (row + 0.5) * _RES


def _expected_codes(array, lat_grid, lon_grid):
    """Independently index ``array`` from a grid's own lat/lon coordinates.

    Deliberately a re-implementation of the index arithmetic rather than a
    copy of the result: it pins the wiring (which coordinate drives which
    axis, and the orientation of the reshape) without being sensitive to the
    floating-point rounding of a grid point that lands exactly on a pixel
    boundary.
    """
    expected = np.empty(lat_grid.shape, dtype=array.dtype)
    for i in range(lat_grid.shape[0]):
        for j in range(lat_grid.shape[1]):
            row = math.floor((lat_grid[i, j] - _NORTH) / -_RES)
            col = math.floor((lon_grid[i, j] - _WEST) / _RES)
            expected[i, j] = array[row, col]
    return expected


# The box used for every grid: about 220 m by 190 m inside the codes raster.
_BOX = dict(minx=32.001, miny=35.001, maxx=32.003, maxy=35.003)


@pytest.fixture()
def lc(unit_toolkit_factory):
    from hera import toolkitHome

    return unit_toolkit_factory(toolkitHome.GIS_LANDCOVER)


@pytest.fixture()
def landcover_grid():
    """A real (i, j) grid from the module's own create_xarray, no data yet."""
    return create_xarray(dxdy=100, **_BOX)


def _with_coords(grid, **coords):
    return grid.assign_coords(
        {name: (("i", "j"), np.full(grid.shape, value, dtype=float if isinstance(value, float) else int))
         for name, value in coords.items()}
    )


# --------------------------------------------------------------------------
# getLandCoverAtPoint
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestGetLandCoverAtPoint:
    """The single-point raster lookup, driven against a synthesised GeoTIFF."""

    @pytest.mark.parametrize("row, col", [(0, 0), (0, 3), (2, 2), (3, 0), (3, 3)])
    def test_it_returns_the_code_of_the_pixel_containing_the_point(self, lc, tmp_path, row, col):
        array, path = _codes_raster(tmp_path)
        lon, lat = _pixel_centre(row, col)
        assert lc.getLandCoverAtPoint(lon=lon, lat=lat, dataSourceName=path) == array[row, col]

    def test_it_returns_a_plain_python_int(self, lc, tmp_path):
        _, path = _codes_raster(tmp_path)
        lon, lat = _pixel_centre(1, 1)
        value = lc.getLandCoverAtPoint(lon=lon, lat=lat, dataSourceName=path)
        assert type(value) is int

    def test_it_reads_a_registered_string_datasource_by_path(self, lc, tmp_path):
        array, path = _codes_raster(tmp_path)
        lc.addDataSource("LANDCOVER", path, "string", version=(0, 0, 1))
        lc.setDataSourceDefaultVersion("LANDCOVER", (0, 0, 1))
        lon, lat = _pixel_centre(2, 1)
        assert lc.getLandCoverAtPoint(lon=lon, lat=lat, dataSourceName="LANDCOVER") == array[2, 1]

    def test_it_falls_back_to_the_configured_default_landcover(self, lc, tmp_path):
        array, path = _codes_raster(tmp_path)
        lc.addDataSource("LANDCOVER", path, "string", version=(0, 0, 1))
        lc.setConfig(defaultLandCover="LANDCOVER")
        lon, lat = _pixel_centre(1, 2)
        assert lc.getLandCoverAtPoint(lon=lon, lat=lat) == array[1, 2]

    def test_without_a_configured_default_the_config_lookup_raises(self, lc):
        """No datasource and no defaultLandCover key: the bare dict lookup surfaces."""
        with pytest.raises(KeyError, match="defaultLandCover"):
            lc.getLandCoverAtPoint(lon=35.0015, lat=32.0015)

    def test_a_point_outside_the_raster_raises_index_error(self, lc, tmp_path):
        _, path = _codes_raster(tmp_path)
        with pytest.raises(IndexError, match="out of raster bounds"):
            lc.getLandCoverAtPoint(lon=35.02, lat=32.0015, dataSourceName=path)

    def test_a_point_north_of_the_raster_raises_index_error(self, lc, tmp_path):
        _, path = _codes_raster(tmp_path)
        with pytest.raises(IndexError, match="out of raster bounds"):
            lc.getLandCoverAtPoint(lon=35.0015, lat=32.05, dataSourceName=path)

    def test_an_unresolvable_datasource_raises_value_error(self, lc):
        with pytest.raises(ValueError, match="Could not load data source: NO_SUCH_SOURCE"):
            lc.getLandCoverAtPoint(lon=35.0015, lat=32.0015, dataSourceName="NO_SUCH_SOURCE")

    def test_itm_input_coordinates_are_reprojected_to_the_raster_crs(self, lc, tmp_path):
        array, path = _codes_raster(tmp_path)
        lon, lat = _pixel_centre(2, 2)
        itm = convertCRS([[lon, lat]], inputCRS=WSG84, outputCRS=ITM)[0]
        assert lc.getLandCoverAtPoint(lon=itm.x, lat=itm.y, inputCRS=ITM,
                                      dataSourceName=path) == array[2, 2]

    def test_a_raster_without_a_crs_is_read_without_reprojection(self, lc, tmp_path):
        array, path = _codes_raster(tmp_path, name="nocrs.tif", crs=None)
        lon, lat = _pixel_centre(2, 2)
        assert lc.getLandCoverAtPoint(lon=lon, lat=lat, dataSourceName=path) == array[2, 2]


# --------------------------------------------------------------------------
# getLandCover
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestGetLandCover:
    """The gridded lookup, and with it the nested vectorizeLandCoverCalc."""

    def test_it_returns_a_grid_with_a_landcover_coordinate(self, lc, tmp_path):
        _, path = _codes_raster(tmp_path)
        result = lc.getLandCover(dxdy=100, dataSourceName=path, **_BOX)
        assert result.dims == ("i", "j")
        assert "landcover" in result.coords
        assert result.landcover.dims == ("i", "j")
        assert result.landcover.shape == result.lat.shape

    def test_the_grid_keeps_the_lat_lon_and_dxdy_coordinates(self, lc, tmp_path):
        _, path = _codes_raster(tmp_path)
        result = lc.getLandCover(dxdy=100, dataSourceName=path, **_BOX)
        assert {"lat", "lon", "dxdy"} <= set(result.coords)
        assert float(result.dxdy.values) == pytest.approx(100.0)

    def test_every_cell_gets_the_code_of_the_pixel_it_falls_in(self, lc, tmp_path):
        """Checked against an independent re-derivation of the pixel index."""
        array, path = _codes_raster(tmp_path)
        result = lc.getLandCover(dxdy=100, dataSourceName=path, **_BOX)
        expected = _expected_codes(array, result.lat.values, result.lon.values)
        np.testing.assert_array_equal(result.landcover.values, expected)

    def test_a_uniform_raster_gives_that_code_everywhere(self, lc, tmp_path):
        path = _write_geotiff(tmp_path / "uniform.tif",
                              np.full((8, 8), 13, dtype="int16"), _WEST, _NORTH, _RES)
        result = lc.getLandCover(dxdy=100, dataSourceName=path, **_BOX)
        assert set(np.unique(result.landcover.values)) == {13}

    def test_it_accepts_a_registered_datasource_name(self, lc, tmp_path):
        array, path = _codes_raster(tmp_path)
        lc.addDataSource("LANDCOVER", path, "string", version=(0, 0, 1))
        lc.setDataSourceDefaultVersion("LANDCOVER", (0, 0, 1))
        by_name = lc.getLandCover(dxdy=100, dataSourceName="LANDCOVER", **_BOX)
        by_path = lc.getLandCover(dxdy=100, dataSourceName=path, **_BOX)
        np.testing.assert_array_equal(by_name.landcover.values, by_path.landcover.values)

    def test_it_falls_back_to_the_configured_default_landcover(self, lc, tmp_path):
        _, path = _codes_raster(tmp_path)
        lc.addDataSource("LANDCOVER", path, "string", version=(0, 0, 1))
        lc.setConfig(defaultLandCover="LANDCOVER")
        assert "landcover" in lc.getLandCover(dxdy=100, **_BOX).coords

    def test_a_coarser_dxdy_gives_a_smaller_grid(self, lc, tmp_path):
        _, path = _codes_raster(tmp_path)
        fine = lc.getLandCover(dxdy=50, dataSourceName=path, **_BOX)
        coarse = lc.getLandCover(dxdy=100, dataSourceName=path, **_BOX)
        assert fine.landcover.size > coarse.landcover.size


@pytest.mark.unit
class TestGetLandCoverIgnoresTheRasterCRS:
    """B185: getLandCover never reprojects, unlike getLandCoverAtPoint."""

    @staticmethod
    def _itm_raster(tmp_path):
        origin = convertCRS([[_WEST, 32.0]], inputCRS=WSG84, outputCRS=ITM)[0]
        array = np.full((20, 20), 13, dtype="int16")
        path = _write_geotiff(tmp_path / "itm.tif", array,
                              origin.x - 500, origin.y + 500, 100.0, crs="EPSG:2039")
        return array, path

    def test_at_point_handles_an_itm_raster(self, lc, tmp_path):
        """The single-point method does reproject, so it is the reference."""
        array, path = self._itm_raster(tmp_path)
        assert lc.getLandCoverAtPoint(lon=35.001, lat=32.001, dataSourceName=path) == 13

    def test_the_grid_lookup_raises_on_an_itm_raster(self, lc, tmp_path):
        """Characterisation of B185."""
        _, path = self._itm_raster(tmp_path)
        with pytest.raises(IndexError, match="out of bounds for axis 0"):
            lc.getLandCover(dxdy=100, dataSourceName=path, **_BOX)

    @pytest.mark.xfail(
        strict=True,
        reason="B185: getLandCover indexes the raster with the grid's WGS84 "
               "degrees directly -- it reads src.transform but never compares "
               "src.crs to the input CRS, unlike getLandCoverAtPoint, which "
               "builds a pyproj Transformer when they differ. On a landcover "
               "raster in any projected CRS (ITM is the project's own second "
               "constant) degrees-as-metres produce indices in the thousands "
               "and the call dies with IndexError; a raster whose projected "
               "extent happened to contain small index values would instead "
               "return silently wrong classes. See the consolidated findings "
               "issue.",
    )
    def test_the_grid_lookup_agrees_with_the_point_lookup_on_an_itm_raster(self, lc, tmp_path):
        _, path = self._itm_raster(tmp_path)
        grid = lc.getLandCover(dxdy=100, dataSourceName=path, **_BOX)
        assert set(np.unique(grid.landcover.values)) == {13}


# --------------------------------------------------------------------------
# getRoughnessAtPoint
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestGetRoughnessAtPoint:
    """Dispatch from a datasource's declared type to a roughness handler."""

    # The table written inline in getRoughnessAtPoint, reached when the
    # datasource declares no type. Pinned as B43 in test_gis_landcover_roughness.
    EXAMPLE_RAMP = {0: 0.01, 5: 0.3, 13: 1.1, 15: 1.3}

    def test_a_datasource_typed_1_routes_to_the_published_table(self, lc, tmp_path):
        array, path = _codes_raster(tmp_path)
        lc.addDataSource("TYPED", path, "string", version=(0, 0, 1), type="1")
        for row, col in [(0, 0), (0, 1), (3, 3)]:
            lon, lat = _pixel_centre(row, col)
            assert lc.getRoughnessAtPoint(lon=lon, lat=lat, dataSourceName="TYPED") == \
                pytest.approx(lc._handleType1(int(array[row, col])))

    def test_it_returns_a_float(self, lc, tmp_path):
        _, path = _codes_raster(tmp_path)
        lc.addDataSource("TYPED", path, "string", version=(0, 0, 1), type="1")
        lon, lat = _pixel_centre(0, 0)
        assert isinstance(lc.getRoughnessAtPoint(lon=lon, lat=lat, dataSourceName="TYPED"), float)

    @pytest.mark.parametrize("code", sorted(EXAMPLE_RAMP))
    def test_an_untyped_datasource_uses_the_inline_example_ramp(self, lc, tmp_path, code):
        """Characterisation of B43 reached through the point API."""
        path = _write_geotiff(tmp_path / f"code{code}.tif",
                              np.full((4, 4), code, dtype="int16"), _WEST, _NORTH, _RES)
        lon, lat = _pixel_centre(1, 1)
        assert lc.getRoughnessAtPoint(lon=lon, lat=lat, dataSourceName=path) == \
            pytest.approx(self.EXAMPLE_RAMP[code])

    def test_an_out_of_table_code_uses_the_documented_fallback(self, lc, tmp_path):
        path = _write_geotiff(tmp_path / "code99.tif",
                              np.full((4, 4), 99, dtype="int16"), _WEST, _NORTH, _RES)
        lon, lat = _pixel_centre(1, 1)
        assert lc.getRoughnessAtPoint(lon=lon, lat=lat, dataSourceName=path) == pytest.approx(0.05)

    def test_an_unknown_declared_type_raises_attribute_error(self, lc, tmp_path):
        _, path = _codes_raster(tmp_path)
        lc.addDataSource("ODD", path, "string", version=(0, 0, 1), type="Zzz")
        lon, lat = _pixel_centre(0, 0)
        with pytest.raises(AttributeError, match="_handleTypeZzz"):
            lc.getRoughnessAtPoint(lon=lon, lat=lat, dataSourceName="ODD")


# --------------------------------------------------------------------------
# getRoughnessFromLandcover
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestGetRoughnessFromLandcover:
    """Adding a z0 coordinate to an existing landcover grid."""

    # The third IGBP table in the module: the inline lambda fallback.
    INLINE_FALLBACK = {
        0: 0.01, 1: 0.02, 2: 0.05, 3: 0.1, 4: 0.15, 5: 0.2, 6: 0.25, 7: 0.3,
        8: 0.35, 9: 0.4, 10: 0.45, 11: 0.5, 12: 0.55, 13: 0.6, 14: 0.01,
        15: 0.001, 16: 0.0001,
    }

    def test_a_typed_datasource_fills_z0_from_the_published_table(self, lc, tmp_path, landcover_grid):
        _, path = _codes_raster(tmp_path)
        lc.addDataSource("TYPED", path, "string", version=(0, 0, 1), type="1")
        grid = _with_coords(landcover_grid, landcover=9)
        result = lc.getRoughnessFromLandcover(grid, dataSourceName="TYPED")
        assert result.z0.dims == ("i", "j")
        np.testing.assert_allclose(result.z0.values, lc._handleType1(9))

    def test_z0_is_added_without_disturbing_the_landcover_coordinate(self, lc, tmp_path, landcover_grid):
        _, path = _codes_raster(tmp_path)
        lc.addDataSource("TYPED", path, "string", version=(0, 0, 1), type="1")
        grid = _with_coords(landcover_grid, landcover=10)
        result = lc.getRoughnessFromLandcover(grid, dataSourceName="TYPED")
        np.testing.assert_array_equal(result.landcover.values, grid.landcover.values)
        np.testing.assert_array_equal(result.lat.values, grid.lat.values)

    def test_a_datasource_document_without_a_type_raises_attribute_error(self, lc, tmp_path, landcover_grid):
        _, path = _codes_raster(tmp_path)
        lc.addDataSource("UNTYPED", path, "string", version=(0, 0, 1))
        grid = _with_coords(landcover_grid, landcover=3)
        with pytest.raises(AttributeError, match="_handleTypeNone"):
            lc.getRoughnessFromLandcover(grid, dataSourceName="UNTYPED")

    def test_an_unknown_declared_type_raises_attribute_error(self, lc, tmp_path, landcover_grid):
        _, path = _codes_raster(tmp_path)
        lc.addDataSource("ODD", path, "string", version=(0, 0, 1), type="Zzz")
        grid = _with_coords(landcover_grid, landcover=3)
        with pytest.raises(AttributeError, match="_handleTypeZzz"):
            lc.getRoughnessFromLandcover(grid, dataSourceName="ODD")

    def test_urban_mode_demands_a_wind_direction_and_a_resolution(self, lc, landcover_grid):
        grid = _with_coords(landcover_grid, landcover=13)
        with pytest.raises(ValueError, match="windMeteorologicalDirection and resolution"):
            lc.getRoughnessFromLandcover(grid, isBuilding=True, dataSourceName="ANY")
        with pytest.raises(ValueError, match="windMeteorologicalDirection and resolution"):
            lc.getRoughnessFromLandcover(grid, windMeteorologicalDirection=270.0,
                                         isBuilding=True, dataSourceName="ANY")

    @pytest.mark.parametrize("code", sorted(INLINE_FALLBACK))
    def test_an_unregistered_datasource_uses_the_inline_fallback_table(self, lc, landcover_grid, code):
        """Characterisation of B184: the third IGBP table in the module."""
        grid = _with_coords(landcover_grid, landcover=code)
        result = lc.getRoughnessFromLandcover(grid, dataSourceName="NO_SUCH_SOURCE")
        np.testing.assert_allclose(result.z0.values, self.INLINE_FALLBACK[code])

    def test_the_inline_fallback_is_used_for_codes_outside_it_too(self, lc, landcover_grid):
        """Characterisation of B184: its own out-of-table default is 0.1."""
        grid = _with_coords(landcover_grid, landcover=99)
        result = lc.getRoughnessFromLandcover(grid, dataSourceName="NO_SUCH_SOURCE")
        np.testing.assert_allclose(result.z0.values, 0.1)

    def test_the_inline_fallback_is_not_monotone_where_the_published_table_is(self, lc):
        """Characterisation of B184: it ramps 0..13 then collapses at 14..16."""
        assert self.INLINE_FALLBACK[13] > self.INLINE_FALLBACK[12] > self.INLINE_FALLBACK[11]
        assert self.INLINE_FALLBACK[14] < self.INLINE_FALLBACK[1]
        assert self.INLINE_FALLBACK[16] == pytest.approx(0.0001)

    def test_snow_and_ice_is_the_only_class_the_two_tables_agree_on(self, lc):
        """Characterisation of B184: 15 out of 17 classes disagree."""
        agree = [code for code in self.INLINE_FALLBACK
                 if self.INLINE_FALLBACK[code] == pytest.approx(lc._handleType1(code))]
        assert agree == [15]

    @pytest.mark.xfail(
        strict=True,
        reason="B184: getRoughnessFromLandcover carries a THIRD IGBP z0 table, "
               "an inline lambda used whenever the datasource document is "
               "missing (the documented 'default type: IGBP' path). It agrees "
               "with neither _handleType1's published Floors et al. (2021) "
               "table nor the ramp in getRoughnessAtPoint already pinned as "
               "B43: water is 0.01 against a published 0.0001 (x100), "
               "croplands 0.55 against 0.15, and barren 0.0001 against 0.01 "
               "(x1/100); only snow and ice agrees. z0 feeds the wind profile "
               "and thus every dispersion result, so which of the three tables "
               "a caller happens to reach materially changes the answer. See "
               "the consolidated findings issue.",
    )
    def test_the_fallback_agrees_with_the_published_table(self, lc):
        for code in sorted(self.INLINE_FALLBACK):
            assert self.INLINE_FALLBACK[code] == pytest.approx(lc._handleType1(code)), code


# --------------------------------------------------------------------------
# getRoughness
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestGetRoughness:
    """The convenience composition of getLandCover and getRoughnessFromLandcover."""

    def test_it_equals_the_two_steps_applied_by_hand(self, lc, tmp_path):
        _, path = _codes_raster(tmp_path)
        lc.addDataSource("TYPED", path, "string", version=(0, 0, 1), type="1")
        combined = lc.getRoughness(dxdy=100, dataSourceName="TYPED", **_BOX)
        stepwise = lc.getRoughnessFromLandcover(
            lc.getLandCover(dxdy=100, dataSourceName="TYPED", **_BOX),
            dataSourceName="TYPED",
        )
        np.testing.assert_array_equal(combined.landcover.values, stepwise.landcover.values)
        np.testing.assert_allclose(combined.z0.values, stepwise.z0.values)

    def test_it_carries_both_the_landcover_and_the_z0_coordinates(self, lc, tmp_path):
        _, path = _codes_raster(tmp_path)
        lc.addDataSource("TYPED", path, "string", version=(0, 0, 1), type="1")
        result = lc.getRoughness(dxdy=100, dataSourceName="TYPED", **_BOX)
        assert {"landcover", "z0"} <= set(result.coords)

    def test_it_covers_the_same_cells_as_getLandCover_alone(self, lc, tmp_path):
        _, path = _codes_raster(tmp_path)
        lc.addDataSource("TYPED", path, "string", version=(0, 0, 1), type="1")
        result = lc.getRoughness(dxdy=100, dataSourceName="TYPED", **_BOX)
        plain = lc.getLandCover(dxdy=100, dataSourceName="TYPED", **_BOX)
        assert result.z0.shape == plain.landcover.shape
        assert result.z0.dims == ("i", "j")

    def test_urban_mode_is_forwarded_through_to_the_urban_branch(self, lc, tmp_path):
        """getRoughness must pass isBuilding on, not swallow it."""
        _, path = _codes_raster(tmp_path)
        lc.addDataSource("TYPED", path, "string", version=(0, 0, 1), type="1")
        with pytest.raises(ValueError, match="windMeteorologicalDirection and resolution"):
            lc.getRoughness(dxdy=100, dataSourceName="TYPED", isBuilding=True, **_BOX)


@pytest.mark.unit
class TestRoughnessDtypeFromTheHandler:
    """B187: np.vectorize with no otypes infers the dtype from cell (0, 0)."""

    # The five forest classes are the only ones _handleType1 spells as the
    # Python int 1 rather than a float; every other class is a float literal.
    FOREST = 1
    GRASSLAND = 10

    @pytest.fixture()
    def typed_source(self, lc, tmp_path):
        _, path = _codes_raster(tmp_path, name="typed.tif")
        lc.addDataSource("TYPED", path, "string", version=(0, 0, 1), type="1")
        return "TYPED"

    @staticmethod
    def _grid_starting_with(grid, first, rest):
        codes = np.full(grid.shape, rest, dtype=int)
        codes[0, 0] = first
        return grid.assign_coords(landcover=(("i", "j"), codes))

    def test_a_grid_of_only_float_valued_classes_keeps_its_precision(self, lc, landcover_grid, typed_source):
        """The control: with no forest cell the handler returns floats throughout."""
        grid = self._grid_starting_with(landcover_grid, self.GRASSLAND, self.GRASSLAND)
        result = lc.getRoughnessFromLandcover(grid, dataSourceName=typed_source)
        assert result.z0.dtype.kind == "f"
        np.testing.assert_allclose(result.z0.values, 0.12)

    def test_a_leading_forest_cell_makes_the_whole_field_integer(self, lc, landcover_grid, typed_source):
        """Characterisation of B187."""
        grid = self._grid_starting_with(landcover_grid, self.FOREST, self.GRASSLAND)
        result = lc.getRoughnessFromLandcover(grid, dataSourceName=typed_source)
        assert result.z0.dtype.kind == "i"

    def test_the_truncation_zeroes_every_sub_metre_roughness(self, lc, landcover_grid, typed_source):
        """Characterisation of B187: 0.12 m grassland becomes exactly 0."""
        grid = self._grid_starting_with(landcover_grid, self.FOREST, self.GRASSLAND)
        result = lc.getRoughnessFromLandcover(grid, dataSourceName=typed_source)
        assert result.z0.values[0, 0] == 1
        assert set(np.unique(result.z0.values[result.z0.values != 1])) == {0}

    def test_it_also_zeroes_urban_roughness(self, lc, landcover_grid, typed_source):
        """Characterisation of B187: urban 0.8 m -> 0, i.e. a smooth city."""
        grid = self._grid_starting_with(landcover_grid, self.FOREST, 13)
        result = lc.getRoughnessFromLandcover(grid, dataSourceName=typed_source)
        assert result.z0.values[0, 1] == 0

    @pytest.mark.xfail(
        strict=True,
        reason="B187: getRoughnessFromLandcover computes z0 with "
               "np.vectorize(handlerFunction) and passes no otypes, so numpy "
               "infers the output dtype by calling the handler on the FIRST "
               "element only. _handleType1's table mixes types -- the five "
               "forest classes (1..5) are the Python int 1, every other class "
               "is a float -- so whenever grid cell (0, 0) is a forest, numpy "
               "picks int64 and silently truncates every other cell: "
               "grassland 0.12 -> 0, urban 0.8 -> 0, water 0.0001 -> 0. A "
               "zero z0 makes the log wind profile singular. Fixing either the "
               "table's literals (1 -> 1.0) or the vectorize call "
               "(otypes=[float]) would do. See the consolidated findings "
               "issue.",
    )
    def test_a_leading_forest_cell_does_not_change_the_other_cells(self, lc, landcover_grid, typed_source):
        grid = self._grid_starting_with(landcover_grid, self.FOREST, self.GRASSLAND)
        result = lc.getRoughnessFromLandcover(grid, dataSourceName=typed_source)
        assert result.z0.values[0, 1] == pytest.approx(lc._handleType1(self.GRASSLAND))


# --------------------------------------------------------------------------
# _getRoughnessFromBuildingsDataFrame
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestRoughnessFromBuildingsDataFrame:
    """The morphometric z0 / displacement height arithmetic, clamps included."""

    @staticmethod
    def _grid(rows):
        import geopandas as gpd
        from shapely.geometry import box

        return gpd.GeoDataFrame(
            {
                "lambdaF": [r[0] for r in rows],
                "lambdaP": [r[1] for r in rows],
                "hc": [r[2] for r in rows],
            },
            geometry=[box(i, i, i + 1, i + 1) for i in range(len(rows))],
            crs=ITM,
        )

    def test_it_clamps_short_canopies_to_a_two_metre_default_cell(self, lc):
        result = lc._getRoughnessFromBuildingsDataFrame(self._grid([(0.1, 0.2, 1.0)]))
        assert result["hc"].iloc[0] == pytest.approx(2.0)
        assert result["lambdaF"].iloc[0] == pytest.approx(0.25)
        assert result["lambdaP"].iloc[0] == pytest.approx(0.25)

    def test_it_caps_the_frontal_and_plan_area_densities(self, lc):
        result = lc._getRoughnessFromBuildingsDataFrame(self._grid([(0.9, 0.95, 10.0)]))
        assert result["lambdaF"].iloc[0] == pytest.approx(0.4)
        assert result["lambdaP"].iloc[0] == pytest.approx(0.6)

    def test_it_leaves_in_range_values_alone(self, lc):
        result = lc._getRoughnessFromBuildingsDataFrame(self._grid([(0.2, 0.3, 12.0)]))
        assert result["lambdaF"].iloc[0] == pytest.approx(0.2)
        assert result["lambdaP"].iloc[0] == pytest.approx(0.3)
        assert result["hc"].iloc[0] == pytest.approx(12.0)

    def test_the_canopy_length_scale_follows_its_definition(self, lc):
        result = lc._getRoughnessFromBuildingsDataFrame(self._grid([(0.2, 0.3, 12.0)]))
        assert result["Lc"].iloc[0] == pytest.approx(12.0 * (1 - 0.3) / 0.2)

    def test_the_mixing_length_scales_with_beta_cubed(self, lc):
        result = lc._getRoughnessFromBuildingsDataFrame(self._grid([(0.2, 0.3, 12.0)]))
        expected_lc = 12.0 * (1 - 0.3) / 0.2
        assert result["ll"].iloc[0] == pytest.approx(2 * BETA ** 3 * expected_lc)

    def test_z0_and_the_displacement_height_follow_the_karman_relations(self, lc):
        result = lc._getRoughnessFromBuildingsDataFrame(self._grid([(0.2, 0.3, 12.0)]))
        mixing = result["ll"].iloc[0]
        assert result["dd"].iloc[0] == pytest.approx(mixing / KARMAN)
        assert result["zz0"].iloc[0] == pytest.approx(mixing / KARMAN * math.exp(-KARMAN / BETA))

    def test_z0_is_always_a_fixed_fraction_of_the_displacement_height(self, lc):
        result = lc._getRoughnessFromBuildingsDataFrame(
            self._grid([(0.2, 0.3, 12.0), (0.5, 0.7, 10.0), (0.1, 0.1, 30.0)])
        )
        ratio = result["zz0"] / result["dd"]
        np.testing.assert_allclose(ratio.values, math.exp(-KARMAN / BETA))

    def test_a_taller_canopy_is_rougher(self, lc):
        result = lc._getRoughnessFromBuildingsDataFrame(
            self._grid([(0.2, 0.3, 5.0), (0.2, 0.3, 20.0)])
        )
        assert result["zz0"].iloc[1] > result["zz0"].iloc[0]

    def test_it_adds_all_four_derived_columns(self, lc):
        result = lc._getRoughnessFromBuildingsDataFrame(self._grid([(0.2, 0.3, 12.0)]))
        assert {"Lc", "ll", "zz0", "dd"} <= set(result.columns)

    def test_it_edits_the_callers_frame_in_place(self, lc):
        """Characterisation: it mutates and returns the same object, so a
        caller's lambdaF/lambdaP/hc are clamped behind its back."""
        given = self._grid([(0.9, 0.95, 1.0)])
        result = lc._getRoughnessFromBuildingsDataFrame(given)
        assert result is given
        assert given["hc"].iloc[0] == pytest.approx(2.0)


# --------------------------------------------------------------------------
# _getUrbanRoughnessFromLandCover
# --------------------------------------------------------------------------

def _fake_buildings_home(lambda_grid, buildings=None, calls=None):
    """A stand-in for the module-level ``toolkitHome`` name.

    Only the two members ``_getUrbanRoughnessFromLandCover`` actually reads
    are provided: ``getBuildingsFromRectangle`` and
    ``analysis.LambdaFromBuildingData``.
    """
    import geopandas as gpd
    from shapely.geometry import box

    if buildings is None:
        buildings = gpd.GeoDataFrame({"id": [0]}, geometry=[box(0, 0, 1, 1)], crs=ITM)
    record = calls if calls is not None else {}

    class FakeAnalysis:
        def LambdaFromBuildingData(self, windMeteorologicalDirection, resolution, buildingsFrame):
            record["lambda"] = (windMeteorologicalDirection, resolution, len(buildingsFrame))
            return lambda_grid.copy()

    class FakeBuildingsToolkit:
        analysis = FakeAnalysis()

        def getBuildingsFromRectangle(self, **kwargs):
            record["rectangle"] = kwargs
            return buildings

    class FakeToolkitHome:
        GIS_BUILDINGS = "GIS_Buildings"

        def getToolkit(self, toolkitName, projectName):
            record["getToolkit"] = (toolkitName, projectName)
            return FakeBuildingsToolkit()

    return FakeToolkitHome(), record


def _covering_lambda_grid(grid, values=(0.2, 0.3, 12.0)):
    """A one-cell lambda GeoDataFrame whose polygon covers the whole grid."""
    import geopandas as gpd
    from shapely.geometry import box

    corners = [
        convertCRS([[float(lon), float(lat)]], inputCRS=WSG84, outputCRS=ITM)[0]
        for lon, lat in zip(grid.lon.values.flatten(), grid.lat.values.flatten())
    ]
    xs = [p.x for p in corners]
    ys = [p.y for p in corners]
    return gpd.GeoDataFrame(
        {"lambdaF": [values[0]], "lambdaP": [values[1]], "hc": [values[2]]},
        geometry=[box(min(xs) - 500, min(ys) - 500, max(xs) + 500, max(ys) + 500)],
        crs=ITM,
    )


@pytest.mark.unit
class TestUrbanRoughnessFromLandCover:
    """The urban branch, with the GIS_Buildings toolkit patched at the module seam."""

    @pytest.fixture()
    def urban_grid(self, lc, tmp_path):
        _, path = _codes_raster(tmp_path, name="urban.tif")
        lc.addDataSource("URBAN", path, "string", version=(0, 0, 1), type="1")
        return lc.getLandCover(dxdy=100, dataSourceName="URBAN", **_BOX)

    def test_it_asks_the_buildings_toolkit_for_the_projects_own_project(self, lc, monkeypatch, urban_grid):
        home, record = _fake_buildings_home(_covering_lambda_grid(urban_grid))
        monkeypatch.setattr(landcover_module, "toolkitHome", home)
        lc.getRoughnessFromLandcover(urban_grid, windMeteorologicalDirection=270.0,
                                     resolution=10.0, isBuilding=True, dataSourceName="URBAN")
        assert record["getToolkit"] == ("GIS_Buildings", lc.projectName)

    def test_it_requests_the_buildings_rectangle_in_itm(self, lc, monkeypatch, urban_grid):
        home, record = _fake_buildings_home(_covering_lambda_grid(urban_grid))
        monkeypatch.setattr(landcover_module, "toolkitHome", home)
        lc.getRoughnessFromLandcover(urban_grid, windMeteorologicalDirection=270.0,
                                     resolution=10.0, isBuilding=True, dataSourceName="URBAN",
                                     GIS_BUILDINGS_dataSourceName="BUILDINGS")
        rectangle = record["rectangle"]
        assert rectangle["inputCRS"] == ITM
        assert rectangle["dataSourceName"] == "BUILDINGS"
        assert rectangle["minx"] < rectangle["maxx"]
        assert rectangle["miny"] < rectangle["maxy"]
        # ITM eastings/northings are six-figure metres, not degrees.
        assert 100000 < rectangle["minx"] < 300000
        assert 400000 < rectangle["miny"] < 800000

    def test_it_forwards_the_wind_direction_and_resolution_unchanged(self, lc, monkeypatch, urban_grid):
        home, record = _fake_buildings_home(_covering_lambda_grid(urban_grid))
        monkeypatch.setattr(landcover_module, "toolkitHome", home)
        lc.getRoughnessFromLandcover(urban_grid, windMeteorologicalDirection=123.5,
                                     resolution=17.0, isBuilding=True, dataSourceName="URBAN")
        assert record["lambda"][:2] == (123.5, 17.0)

    def test_an_empty_buildings_frame_raises_value_error(self, lc, monkeypatch, urban_grid):
        import geopandas as gpd

        empty = gpd.GeoDataFrame({"id": []}, geometry=[], crs=ITM)
        home, _ = _fake_buildings_home(_covering_lambda_grid(urban_grid), buildings=empty)
        monkeypatch.setattr(landcover_module, "toolkitHome", home)
        with pytest.raises(ValueError, match="Buildings DataFrame .* is empty"):
            lc.getRoughnessFromLandcover(urban_grid, windMeteorologicalDirection=270.0,
                                         resolution=10.0, isBuilding=True, dataSourceName="URBAN")

    def test_it_adds_the_six_urban_fields(self, lc, monkeypatch, urban_grid):
        home, _ = _fake_buildings_home(_covering_lambda_grid(urban_grid))
        monkeypatch.setattr(landcover_module, "toolkitHome", home)
        result = lc.getRoughnessFromLandcover(urban_grid, windMeteorologicalDirection=270.0,
                                              resolution=10.0, isBuilding=True,
                                              dataSourceName="URBAN")
        assert {"z0", "dd", "lambdaF", "lambdaP", "hc", "ll"} <= set(result.coords)

    def test_a_cell_inside_a_lambda_polygon_takes_that_polygons_values(self, lc, monkeypatch, urban_grid):
        lambda_grid = _covering_lambda_grid(urban_grid, values=(0.2, 0.3, 12.0))
        home, _ = _fake_buildings_home(lambda_grid)
        monkeypatch.setattr(landcover_module, "toolkitHome", home)
        result = lc.getRoughnessFromLandcover(urban_grid, windMeteorologicalDirection=270.0,
                                              resolution=10.0, isBuilding=True,
                                              dataSourceName="URBAN")
        expected = lc._getRoughnessFromBuildingsDataFrame(lambda_grid.copy())
        np.testing.assert_allclose(result["z0"].values, expected["zz0"].iloc[0])
        np.testing.assert_allclose(result["dd"].values, expected["dd"].iloc[0])
        np.testing.assert_allclose(result["ll"].values, expected["ll"].iloc[0])
        np.testing.assert_allclose(result["hc"].values, 12.0)
        np.testing.assert_allclose(result["lambdaF"].values, 0.2)
        np.testing.assert_allclose(result["lambdaP"].values, 0.3)


@pytest.mark.unit
class TestUrbanRoughnessOutsideTheBuildings:
    """B183: cells with no buildings get a land-cover CODE in z0."""

    @pytest.fixture()
    def setup(self, lc, tmp_path, monkeypatch):
        """A grid whose cells fall in no lambda polygon at all."""
        import geopandas as gpd
        from shapely.geometry import box

        _, path = _codes_raster(tmp_path, name="urban.tif")
        lc.addDataSource("URBAN", path, "string", version=(0, 0, 1), type="1")
        grid = lc.getLandCover(dxdy=100, dataSourceName="URBAN", **_BOX)
        far_away = gpd.GeoDataFrame(
            {"lambdaF": [0.2], "lambdaP": [0.3], "hc": [12.0]},
            geometry=[box(0.0, 0.0, 1.0, 1.0)],
            crs=ITM,
        )
        home, _ = _fake_buildings_home(far_away)
        monkeypatch.setattr(landcover_module, "toolkitHome", home)
        result = lc.getRoughnessFromLandcover(grid, windMeteorologicalDirection=270.0,
                                              resolution=10.0, isBuilding=True,
                                              dataSourceName="URBAN")
        return grid, result

    def test_the_other_urban_fields_are_zeroed(self, setup):
        _, result = setup
        for field in ("dd", "lambdaF", "lambdaP", "hc", "ll"):
            np.testing.assert_allclose(result[field].values, 0.0)

    def test_z0_is_the_land_cover_class_code(self, setup):
        """Characterisation of B183."""
        grid, result = setup
        np.testing.assert_allclose(result["z0"].values, grid.landcover.values.astype(float))

    def test_z0_exceeds_any_physical_roughness_length(self, setup):
        """Characterisation of B183: class codes run to 16, z0 does not."""
        _, result = setup
        assert result["z0"].values.max() > 3.0

    @pytest.mark.xfail(
        strict=True,
        reason="B183: in _getUrbanRoughnessFromLandCover, a grid cell that "
               "intersects no lambda polygon is filled with "
               "self.getLandCoverAtPoint(...) -- the IGBP CLASS CODE (0..16) -- "
               "assigned straight into the z0 field, with no roughness lookup. "
               "Every other assignment in the same block writes a real z0 "
               "(lambdas['zz0']). So the rural cells of an urban domain get z0 "
               "= 12 m or 14 m instead of 0.15 m, an error of two orders of "
               "magnitude in the field that sets the wind profile. The fix is "
               "presumably getRoughnessAtPoint / _handleType1 rather than "
               "getLandCoverAtPoint. See the consolidated findings issue.",
    )
    def test_z0_outside_the_buildings_is_a_roughness_length(self, setup):
        _, result = setup
        assert result["z0"].values.max() <= 3.0
        assert result["z0"].values.min() > 0.0


# --------------------------------------------------------------------------
# presentation
# --------------------------------------------------------------------------

@pytest.fixture()
def axes_image():
    """A real matplotlib AxesImage, the 'plot' argument the layer expects."""
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots()
    return axes, axes.imshow(np.zeros((4, 4)), extent=(10.0, 20.0, 30.0, 40.0))


@pytest.mark.unit
class TestPresentationBasics:
    def test_the_datalayer_property_returns_the_toolkit(self, lc):
        assert lc.presentation.datalayer is lc

    def test_the_colour_map_covers_every_igbp_class(self, lc):
        assert set(lc.presentation.landcover_colors_map) == set(range(17))

    def test_water_and_urban_get_distinguishable_colours(self, lc):
        colors = lc.presentation.landcover_colors_map
        assert colors[0] == "blue"
        assert colors[13] == "#2f2f2f"

    def test_a_freshly_built_presentation_wraps_the_datalayer_it_is_given(self, lc):
        other = landcover_module.presentation(dataLayer=lc)
        assert other.datalayer is lc
        assert other.landcover_colors_map == lc.presentation.landcover_colors_map


@pytest.mark.unit
class TestPresentationRectangles:
    def test_adding_rectangles_puts_one_patch_per_spec_on_the_axes(self, lc, axes_image):
        axes, _ = axes_image
        rectangles = [(0.0, 0.0, 1.0, 2.0, "red"), (5.0, 5.0, 1.0, 1.0, "blue")]
        returned = lc.presentation._adddRectanglesToPlot(axes, rectangles, 0.3)
        assert returned is axes
        assert len(axes.patches) == len(rectangles)

    def test_each_patch_keeps_its_geometry_colour_and_alpha(self, lc, axes_image):
        axes, _ = axes_image
        lc.presentation._adddRectanglesToPlot(axes, [(1.0, 2.0, 3.0, 4.0, "red")], 0.25)
        patch = axes.patches[0]
        assert patch.get_xy() == (1.0, 2.0)
        assert patch.get_width() == pytest.approx(3.0)
        assert patch.get_height() == pytest.approx(4.0)
        assert patch.get_alpha() == pytest.approx(0.25)

    def test_an_empty_rectangle_list_adds_nothing(self, lc, axes_image):
        axes, _ = axes_image
        lc.presentation._adddRectanglesToPlot(axes, [], 0.5)
        assert len(axes.patches) == 0

    def test_plotting_with_rectangles_adopts_the_images_extent(self, lc, axes_image):
        axes, image = axes_image
        returned = lc.presentation._plotWithRectangles(
            axes, image, [(11.0, 31.0, 1.0, 1.0, "red")], 0.5
        )
        assert returned is axes
        assert axes.get_xlim() == pytest.approx((10.0, 20.0))
        assert axes.get_ylim() == pytest.approx((30.0, 40.0))
        assert len(axes.patches) == 1


@pytest.mark.unit
class TestPresentationRectangleBuilders:
    """The three grid -> rectangle converters, on a real (i, j) grid."""

    @pytest.fixture()
    def grid(self, landcover_grid):
        return landcover_grid.assign_coords(
            landcover=(("i", "j"), np.full(landcover_grid.shape, 13, dtype=int)),
            z0=(("i", "j"), np.linspace(0.1, 0.9, landcover_grid.size).reshape(landcover_grid.shape)),
            lambdaF=(("i", "j"), np.linspace(0.05, 0.4, landcover_grid.size).reshape(landcover_grid.shape)),
        )

    def test_landcover_rectangles_are_itm_metres_sized_by_dxdy(self, lc, grid):
        rectangles = lc.presentation._getLandcoverRectangles(grid)
        assert len(rectangles) == grid.size
        for x, y, width, height, color in rectangles:
            assert 100000 < x < 300000
            assert 400000 < y < 800000
            assert (width, height) == (pytest.approx(100.0), pytest.approx(100.0))
            assert color == lc.presentation.landcover_colors_map[13]

    def test_landcover_rectangles_use_the_class_colour(self, lc, landcover_grid):
        grid = landcover_grid.assign_coords(
            landcover=(("i", "j"), np.full(landcover_grid.shape, 0, dtype=int))
        )
        assert {r[4] for r in lc.presentation._getLandcoverRectangles(grid)} == {"blue"}

    def test_roughness_rectangles_are_coloured_by_z0(self, lc, grid):
        rectangles = lc.presentation._getRoughnessRectangles(grid)
        assert len(rectangles) == grid.size
        # z0 spans a range, so the viridis colours must not all collapse to one.
        assert len({r[4] for r in rectangles}) > 1
        for _, _, _, _, color in rectangles:
            assert len(color) == 4  # RGBA

    def test_lambda_rectangles_are_coloured_by_the_named_field(self, lc, grid):
        rectangles = lc.presentation._getLambdasRectangles("lambdaF", grid)
        assert len(rectangles) == grid.size
        assert len({r[4] for r in rectangles}) > 1

    def test_a_constant_field_gives_one_colour(self, lc, landcover_grid):
        grid = landcover_grid.assign_coords(
            lambdaP=(("i", "j"), np.full(landcover_grid.shape, 0.3))
        )
        assert len({r[4] for r in lc.presentation._getLambdasRectangles("lambdaP", grid)}) == 1

    def test_an_unknown_field_raises_a_key_error(self, lc, grid):
        with pytest.raises(KeyError, match="noSuchField"):
            lc.presentation._getLambdasRectangles("noSuchField", grid)


@pytest.mark.unit
class TestPresentationPlots:
    """The three public plot entry points, headless (matplotlib Agg)."""

    @pytest.fixture()
    def grid(self, landcover_grid):
        return landcover_grid.assign_coords(
            landcover=(("i", "j"), np.full(landcover_grid.shape, 13, dtype=int)),
            z0=(("i", "j"), np.linspace(0.1, 0.9, landcover_grid.size).reshape(landcover_grid.shape)),
            lambdaF=(("i", "j"), np.linspace(0.05, 0.4, landcover_grid.size).reshape(landcover_grid.shape)),
        )

    def test_plotting_landcover_opens_a_figure_with_one_patch_per_cell(self, lc, grid, axes_image):
        import matplotlib.pyplot as plt

        _, image = axes_image
        before = len(plt.get_fignums())
        lc.presentation.plotLandcover(image, grid, figsize=(2, 2))
        assert len(plt.get_fignums()) == before + 1
        assert len(plt.gcf().axes[0].patches) == grid.size

    def test_plotting_roughness_adds_a_labelled_colour_bar(self, lc, grid, axes_image):
        import matplotlib.pyplot as plt

        _, image = axes_image
        lc.presentation.plotRoughness(image, grid, figsize=(2, 2))
        figure = plt.gcf()
        # main axes plus the colour bar's own axes
        assert len(figure.axes) == 2
        assert figure.axes[1].get_ylabel() == "Roughness Value (z0)"

    def test_plotting_a_lambda_field_labels_the_colour_bar_with_that_field(self, lc, grid, axes_image):
        import matplotlib.pyplot as plt

        _, image = axes_image
        lc.presentation.plotLambdas("lambdaF", image, grid, figsize=(2, 2))
        assert plt.gcf().axes[1].get_ylabel() == "lambdaF Value"

    def test_the_alpha_argument_reaches_the_patches(self, lc, grid, axes_image):
        import matplotlib.pyplot as plt

        _, image = axes_image
        lc.presentation.plotLandcover(image, grid, alpha=0.75, figsize=(2, 2))
        assert plt.gcf().axes[0].patches[0].get_alpha() == pytest.approx(0.75)


@pytest.mark.unit
class TestPlotLambdasNormalisesOnZ0:
    """B186: plotLambdas' colour bar is built from z0, not from `field`."""

    @pytest.fixture()
    def grid_without_z0(self, landcover_grid):
        return landcover_grid.assign_coords(
            lambdaF=(("i", "j"), np.linspace(0.05, 0.4, landcover_grid.size).reshape(landcover_grid.shape))
        )

    def test_the_rectangles_alone_need_only_the_named_field(self, lc, grid_without_z0):
        """_getLambdasRectangles normalises on `field`, as intended."""
        rectangles = lc.presentation._getLambdasRectangles("lambdaF", grid_without_z0)
        assert len(rectangles) == grid_without_z0.size

    def test_plotting_a_lambda_field_without_z0_raises(self, lc, grid_without_z0, axes_image):
        """Characterisation of B186."""
        _, image = axes_image
        with pytest.raises(AttributeError, match="z0"):
            lc.presentation.plotLambdas("lambdaF", image, grid_without_z0, figsize=(2, 2))

    def test_the_colour_bar_range_comes_from_z0_not_from_the_field(self, lc, landcover_grid, axes_image):
        """Characterisation of B186: the bar's limits track z0's values."""
        import matplotlib.pyplot as plt

        grid = landcover_grid.assign_coords(
            lambdaF=(("i", "j"), np.full(landcover_grid.shape, 0.3)),
            z0=(("i", "j"), np.linspace(1.0, 5.0, landcover_grid.size).reshape(landcover_grid.shape)),
        )
        _, image = axes_image
        lc.presentation.plotLambdas("lambdaF", image, grid, figsize=(2, 2))
        colorbar_axes = plt.gcf().axes[1]
        assert colorbar_axes.get_ylim() == pytest.approx((1.0, 5.0))

    @pytest.mark.xfail(
        strict=True,
        reason="B186: presentation.plotLambdas colours its rectangles by "
               "normalising landcover[field] (in _getLambdasRectangles) but "
               "builds its colour bar by normalising landcover.z0, then labels "
               "that bar f'{field} Value'. The legend therefore describes the "
               "wrong variable whenever field and z0 have different ranges, and "
               "the call fails outright with AttributeError on a landcover that "
               "carries lambda fields but no z0 -- which is exactly the "
               "non-urban case. Both norms should use landcover[field]. See the "
               "consolidated findings issue.",
    )
    def test_the_colour_bar_range_matches_the_field_being_plotted(self, lc, landcover_grid, axes_image):
        import matplotlib.pyplot as plt

        grid = landcover_grid.assign_coords(
            lambdaF=(("i", "j"), np.linspace(0.05, 0.4, landcover_grid.size).reshape(landcover_grid.shape)),
            z0=(("i", "j"), np.linspace(1.0, 5.0, landcover_grid.size).reshape(landcover_grid.shape)),
        )
        _, image = axes_image
        lc.presentation.plotLambdas("lambdaF", image, grid, figsize=(2, 2))
        assert plt.gcf().axes[1].get_ylim() == pytest.approx((0.05, 0.4))


# --------------------------------------------------------------------------
# source-level observations
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestSourceLevelObservations:
    def test_roughnesslength2sandgrainroughness_is_defined_twice(self, lc):
        """The class body defines the same staticmethod twice; the second wins.

        Both bodies compute z0 * 30, so this is a smell rather than a bug --
        recorded here so that a future edit to only one of them is caught.
        """
        source = inspect.getsource(landcover_module.LandCoverToolkit)
        assert source.count("def roughnesslength2sandgrainroughness") == 2
        assert lc.roughnesslength2sandgrainroughness(2.0) == pytest.approx(60.0)

    def test_the_module_carries_three_separate_igbp_roughness_tables(self):
        """B43 plus B184: _handleType1, getRoughnessAtPoint, getRoughnessFromLandcover."""
        source = inspect.getsource(landcover_module.LandCoverToolkit)
        assert source.count("0.0001,  # Water") == 1          # _handleType1
        assert source.count("0.01,  # Example values") == 1    # getRoughnessAtPoint
        assert source.count("0: 0.01, 1: 0.02, 2: 0.05") == 1  # the inline fallback
