"""GIS coordinate handling.

CLAUDE.md requires the named constants rather than raw EPSG integers, and
explicit conversion rather than assumed CRS.  The assertions below use
independently known coordinates: Tel Aviv sits near 34.78E 32.08N in WGS84,
and the Israeli Transverse Mercator grid puts it near easting 178000,
northing 663000.
"""
import pathlib

import pytest

from hera.measurements.GIS.utils import ITM, WSG84, convertCRS

TEL_AVIV_WGS84 = (34.78, 32.08)


@pytest.mark.unit
class TestCoordinateSystemConstants:
    def test_wgs84_is_epsg_4326(self):
        assert WSG84 == 4326

    def test_itm_is_epsg_2039(self):
        """The Israeli Transverse Mercator grid."""
        assert ITM == 2039

    def test_the_two_are_distinct(self):
        assert WSG84 != ITM


@pytest.mark.unit
class TestConvertCRS:
    def test_tel_aviv_lands_in_the_israeli_grid(self):
        """ITM eastings run roughly 100-300 km, northings 400-800 km."""
        point = convertCRS([TEL_AVIV_WGS84], WSG84, ITM)[0]
        assert 100_000 < point.x < 300_000
        assert 400_000 < point.y < 800_000

    def test_the_conversion_round_trips(self):
        forward = convertCRS([TEL_AVIV_WGS84], WSG84, ITM)[0]
        recovered = convertCRS([(forward.x, forward.y)], ITM, WSG84)[0]
        assert recovered.x == pytest.approx(TEL_AVIV_WGS84[0], abs=1e-6)
        assert recovered.y == pytest.approx(TEL_AVIV_WGS84[1], abs=1e-6)

    def test_the_result_is_a_list_of_points(self):
        """Documented: "list of shapely.geometry.Point in the correct CRS".

        Its sibling TopographyToolkit.convertPointsCRS returns a GeoDataFrame
        instead, so the two conversion entry points differ in return type.
        """
        from shapely.geometry import Point

        converted = convertCRS([TEL_AVIV_WGS84], WSG84, ITM)
        assert isinstance(converted, list)
        assert isinstance(converted[0], Point)

    def test_every_point_survives(self):
        points = [(34.0, 31.0), (35.0, 32.0), (36.0, 33.0)]
        assert len(convertCRS(points, WSG84, ITM)) == len(points)

    def test_a_numpy_array_is_accepted(self):
        import numpy as np

        converted = convertCRS(np.array([[34.78, 32.08], [35.0, 32.5]]), WSG84, ITM)
        assert len(converted) == 2

    def test_a_dataframe_is_accepted(self):
        import pandas as pd

        frame = pd.DataFrame({"x": [34.78], "y": [32.08]})
        assert len(convertCRS(frame, WSG84, ITM)) == 1

    def test_an_unsupported_container_names_its_type(self):
        with pytest.raises(ValueError, match="must be numpy.array, pandas.DataFrame, or list"):
            convertCRS("not points", WSG84, ITM)

    def test_converting_to_the_same_system_is_the_identity(self):
        point = convertCRS([TEL_AVIV_WGS84], WSG84, WSG84)[0]
        assert point.x == pytest.approx(TEL_AVIV_WGS84[0])
        assert point.y == pytest.approx(TEL_AVIV_WGS84[1])

    def test_northward_stays_northward(self):
        """A point further north must have a larger ITM northing."""
        south = convertCRS([(34.78, 31.0)], WSG84, ITM)[0]
        north = convertCRS([(34.78, 33.0)], WSG84, ITM)[0]
        assert north.y > south.y

    def test_eastward_stays_eastward(self):
        west = convertCRS([(34.2, 32.08)], WSG84, ITM)[0]
        east = convertCRS([(35.4, 32.08)], WSG84, ITM)[0]
        assert east.x > west.x

    def test_a_metre_grid_gives_metre_scale_distances(self):
        """One degree of latitude is about 111 km, whatever the projection."""
        south = convertCRS([(34.78, 32.00)], WSG84, ITM)[0]
        north = convertCRS([(34.78, 33.00)], WSG84, ITM)[0]
        assert 100_000 < (north.y - south.y) < 120_000


@pytest.mark.unit
class TestImportSideEffects:
    """B42: one module runs a demo at import time and writes 8 MB to the CWD."""

    SOURCE = pathlib.Path("hera/measurements/GIS/raster/hill2stl.py")

    def test_the_driver_code_sits_at_module_level(self):
        source = self.SOURCE.read_text(encoding="utf-8")
        assert "# Run the function" in source
        assert "generate_solid_stl(function" in source

    def test_the_output_filename_is_hard_coded(self):
        assert "filename='test1.stl'" in self.SOURCE.read_text(encoding="utf-8")

    def test_the_package_import_does_not_reach_it(self):
        """Mitigating detail: only a direct import of hill2stl triggers it."""
        import hera.measurements.GIS as gis

        assert gis is not None

    @pytest.mark.xfail(
        strict=True,
        reason="B42: importing hera.measurements.GIS.raster.hill2stl executes a "
               "demo at module level -- it prints to stdout and writes an 8.1 MB "
               "test1.stl into the current working directory. Anything that walks "
               "the package tree (docs builds, IDE indexing, collection) leaves "
               "the file behind. See the consolidated findings issue.",
    )
    def test_importing_it_writes_nothing(self, tmp_path, monkeypatch):
        import importlib
        import sys

        monkeypatch.chdir(tmp_path)
        sys.modules.pop("hera.measurements.GIS.raster.hill2stl", None)
        importlib.import_module("hera.measurements.GIS.raster.hill2stl")

        assert list(tmp_path.iterdir()) == []
