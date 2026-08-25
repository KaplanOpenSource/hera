"""stlFactory: converting height-contour geopandas/pandas data to an STL raster.

One defect surfaced while probing ``vectorToSTL``:

* B74: it dispatches with
  ``isinstance(gpandas, pandas.DataFrame) or isinstance(gpandas, dask.dataframe.DataFrame)``,
  but ``dask`` is never imported anywhere in this module. Passing a
  ``geopandas.GeoDataFrame`` or a plain ``pandas.DataFrame`` never reaches
  the second operand (short-circuit ``or``), so those work -- but anything
  that is neither raises ``NameError: name 'dask' is not defined`` instead
  of the intended ``TypeError``-shaped rejection.
"""
import geopandas
import numpy
import pytest
from shapely.geometry import LineString

from hera.measurements.GIS.utils import stlFactory


@pytest.mark.unit
class TestHeightColumnsNames:
    def test_defaults_to_height(self):
        assert stlFactory().heightColumnsNames == "HEIGHT"

    def test_can_be_reassigned(self):
        factory = stlFactory()
        factory.heightColumnsNames = "Z"
        assert factory.heightColumnsNames == "Z"


@pytest.mark.unit
class TestRasterizeGeopandas:
    def test_it_returns_x_y_and_height_grids(self):
        gdf = geopandas.GeoDataFrame({
            "HEIGHT": [0.0, 10.0],
            "geometry": [LineString([(0, 0), (100, 0)]), LineString([(0, 100), (100, 100)])],
        })
        result = stlFactory().rasterizeGeopandas(gdf, dxdy=20)
        assert set(result) == {"x", "y", "height"}
        assert result["x"].shape == result["y"].shape == result["height"].shape

    def test_the_interpolated_height_increases_from_the_low_to_the_high_contour(self):
        gdf = geopandas.GeoDataFrame({
            "HEIGHT": [0.0, 10.0],
            "geometry": [LineString([(0, 0), (100, 0)]), LineString([(0, 100), (100, 100)])],
        })
        result = stlFactory().rasterizeGeopandas(gdf, dxdy=20)
        # each column of the grid runs along y: low contour first, high contour last
        assert result["height"][0][0] < result["height"][0][-1]

    def test_a_custom_height_column_name_is_respected(self):
        gdf = geopandas.GeoDataFrame({
            "Z": [5.0, 5.0],
            "geometry": [LineString([(0, 0), (50, 0)]), LineString([(0, 50), (50, 50)])],
        })
        factory = stlFactory()
        factory.heightColumnsNames = "Z"
        result = factory.rasterizeGeopandas(gdf, dxdy=20)
        assert numpy.allclose(result["height"], 5.0)


@pytest.mark.unit
class TestVectorToSTLDispatch:
    @pytest.mark.xfail(
        strict=True,
        reason="B74: the dask branch of vectorToSTL's isinstance check "
               "references the name `dask`, which this module never "
               "imports. Anything that is neither a GeoDataFrame nor a "
               "plain pandas.DataFrame raises NameError instead of a "
               "clean type rejection. See the consolidated findings issue.",
    )
    def test_an_unsupported_type_should_be_rejected_cleanly(self):
        stlFactory().vectorToSTL(object())

    def test_an_unsupported_type_currently_raises_nameerror(self):
        """Characterisation of B74."""
        with pytest.raises(NameError, match="dask"):
            stlFactory().vectorToSTL(object())

    def test_a_geodataframe_produces_an_ascii_stl_string(self):
        gdf = geopandas.GeoDataFrame({
            "HEIGHT": [0.0, 10.0],
            "geometry": [LineString([(0, 0), (100, 0)]), LineString([(0, 100), (100, 100)])],
        })
        stl_text = stlFactory().vectorToSTL(gdf, dxdy=20, solidName="Test")
        assert "solid Test" in stl_text
        assert "endsolid Test" in stl_text
