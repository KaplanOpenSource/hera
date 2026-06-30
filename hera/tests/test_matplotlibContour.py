"""
Tests for hera.utils.matplotlibCountour.

Covers matplotlib >= 3.10 compatibility (ContourSet.collections removed),
correct output structure, and valid geometry generation.
"""

import warnings

import geopandas
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
from shapely.geometry import Polygon

from hera.utils.matplotlibCountour import toGeopandas, standardize_polygon


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_circular_contourf(levels=None):
    """Return a contourf ContourSet for Z = X² + Y² on a 100×100 grid."""
    if levels is None:
        levels = [1, 2, 4, 8]
    x = np.linspace(-3, 3, 100)
    y = np.linspace(-3, 3, 100)
    X, Y = np.meshgrid(x, y)
    Z = X**2 + Y**2
    fig, ax = plt.subplots()
    cs = ax.contourf(X, Y, Z, levels=levels)
    return cs, fig


# ---------------------------------------------------------------------------
# standardize_polygon
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestStandardizePolygon:
    def test_array_input_identity(self):
        pts = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
        result = standardize_polygon(pts, 1.0)
        assert len(result) == 4
        assert result[0] == (0.0, 0.0)
        assert result[2] == (1.0, 1.0)

    def test_array_input_scaling(self):
        pts = np.array([[1.0, 2.0], [3.0, 4.0]])
        result = standardize_polygon(pts, 2.0)
        assert result[0] == (2.0, 4.0)
        assert result[1] == (6.0, 8.0)

    def test_list_of_tuples_input(self):
        pts = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
        result = standardize_polygon(pts, 1.0)
        assert len(result) == 3
        assert result[0] == (0.0, 0.0)


# ---------------------------------------------------------------------------
# toGeopandas — output structure
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestToGeopandasStructure:
    def setup_method(self):
        plt.switch_backend("Agg")

    def teardown_method(self):
        plt.close("all")

    def test_returns_geodataframe(self):
        cs, _ = _make_circular_contourf()
        result = toGeopandas(cs)
        assert isinstance(result, geopandas.GeoDataFrame)

    def test_has_level_and_contour_columns(self):
        cs, _ = _make_circular_contourf()
        result = toGeopandas(cs)
        assert "Level" in result.columns
        assert result.geometry.name == "contour"

    def test_nonempty_result(self):
        cs, _ = _make_circular_contourf()
        result = toGeopandas(cs)
        assert len(result) > 0

    def test_geometry_type_is_polygon(self):
        cs, _ = _make_circular_contourf()
        result = toGeopandas(cs)
        assert set(result.geometry.geom_type.unique()) == {"Polygon"}

    def test_all_geometries_valid(self):
        cs, _ = _make_circular_contourf()
        result = toGeopandas(cs)
        assert result.geometry.is_valid.all()

    def test_level_count_matches_intervals(self):
        levels = [1, 2, 4, 8]   # 4 levels → 3 filled bands
        cs, _ = _make_circular_contourf(levels=levels)
        result = toGeopandas(cs)
        unique_levels = sorted(result["Level"].unique())
        # Each band's lower-bound level should appear
        assert len(unique_levels) == len(levels) - 1

    def test_level_values_are_lower_bounds(self):
        levels = [1.0, 2.0, 4.0, 8.0]
        cs, _ = _make_circular_contourf(levels=levels)
        result = toGeopandas(cs)
        for lv in result["Level"]:
            assert lv in levels


# ---------------------------------------------------------------------------
# toGeopandas — geometry correctness
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestToGeopandasGeometry:
    def setup_method(self):
        plt.switch_backend("Agg")

    def teardown_method(self):
        plt.close("all")

    def test_polygons_have_nonzero_area(self):
        cs, _ = _make_circular_contourf()
        result = toGeopandas(cs)
        assert (result.geometry.area > 0).all()

    def test_total_area_is_positive(self):
        """Sum of polygon areas must be positive."""
        cs, _ = _make_circular_contourf(levels=[1, 2, 4, 8])
        result = toGeopandas(cs)
        assert result.geometry.area.sum() > 0

    def test_unit_conversion_scales_area(self):
        from hera.utils.unitHandler import ureg
        cs, _ = _make_circular_contourf()
        result_m = toGeopandas(cs)          # default: inunits=None → 1 m
        result_km = toGeopandas(cs, inunits=1 * ureg.km)
        # km → m scaling: areas should differ by 1e6
        ratio = result_km.geometry.area.sum() / result_m.geometry.area.sum()
        assert abs(ratio - 1e6) / 1e6 < 1e-3, f"Expected ratio ~1e6, got {ratio}"


# ---------------------------------------------------------------------------
# matplotlib version compatibility
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestMatplotlibCompatibility:
    def setup_method(self):
        plt.switch_backend("Agg")

    def teardown_method(self):
        plt.close("all")

    def test_no_deprecation_warning(self):
        cs, _ = _make_circular_contourf()
        with warnings.catch_warnings():
            warnings.simplefilter("error", DeprecationWarning)
            warnings.simplefilter("error", PendingDeprecationWarning)
            # Should not raise any deprecation warning
            toGeopandas(cs)

    def test_collections_attribute_absent_in_310_plus(self):
        """Verify the test environment uses the new API path."""
        cs, _ = _make_circular_contourf()
        mpl_major, mpl_minor = (int(x) for x in matplotlib.__version__.split(".")[:2])
        if mpl_major > 3 or (mpl_major == 3 and mpl_minor >= 10):
            assert not hasattr(cs, "collections"), (
                "Expected ContourSet.collections to be absent in matplotlib >= 3.10"
            )

    def test_works_with_two_levels(self):
        """Edge case: minimum meaningful contourf (2 levels = 1 band)."""
        cs, _ = _make_circular_contourf(levels=[1, 4])
        result = toGeopandas(cs)
        assert len(result) >= 1
        assert isinstance(result.geometry.iloc[0], Polygon)
