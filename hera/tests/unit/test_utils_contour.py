"""Converting matplotlib contours into geopandas polygons.

The conversion is geometric, so the assertions are geometric: a known contour
of a known field must come back as a polygon enclosing a known area.
"""
import numpy as np
import pytest

from hera.utils.matplotlibCountour import standardize_polygon, toGeopandas


@pytest.mark.unit
class TestStandardizePolygon:
    """Scales coordinates by a conversion factor, accepting a list or an array."""

    def test_list_of_pairs_is_scaled(self):
        assert standardize_polygon([(1.0, 2.0), (3.0, 4.0)], 2.0) == [
            (2.0, 4.0),
            (6.0, 8.0),
        ]

    def test_numpy_array_is_scaled(self):
        scaled = standardize_polygon(np.array([[1.0, 2.0], [3.0, 4.0]]), 2.0)
        assert scaled == [(2.0, 4.0), (6.0, 8.0)]

    def test_the_two_input_forms_agree(self):
        pairs = [(1.5, 2.5), (3.5, 4.5)]
        assert standardize_polygon(pairs, 3.0) == standardize_polygon(
            np.array(pairs), 3.0
        )

    def test_unit_factor_is_the_identity(self):
        pairs = [(1.0, 2.0), (3.0, 4.0)]
        assert standardize_polygon(pairs, 1.0) == pairs

    def test_output_is_a_list_of_tuples(self):
        result = standardize_polygon(np.array([[1.0, 2.0]]), 1.0)
        assert isinstance(result, list)
        assert isinstance(result[0], tuple)

    def test_empty_input_gives_empty_output(self):
        assert standardize_polygon([], 2.0) == []

    def test_vertex_count_is_preserved(self):
        pairs = [(float(i), float(i)) for i in range(7)]
        assert len(standardize_polygon(pairs, 0.5)) == 7


@pytest.fixture()
def circular_contour():
    """A single contour of r^2 = x^2 + y^2 at level 1, i.e. the unit circle.

    Chosen because its enclosed area is pi, which makes the geometry
    assertions statements about correctness rather than snapshots.
    """
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    grid = np.linspace(-2, 2, 400)
    x, y = np.meshgrid(grid, grid)
    figure, axes = plt.subplots()
    contour = axes.contour(x, y, x**2 + y**2, levels=[1.0])
    yield contour
    plt.close(figure)


@pytest.mark.unit
class TestToGeopandas:
    def test_returns_a_geodataframe_with_level_and_contour(self, circular_contour):
        import geopandas as gpd

        result = toGeopandas(circular_contour)
        assert isinstance(result, gpd.GeoDataFrame)
        assert set(result.columns) == {"Level", "contour"}
        assert result.geometry.name == "contour"

    def test_the_level_is_carried_over(self, circular_contour):
        result = toGeopandas(circular_contour)
        assert result["Level"].tolist() == [pytest.approx(1.0)]

    def test_the_polygon_encloses_the_right_area(self, circular_contour):
        """The unit circle: area pi, to within the contour's discretisation."""
        result = toGeopandas(circular_contour)
        assert result.geometry.iloc[0].area == pytest.approx(np.pi, rel=1e-3)

    def test_the_polygon_is_valid(self, circular_contour):
        result = toGeopandas(circular_contour)
        assert result.geometry.iloc[0].is_valid

    def test_input_units_scale_the_geometry(self, circular_contour):
        """Output is metres; a kilometre input must scale linearly, so area x1e6."""
        from hera.utils.unitHandler import ureg

        in_metres = toGeopandas(circular_contour)
        in_kilometres = toGeopandas(circular_contour, inunits=1 * ureg.km)

        ratio = in_kilometres.geometry.iloc[0].area / in_metres.geometry.iloc[0].area
        assert ratio == pytest.approx(1e6, rel=1e-6)

    def test_default_units_are_metres(self, circular_contour):
        from hera.utils.unitHandler import ureg

        assert toGeopandas(circular_contour).geometry.iloc[0].area == pytest.approx(
            toGeopandas(circular_contour, inunits=1 * ureg.m).geometry.iloc[0].area
        )

    def test_multiple_levels_produce_multiple_rows(self):
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        grid = np.linspace(-3, 3, 300)
        x, y = np.meshgrid(grid, grid)
        figure, axes = plt.subplots()
        contour = axes.contour(x, y, x**2 + y**2, levels=[1.0, 4.0])
        try:
            result = toGeopandas(contour)
            assert sorted(result["Level"].tolist()) == [
                pytest.approx(1.0),
                pytest.approx(4.0),
            ]
            areas = sorted(result.geometry.area.tolist())
            assert areas[0] == pytest.approx(np.pi, rel=1e-2)
            assert areas[1] == pytest.approx(4 * np.pi, rel=1e-2)
        finally:
            plt.close(figure)
