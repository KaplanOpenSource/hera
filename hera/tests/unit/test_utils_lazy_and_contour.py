"""utils/lazy.py: _LazyModule.__call__ (not covered by attribute-access
tests elsewhere). utils/matplotlibCountour.py: plot_polygons, exercised
against a real (Agg-backed) matplotlib Axes."""
import pytest

from hera.utils.lazy import _LazyModule


@pytest.mark.unit
class TestLazyModuleCall:
    def test_calling_a_non_callable_module_raises_typeerror(self):
        os_proxy = _LazyModule("os")
        with pytest.raises(TypeError):
            os_proxy()

    def test_a_second_call_reuses_the_already_loaded_module(self):
        os_proxy = _LazyModule("os")
        with pytest.raises(TypeError):
            os_proxy()
        assert object.__getattribute__(os_proxy, "_mod") is not None
        with pytest.raises(TypeError):
            os_proxy()


@pytest.mark.unit
class TestPlotPolygons:
    def test_an_empty_geometry_is_a_no_op(self):
        from shapely.geometry import Polygon
        from hera.utils.matplotlibCountour import plot_polygons

        ax = _FakeAx()
        plot_polygons(ax, Polygon())
        assert ax.patches == [] and ax.lines == []

    def test_a_polygon_adds_a_patch(self):
        from shapely.geometry import Polygon
        from hera.utils.matplotlibCountour import plot_polygons

        ax = _FakeAx()
        plot_polygons(ax, Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]))
        assert len(ax.patches) == 1

    def test_a_point_plots_a_marker(self):
        from shapely.geometry import Point
        from hera.utils.matplotlibCountour import plot_polygons

        ax = _FakeAx()
        plot_polygons(ax, Point(1, 2))
        assert ax.lines == [((1.0, 2.0), "o")]

    def test_an_unsupported_geometry_type_raises(self):
        from hera.utils.matplotlibCountour import plot_polygons

        class _NotAGeom:
            is_empty = False
            geom_type = "Weird"

        with pytest.raises(TypeError, match="Unsupported geometry type"):
            plot_polygons(_FakeAx(), _NotAGeom())


class _FakeAx:
    def __init__(self):
        self.patches = []
        self.lines = []

    def add_patch(self, patch):
        self.patches.append(patch)

    def plot(self, x, y, *args, **kwargs):
        self.lines.append(((x, y), args[0] if args else None))
