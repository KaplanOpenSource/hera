"""``riskassessment/presentation/casualtiesFigs.py``: the casualty rose, the
casualty projection map, and the shapely-to-matplotlib ``plot_geom`` helper.

These are the presentation layer of ``RiskToolkit`` (``risk.presentation``).
They take a ``thresholdGeoDataFrame`` of injury isopleths and draw them, so
the tests here supply a stand-in for that object exposing exactly the two
methods the plots call -- ``project(area, loc, mathematical_angle=...)`` and
``shiftLocationAndAngle(loc, mathematical_angle=..., geometry=...)``.  That
keeps the tests about the plotting arithmetic rather than about
``thresholdGeoDataFrame``, which is covered in its own files (and whose
``_project`` is blocked behind the already-reported B108).

The assertions are on the artists the methods leave on the axes -- bar
widths, heights, patch counts, line counts -- and on the pivoted table
``plotCasualtiesRose`` returns, whose numbers are derived from the
stand-in's own input rather than from hera's output.  The conftest forces
the Agg backend, so nothing is displayed.

Two conventions matter throughout:

* ``meteorological_angles`` are converted with ``toMathematicalAngle``
  before use, and the returned table reports the *mathematical* angle in
  degrees.  ``toMathematicalAngle(0) == 270``, so a north wind appears at
  270 in the result.
* the bars are *stacked*: each severity's bar starts at the running total of
  the severities drawn before it, which is what stops a casualty from being
  counted twice across severities.

One bug is pinned below as B265: asking ``plotCasualtiesProjection`` for
more than one ``plumSeverity`` silently draws no plume outline at all,
because the pandas ``query`` expression is built by ``%``-formatting a numpy
array.

Robustness note (not pinned as a bug): ``plotCasualtiesRose`` builds nine
tick labels from ``numpy.linspace(0, 360, 9)`` and hands them to
``ax.set_xticklabels`` on a polar axis that has eight ticks, which
matplotlib 3.9 answers with a UserWarning.  The eight labels that are used
are the correct ones (the ninth duplicates the first, 360 deg == 0 deg), so
the figure is right and only the call is untidy.

Deliberately not covered here: the wind-distribution overlay's *numeric*
placement.  ``plotCasualtiesRose(windDistribution=...)`` builds a second,
floating polar axis through ``mpl_toolkits.axisartist``; the test below
asserts that the auxiliary axis is created and carries the series, but the
grid-helper tick geometry is matplotlib's business, not hera's.
"""
import numpy
import pandas
import pytest
from shapely.geometry import (GeometryCollection, LineString, MultiPolygon,
                              Point, Polygon)

from hera.riskassessment.presentation.casualtiesFigs import casualtiesPlot, plot_geom
from hera.utils import toMathematicalAngle

POPULATION = "effectedtotal_pop"


def _square(x0, y0, side):
    return Polygon([(x0, y0), (x0 + side, y0), (x0 + side, y0 + side), (x0, y0 + side)])


class _Isopleths:
    """The two methods the plots need from a thresholdGeoDataFrame."""

    def __init__(self, table=None, geoTable=None, blindAngles=()):
        self._table = table
        self._geoTable = geoTable
        self._blindAngles = set(blindAngles)
        self.projectCalls = []

    def project(self, area, loc=None, mathematical_angle=None):
        self.projectCalls.append((area, loc, mathematical_angle))
        if mathematical_angle in self._blindAngles:
            return None
        source = self._table if self._table is not None else self._geoTable
        return source.copy()

    def shiftLocationAndAngle(self, loc, mathematical_angle=None,
                              geometry="TotalPolygon"):
        return self._geoTable.copy()


def _casualties(**perSeverity):
    """A flat projection result: one row per severity, one population column."""
    rows = []
    for severity, values in perSeverity.items():
        for value in numpy.atleast_1d(values):
            rows.append({"severity": severity, POPULATION: float(value)})
    return pandas.DataFrame(rows)


@pytest.fixture()
def plotter():
    return casualtiesPlot()


# ---------------------------------------------------------------------------
# plot_geom
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPlotGeom:
    @pytest.fixture()
    def axes(self):
        import matplotlib.pyplot as plt

        figure, axes = plt.subplots()
        return axes

    def test_a_polygon_becomes_one_patch(self, axes):
        plot_geom(axes, _square(0, 0, 1), dict(facecolor="red"), dict(color="black"))

        assert len(axes.patches) == 1

    def test_the_patch_keeps_the_polygon_vertices(self, axes):
        plot_geom(axes, _square(2, 3, 4), dict(facecolor="red"), dict(color="black"))

        assert axes.patches[0].get_xy().min(axis=0).tolist() == [2.0, 3.0]
        assert axes.patches[0].get_xy().max(axis=0).tolist() == [6.0, 7.0]

    def test_the_patch_properties_are_forwarded(self, axes):
        plot_geom(axes, _square(0, 0, 1), dict(facecolor="red"), dict(color="black"))

        assert axes.patches[0].get_facecolor() == (1.0, 0.0, 0.0, 1.0)

    def test_an_interior_ring_is_drawn_as_a_line_on_top_of_the_patch(self, axes):
        holey = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)],
                        [[(2, 2), (4, 2), (4, 4), (2, 4)]])
        plot_geom(axes, holey, dict(facecolor="red"), dict(color="black"))

        assert len(axes.patches) == 1
        assert len(axes.lines) == 1

    def test_a_point_becomes_a_marker(self, axes):
        plot_geom(axes, Point(3, 4), {}, dict(color="blue"))

        assert len(axes.lines) == 1
        assert axes.lines[0].get_xydata().tolist() == [[3.0, 4.0]]

    def test_a_linestring_becomes_a_line(self, axes):
        plot_geom(axes, LineString([(0, 0), (1, 2)]), {}, dict(color="green"))

        assert len(axes.lines) == 1
        assert axes.lines[0].get_xydata().tolist() == [[0.0, 0.0], [1.0, 2.0]]

    def test_a_multipolygon_is_recursed_into(self, axes):
        multi = MultiPolygon([_square(0, 0, 1), _square(5, 5, 1), _square(9, 9, 1)])
        plot_geom(axes, multi, dict(facecolor="red"), dict(color="black"))

        assert len(axes.patches) == 3

    def test_a_mixed_collection_dispatches_per_member(self, axes):
        collection = GeometryCollection([_square(0, 0, 1), Point(5, 5)])
        plot_geom(axes, collection, dict(facecolor="red"), dict(color="black"))

        assert len(axes.patches) == 1
        assert len(axes.lines) == 1

    def test_an_empty_geometry_draws_nothing(self, axes):
        plot_geom(axes, Polygon(), dict(facecolor="red"), dict(color="black"))

        assert len(axes.patches) == 0
        assert len(axes.lines) == 0

    def test_an_empty_member_of_a_collection_is_skipped(self, axes):
        multi = MultiPolygon([_square(0, 0, 1)])
        plot_geom(axes, GeometryCollection([Polygon(), multi]),
                  dict(facecolor="red"), dict(color="black"))

        assert len(axes.patches) == 1


# ---------------------------------------------------------------------------
# plotCasualtiesRose
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestCasualtiesRoseAngles:
    def test_neither_angle_convention_is_rejected(self, plotter):
        with pytest.raises(ValueError, match="Must supply meteorology or mathematical angles"):
            plotter.plotCasualtiesRose(_Isopleths(_casualties(Severe=1.0)),
                                       area={}, severityList=["Severe"], loc=(0, 0))

    def test_one_projection_is_requested_per_angle(self, plotter):
        isopleths = _Isopleths(_casualties(Severe=1.0))
        plotter.plotCasualtiesRose(isopleths, area={}, severityList=["Severe"],
                                   loc=(0, 0), mathematical_angles=[0.0, 90.0, 180.0],
                                   legend=False)

        assert [call[2] for call in isopleths.projectCalls] == [0.0, 90.0, 180.0]

    def test_the_area_and_location_reach_the_projection(self, plotter):
        isopleths = _Isopleths(_casualties(Severe=1.0))
        plotter.plotCasualtiesRose(isopleths, area={"radius": 5}, severityList=["Severe"],
                                   loc=(100, 200), mathematical_angles=[0.0], legend=False)

        assert isopleths.projectCalls[0][0] == {"radius": 5}
        assert isopleths.projectCalls[0][1] == (100, 200)

    def test_meteorological_angles_are_converted_before_projecting(self, plotter):
        isopleths = _Isopleths(_casualties(Severe=1.0))
        plotter.plotCasualtiesRose(isopleths, area={}, severityList=["Severe"],
                                   loc=(0, 0), meteorological_angles=[0.0, 90.0],
                                   legend=False)

        assert [call[2] for call in isopleths.projectCalls] == [
            toMathematicalAngle(0.0), toMathematicalAngle(90.0)]

    def test_mathematical_angles_win_only_when_no_met_angle_is_given(self, plotter):
        """Documented precedence: meteorological angles take priority."""
        isopleths = _Isopleths(_casualties(Severe=1.0))
        plotter.plotCasualtiesRose(isopleths, area={}, severityList=["Severe"],
                                   loc=(0, 0), meteorological_angles=[0.0],
                                   mathematical_angles=[42.0], legend=False)

        assert [call[2] for call in isopleths.projectCalls] == [toMathematicalAngle(0.0)]

    def test_the_returned_table_reports_angles_in_degrees(self, plotter):
        _, table = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[0.0, 90.0], legend=False)

        assert sorted(table["angle"]) == [0.0, 90.0]

    def test_a_north_wind_lands_at_the_mathematical_angle(self, plotter):
        _, table = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), meteorological_angles=[0.0], legend=False)

        assert table["angle"].tolist() == [toMathematicalAngle(0.0)]

    def test_an_angle_without_casualties_is_dropped(self, plotter, capsys):
        _, table = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=5.0), blindAngles=[90.0]),
            area={}, severityList=["Severe"], loc=(0, 0),
            mathematical_angles=[0.0, 90.0], legend=False)

        assert table["angle"].tolist() == [0.0]
        assert "90.0 does not have casualties" in capsys.readouterr().out


@pytest.mark.unit
class TestCasualtiesRoseTotals:
    def test_rows_of_the_same_severity_are_summed(self, plotter):
        _, table = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=[1.0, 2.0, 4.0])), area={},
            severityList=["Severe"], loc=(0, 0), mathematical_angles=[0.0],
            legend=False)

        assert table["Severe"].tolist() == [7.0]

    def test_each_severity_gets_its_own_column(self, plotter):
        _, table = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=[1.0, 2.0], Light=4.0)), area={},
            severityList=["Severe", "Light"], loc=(0, 0),
            mathematical_angles=[0.0], legend=False)

        assert table["Severe"].tolist() == [3.0]
        assert table["Light"].tolist() == [4.0]

    def test_the_total_column_sums_the_requested_severities(self, plotter):
        _, table = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=3.0, Light=4.0)), area={},
            severityList=["Severe", "Light"], loc=(0, 0),
            mathematical_angles=[0.0], legend=False)

        assert table["total"].tolist() == [7.0]

    def test_a_severity_left_out_of_the_list_is_not_totalled(self, plotter):
        _, table = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=3.0, Light=4.0)), area={},
            severityList=["Severe"], loc=(0, 0), mathematical_angles=[0.0],
            legend=False)

        assert table["total"].tolist() == [3.0]

    def test_a_severity_absent_from_the_data_is_skipped_without_error(self, plotter):
        axes, table = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=3.0)), area={},
            severityList=["Severe", "NeverHappens"], loc=(0, 0),
            mathematical_angles=[0.0], legend=False)

        assert "NeverHappens" not in table.columns
        assert len(axes["innerax"].patches) == 1

    def test_a_missing_severity_at_one_angle_becomes_zero(self, plotter):
        class Uneven(_Isopleths):
            def project(self, area, loc=None, mathematical_angle=None):
                self.projectCalls.append((area, loc, mathematical_angle))
                if mathematical_angle == 0.0:
                    return _casualties(Severe=3.0, Light=4.0)
                return _casualties(Severe=1.0)

        _, table = plotter.plotCasualtiesRose(
            Uneven(), area={}, severityList=["Severe", "Light"], loc=(0, 0),
            mathematical_angles=[0.0, 90.0], legend=False)

        byAngle = table.set_index("angle")
        assert byAngle.loc[90.0, "Light"] == 0.0
        assert byAngle.loc[90.0, "total"] == 1.0

    def test_a_different_population_column_is_honoured(self, plotter):
        table = pandas.DataFrame({"severity": ["Severe"], "effectedchildren": [9.0]})
        _, pivoted = plotter.plotCasualtiesRose(
            _Isopleths(table), area={}, severityList=["Severe"], loc=(0, 0),
            mathematical_angles=[0.0], effectedPopulation="effectedchildren",
            legend=False)

        assert pivoted["Severe"].tolist() == [9.0]


@pytest.mark.unit
class TestCasualtiesRoseArtists:
    def test_the_inner_axes_is_polar(self, plotter):
        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[0.0], legend=False)

        assert axes["innerax"].name == "polar"

    def test_without_a_wind_distribution_there_is_no_outer_axes(self, plotter):
        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[0.0], legend=False)

        assert axes["outerax"] is None

    def test_one_bar_per_severity_per_angle(self, plotter):
        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0, Light=2.0)), area={},
            severityList=["Severe", "Light"], loc=(0, 0),
            mathematical_angles=[0.0, 90.0, 180.0], legend=False)

        assert len(axes["innerax"].patches) == 6

    def test_each_bar_is_as_tall_as_its_casualty_count(self, plotter):
        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=3.0, Light=4.0)), area={},
            severityList=["Severe", "Light"], loc=(0, 0),
            mathematical_angles=[0.0], legend=False)

        assert sorted(p.get_height() for p in axes["innerax"].patches) == [3.0, 4.0]

    def test_the_bars_are_stacked_so_the_second_starts_at_the_first(self, plotter):
        """The stacking is the visual form of "do not double count"."""
        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=3.0, Light=4.0)), area={},
            severityList=["Severe", "Light"], loc=(0, 0),
            mathematical_angles=[0.0], legend=False)

        bottoms = sorted(p.get_y() for p in axes["innerax"].patches)
        assert bottoms == [0.0, 3.0]

    def test_the_bars_are_centred_on_the_angle_in_radians(self, plotter):
        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[90.0], legend=False)

        bar = axes["innerax"].patches[0]
        assert bar.get_x() + bar.get_width() / 2 == pytest.approx(numpy.pi / 2)

    def test_each_bar_series_is_labelled_with_its_severity(self, plotter):
        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0, Light=2.0)), area={},
            severityList=["Severe", "Light"], loc=(0, 0),
            mathematical_angles=[0.0], legend=False)

        assert {c.get_label() for c in axes["innerax"].containers} == {"Severe", "Light"}

    def test_the_default_bar_width_is_the_documented_default(self, plotter):
        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[0.0], legend=False)

        assert axes["innerax"].patches[0].get_width() == pytest.approx(0.18)

    def test_a_weight_overrides_the_bar_width(self, plotter):
        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[0.0], weights=0.5, legend=False)

        assert axes["innerax"].patches[0].get_width() == pytest.approx(0.5)

    def test_a_supplied_cycler_is_used_and_still_gets_a_width(self, plotter):
        """A cycler without a 'width' key has the default width added to it,
        so the caller's properties and the bar geometry both apply."""
        import matplotlib.pyplot as plt

        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[0.0],
            cycler=plt.cycler(color=["red"]), legend=False)

        patch = axes["innerax"].patches[0]
        assert patch.get_facecolor() == (1.0, 0.0, 0.0, 1.0)
        assert patch.get_width() == pytest.approx(0.18)

    def test_a_supplied_cycler_with_weights_takes_the_weight(self, plotter):
        import matplotlib.pyplot as plt

        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[0.0], weights=0.9,
            cycler=plt.cycler(color=["red"]), legend=False)

        assert axes["innerax"].patches[0].get_width() == pytest.approx(0.9)

    def test_a_cycler_that_already_sets_the_width_is_left_alone(self, plotter):
        import matplotlib.pyplot as plt

        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[0.0],
            cycler=plt.cycler(width=[0.4]), legend=False)

        assert axes["innerax"].patches[0].get_width() == pytest.approx(0.4)

    def test_a_subplot_spec_places_the_axes(self, plotter):
        """Passing a list is documented as the add_subplot signature."""
        import matplotlib.pyplot as plt

        figure = plt.figure()
        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[0.0], ax=[2, 1, 2], legend=False)

        assert axes["innerax"] in figure.axes

    def test_the_meteorological_tick_labels_are_written(self, plotter):
        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[0.0], legend=False)

        labels = [t.get_text() for t in axes["innerax"].get_xticklabels()]
        assert labels[0] == "$90^o$"
        assert labels[2] == "$0^o$"

    def test_a_legend_is_added_when_asked_for(self, plotter):
        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[0.0], legend=True)

        assert axes["innerax"].get_legend() is not None

    def test_no_legend_is_added_when_not_asked_for(self, plotter):
        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[0.0], legend=False)

        assert axes["innerax"].get_legend() is None


@pytest.mark.unit
class TestCasualtiesRoseWindOverlay:
    @staticmethod
    def _distribution():
        return pandas.DataFrame({"angle": [0.0, 90.0, 180.0, 270.0],
                                 "distribution": [40.0, 30.0, 20.0, 10.0]})

    def test_an_outer_axes_is_created_for_the_overlay(self, plotter):
        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[0.0],
            windDistribution=self._distribution(), legend=False)

        assert axes["outerax"] is not None

    def test_the_distribution_is_drawn_as_a_line_by_default(self, plotter):
        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[0.0],
            windDistribution=self._distribution(), legend=False)

        assert len(axes["outerax"].lines) == 1
        assert len(axes["outerax"].lines[0].get_xydata()) == 4

    def test_the_plot_type_can_be_switched_to_a_scatter(self, plotter):
        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[0.0],
            windDistribution=self._distribution(), plotType="scatter", legend=False)

        assert len(axes["outerax"].lines) == 0
        assert len(axes["outerax"].collections) == 1

    def test_the_supplied_distribution_frame_is_not_mutated(self, plotter):
        """The overlay rescales angle and distribution; it copies first."""
        distribution = self._distribution()
        before = distribution.copy(deep=True)

        plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[0.0],
            windDistribution=distribution, legend=False)

        pandas.testing.assert_frame_equal(distribution, before)

    def test_the_radial_limit_leaves_room_for_the_overlay(self, plotter):
        """The casualty bars are squeezed into the lower half so the wind
        rose can occupy the outer ring: the limit is twice the largest total.
        """
        axes, table = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=3.0, Light=4.0)), area={},
            severityList=["Severe", "Light"], loc=(0, 0),
            mathematical_angles=[0.0], windDistribution=self._distribution(),
            legend=False)

        assert axes["innerax"].get_ylim()[1] == pytest.approx(2 * table["total"].max())

    def test_custom_wind_ticks_are_accepted(self, plotter):
        axes, _ = plotter.plotCasualtiesRose(
            _Isopleths(_casualties(Severe=1.0)), area={}, severityList=["Severe"],
            loc=(0, 0), mathematical_angles=[0.0],
            windDistribution=self._distribution(), windTicks=[10, 20, 40],
            legend=False)

        assert axes["outerax"] is not None


# ---------------------------------------------------------------------------
# plotCasualtiesProjection
# ---------------------------------------------------------------------------

def _projected():
    import geopandas

    return geopandas.GeoDataFrame(
        {"severity": ["Severe", "Light"],
         "geometry": [_square(0, 0, 1), _square(0, 0, 3)]},
        geometry="geometry")


def _totals():
    import geopandas

    return geopandas.GeoDataFrame(
        {"severity": ["Severe", "Light"],
         "TotalPolygon": [_square(0, 0, 1), _square(0, 0, 3)]},
        geometry="TotalPolygon")


@pytest.mark.unit
class TestCasualtiesProjection:
    @pytest.fixture()
    def isopleths(self):
        return _Isopleths(geoTable=_projected())

    @staticmethod
    def _both():
        source = _Isopleths(geoTable=_projected())
        source._geoTable = _projected()
        return source

    def test_neither_angle_convention_is_rejected(self, plotter, isopleths):
        with pytest.raises(ValueError, match="Must supply meteorology or mathematical angle"):
            plotter.plotCasualtiesProjection(isopleths, area={},
                                             severityList=["Severe"], loc=(0, 0))

    def test_no_casualties_at_all_is_rejected(self, plotter):
        source = _Isopleths(geoTable=_projected(), blindAngles=[0.0])
        with pytest.raises(ValueError, match="no casualties in this parameters set"):
            plotter.plotCasualtiesProjection(source, area={}, severityList=["Severe"],
                                             loc=(0, 0), mathematical_angle=0.0)

    def test_a_meteorological_angle_is_converted_before_projecting(self, plotter,
                                                                   isopleths):
        plotter.plotCasualtiesProjection(isopleths, area={}, severityList=["Severe"],
                                         loc=(0, 0), meteorological_angle=90.0)

        assert isopleths.projectCalls[0][2] == toMathematicalAngle(90.0)

    def test_a_mathematical_angle_is_passed_through(self, plotter, isopleths):
        plotter.plotCasualtiesProjection(isopleths, area={}, severityList=["Severe"],
                                         loc=(0, 0), mathematical_angle=42.0)

        assert isopleths.projectCalls[0][2] == 42.0

    def test_the_projection_is_returned_alongside_the_axes(self, plotter, isopleths):
        axes, projection = plotter.plotCasualtiesProjection(
            isopleths, area={}, severityList=["Severe"], loc=(0, 0),
            mathematical_angle=0.0)

        assert sorted(projection["severity"]) == ["Light", "Severe"]

    def test_the_axes_it_was_given_is_the_axes_it_returns(self, plotter, isopleths):
        import matplotlib.pyplot as plt

        figure, given = plt.subplots()
        axes, _ = plotter.plotCasualtiesProjection(
            isopleths, area={}, severityList=["Severe"], loc=(0, 0),
            mathematical_angle=0.0, ax=given)

        assert axes is given

    def test_one_polygon_is_drawn_per_requested_severity(self, plotter, isopleths):
        axes, _ = plotter.plotCasualtiesProjection(
            isopleths, area={}, severityList=["Severe", "Light"], loc=(0, 0),
            mathematical_angle=0.0)

        assert len(axes.patches) == 2

    def test_a_severity_left_out_of_the_list_is_not_drawn(self, plotter, isopleths):
        axes, _ = plotter.plotCasualtiesProjection(
            isopleths, area={}, severityList=["Severe"], loc=(0, 0),
            mathematical_angle=0.0)

        assert len(axes.patches) == 1

    def test_a_severity_absent_from_the_projection_is_skipped(self, plotter, isopleths):
        axes, _ = plotter.plotCasualtiesProjection(
            isopleths, area={}, severityList=["Severe", "NeverHappens"], loc=(0, 0),
            mathematical_angle=0.0)

        assert len(axes.patches) == 1

    def test_without_an_axes_one_is_created(self, plotter, isopleths):
        axes, _ = plotter.plotCasualtiesProjection(
            isopleths, area={}, severityList=["Severe"], loc=(0, 0),
            mathematical_angle=0.0)

        assert axes.name == "rectilinear"

    def test_a_subplot_spec_places_the_axes(self, plotter, isopleths):
        import matplotlib.pyplot as plt

        figure = plt.figure()
        axes, _ = plotter.plotCasualtiesProjection(
            isopleths, area={}, severityList=["Severe"], loc=(0, 0),
            mathematical_angle=0.0, ax=[2, 1, 1])

        assert axes in figure.axes

    def test_no_plume_outline_is_drawn_by_default(self, plotter, isopleths):
        axes, _ = plotter.plotCasualtiesProjection(
            isopleths, area={}, severityList=["Severe"], loc=(0, 0),
            mathematical_angle=0.0)

        assert len(axes.lines) == 0

    def test_a_plume_severity_adds_its_outline(self, plotter):
        source = _Isopleths(geoTable=_projected())
        source._geoTable = _projected()
        source.shiftLocationAndAngle = lambda loc, mathematical_angle=None, \
            geometry="TotalPolygon": _totals()

        axes, _ = plotter.plotCasualtiesProjection(
            source, area={}, severityList=["Severe"], loc=(0, 0),
            mathematical_angle=0.0, plumSeverity=["Severe"])

        assert len(axes.lines) == 1

    def test_the_outline_traces_the_convex_hull_of_the_total_polygon(self, plotter):
        source = _Isopleths(geoTable=_projected())
        source.shiftLocationAndAngle = lambda loc, mathematical_angle=None, \
            geometry="TotalPolygon": _totals()

        axes, _ = plotter.plotCasualtiesProjection(
            source, area={}, severityList=["Severe"], loc=(0, 0),
            mathematical_angle=0.0, plumSeverity=["Light"])

        drawn = axes.lines[0].get_xydata()
        assert drawn.min(axis=0).tolist() == [0.0, 0.0]
        assert drawn.max(axis=0).tolist() == [3.0, 3.0]

    def test_a_supplied_cycler_colours_the_polygons(self, plotter, isopleths):
        import matplotlib.pyplot as plt

        axes, _ = plotter.plotCasualtiesProjection(
            isopleths, area={}, severityList=["Severe"], loc=(0, 0),
            mathematical_angle=0.0,
            cycler=plt.cycler(facecolor=["red"]) * plt.cycler(edgecolor=["blue"]))

        assert axes.patches[0].get_facecolor() == (1.0, 0.0, 0.0, 1.0)

    def test_a_supplied_boundary_cycler_colours_the_outline(self, plotter):
        import matplotlib.pyplot as plt

        source = _Isopleths(geoTable=_projected())
        source.shiftLocationAndAngle = lambda loc, mathematical_angle=None, \
            geometry="TotalPolygon": _totals()

        axes, _ = plotter.plotCasualtiesProjection(
            source, area={}, severityList=["Severe"], loc=(0, 0),
            mathematical_angle=0.0, plumSeverity=["Severe"],
            boundarycycler=plt.cycler(color=["magenta"]))

        assert axes.lines[0].get_color() == "magenta"


@pytest.mark.unit
class TestSeveralPlumeSeveritiesDrawNothing:
    """B265: see the reasons below."""

    @staticmethod
    def _source():
        source = _Isopleths(geoTable=_projected())
        source.shiftLocationAndAngle = lambda loc, mathematical_angle=None, \
            geometry="TotalPolygon": _totals()
        return source

    @pytest.mark.xfail(
        strict=True,
        reason="B265: plotCasualtiesProjection selects the plume outlines "
               "with .query(\"severity in %s\" % "
               "numpy.atleast_1d(plumSeverity)).  str() of a numpy array "
               "separates elements with spaces and no commas, so two "
               "severities render as \"severity in ['Severe' 'Light']\"; "
               "Python's adjacent-string-literal concatenation collapses that "
               "to the single name ['SevereLight'], the query matches no row, "
               "and every requested outline is silently dropped -- no error, "
               "no warning, just a missing plume boundary.  One severity "
               "happens to work because a single-element array has nothing to "
               "concatenate with.  The list should be interpolated as a "
               "Python list (e.g. list(numpy.atleast_1d(plumSeverity))) or "
               "bound through @-variables.  See the consolidated findings "
               "issue.",
    )
    def test_two_plume_severities_should_draw_two_outlines(self, plotter):
        axes, _ = plotter.plotCasualtiesProjection(
            self._source(), area={}, severityList=["Severe"], loc=(0, 0),
            mathematical_angle=0.0, plumSeverity=["Severe", "Light"])

        assert len(axes.lines) == 2

    def test_two_plume_severities_currently_draw_none(self, plotter):
        """Characterisation of B265."""
        axes, _ = plotter.plotCasualtiesProjection(
            self._source(), area={}, severityList=["Severe"], loc=(0, 0),
            mathematical_angle=0.0, plumSeverity=["Severe", "Light"])

        assert len(axes.lines) == 0

    def test_one_plume_severity_still_works(self, plotter):
        """Characterisation of B265: isolates the element count as the
        trigger, so the diagnosis is unambiguous."""
        axes, _ = plotter.plotCasualtiesProjection(
            self._source(), area={}, severityList=["Severe"], loc=(0, 0),
            mathematical_angle=0.0, plumSeverity=["Severe"])

        assert len(axes.lines) == 1

    def test_the_query_string_it_builds_selects_a_concatenated_name(self):
        """Characterisation of B265's mechanism, derived from the standard
        library rather than from hera's output."""
        expression = "severity in %s" % numpy.atleast_1d(["Severe", "Light"])

        assert expression == "severity in ['Severe' 'Light']"
        assert _totals().query(expression).empty
