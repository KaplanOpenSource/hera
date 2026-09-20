"""measurements/meteorology/lowfreqdata/presentationLayer.py: the four
drawing methods batch-31's file left uncovered --
``DailyPlots.plotScatter``, ``DailyPlots.dateLinePlot``,
``DailyPlots.plotProbContourf`` and
``SeasonalPlots.plotProbContourf_bySeason``.  The dict/colormap builders and
the wiring classes are already covered by
test_meteorology_lowfreq_presentation.py.

Everything is asserted against the resulting matplotlib artists under the
conftest's Agg backend: scatter offsets, line data, contour levels, axis
limits, labels and subplot titles.  The expected numbers are derived from
the documented conventions rather than from the implementation:

* the daily axis is the full 24-hour day, so the hour coordinate of a
  sample is ``hour + minute/60`` and the axis runs 0..24;
* ``_contourvalsdict`` declares ``under_value=0.1``, ``max_value=1.0``,
  ``contourfnum=10`` and ``contourskip=2``, so the filled levels are the
  ten equally spaced probabilities 0.1..1.0 and the line levels are every
  second one of those;
* the seasons and their bracket codes come from ``seasonsdict``, so the
  four subplot titles are fixed by the module's own declaration.

Five defects surfaced:

* B296: ``Plots._scatterdict`` passes ``size=0.5`` to
  ``seaborn.scatterplot``.  seaborn's marker-area parameter is ``s``;
  ``size`` is the name of a *semantic grouping variable*.  A scalar is
  accepted and silently ignored, so every marker is drawn at matplotlib's
  default area instead of the requested 0.5 -- and because the stale
  ``size`` key stays in the dict, an explicit ``s`` from the caller is
  masked too.
* B297: ``dateLinePlot`` applies the field-specific ``set_ylabel`` from
  ``_plotfieldaxfuncdict`` and then unconditionally overwrites it with the
  bare column name.  Its two sibling methods (``plotScatter`` and
  ``plotProbContourf``) keep the readable label, so the same field is
  labelled ``'Wind Speed [m/s]'`` on a scatter plot and ``'WS'`` on a line
  plot.  A caller-supplied ``set_ylabel`` is discarded the same way.
* B298: ``_labelsdict['levels']`` is computed once in ``Plots.__init__``
  from the *default* ``_contourvalsdict`` and is never rebuilt from the
  levels actually drawn, so any ``contour_values`` override leaves the
  clabel levels outside ``CS.levels`` and ``ax.clabel`` raises
  ``ValueError``.  ``withLabels`` defaults to ``True``, so passing
  ``contour_values`` at all is enough to break the call.
* B299: in ``plotProbContourf``'s ``y_normalized`` branch, the
  "Contour levels must be increasing" stability clamp is applied to the
  contour and contourf levels but *not* to the clabel levels, which are
  rebuilt from the un-clamped ``M_hist.max()``.  For a histogram whose
  maximum falls below ``under_value`` the label levels come out decreasing
  and are not a subset of the contour levels, so ``ax.clabel`` raises.
* B300: ``plotProbContourf_bySeason`` documents an optional ``ax``, but
  it both indexes ``ax.shape``/``ax[i, j]`` (so it needs a 2x2 array) and
  hands the same object to ``plt.sca`` (which takes a single Axes).  An
  array fails in ``plt.sca``, a lone Axes fails on ``.shape``: the
  parameter cannot be used at all.
"""
import matplotlib.pyplot as plt
import numpy
import pandas
import pytest

from hera.measurements.meteorology.lowfreqdata.presentationLayer import (
    DailyPlots,
    SeasonalPlots,
    presenation,
)
from hera.measurements.meteorology.lowfreqdata.analysis import seasonsdict

FILLED_LEVELS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
LINE_LEVELS = [0.1, 0.3, 0.5, 0.7, 0.9]
MATPLOTLIB_DEFAULT_MARKER_AREA = 18.0


def _frame(hours, values, field="WS", day="2020-06-15"):
    """A frame with a ``datetime`` column at the given hours of one day."""
    base = pandas.Timestamp(day)
    stamps = pandas.DatetimeIndex(
        [base + pandas.Timedelta(seconds=round(h * 3600)) for h in hours])
    return pandas.DataFrame({"datetime": stamps, field: list(values)})


def _scattered(n=400, seed=1, field="WS"):
    """Enough well-spread samples, at distinct timestamps, for a histogram."""
    rng = numpy.random.default_rng(seed)
    seconds = rng.choice(24 * 3600, size=n, replace=False)
    return _frame(seconds / 3600.0, rng.uniform(0, 10, n), field=field)


@pytest.fixture()
def daily():
    return DailyPlots(presentation=None)


@pytest.fixture()
def seasonal():
    return SeasonalPlots(presentation=None)


def _offsets(ax):
    scatters = [c for c in ax.collections if type(c).__name__ == "PathCollection"]
    assert len(scatters) == 1, "expected exactly one scatter layer"
    return scatters[0].get_offsets()


# ---------------------------------------------------------------------------
# plotScatter
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPlotScatterGeometry:
    def test_each_sample_is_placed_at_its_fractional_hour_of_the_day(self, daily):
        ax = daily.plotScatter(_frame([1.0, 6.5, 12.25], [2.0, 4.0, 6.0]), "WS")
        assert _offsets(ax).tolist() == [[1.0, 2.0], [6.5, 4.0], [12.25, 6.0]]

    def test_the_hour_axis_covers_the_whole_day_with_one_tick_per_hour(self, daily):
        ax = daily.plotScatter(_frame([3.0, 4.0], [1.0, 2.0]), "WS")
        assert ax.get_xlim() == (0.0, 24.0)
        assert list(ax.get_xticks()) == list(range(25))
        assert ax.get_xlabel() == "Time [Hours]"

    def test_only_every_second_hour_is_labelled(self, daily):
        ax = daily.plotScatter(_frame([3.0, 4.0], [1.0, 2.0]), "WS")
        labels = [t.get_text() for t in ax.get_xticklabels()]
        assert labels[:5] == ["0", "", "2", "", "4"]

    @pytest.mark.parametrize("field, label", [
        ("WS", "Wind Speed [m/s]"),
        ("WD", "Wind Direction [° from N]"),
        ("RH", "Relative Humidity [%]"),
    ])
    def test_the_known_fields_get_their_readable_axis_label(self, daily, field, label):
        ax = daily.plotScatter(_frame([3.0, 4.0], [1.0, 2.0], field=field), field)
        assert ax.get_ylabel() == label

    def test_an_unknown_field_is_left_unlabelled(self, daily):
        ax = daily.plotScatter(_frame([3.0, 4.0], [1.0, 2.0], field="PRESSURE"), "PRESSURE")
        assert ax.get_ylabel() == ""

    def test_it_draws_into_a_supplied_axes_and_returns_it(self, daily):
        _, target = plt.subplots()
        ax = daily.plotScatter(_frame([3.0], [1.0]), "WS", ax=target)
        assert ax is target


@pytest.mark.unit
class TestPlotScatterCleaning:
    def test_the_minus_9999_sentinel_is_not_plotted(self, daily):
        ax = daily.plotScatter(_frame([1.0, 6.5, 12.0], [2.0, -9999.0, 6.0]), "WS")
        assert _offsets(ax).tolist() == [[1.0, 2.0], [12.0, 6.0]]

    def test_the_sentinel_threshold_is_exclusive_at_minus_5000(self, daily):
        kept = daily.plotScatter(_frame([1.0, 6.5], [2.0, -4999.0]), "WS")
        assert len(_offsets(kept)) == 2
        plt.close("all")
        dropped = daily.plotScatter(_frame([1.0, 6.5], [2.0, -5000.0]), "WS")
        assert len(_offsets(dropped)) == 1

    def test_a_repeated_timestamp_keeps_only_its_first_sample(self, daily):
        ax = daily.plotScatter(_frame([6.5, 6.5, 12.0], [1.0, 2.0, 3.0]), "WS")
        assert _offsets(ax).tolist() == [[6.5, 1.0], [12.0, 3.0]]

    def test_an_empty_frame_yields_bare_axes_rather_than_an_error(self, daily):
        ax = daily.plotScatter(_frame([], [], field="WS"), "WS")
        assert len(ax.collections) == 0
        assert ax.get_xlim() != (0.0, 24.0)

    def test_an_index_already_called_datetime_is_accepted(self, daily):
        data = _frame([1.0, 2.0], [3.0, 4.0]).set_index("datetime")
        ax = daily.plotScatter(data, "WS")
        assert _offsets(ax).tolist() == [[1.0, 3.0], [2.0, 4.0]]

    def test_a_frame_with_no_datetime_anywhere_is_rejected(self, daily):
        with pytest.raises(Exception, match="No column 'datetime'"):
            daily.plotScatter(pandas.DataFrame({"WS": [1.0, 2.0]}), "WS")


@pytest.mark.unit
class TestPlotScatterProperties:
    def test_the_caller_can_override_the_default_marker_colour(self, daily):
        ax = daily.plotScatter(_frame([1.0, 2.0], [3.0, 4.0]), "WS",
                               scatter_properties=dict(color="red"))
        assert ax.collections[0].get_facecolor()[0][:3] == pytest.approx((1.0, 0.0, 0.0))

    def test_the_caller_can_override_the_axis_functions(self, daily):
        ax = daily.plotScatter(_frame([1.0, 2.0], [3.0, 4.0]), "WS",
                               ax_functions_properties=dict(set_xlabel="hour of day"))
        assert ax.get_xlabel() == "hour of day"

    def test_the_default_scatter_dict_is_not_mutated_by_an_override(self, daily):
        daily.plotScatter(_frame([1.0], [3.0]), "WS", scatter_properties=dict(color="red"))
        assert daily._scatterdict["color"] == "k"

    @pytest.mark.xfail(
        strict=True,
        reason="B296: Plots._scatterdict asks seaborn.scatterplot for "
               "size=0.5, but seaborn's marker-area parameter is `s` -- "
               "`size` names a semantic grouping variable, so the scalar is "
               "accepted and ignored and the markers keep matplotlib's "
               "default area. See the consolidated findings issue.",
    )
    def test_the_markers_should_be_drawn_at_the_requested_half_point_size(self, daily):
        ax = daily.plotScatter(_frame([1.0, 2.0], [3.0, 4.0]), "WS")
        assert ax.collections[0].get_sizes().tolist() == [0.5, 0.5]

    def test_the_markers_currently_keep_the_matplotlib_default_size(self, daily):
        """Characterisation of B296."""
        ax = daily.plotScatter(_frame([1.0, 2.0], [3.0, 4.0]), "WS")
        sizes = ax.collections[0].get_sizes().tolist()
        assert sizes == [MATPLOTLIB_DEFAULT_MARKER_AREA] * 2

    def test_an_explicit_marker_area_is_masked_by_the_default_size_semantic(self, daily):
        """Characterisation of B296: `s` loses to the stale `size` key."""
        ax = daily.plotScatter(_frame([1.0, 2.0], [3.0, 4.0]), "WS",
                               scatter_properties=dict(s=0.5))
        assert ax.collections[0].get_sizes().tolist() == [
            MATPLOTLIB_DEFAULT_MARKER_AREA] * 2

    def test_clearing_the_size_semantic_lets_the_marker_area_through(self, daily):
        """Characterisation of B296: the fix would be `s` and no `size`."""
        ax = daily.plotScatter(_frame([1.0, 2.0], [3.0, 4.0]), "WS",
                               scatter_properties=dict(size=None, s=0.5))
        assert ax.collections[0].get_sizes().tolist() == [0.5]


# ---------------------------------------------------------------------------
# dateLinePlot
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestDateLinePlot:
    @pytest.fixture()
    def two_days(self):
        wanted = _frame([1.0, 6.5, 12.0], [2.0, 4.0, 6.0])
        other = _frame([3.0, 9.0], [90.0, 91.0], day="2020-06-16")
        return pandas.concat([wanted, other], ignore_index=True)

    def test_it_returns_the_axes_and_the_line_artists(self, daily, two_days):
        ax, line = daily.dateLinePlot(two_days, "WS", "2020-06-15")
        assert isinstance(line, list) and len(line) == 1
        assert line[0] in ax.lines

    def test_only_the_requested_day_is_drawn(self, daily, two_days):
        _, line = daily.dateLinePlot(two_days, "WS", "2020-06-15")
        assert line[0].get_xdata().tolist() == [1.0, 6.5, 12.0]
        assert line[0].get_ydata().tolist() == [2.0, 4.0, 6.0]

    def test_a_day_with_no_samples_draws_an_empty_line(self, daily, two_days):
        _, line = daily.dateLinePlot(two_days, "WS", "2020-06-17")
        assert line[0].get_xdata().tolist() == []

    def test_the_line_uses_the_declared_default_style(self, daily, two_days):
        _, line = daily.dateLinePlot(two_days, "WS", "2020-06-15")
        assert line[0].get_linewidth() == 3
        assert line[0].get_linestyle() == "-"
        assert line[0].get_zorder() == 4

    def test_line_properties_override_the_defaults(self, daily, two_days):
        _, line = daily.dateLinePlot(two_days, "WS", "2020-06-15",
                                     line_properties=dict(linewidth=1, linestyle=":"))
        assert line[0].get_linewidth() == 1
        assert line[0].get_linestyle() == ":"

    def test_the_title_names_the_field_and_the_day(self, daily, two_days):
        ax, _ = daily.dateLinePlot(two_days, "WS", "2020-06-15")
        assert ax.get_title() == "WS for 2020-06-15"

    def test_the_legend_is_labelled_with_the_field(self, daily, two_days):
        ax, _ = daily.dateLinePlot(two_days, "WS", "2020-06-15")
        assert [t.get_text() for t in ax.get_legend().get_texts()] == ["WS"]

    def test_the_legend_can_be_switched_off(self, daily, two_days):
        ax, _ = daily.dateLinePlot(two_days, "WS", "2020-06-15", legend=False)
        assert ax.get_legend() is None

    def test_the_hour_axis_covers_the_whole_day(self, daily, two_days):
        ax, _ = daily.dateLinePlot(two_days, "WS", "2020-06-15")
        assert ax.get_xlim() == (0.0, 24.0)
        assert ax.get_xlabel() == "Time [Hours]"

    def test_the_minus_9999_sentinel_becomes_a_gap(self, daily):
        data = _frame([1.0, 6.5, 12.0], [2.0, -9999.0, 6.0])
        _, line = daily.dateLinePlot(data, "WS", "2020-06-15")
        ydata = line[0].get_ydata()
        assert numpy.isnan(ydata[1])
        assert ydata[[0, 2]].tolist() == [2.0, 6.0]

    def test_a_frame_with_no_datetime_anywhere_is_rejected(self, daily):
        with pytest.raises(Exception, match="No column 'datetime'"):
            daily.dateLinePlot(pandas.DataFrame({"WS": [1.0]}), "WS", "2020-06-15")

    @pytest.mark.xfail(
        strict=True,
        reason="B297: dateLinePlot applies the field-specific set_ylabel "
               "from _plotfieldaxfuncdict and then unconditionally "
               "overwrites it with the bare column name, so the same field "
               "is labelled 'Wind Speed [m/s]' by plotScatter and "
               "plotProbContourf but 'WS' here. See the consolidated "
               "findings issue.",
    )
    def test_the_value_axis_should_keep_the_readable_field_label(self, daily, two_days):
        ax, _ = daily.dateLinePlot(two_days, "WS", "2020-06-15")
        assert ax.get_ylabel() == "Wind Speed [m/s]"

    def test_the_value_axis_is_currently_labelled_with_the_raw_column_name(self, daily,
                                                                          two_days):
        """Characterisation of B297."""
        ax, _ = daily.dateLinePlot(two_days, "WS", "2020-06-15")
        assert ax.get_ylabel() == "WS"

    def test_an_explicit_ylabel_is_still_overwritten(self, daily, two_days):
        """Characterisation of B297: even ax_functions_properties loses."""
        ax, _ = daily.dateLinePlot(two_days, "WS", "2020-06-15",
                                   ax_functions_properties=dict(set_ylabel="mine"))
        assert ax.get_ylabel() == "WS"


# ---------------------------------------------------------------------------
# plotProbContourf
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPlotProbContourf:
    def test_it_returns_the_line_set_the_filled_set_and_the_axes(self, daily):
        CS, CFS, ax = daily.plotProbContourf(_scattered(), "WS")
        assert CS is not CFS
        assert CS in ax.collections or CS.axes is ax

    def test_the_filled_levels_are_the_ten_declared_probabilities(self, daily):
        _, CFS, _ = daily.plotProbContourf(_scattered(), "WS")
        assert numpy.allclose(CFS.levels, FILLED_LEVELS)

    def test_the_line_levels_are_every_second_filled_level(self, daily):
        CS, _, _ = daily.plotProbContourf(_scattered(), "WS")
        assert numpy.allclose(CS.levels, LINE_LEVELS)

    def test_explicit_levels_replace_both_sets(self, daily):
        CS, CFS, _ = daily.plotProbContourf(_scattered(), "WS", levels=[0.2, 0.4, 0.6])
        assert numpy.allclose(CS.levels, [0.2, 0.4, 0.6])
        assert numpy.allclose(CFS.levels, [0.2, 0.4, 0.6])

    def test_contour_values_shift_the_level_range(self, daily):
        _, CFS, _ = daily.plotProbContourf(
            _scattered(), "WS", withLabels=False,
            contour_values=dict(under_value=0.0, max_value=0.5, contourfnum=6))
        assert numpy.allclose(CFS.levels, [0.0, 0.1, 0.2, 0.3, 0.4, 0.5])

    @pytest.mark.xfail(
        strict=True,
        reason="B298: _labelsdict['levels'] is frozen in Plots.__init__ "
               "from the *default* _contourvalsdict and is never rebuilt "
               "from the effective levels, so any contour_values override "
               "leaves the clabel levels outside CS.levels and ax.clabel "
               "raises. withLabels defaults to True, so overriding "
               "contour_values at all is enough to break the call. See the "
               "consolidated findings issue.",
    )
    def test_overriding_the_contour_values_should_still_label_the_lines(self, daily):
        daily.plotProbContourf(
            _scattered(), "WS",
            contour_values=dict(under_value=0.0, max_value=0.5, contourfnum=6))

    def test_overriding_the_contour_values_currently_breaks_the_labels(self, daily):
        """Characterisation of B298."""
        with pytest.raises(ValueError, match="don't match available levels"):
            daily.plotProbContourf(
                _scattered(), "WS",
                contour_values=dict(under_value=0.0, max_value=0.5, contourfnum=6))

    def test_labels_properties_can_repair_the_mismatch(self, daily):
        CS, _, _ = daily.plotProbContourf(
            _scattered(), "WS",
            contour_values=dict(under_value=0.0, max_value=0.5, contourfnum=6),
            labels_properties=dict(levels=[0.2]))
        assert numpy.allclose(CS.levels, [0.0, 0.2, 0.4])

    def test_the_colorbar_adds_a_second_axes_to_the_figure(self, daily):
        _, _, ax = daily.plotProbContourf(_scattered(), "WS", colorbar=True)
        assert len(ax.figure.axes) == 2

    def test_without_a_colorbar_the_figure_holds_only_the_plot(self, daily):
        _, _, ax = daily.plotProbContourf(_scattered(), "WS", colorbar=False)
        assert len(ax.figure.axes) == 1

    def test_the_scatter_layer_is_added_by_default(self, daily):
        _, _, ax = daily.plotProbContourf(_scattered(), "WS", colorbar=False)
        assert len(_offsets(ax)) == 400

    def test_the_scatter_layer_can_be_switched_off(self, daily):
        _, _, ax = daily.plotProbContourf(_scattered(), "WS", scatter=False, colorbar=False)
        assert [c for c in ax.collections if type(c).__name__ == "PathCollection"] == []

    def test_the_axes_get_the_daily_hour_axis_and_the_field_label(self, daily):
        _, _, ax = daily.plotProbContourf(_scattered(), "WS", colorbar=False)
        assert ax.get_xlim() == (0.0, 24.0)
        assert ax.get_ylabel() == "Wind Speed [m/s]"

    def test_a_frame_with_no_datetime_anywhere_is_rejected(self, daily):
        with pytest.raises(Exception, match="No column 'datetime'"):
            daily.plotProbContourf(pandas.DataFrame({"WS": [1.0, 2.0]}), "WS")

    def test_the_presentations_own_analysis_object_is_never_consulted(self):
        """It calls analysis.calcHourlyDist(analysis, ...) on the class."""

        class Exploding:
            def calcHourlyDist(self, *args, **kwargs):
                raise AssertionError("the wired analysis object was used")

        plots = DailyPlots(presentation=presenation(dataLayer=None, analysis=Exploding()))
        _, CFS, _ = plots.plotProbContourf(_scattered(), "WS", colorbar=False)
        assert numpy.allclose(CFS.levels, FILLED_LEVELS)


@pytest.mark.unit
class TestPlotProbContourfYNormalized:
    @staticmethod
    def _spread_thinly():
        """One value seen once in each of twenty different hours.

        ``calcDist2d``'s ``y_normalized`` divides each value-bin row by its
        own sum, so this histogram's largest cell is 1/20 = 0.05 -- below
        the declared ``under_value`` of 0.1, which is exactly the case the
        stability clamp exists for.
        """
        return _frame(list(range(2, 22)), [5.0] * 20)

    def test_the_levels_follow_the_histogram_maximum(self, daily):
        _, CFS, _ = daily.plotProbContourf(_scattered(), "WS",
                                           normalization="y_normalized",
                                           withLabels=False, colorbar=False)
        assert CFS.levels[0] == pytest.approx(0.1)
        assert CFS.levels[-1] < 1.0

    def test_the_clamp_keeps_the_levels_strictly_increasing(self, daily):
        CS, CFS, _ = daily.plotProbContourf(self._spread_thinly(), "WS",
                                            normalization="y_normalized",
                                            withLabels=False, colorbar=False)
        assert (numpy.diff(CS.levels) > 0).all()
        assert (numpy.diff(CFS.levels) > 0).all()

    def test_the_clamp_lifts_the_top_level_a_tenth_above_the_floor(self, daily):
        CS, _, _ = daily.plotProbContourf(self._spread_thinly(), "WS",
                                          normalization="y_normalized",
                                          withLabels=False, colorbar=False)
        expected = numpy.linspace(0.1, 0.2, 10)[::2]
        assert numpy.allclose(CS.levels, expected)

    @pytest.mark.xfail(
        strict=True,
        reason="B299: plotProbContourf's y_normalized stability clamp is "
               "applied to the contour and contourf levels but not to the "
               "clabel levels, which are rebuilt from the un-clamped "
               "M_hist.max(). When the histogram maximum falls below "
               "under_value the label levels come out decreasing and are "
               "not a subset of the contour levels, so ax.clabel raises -- "
               "and withLabels defaults to True. See the consolidated "
               "findings issue.",
    )
    def test_a_low_variance_histogram_should_still_be_labelled(self, daily):
        daily.plotProbContourf(self._spread_thinly(), "WS",
                               normalization="y_normalized", colorbar=False)

    def test_a_low_variance_histogram_currently_fails_to_label(self, daily):
        """Characterisation of B299."""
        with pytest.raises(ValueError, match="don't match available levels"):
            daily.plotProbContourf(self._spread_thinly(), "WS",
                                   normalization="y_normalized", colorbar=False)

    def test_supplying_the_label_levels_explicitly_works_around_it(self, daily):
        CS, _, _ = daily.plotProbContourf(self._spread_thinly(), "WS",
                                          normalization="y_normalized",
                                          labels_properties=dict(levels=[0.1]),
                                          colorbar=False)
        assert CS.levels[0] == pytest.approx(0.1)


# ---------------------------------------------------------------------------
# plotProbContourf_bySeason
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPlotProbContourfBySeason:
    @staticmethod
    def _year(seed=2, n=400):
        rng = numpy.random.default_rng(seed)
        stamps = pandas.date_range("2020-01-01", periods=n, freq="17h")
        return pandas.DataFrame({"datetime": stamps, "WS": rng.uniform(0, 10, n)})

    def test_it_returns_a_two_by_two_grid_of_axes(self, seasonal):
        ax = seasonal.plotProbContourf_bySeason(self._year(), "WS")
        assert ax.shape == (2, 2)

    def test_each_panel_is_titled_with_its_season_and_month_code(self, seasonal):
        ax = seasonal.plotProbContourf_bySeason(self._year(), "WS")
        expected = ["%s %s" % (season, seasonsdict[season]["strmonths"])
                    for season in seasonsdict]
        assert [a.get_title() for a in ax.ravel()] == expected

    def test_the_panels_are_filled_in_reading_order(self, seasonal):
        ax = seasonal.plotProbContourf_bySeason(self._year(), "WS")
        assert ax[0, 0].get_title().startswith("Winter")
        assert ax[1, 1].get_title().startswith("Autumn")

    def test_every_panel_carries_the_daily_hour_axis(self, seasonal):
        ax = seasonal.plotProbContourf_bySeason(self._year(), "WS")
        for panel in ax.ravel():
            assert panel.get_xlim() == (0.0, 24.0)
            assert panel.get_xlabel() == "Time [Hours]"

    def test_a_single_shared_colorbar_is_added_for_the_grid(self, seasonal):
        ax = seasonal.plotProbContourf_bySeason(self._year(), "WS")
        assert len(ax[0, 0].figure.axes) == 5

    def test_no_colorbar_leaves_only_the_four_panels(self, seasonal):
        ax = seasonal.plotProbContourf_bySeason(self._year(), "WS", colorbar=False)
        assert len(ax[0, 0].figure.axes) == 4

    @pytest.mark.xfail(
        strict=True,
        reason="B300: plotProbContourf_bySeason documents an optional `ax` "
               "but needs a 2x2 grid (it indexes ax.shape and "
               "ax[i, j]) while feeding that same argument to plt.sca, "
               "which only accepts a single Axes. An array raises "
               "AttributeError in plt.sca and a lone Axes raises "
               "AttributeError on .shape, so the parameter cannot be used "
               "at all. See the consolidated findings issue.",
    )
    def test_it_should_draw_into_a_supplied_grid_of_axes(self, seasonal):
        _, target = plt.subplots(2, 2)
        ax = seasonal.plotProbContourf_bySeason(self._year(), "WS", ax=target,
                                                colorbar=False)
        assert ax is target

    def test_a_supplied_grid_currently_breaks_in_plt_sca(self, seasonal):
        """Characterisation of B300."""
        _, target = plt.subplots(2, 2)
        with pytest.raises(AttributeError, match="no attribute 'figure'"):
            seasonal.plotProbContourf_bySeason(self._year(), "WS", ax=target,
                                               colorbar=False)

    def test_a_supplied_single_axes_currently_breaks_on_shape(self, seasonal):
        """Characterisation of B300."""
        _, target = plt.subplots()
        with pytest.raises(AttributeError, match="shape"):
            seasonal.plotProbContourf_bySeason(self._year(), "WS", ax=target,
                                               colorbar=False)

    def test_each_panel_scatters_only_its_own_seasons_samples(self, seasonal):
        data = self._year()
        ax = seasonal.plotProbContourf_bySeason(data, "WS", colorbar=False)
        months = data["datetime"].dt.month
        for panel, season in zip(ax.ravel(), seasonsdict):
            wanted = months.isin(seasonsdict[season]["months"]).sum()
            assert len(_offsets(panel)) == wanted

    def test_the_samples_are_partitioned_across_the_four_panels(self, seasonal):
        data = self._year()
        ax = seasonal.plotProbContourf_bySeason(data, "WS", colorbar=False)
        assert sum(len(_offsets(panel)) for panel in ax.ravel()) == len(data)

    def test_an_index_already_called_datetime_is_accepted(self, seasonal):
        data = self._year().set_index("datetime")
        ax = seasonal.plotProbContourf_bySeason(data, "WS", colorbar=False)
        assert ax.shape == (2, 2)

    def test_a_frame_with_no_datetime_anywhere_is_rejected(self, seasonal):
        with pytest.raises(Exception, match="No column 'datetime'"):
            seasonal.plotProbContourf_bySeason(pandas.DataFrame({"WS": [1.0, 2.0]}), "WS")
