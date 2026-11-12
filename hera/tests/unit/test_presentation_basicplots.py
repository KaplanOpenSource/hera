"""presentation/basicplots.py: PlotFields, the generic time-series /
scatter / stacked-bar helpers.

These methods draw and return nothing, so the tests assert on the state
of the resulting matplotlib axes (lines drawn, labels, limits, legend)
under the Agg backend the unit conftest installs.

Two defects surfaced:

* B143: ``plotFields`` cannot plot a single field. It always calls
  ``plt.subplots(1, plotsNum)``, which for ``plotsNum == 1`` returns a
  bare ``Axes`` rather than an array, and then indexes it as
  ``axes[i]`` -- ``TypeError: 'Axes' object is not subscriptable``. Only
  two or more fields work.
* B144: ``plotFields`` accepts an ``axes`` argument and immediately
  discards it -- the very next statement rebinds ``axes`` to the result
  of ``plt.subplots(...)``. A caller handing in prepared axes is silently
  ignored and gets a brand-new figure instead. (Same family as B113 in
  the experiment presentation layer.)
"""
import matplotlib.pyplot as plt
import numpy
import pandas
import pytest

from hera.presentation.basicplots import PlotFields


@pytest.fixture(autouse=True)
def _fresh_figure():
    plt.close("all")
    yield
    plt.close("all")


@pytest.fixture()
def frame():
    index = pandas.date_range("2020-01-01", periods=60, freq="1min")
    return pandas.DataFrame(
        {"u": numpy.arange(60.0), "w": numpy.arange(60.0) * 2}, index=index,
    )


@pytest.fixture()
def plotter(frame):
    return PlotFields(data=frame)


@pytest.mark.unit
class TestConstruction:
    def test_the_data_is_stored(self, frame):
        assert PlotFields(data=frame)._Data is frame

    def test_it_can_be_built_without_data(self):
        assert PlotFields()._Data is None


@pytest.mark.unit
class TestPlotField:
    def test_the_raw_series_is_drawn_when_asked(self, plotter, frame):
        plotter.plotField("u", frame.index[0], frame.index[-1], plotRaw=True)
        lines = plt.gca().get_lines()
        assert len(lines) == 1
        assert list(lines[0].get_ydata()) == list(frame["u"])

    def test_nothing_is_drawn_without_raw_or_averages(self, plotter, frame):
        plotter.plotField("u", frame.index[0], frame.index[-1])
        assert plt.gca().get_lines() == []

    def test_each_resample_period_adds_its_own_line(self, plotter, frame):
        plotter.plotField(
            "u", frame.index[0], frame.index[-1], resampleList=["5min", "10min"],
        )
        assert len(plt.gca().get_lines()) == 2

    def test_a_resampled_line_is_shorter_than_the_raw_series(self, plotter, frame):
        plotter.plotField("u", frame.index[0], frame.index[-1], resampleList=["10min"])
        assert len(plt.gca().get_lines()[0].get_ydata()) == 6

    def test_each_moving_window_adds_its_own_line(self, plotter, frame):
        plotter.plotField("u", frame.index[0], frame.index[-1], movingList=[5])
        assert len(plt.gca().get_lines()) == 1

    def test_the_date_range_is_honoured(self, plotter, frame):
        plotter.plotField("u", frame.index[10], frame.index[20], plotRaw=True)
        assert len(plt.gca().get_lines()[0].get_ydata()) == 11

    def test_the_y_label_is_applied(self, plotter, frame):
        plotter.plotField("u", frame.index[0], frame.index[-1], yLabel="u [m/s]", plotRaw=True)
        assert plt.gca().get_ylabel() == "u [m/s]"

    def test_a_legend_is_added_only_when_requested(self, plotter, frame):
        plotter.plotField("u", frame.index[0], frame.index[-1], plotRaw=True, legend=True)
        assert plt.gca().get_legend() is not None

    def test_it_draws_on_the_axes_it_is_given(self, plotter, frame):
        _, ax = plt.subplots()
        plotter.plotField("u", frame.index[0], frame.index[-1], ax=ax, plotRaw=True)
        assert len(ax.get_lines()) == 1

    def test_saving_writes_a_file(self, plotter, frame, tmp_path):
        target = tmp_path / "field.png"
        plotter.plotField(
            "u", frame.index[0], frame.index[-1], plotRaw=True,
            save=True, saveTo=str(target),
        )
        assert target.is_file()


@pytest.mark.unit
class TestPlotFieldsMultiple:
    def test_two_fields_produce_two_subplots(self, plotter, frame):
        plotter.plotFields(["u", "w"], frame.index[0], frame.index[-1], plotRaw=True)
        assert len(plt.gcf().axes) == 2

    def test_each_subplot_gets_its_own_field(self, plotter, frame):
        plotter.plotFields(["u", "w"], frame.index[0], frame.index[-1], plotRaw=True)
        first, second = plt.gcf().axes
        assert list(first.get_lines()[0].get_ydata()) == list(frame["u"])
        assert list(second.get_lines()[0].get_ydata()) == list(frame["w"])

    def test_the_y_labels_are_distributed_per_field(self, plotter, frame):
        plotter.plotFields(
            ["u", "w"], frame.index[0], frame.index[-1], plotRaw=True, yLabels=["U", "W"],
        )
        assert [ax.get_ylabel() for ax in plt.gcf().axes] == ["U", "W"]

    def test_saving_writes_a_file(self, plotter, frame, tmp_path):
        target = tmp_path / "fields.png"
        plotter.plotFields(
            ["u", "w"], frame.index[0], frame.index[-1], plotRaw=True,
            save=True, saveTo=str(target),
        )
        assert target.is_file()


@pytest.mark.unit
class TestPlotFieldsSingleFieldIsBroken:
    """B143: see the module docstring."""

    @pytest.mark.xfail(
        strict=True,
        reason="B143: plotFields calls plt.subplots(1, 1) for a single "
               "field, which returns a bare Axes, then indexes it as "
               "axes[i]. See the consolidated findings issue.",
    )
    def test_a_single_field_should_plot(self, plotter, frame):
        plotter.plotFields(["u"], frame.index[0], frame.index[-1], plotRaw=True)

    def test_a_single_field_currently_raises(self, plotter, frame):
        """Characterisation of B143."""
        with pytest.raises(TypeError, match="not subscriptable"):
            plotter.plotFields(["u"], frame.index[0], frame.index[-1], plotRaw=True)

    def test_two_fields_are_the_smallest_working_call(self, plotter, frame):
        """Characterisation of B143: the bug is specific to plotsNum == 1."""
        plotter.plotFields(["u", "w"], frame.index[0], frame.index[-1], plotRaw=True)
        assert len(plt.gcf().axes) == 2


@pytest.mark.unit
class TestPlotFieldsIgnoresItsAxes:
    """B144: see the module docstring."""

    @pytest.mark.xfail(
        strict=True,
        reason="B144: plotFields rebinds its `axes` parameter to the "
               "result of plt.subplots on the very next line, so caller "
               "supplied axes are silently discarded. See the "
               "consolidated findings issue.",
    )
    def test_the_axes_that_were_passed_in_should_be_used(self, plotter, frame):
        fig, axes = plt.subplots(1, 2)
        plotter.plotFields(
            ["u", "w"], frame.index[0], frame.index[-1], axes=axes, plotRaw=True,
        )
        assert len(axes[0].get_lines()) == 1

    def test_the_axes_that_were_passed_in_are_currently_left_empty(self, plotter, frame):
        """Characterisation of B144."""
        fig, axes = plt.subplots(1, 2)
        plotter.plotFields(
            ["u", "w"], frame.index[0], frame.index[-1], axes=axes, plotRaw=True,
        )
        assert axes[0].get_lines() == []
        assert axes[1].get_lines() == []

    def test_a_brand_new_figure_is_created_instead(self, plotter, frame):
        """Characterisation of B144."""
        first = plt.subplots(1, 2)[0]
        plotter.plotFields(
            ["u", "w"], frame.index[0], frame.index[-1], axes=first.axes, plotRaw=True,
        )
        assert plt.gcf() is not first


@pytest.mark.unit
class TestTest:
    """`PlotFields.test` is a debug helper that just echoes its arguments."""

    def test_it_echoes_every_argument_it_was_given(self, plotter, capsys):
        plotter.test(["u"], "2020-01-01", "2020-01-02", resampleList=["5min"], yLabels=["U"])
        out = capsys.readouterr().out
        assert "fieldNames: ['u']" in out
        assert "resampleList: ['5min']" in out
        assert "yLabels: ['U']" in out

    def test_it_returns_nothing(self, plotter):
        assert plotter.test(["u"], "2020-01-01", "2020-01-02") is None


@pytest.mark.unit
class TestPlotFrequencyBar:
    @pytest.fixture()
    def categorical(self):
        index = pandas.date_range("2020-01-01", periods=24, freq="1h")
        return pandas.DataFrame(
            {"hour": index.hour, "sector": ["N", "S"] * 12}, index=index,
        )

    def test_it_draws_one_stacked_bar_group_per_category(self, plotter, categorical):
        plotter.plot_frequency_bar(categorical, category="sector", grouped_by="hour")
        assert len(plt.gca().containers) == 2

    def test_the_bars_are_normalised_to_fractions(self, plotter, categorical):
        """Each group is divided by its own row sum, so every bar totals 1."""
        plotter.plot_frequency_bar(categorical, category="sector", grouped_by="hour")
        heights = [
            sum(bar.get_height() for bar in container) for container in plt.gca().containers
        ]
        assert sum(heights) == pytest.approx(len(categorical["hour"].unique()))

    def test_an_explicit_category_list_adds_missing_categories_as_zero(self, plotter, categorical):
        plotter.plot_frequency_bar(
            categorical, category="sector", grouped_by="hour", cat_list=["N", "S", "E"],
        )
        assert len(plt.gca().containers) == 3

    def test_it_draws_on_the_axes_it_is_given(self, plotter, categorical):
        _, ax = plt.subplots()
        plotter.plot_frequency_bar(categorical, category="sector", grouped_by="hour", axis=ax)
        assert len(ax.containers) == 2

    def test_axis_properties_are_applied_by_name(self, plotter, categorical):
        plotter.plot_frequency_bar(
            categorical, category="sector", grouped_by="hour",
            ax_props={"set_title": {"label": "sectors"}},
        )
        assert plt.gca().get_title() == "sectors"


@pytest.mark.unit
class TestPlotScatterFunc:
    def test_vectors_are_scattered_directly_when_no_frame_is_given(self, plotter):
        plotter.plot_scatter_func(x=[1.0, 2.0, 3.0], y=[2.0, 4.0, 6.0])
        assert len(plt.gca().collections) == 1

    def test_column_names_are_resolved_against_the_frame(self, plotter):
        data = pandas.DataFrame({"a": [1.0, 2.0, 3.0], "b": [1.0, 4.0, 9.0]})
        plotter.plot_scatter_func(x="a", y="b", data=data)
        assert len(plt.gca().collections) == 1

    def test_the_self_correlation_overlay_adds_a_second_scatter(self, plotter):
        plotter.plot_scatter_func(
            x=[1.0, 2.0, 3.0], y=[2.0, 4.0, 6.0], SelfCorrelation=True,
        )
        assert len(plt.gca().collections) == 2

    def test_a_function_is_drawn_either_side_of_the_origin(self, plotter):
        """The x range is split at zero into two arrays, so one function
        yields two lines sharing a colour."""
        plotter.plot_scatter_func(
            x=[-1.0, 1.0], y=[1.0, 1.0], func_dict={0: (lambda z: z ** 2, "z^2")},
        )
        lines = plt.gca().get_lines()
        assert len(lines) == 2
        assert lines[0].get_color() == lines[1].get_color()

    def test_axis_properties_are_applied_by_name(self, plotter):
        plotter.plot_scatter_func(
            x=[1.0, 2.0], y=[1.0, 2.0], ax_props={"set_xlabel": {"xlabel": "a"}},
        )
        assert plt.gca().get_xlabel() == "a"
