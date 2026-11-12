"""simulations/gaussian/toolkit.py: the presentation layer's concentration,
dosage and TIAC plots.

Each method selects a slice out of an xarray concentration field and draws
one line, so the tests build a small synthetic field with the real dim
names (x, y, z, time) and assert on the resulting matplotlib axes -- the
labels, title and the plotted data -- rather than on the return value
(these methods return nothing and end in ``plt.show()``, which is inert
under the Agg backend the unit conftest installs).

Note the two families disagree about what ``attrs['Q']`` is: the plain
methods do ``str(C.attrs['Q']).split(' ')[1]`` (a string form), while the
``_noQ`` ones do ``C.attrs['Q'].units`` (a live pint Quantity). A pint
Quantity satisfies both, so that is what the fixture supplies.

``getSpaceTime``/``getGasCloud`` are not covered here -- they drive the
full dispersion solve and belong with the physics tests.
"""
import matplotlib.pyplot as plt
import numpy
import pytest
import xarray

from hera.utils.unitHandler import ureg


@pytest.fixture(autouse=True)
def _fresh_figure():
    plt.close("all")
    yield
    plt.close("all")


@pytest.fixture()
def presentation(unit_toolkit_factory):
    from hera import toolkitHome

    return unit_toolkit_factory(toolkitHome.GAUSSIANDISPERSION).presentation


@pytest.fixture()
def field():
    """A 4x3x2x5 (x, y, z, time) field whose values rise with x and fall
    with time, so the plots have something monotone to show."""
    x = numpy.array([0.0, 100.0, 200.0, 300.0])
    y = numpy.array([-50.0, 0.0, 50.0])
    z = numpy.array([0.0, 10.0])
    time = numpy.array([0.0, 1.0, 2.0, 3.0, 4.0])
    values = (
        x[:, None, None, None]
        + numpy.zeros((1, len(y), len(z), len(time)))
        + (len(time) - time)[None, None, None, :]
    )
    return xarray.DataArray(
        values, dims=("x", "y", "z", "time"),
        coords={"x": x, "y": y, "z": z, "time": time},
        attrs={"Q": 1.0 * ureg.kg},
    )


def _line():
    """The single line drawn on the current axes, as (xdata, ydata)."""
    lines = plt.gca().get_lines()
    assert len(lines) == 1, f"expected exactly one line, got {len(lines)}"
    return lines[0].get_xdata(), lines[0].get_ydata()


@pytest.mark.unit
class TestPlotMaxConcentrationOverTime:
    def test_it_plots_the_time_maximum_against_distance(self, presentation, field):
        presentation.plotMaxConcentrationOverTime(field, y=0 * ureg.m, z=0 * ureg.m)
        xdata, ydata = _line()
        assert list(xdata) == [0.0, 100.0, 200.0, 300.0]
        expected = field.sel(y=0.0, z=0.0).max(dim="time").values
        assert ydata == pytest.approx(expected)

    def test_it_labels_the_axes_with_distance_and_the_mass_unit(self, presentation, field):
        presentation.plotMaxConcentrationOverTime(field, y=0 * ureg.m, z=0 * ureg.m)
        ax = plt.gca()
        assert "Distance from source" in ax.get_xlabel()
        assert "kilogram" in ax.get_ylabel()

    def test_the_title_records_the_selected_line(self, presentation, field):
        presentation.plotMaxConcentrationOverTime(field, y=50 * ureg.m, z=10 * ureg.m)
        assert "y=50[m], z=10[m]" in plt.gca().get_title()

    def test_an_explicit_x_window_clamps_the_axis(self, presentation, field):
        presentation.plotMaxConcentrationOverTime(
            field, y=0 * ureg.m, z=0 * ureg.m, x_min=50 * ureg.m, x_max=250 * ureg.m,
        )
        assert plt.gca().get_xlim() == pytest.approx((50.0, 250.0))

    def test_a_height_in_other_units_is_converted(self, presentation, field):
        """z given in centimetres must land on the nearest metre-based
        grid level, not be treated as 1000 m."""
        presentation.plotMaxConcentrationOverTime(field, y=0 * ureg.m, z=1000 * ureg.cm)
        assert "z=10.0[m]" in plt.gca().get_title()


@pytest.mark.unit
class TestPlotFixedPointConcentrationOverTime:
    def test_it_plots_the_series_at_one_point_against_time(self, presentation, field):
        presentation.plotFixedPointConcentrationOverTime(
            field, x=100 * ureg.m, y=0 * ureg.m, z=0 * ureg.m,
        )
        xdata, ydata = _line()
        assert list(xdata) == [0.0, 1.0, 2.0, 3.0, 4.0]
        expected = field.sel(x=100.0, y=0.0, z=0.0).values
        assert ydata == pytest.approx(expected)

    def test_it_labels_the_time_axis(self, presentation, field):
        presentation.plotFixedPointConcentrationOverTime(
            field, x=100 * ureg.m, y=0 * ureg.m, z=0 * ureg.m,
        )
        assert "Time from release" in plt.gca().get_xlabel()

    def test_the_title_records_the_detector_position(self, presentation, field):
        presentation.plotFixedPointConcentrationOverTime(
            field, x=200 * ureg.m, y=50 * ureg.m, z=10 * ureg.m,
        )
        assert "x=200[m], y=50[m], z=10[m]" in plt.gca().get_title()

    def test_an_explicit_time_window_clamps_the_axis(self, presentation, field):
        presentation.plotFixedPointConcentrationOverTime(
            field, x=100 * ureg.m, y=0 * ureg.m, z=0 * ureg.m,
            t_min=1 * ureg.min, t_max=3 * ureg.min,
        )
        assert plt.gca().get_xlim() == pytest.approx((1.0, 3.0))


@pytest.mark.unit
class TestPlotDosagePerDistance:
    def test_it_plots_one_time_slice_against_distance(self, presentation, field):
        presentation.plotDosagePerDistance(
            field, y=0 * ureg.m, z=0 * ureg.m, time=2 * ureg.min,
        )
        xdata, ydata = _line()
        assert list(xdata) == [0.0, 100.0, 200.0, 300.0]
        expected = field.sel(y=0.0, z=0.0, time=2.0).values
        assert ydata == pytest.approx(expected)

    def test_it_labels_the_dosage_axis(self, presentation, field):
        presentation.plotDosagePerDistance(
            field, y=0 * ureg.m, z=0 * ureg.m, time=2 * ureg.min,
        )
        assert "Dosage" in plt.gca().get_ylabel()

    def test_the_title_records_the_time_slice(self, presentation, field):
        presentation.plotDosagePerDistance(
            field, y=0 * ureg.m, z=0 * ureg.m, time=3 * ureg.min,
        )
        assert "time=3[min]" in plt.gca().get_title()


@pytest.mark.unit
class TestPlotFixedPointDosageOverTime:
    def test_it_plots_the_dosage_series_at_one_point(self, presentation, field):
        presentation.plotFixedPointDosageOverTime(
            field, x=100 * ureg.m, y=0 * ureg.m, z=0 * ureg.m,
        )
        xdata, ydata = _line()
        assert list(xdata) == [0.0, 1.0, 2.0, 3.0, 4.0]
        expected = field.sel(x=100.0, y=0.0, z=0.0).values
        assert ydata == pytest.approx(expected)

    def test_it_labels_the_axis_as_tiac_over_time(self, presentation, field):
        presentation.plotFixedPointDosageOverTime(
            field, x=100 * ureg.m, y=0 * ureg.m, z=0 * ureg.m,
        )
        assert "TIAC over time" in plt.gca().get_ylabel()


@pytest.mark.unit
class TestNoQVariants:
    """The ``_noQ`` family renders the unit from the live pint Quantity in
    ``attrs['Q']`` instead of a string split, and defaults its x/time
    window to the extent of the data rather than leaving it unclamped."""

    def test_max_concentration_plots_the_time_maximum(self, presentation, field):
        presentation.plotMaxConcentrationOverTime_noQ(field, y=0 * ureg.m, z=0 * ureg.m)
        xdata, ydata = _line()
        expected = field.sel(y=0.0, z=0.0).max(dim="time").values
        assert ydata == pytest.approx(expected)

    def test_max_concentration_defaults_the_window_to_the_data_extent(self, presentation, field):
        presentation.plotMaxConcentrationOverTime_noQ(field, y=0 * ureg.m, z=0 * ureg.m)
        assert plt.gca().get_xlim() == pytest.approx((0.0, 300.0))

    def test_max_concentration_renders_the_unit_from_the_quantity(self, presentation, field):
        presentation.plotMaxConcentrationOverTime_noQ(field, y=0 * ureg.m, z=0 * ureg.m)
        assert "Concentration" in plt.gca().get_ylabel()

    def test_fixed_point_concentration_plots_the_series_at_one_point(self, presentation, field):
        presentation.plotFixedPointConcentrationOverTime_noQ(
            field, x=100 * ureg.m, y=0 * ureg.m, z=0 * ureg.m,
        )
        xdata, ydata = _line()
        expected = field.sel(x=100.0, y=0.0, z=0.0).values
        assert ydata == pytest.approx(expected)

    def test_tiac_per_distance_plots_one_time_slice(self, presentation, field):
        presentation.plotTIACPerDistance_noQ(
            field, y=0 * ureg.m, z=0 * ureg.m, time=2 * ureg.min,
        )
        xdata, ydata = _line()
        expected = field.sel(y=0.0, z=0.0, time=2.0).values
        assert ydata == pytest.approx(expected)

    def test_fixed_point_tiac_plots_the_series_at_one_point(self, presentation, field):
        presentation.plotFixedPointTIACOverTime_noQ(
            field, x=200 * ureg.m, y=0 * ureg.m, z=0 * ureg.m,
        )
        xdata, ydata = _line()
        expected = field.sel(x=200.0, y=0.0, z=0.0).values
        assert ydata == pytest.approx(expected)

    def test_every_noq_plot_selects_the_nearest_grid_point(self, presentation, field):
        """An off-grid request must snap, not interpolate or fail."""
        presentation.plotFixedPointConcentrationOverTime_noQ(
            field, x=97 * ureg.m, y=1 * ureg.m, z=1 * ureg.m,
        )
        _, ydata = _line()
        expected = field.sel(x=100.0, y=0.0, z=0.0).values
        assert ydata == pytest.approx(expected)
