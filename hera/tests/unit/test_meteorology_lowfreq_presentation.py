"""measurements/meteorology/lowfreqdata/presentationLayer.py: the plot-dict
and colormap builders, plus the thin wiring classes (`presenation`, `Plots`).
The heavier SeasonalPlots/DailyPlots contour-plotting methods need real wind
data shaped a specific way and are left uncovered."""
import matplotlib
import numpy
import pytest

from hera.measurements.meteorology.lowfreqdata.presentationLayer import Plots, presenation


@pytest.mark.unit
class TestPresenationWiring:
    def test_it_stores_the_datalayer_and_analysis(self):
        p = presenation(dataLayer="myDataLayer", analysis="myAnalysis")
        assert p.datalayer == "myDataLayer"
        assert p.analysis == "myAnalysis"

    def test_it_builds_a_seasonal_and_a_daily_plots_object(self):
        p = presenation(dataLayer="myDataLayer", analysis="myAnalysis")
        assert p.seasonalPlots is not None
        assert p.dailyPlots is not None
        assert p.seasonalPlots.presenatation is p
        assert p.dailyPlots.presenatation is p


@pytest.mark.unit
class TestPlotsConstruction:
    def test_it_stores_a_back_reference_to_the_presentation(self):
        plots = Plots(presentation="fakePresentation")
        assert plots.presenatation == "fakePresentation"

    def test_the_default_contour_value_dict_has_the_expected_keys(self):
        plots = Plots(presentation=None)
        assert plots._contourvalsdict == dict(
            under_value=0.1, contourskip=2, contourfnum=10, max_value=1.0
        )


@pytest.mark.unit
class TestGetCountourDict:
    def test_it_builds_a_black_line_contour_spec(self):
        plots = Plots(presentation=None)
        params = dict(under_value=0.0, max_value=1.0, contourfnum=5)
        spec = plots._getCountourDict(params)
        assert spec["colors"] == "k"
        assert spec["linewidths"] == 0.5
        assert spec["zorder"] == 3

    def test_levels_are_evenly_spaced_between_under_and_max(self):
        plots = Plots(presentation=None)
        params = dict(under_value=0.0, max_value=1.0, contourfnum=5)
        spec = plots._getCountourDict(params)
        assert list(spec["levels"]) == pytest.approx([0.0, 0.5, 1.0])  # every 2nd of 5 levels


@pytest.mark.unit
class TestGetContourfDict:
    def test_it_builds_a_filled_contour_spec_with_a_colormap(self):
        plots = Plots(presentation=None)
        params = dict(under_value=0.0, max_value=1.0, contourfnum=5)
        spec = plots._getContourfDict(params)
        assert spec["extend"] == "min"
        assert spec["zorder"] == 2
        assert isinstance(spec["cmap"], matplotlib.colors.Colormap)

    def test_all_five_levels_are_present_unlike_the_line_contour(self):
        plots = Plots(presentation=None)
        params = dict(under_value=0.0, max_value=1.0, contourfnum=5)
        spec = plots._getContourfDict(params)
        assert len(spec["levels"]) == 5


@pytest.mark.unit
class TestGetCmap:
    def test_the_default_cmap_has_no_special_under_or_over_color(self):
        plots = Plots(presentation=None)
        cmap = plots._getcmap()
        default = matplotlib.colormaps["jet"]
        assert cmap(-1.0) == default(-1.0)

    def test_setting_under_changes_the_below_range_color(self):
        plots = Plots(presentation=None)
        cmap = plots._getcmap(under=True, undercolor="0.9")
        under_rgba = cmap.get_under()
        assert under_rgba[:3] == pytest.approx((0.9, 0.9, 0.9))

    def test_setting_over_changes_the_above_range_color(self):
        plots = Plots(presentation=None)
        cmap = plots._getcmap(over=True, overcolor="0.1")
        over_rgba = cmap.get_over()
        assert over_rgba[:3] == pytest.approx((0.1, 0.1, 0.1))

    def test_it_does_not_mutate_the_shared_registered_jet_colormap(self):
        plots = Plots(presentation=None)
        before = matplotlib.colormaps["jet"].get_under()
        plots._getcmap(under=True, undercolor="0.1")
        after = matplotlib.colormaps["jet"].get_under()
        assert numpy.array_equal(before, after)
