"""measurements/experiment/presentation.py: experimentPresentation's
construction, properties and pure helpers.

``experimentPresentation.__init__`` takes the datalayer and analysis
objects and does nothing with them but store them, so the tests below
drive it with a small stand-in exposing only the two attributes the
methods under test actually read. The plotting methods that need a real
loaded experiment (``plotImage``, ``plotMap``, ``plotDevices``,
``plotDevicesOnImage``, ``_plotEntityLocation*``,
``plotDeviceTypeFunctionality``) need a populated argos experiment
directory and map tiles, which is integration-shaped and out of scope.

Two defects surfaced:

* B113: ``plotOrigin`` accepts an ``ax`` argument and returns it, but
  draws with the module-level ``plt.scatter(...)`` instead of
  ``ax.scatter(...)``. When the caller passes an axes that is not the
  current one, the marker silently lands on whichever axes matplotlib
  happens to consider current, and the returned axes stays empty.
* B114: ``generateLatexTable`` calls
  ``self.plotDevices(trialSetName=..., trialName=..., device=...,
  display=False)``, but ``plotDevices``'s signature is
  ``(trialSetName, trialName, deviceType, mapName, ax=None,
  plotkwargs=None)`` -- there is no ``device`` or ``display`` parameter
  and ``mapName`` is mandatory. Every call raises ``TypeError`` before
  producing anything, so the whole LaTeX export path is dead.
"""
import matplotlib.pyplot as plt
import pandas
import pytest

from hera.measurements.experiment.presentation import experimentPresentation


class _FakeEntities:
    def __init__(self, table):
        self.entitiesTable = table


class _FakeDatalayer:
    """Exposes only what the methods under test read: `setup` (trial sets)
    and `trialSet` (the entities table)."""

    def __init__(self):
        table = pandas.DataFrame({
            "deviceTypeName": ["sonic", "sonic"],
            "deviceItemName": ["sonic_1", "sonic_2"],
            "stationName": ["st_1", "st_2"],
            "Latitude": [32.0, 32.1],
            "Longitude": [34.0, 34.1],
            "height": [3.0, 10.0],
        })
        self.setup = {"trialSets": [{"name": "TS", "trials": [{"name": "T1"}]}]}
        self.trialSet = {"Measurements": {"Measurements": _FakeEntities(table)}}


@pytest.fixture(autouse=True)
def _fresh_figure():
    plt.close("all")
    yield
    plt.close("all")


@pytest.fixture()
def datalayer():
    return _FakeDatalayer()


@pytest.fixture()
def presentation(datalayer):
    return experimentPresentation(datalayer=datalayer, analysis=object())


@pytest.mark.unit
class TestConstruction:
    def test_the_datalayer_and_analysis_are_stored(self, datalayer):
        analysis = object()
        p = experimentPresentation(datalayer=datalayer, analysis=analysis)
        assert p.datalayer is datalayer
        assert p.analysis is analysis

    def test_a_colormap_is_built_during_construction(self, presentation):
        from matplotlib.colors import LinearSegmentedColormap

        assert isinstance(presentation.cmap, LinearSegmentedColormap)

    def test_the_colormap_spans_red_through_green(self, presentation):
        """The level list runs 0 -> 1 over red, orange, green, so the two
        ends of the map must be the red and the green."""
        red = presentation.cmap(0.0)
        green = presentation.cmap(1.0)
        assert red[0] > red[1] and red[0] > red[2]
        assert green[1] > green[0] and green[1] > green[2]


@pytest.mark.unit
class TestSaveProperties:
    def test_savefigures_round_trips(self, presentation):
        assert presentation.saveFigures is None
        presentation.saveFigures = True
        assert presentation.saveFigures is True

    def test_savepath_is_returned_as_an_absolute_path(self, presentation, tmp_path):
        presentation.savePath = str(tmp_path)
        assert presentation.savePath == str(tmp_path)

    def test_a_relative_savepath_is_absolutised(self, presentation):
        presentation.savePath = "figures"
        import os

        assert presentation.savePath == os.path.abspath("figures")


@pytest.mark.unit
class TestSplitName:
    def test_a_two_token_name_yields_the_numeric_suffix(self, presentation):
        assert presentation._splitName("sonic 12") == 12

    def test_a_single_token_name_is_returned_whole(self, presentation):
        assert presentation._splitName("sonic") == "sonic"

    def test_a_non_numeric_second_token_raises(self, presentation):
        with pytest.raises(ValueError):
            presentation._splitName("sonic left")


@pytest.mark.unit
class TestScatterHeightColour:
    def test_it_sets_the_colour_range_from_the_heights(self, presentation, datalayer):
        table = datalayer.trialSet["Measurements"]["Measurements"].entitiesTable
        attrs = presentation._scatter_height_color(table, {})
        assert attrs["vmin"] == pytest.approx(3.0)
        assert attrs["vmax"] == pytest.approx(10.0)

    def test_an_empty_table_collapses_the_range_to_zero(self, presentation):
        empty = pandas.DataFrame({"height": []})
        attrs = presentation._scatter_height_color(empty, {})
        assert attrs["vmin"] == 0
        assert attrs["vmax"] == 0

    def test_existing_keys_are_preserved(self, presentation, datalayer):
        table = datalayer.trialSet["Measurements"]["Measurements"].entitiesTable
        attrs = presentation._scatter_height_color(table, {"s": 50})
        assert attrs["s"] == 50


@pytest.mark.unit
class TestGetContinuousCmap:
    def test_it_builds_a_colormap_from_a_colour_list(self, presentation):
        from matplotlib.colors import LinearSegmentedColormap

        cmap = presentation._get_continuous_cmap([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
        assert isinstance(cmap, LinearSegmentedColormap)

    def test_the_ends_of_the_map_are_the_ends_of_the_list(self, presentation):
        cmap = presentation._get_continuous_cmap([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
        assert cmap(0.0)[:3] == pytest.approx((0.0, 0.0, 0.0))
        assert cmap(1.0)[:3] == pytest.approx((1.0, 1.0, 1.0))

    def test_explicit_stops_place_the_colours(self, presentation):
        cmap = presentation._get_continuous_cmap(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 1.0]], [0.0, 0.25, 1.0],
        )
        assert cmap(0.25)[:3] == pytest.approx((1.0, 0.0, 0.0), abs=0.01)


@pytest.mark.unit
class TestPlotOriginIgnoresItsAxes:
    """B113: see the module docstring."""

    def test_it_returns_the_axes_it_was_given(self, presentation):
        _, ax = plt.subplots()
        assert presentation.plotOrigin(ax=ax) is ax

    def test_it_creates_axes_when_given_none(self, presentation):
        assert presentation.plotOrigin() is not None

    @pytest.mark.xfail(
        strict=True,
        reason="B113: plotOrigin draws with plt.scatter instead of "
               "ax.scatter, so the marker lands on the current axes rather "
               "than the one passed in. See the consolidated findings issue.",
    )
    def test_the_marker_should_land_on_the_axes_that_was_passed_in(self, presentation):
        _, target = plt.subplots()
        plt.subplots()  # a second figure becomes the current axes
        presentation.plotOrigin(ax=target)
        assert len(target.collections) == 1

    def test_the_marker_currently_lands_on_the_current_axes_instead(self, presentation):
        """Characterisation of B113."""
        _, target = plt.subplots()
        _, current = plt.subplots()
        presentation.plotOrigin(ax=target)
        assert len(target.collections) == 0
        assert len(current.collections) == 1


@pytest.mark.unit
class TestGenerateLatexTableIsDead:
    """B114: see the module docstring."""

    @pytest.mark.xfail(
        strict=True,
        reason="B114: generateLatexTable calls plotDevices with `device=` "
               "and `display=`, neither of which exists in its signature "
               "(the real parameters are deviceType and a mandatory "
               "mapName), so it always raises TypeError. See the "
               "consolidated findings issue.",
    )
    def test_it_should_write_a_tex_document(self, presentation, tmp_path):
        presentation.generateLatexTable("{{ trialSets }}", str(tmp_path / "out"))
        assert (tmp_path / "out" / "latex_document.tex").is_file()

    def test_it_currently_raises_typeerror_on_the_plotdevices_call(self, presentation, tmp_path):
        """Characterisation of B114."""
        with pytest.raises(TypeError, match="device"):
            presentation.generateLatexTable("{{ trialSets }}", str(tmp_path / "out"))

    def test_plotdevices_really_has_no_device_or_display_parameter(self):
        """So B114 is a call-signature mismatch, not a missing feature."""
        import inspect

        params = inspect.signature(experimentPresentation.plotDevices).parameters
        assert "device" not in params
        assert "display" not in params
        assert "deviceType" in params and "mapName" in params
