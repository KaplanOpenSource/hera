"""measurements/experiment/presentation.py: the plotting methods that need
a loaded experiment.

``test_experiment_presentation.py`` covers ``experimentPresentation``'s
construction, its properties and its pure helpers (and pins B113/B114)
against a hand-built stand-in datalayer, and deliberately leaves out the
seven methods that need a real experiment. This file covers exactly those
seven: ``plotImage``, ``plotMap``, ``_plotEntityLocationScatter``,
``_plotEntityLocationNames``, ``plotDevicesOnImage``, ``plotDevices`` and
``plotDeviceTypeFunctionality``.

Construction strategy
---------------------
No stubbing. Following ``test_experiment_experiment.py``, a genuine argos
**v3.0.0** setup dict is synthesised and written as
``<root>/runtimeExperimentData/Datasources_Configurations.json`` plus
``<root>/runtimeExperimentData/UNITEXP.zip`` under ``tmp_path``; the zip
also carries a real ``images/map1.png`` (a 4x4 array written with
``plt.imsave``) so that ``argos``'s ``ExperimentZipFile.getImage`` -- which
opens the member from inside the archive -- succeeds. The experiment is
then loaded through the real ``experimentSetupWithData`` constructor, and
the presentation layer under test is the one that constructor wires up
(``experiment.presentation``). The conftest forces the Agg backend and
closes every figure after each test, so the assertions below are made on
the resulting matplotlib artists (collections, texts, images, figure
size).

Three device layouts are used, because ``plotDevicesOnImage`` behaves
differently in each: devices positioned directly on the map
("uncontained"), a mix of one contained in a pole and one not, and every
device contained in a pole. Containment matters because argos only adds
the ``containedIn`` column to a trial's entities table when some device
on the trial declares it.

Covered
-------
Only one of the seven methods can run at all: ``plotDevicesOnImage``, and
only in the all-contained layout. Its working behaviour is covered in
full (return value, one scatter per device, marker defaults and
overrides, label source, text offsets, a caller-supplied axes,
``plotkwargs`` forwarding and the empty selection). The other six are dead
on their first or second statement and are pinned as B218-B221 and B226;
``plotDevicesOnImage``'s own three defects are pinned as B222-B225.

Deliberately not covered
------------------------
``plotOrigin``, ``generateLatexTable``, ``_splitName``,
``_scatter_height_color`` and ``_get_continuous_cmap`` -- all already
covered in ``test_experiment_presentation.py``.
``plotImage``'s ``withGrid`` / ``plt_kwargs`` / ``majorLocator``
arguments, and ``plotDevicesOnImage``'s never-plotted background image:
unreachable while B218 stands.

Bugs pinned here
----------------
* B218: ``plotImage`` builds the ``imshow`` extent from
  ``metadata['xleft']``, ``['xright']``, ``['ybottom']``, ``['ytop']``,
  but ``argos.Experiment.getImageMetadata`` returns the raw
  ``imageStandalone`` entry, whose bounding-box keys are ``left``,
  ``right``, ``lower`` and ``upper``. Every call raises
  ``KeyError: 'xleft'`` -- after the image itself has been decoded, so the
  method is dead for every experiment.
* B219: ``plotMap`` reads ``self.trialSet[trialSetName][trialName]``, but
  ``trialSet`` lives on the datalayer, not on the presentation layer
  (``self.datalayer.trialSet``), so the first statement raises
  ``AttributeError``. Three further defects sit behind it: ``ax`` and
  ``plot_kwargs`` are read although neither is a parameter of the method
  (its only parameters are ``trialSetName`` and ``trialName``),
  ``deviceType`` and ``toolkitDataSource`` are undefined names, and the
  method returns nothing at all.
* B220: ``_plotEntityLocationScatter`` and ``_plotEntityLocationNames``
  both start from ``self.datalayer.experimentSetup.trialSet[...]``, but
  ``experimentSetupWithData`` has no ``experimentSetup`` attribute -- it
  *is* the experiment setup and exposes ``.trialSet`` directly (the same
  mistake as B165 in ``dataEngine.py``), so both raise ``AttributeError``
  immediately. Behind that: ``entitiesTable(status)`` calls an argos
  property, ``self._entityMarkers`` is never defined on the class, and
  ``FLOOR_PLATFORM``/``FLOOR_CONCOURSE`` are undefined names in the
  module.
* B221: ``plotDevices`` opens with ``plot_kwargs = plot_kwargs or {}``,
  reading a local before it is assigned -- the parameter is spelled
  ``plotkwargs`` -- so with the default ``ax=None`` it raises
  ``UnboundLocalError``; passing an ``ax`` skips that line only to hit
  ``self.trialSet`` (B219's mistake again) and, behind it,
  ``self.datalayer._process_row`` (dead per B160) and ``row.stationName``,
  a column argos does not produce.
* B222: ``plotDevicesOnImage`` reads ``row.containedIn`` for every device,
  but argos's ``fillContained`` only creates that column when at least one
  device on the trial declares containment. For the ordinary experiment
  whose devices sit directly on the map the loop raises
  ``AttributeError: 'Pandas' object has no attribute 'containedIn'``.
* B223: ``plotDevicesOnImage`` computes its bounding box as
  ``devices_df.max().longitude`` / ``devices_df.min().latitude`` -- an
  aggregate over the *whole* mixed-dtype frame rather than over the two
  coordinate columns. As soon as any string column holds a NaN (which is
  exactly what a mix of contained and uncontained devices produces in
  ``containedIn``) the aggregate raises ``TypeError: '>=' not supported
  between instances of 'str' and 'float'``. Together with B222 this leaves
  a single layout -- every device of the type contained in another device
  -- in which the method runs at all.
* B224: ``plotDevicesOnImage`` scatters ``x=row.latitude, y=row.longitude``,
  transposing the two axes: the image these devices are meant to sit on is
  drawn by ``plotImage`` with the longitudes (``left``/``right``) on x and
  the latitudes (``lower``/``upper``) on y. The text offset compounds it,
  shifting a longitude value by ``(maxy - miny) * textDeltaY`` where
  ``maxy``/``miny`` are latitudes.
* B225: ``plotDevicesOnImage`` mutates its own mutable default argument --
  ``scatterkwargs.setdefault("s", 50)`` / ``setdefault("c", "r")`` on a
  ``scatterkwargs={}`` default -- so after the first call that omits the
  argument the shared default dict permanently carries ``s``/``c`` for the
  rest of the process, and a later caller can no longer get an unstyled
  scatter by omitting it.
* B226: ``plotDeviceTypeFunctionality`` cannot run, because it calls
  ``analysis.getDeviceTypeTransmissionFrequencyOfTrial`` with
  ``normalize=True`` and without ``recalculate=False``, so it always
  enters the recompute branch killed by B216 (and would then hit B215's
  missing ``getOptimalFrequencyHz``). Both are pinned in
  ``test_experiment_analysis.py``; the tests here record that the
  presentation entry point inherits their deadness.
"""
import inspect
import json
import os
import zipfile

import matplotlib.pyplot as plt
import pandas
import pytest

from hera.measurements.experiment.experiment import experimentSetupWithData
from hera.measurements.experiment.presentation import experimentPresentation
from hera.tests.unit.conftest import UNIT_PROJECT_NAME

TRIAL_START = "2020-01-01 00:00:00"
TRIAL_END = "2020-01-01 00:09:00"
TZ = "Asia/Jerusalem"

# The map the devices are placed on, in the units argos stores them in.
MAP_LEFT, MAP_RIGHT, MAP_LOWER, MAP_UPPER = 34.0, 35.0, 32.0, 33.0

# Device positions, longitude first (that is the order argos's
# ``location.coordinates`` uses).
SONIC_1 = (34.0, 32.0)
SONIC_2 = (34.5, 32.5)


# ---------------------------------------------------------------------------
# A genuine argos v3.0.0 experiment, laid out on disk with a real image
# ---------------------------------------------------------------------------

def _sonic(name, longitude=None, latitude=None, containedIn=None):
    """One entry of a v3.0.0 ``devicesOnTrial`` list.

    A device with ``containedIn`` carries no location of its own: argos
    makes it inherit its container's.
    """
    device = {
        "deviceTypeName": "sonic",
        "deviceItemName": name,
        "attributes": [{"name": "height", "value": "3"}],
    }
    if containedIn is None:
        device["location"] = {"name": "map1", "coordinates": [longitude, latitude]}
    else:
        device["containedIn"] = {"deviceTypeName": "pole", "deviceItemName": containedIn}
    return device


def _pole(name, longitude, latitude):
    return {
        "deviceTypeName": "pole",
        "deviceItemName": name,
        "location": {"name": "map1", "coordinates": [longitude, latitude]},
        "attributes": [],
    }


SONIC_TYPE = {
    "name": "sonic",
    "attributeTypes": [
        {"name": "height", "type": "Number", "label": "height",
         "description": "", "scope": "Device"},
        {"name": "StoreDataPerDevice", "type": "String",
         "label": "StoreDataPerDevice", "description": "",
         "scope": "Constant", "defaultValue": False},
    ],
    "devices": [
        {"name": "sonic 1", "attributes": [{"name": "height", "value": "3"}]},
        {"name": "sonic 2", "attributes": [{"name": "height", "value": "6"}]},
    ],
}


def _poleType(names):
    return {
        "name": "pole",
        "attributeTypes": [],
        "devices": [{"name": name, "attributes": []} for name in names],
    }


def _experimentSetup(layout="uncontained"):
    """A minimal but genuine argos v3.0.0 setup dict.

    layout
        ``uncontained``  -- both sonics sit directly on the map;
        ``mixed``        -- sonic 1 is contained in a pole, sonic 2 is not;
        ``allContained`` -- each sonic is contained in its own pole.
    """
    deviceTypes = [dict(SONIC_TYPE)]
    if layout == "uncontained":
        devices = [_sonic("sonic 1", *SONIC_1), _sonic("sonic 2", *SONIC_2)]
    elif layout == "mixed":
        devices = [
            _sonic("sonic 1", containedIn="pole 1"),
            _sonic("sonic 2", *SONIC_2),
            _pole("pole 1", *SONIC_1),
        ]
        deviceTypes.append(_poleType(["pole 1"]))
    elif layout == "allContained":
        devices = [
            _sonic("sonic 1", containedIn="pole 1"),
            _sonic("sonic 2", containedIn="pole 2"),
            _pole("pole 1", *SONIC_1),
            _pole("pole 2", *SONIC_2),
        ]
        deviceTypes.append(_poleType(["pole 1", "pole 2"]))
    else:
        raise AssertionError(f"unknown layout {layout}")

    return {
        "version": "3.0.0",
        "name": "UNITEXP",
        "description": "a synthetic experiment",
        "startDate": "2020-01-01T00:00:00.000Z",
        "endDate": "2020-01-02T00:00:00.000Z",
        "trialTypes": [{
            "name": "Measurements",
            "description": "the measurement trials",
            "attributeTypes": [
                {"key": "TrialStart", "name": "TrialStart", "type": "datetime-local",
                 "label": "TrialStart", "description": ""},
                {"key": "TrialEnd", "name": "TrialEnd", "type": "datetime-local",
                 "label": "TrialEnd", "description": ""},
            ],
            "trials": [{
                "name": "T1",
                "createdDate": "2020-01-01T00:00:00.000Z",
                "properties": [
                    {"key": "TrialStart", "val": TRIAL_START},
                    {"key": "TrialEnd", "val": TRIAL_END},
                ],
                "devicesOnTrial": devices,
            }],
        }],
        "deviceTypes": deviceTypes,
        "imageStandalone": [{
            "name": "map1", "filename": "map1.png",
            "left": MAP_LEFT, "right": MAP_RIGHT,
            "lower": MAP_LOWER, "upper": MAP_UPPER,
            "width": 100, "height": 100,
        }],
    }


@pytest.fixture()
def experimentFactory(tmp_path, unit_files_directory, unit_project):
    """Build a real ``experimentSetupWithData``, image included.

    ``unit_project`` is requested first on purpose: it persists the tmp
    files directory in the project config, so the ``Project`` that
    ``parquetDataEngineHera`` builds for itself (which forwards no
    ``filesDirectory``) does not fall back to ``~/.hera``.
    """
    import numpy

    imagePath = str(tmp_path / "map1.png")
    plt.imsave(imagePath, numpy.zeros((4, 4, 3)))

    counter = {"n": 0}

    def _build(layout="uncontained", defaultTrialSetName="Measurements"):
        counter["n"] += 1
        root = str(tmp_path / f"exp{counter['n']}")
        runtime = os.path.join(root, "runtimeExperimentData")
        os.makedirs(runtime, exist_ok=True)
        with open(os.path.join(runtime, "Datasources_Configurations.json"), "w") as handle:
            json.dump({"experimentName": "UNITEXP"}, handle)
        with zipfile.ZipFile(os.path.join(runtime, "UNITEXP.zip"), "w") as archive:
            archive.writestr("data.json", json.dumps(_experimentSetup(layout)))
            archive.write(imagePath, "images/map1.png")
        return experimentSetupWithData(
            projectName=UNIT_PROJECT_NAME,
            pathToExperiment=root,
            filesDirectory=unit_files_directory,
            defaultTrialSetName=defaultTrialSetName,
        )

    return _build


@pytest.fixture()
def experiment(experimentFactory):
    """Devices placed directly on the map -- the ordinary layout."""
    return experimentFactory("uncontained")


@pytest.fixture()
def presentation(experiment):
    """The presentation layer the experiment constructor wired up."""
    return experiment.presentation


@pytest.fixture()
def containedExperiment(experimentFactory):
    """Every sonic contained in its own pole -- the only layout in which
    ``plotDevicesOnImage`` runs (see B222/B223)."""
    return experimentFactory("allContained")


@pytest.fixture()
def containedPresentation(containedExperiment):
    return containedExperiment.presentation


@pytest.fixture()
def mixedPresentation(experimentFactory):
    """One contained device and one not."""
    return experimentFactory("mixed").presentation


def _deviceData(project, tmp_path):
    """Register the shared-file raw-data document for the sonic type."""
    index = pandas.date_range("2020-01-01 00:00:00", periods=20, freq="1min", tz=TZ)
    frame = pandas.DataFrame(
        {
            "deviceName": ["sonic 1" if i % 2 == 0 else "sonic 2" for i in range(20)],
            "TC": [20.0 + i for i in range(20)],
        },
        index=index,
    ).rename_axis("timestamp")
    path = tmp_path / "rawData_sonic.pkl"
    frame.to_pickle(str(path))
    project.addMeasurementsDocument(
        resource=str(path),
        dataFormat="pickle",
        type="Experiment_rawData",
        desc=dict(experimentName="UNITEXP", deviceType="sonic"),
    )


@pytest.fixture()
def deviceData(unit_project, tmp_path):
    _deviceData(unit_project, tmp_path)


# ---------------------------------------------------------------------------
# plotDevicesOnImage -- the one method that runs
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPlotDevicesOnImage:
    """The all-contained layout, which is the only one that survives."""

    def test_it_returns_the_figure_and_the_axes_it_drew_on(self, containedPresentation):
        fig, ax = containedPresentation.plotDevicesOnImage(
            "Measurements", "T1", "sonic", "map1", scatterkwargs=dict())
        assert ax.figure is fig

    def test_it_draws_one_scatter_per_device_of_the_type(self, containedPresentation):
        _, ax = containedPresentation.plotDevicesOnImage(
            "Measurements", "T1", "sonic", "map1", scatterkwargs=dict())
        assert len(ax.collections) == 2

    def test_the_marker_defaults_are_a_red_dot_of_size_fifty(self, containedPresentation):
        _, ax = containedPresentation.plotDevicesOnImage(
            "Measurements", "T1", "sonic", "map1", scatterkwargs=dict())
        from matplotlib.colors import to_rgba

        assert ax.collections[0].get_sizes().tolist() == [50]
        # Compare against whatever "r" currently resolves to rather than a
        # literal RGBA: importing hera.presentation.basicplots runs
        # seaborn.set() at import time, which rebinds matplotlib's
        # single-letter colour shorthands process-wide (B309), so a literal
        # [1,0,0,1] here passes alone and fails in a full run.
        assert ax.collections[0].get_facecolor()[0].tolist() == list(to_rgba("r"))

    def test_explicit_scatter_keywords_override_the_defaults(self, containedPresentation):
        _, ax = containedPresentation.plotDevicesOnImage(
            "Measurements", "T1", "sonic", "map1", scatterkwargs=dict(s=7, c="blue"))
        from matplotlib.colors import to_rgba

        assert ax.collections[0].get_sizes().tolist() == [7]
        assert ax.collections[0].get_facecolor()[0].tolist() == list(to_rgba("blue"))

    def test_each_device_is_labelled_by_the_item_that_contains_it(self, containedPresentation):
        _, ax = containedPresentation.plotDevicesOnImage(
            "Measurements", "T1", "sonic", "map1", scatterkwargs=dict())
        assert sorted(text.get_text() for text in ax.texts) == ["pole 1", "pole 2"]

    def test_a_zero_text_offset_puts_the_label_on_the_marker(self, containedPresentation):
        _, ax = containedPresentation.plotDevicesOnImage(
            "Measurements", "T1", "sonic", "map1", scatterkwargs=dict(),
            textDeltaX=0, textDeltaY=0)
        offsets = sorted(tuple(c.get_offsets().tolist()[0]) for c in ax.collections)
        positions = sorted(tuple(text.get_position()) for text in ax.texts)
        assert offsets == positions

    def test_the_text_offset_shifts_the_label_off_the_marker(self, containedPresentation):
        _, ax = containedPresentation.plotDevicesOnImage(
            "Measurements", "T1", "sonic", "map1", scatterkwargs=dict(),
            textDeltaX=2, textDeltaY=0)
        marker = sorted(c.get_offsets().tolist()[0] for c in ax.collections)[0]
        label = sorted(text.get_position() for text in ax.texts)[0]
        assert label[0] == pytest.approx(marker[0] + 2)
        assert label[1] == pytest.approx(marker[1])

    def test_a_caller_supplied_axes_is_drawn_on_and_its_figure_returned(
            self, containedPresentation):
        figure, axes = plt.subplots()
        fig, ax = containedPresentation.plotDevicesOnImage(
            "Measurements", "T1", "sonic", "map1", ax=axes, scatterkwargs=dict())
        assert ax is axes
        assert fig is figure
        assert len(axes.collections) == 2

    def test_plot_keywords_are_forwarded_to_the_new_figure(self, containedPresentation):
        fig, _ = containedPresentation.plotDevicesOnImage(
            "Measurements", "T1", "sonic", "map1",
            plotkwargs=dict(figsize=(3, 4)), scatterkwargs=dict())
        assert tuple(fig.get_size_inches()) == (3.0, 4.0)

    def test_an_unknown_device_type_draws_nothing(self, containedPresentation):
        _, ax = containedPresentation.plotDevicesOnImage(
            "Measurements", "T1", "noSuchType", "map1", scatterkwargs=dict())
        assert len(ax.collections) == 0
        assert len(ax.texts) == 0

    def test_an_unknown_map_name_draws_nothing(self, containedPresentation):
        _, ax = containedPresentation.plotDevicesOnImage(
            "Measurements", "T1", "sonic", "noSuchMap", scatterkwargs=dict())
        assert len(ax.collections) == 0

    def test_an_unknown_trial_name_raises_key_error(self, containedPresentation):
        with pytest.raises(KeyError):
            containedPresentation.plotDevicesOnImage(
                "Measurements", "noSuchTrial", "sonic", "map1", scatterkwargs=dict())


@pytest.mark.unit
class TestPlotDevicesOnImageNeedsContainment:
    """B222: an optional argos column is read unconditionally."""

    def test_an_uncontained_trial_has_no_containedin_column(self, experiment):
        """Characterisation of B222: the table the loop is given."""
        table = experiment.trialSet["Measurements"]["T1"].entitiesTable
        assert "deviceItemName" in table.columns
        assert "containedIn" not in table.columns

    @pytest.mark.xfail(
        strict=True,
        reason="B222: plotDevicesOnImage reads row.containedIn for every "
               "device, but argos's fillContained only adds that column "
               "when some device on the trial declares containment, so the "
               "ordinary map-positioned experiment raises AttributeError -- "
               "even though the method already has a fallback branch meant "
               "for exactly that case. See the consolidated findings issue.",
    )
    def test_it_should_label_uncontained_devices_by_their_own_name(self, presentation):
        _, ax = presentation.plotDevicesOnImage(
            "Measurements", "T1", "sonic", "map1", scatterkwargs=dict())
        assert sorted(text.get_text() for text in ax.texts) == ["sonic 1", "sonic 2"]

    def test_it_currently_raises_an_attribute_error(self, presentation):
        """Characterisation of B222."""
        with pytest.raises(AttributeError, match="containedIn"):
            presentation.plotDevicesOnImage(
                "Measurements", "T1", "sonic", "map1", scatterkwargs=dict())


@pytest.mark.unit
class TestPlotDevicesOnImageAggregatesTheWholeFrame:
    """B223: min/max over every column, not just the coordinates."""

    def test_a_mixed_trial_leaves_a_nan_in_a_string_column(self, experimentFactory):
        """Characterisation of B223: the frame the aggregate is given."""
        experiment = experimentFactory("mixed")
        table = experiment.trialSet["Measurements"]["T1"].entitiesTable
        sonics = table.query("deviceTypeName=='sonic'")
        assert sonics.containedIn.isna().sum() == 1
        with pytest.raises(TypeError, match="'>=' not supported"):
            sonics.max()

    def test_aggregating_only_the_coordinates_works_on_the_same_frame(self, experimentFactory):
        """Characterisation of B223: what the method actually needs."""
        experiment = experimentFactory("mixed")
        sonics = experiment.trialSet["Measurements"]["T1"].entitiesTable.query(
            "deviceTypeName=='sonic'")
        assert sonics.latitude.max() == pytest.approx(SONIC_2[1])
        assert sonics.longitude.min() == pytest.approx(SONIC_1[0])

    @pytest.mark.xfail(
        strict=True,
        reason="B223: plotDevicesOnImage takes its bounding box from "
               "devices_df.max()/min() over the whole mixed-dtype frame "
               "rather than the latitude/longitude columns, so a trial "
               "mixing contained and uncontained devices (a NaN in the "
               "containedIn string column) raises TypeError. "
               "See the consolidated findings issue.",
    )
    def test_it_should_plot_a_trial_mixing_contained_and_uncontained_devices(
            self, mixedPresentation):
        _, ax = mixedPresentation.plotDevicesOnImage(
            "Measurements", "T1", "sonic", "map1", scatterkwargs=dict())
        assert len(ax.collections) == 2

    def test_it_currently_raises_a_type_error(self, mixedPresentation):
        """Characterisation of B223."""
        with pytest.raises(TypeError, match="'>=' not supported"):
            mixedPresentation.plotDevicesOnImage(
                "Measurements", "T1", "sonic", "map1", scatterkwargs=dict())


@pytest.mark.unit
class TestPlotDevicesOnImageSwapsTheAxes:
    """B224: x gets the latitude and y gets the longitude."""

    def test_the_image_these_devices_sit_on_is_longitude_by_latitude(self, experiment):
        """Characterisation of B224: the map's own orientation.

        ``left``/``right`` bound the same values as the device longitudes
        and ``lower``/``upper`` the latitudes, so x is longitude.
        """
        metadata = experiment.getImageMetadata("map1")
        assert (metadata["left"], metadata["right"]) == (MAP_LEFT, MAP_RIGHT)
        assert (metadata["lower"], metadata["upper"]) == (MAP_LOWER, MAP_UPPER)
        table = experiment.trialSet["Measurements"]["T1"].entitiesTable
        assert metadata["left"] <= table.longitude.min()
        assert metadata["lower"] <= table.latitude.min()

    @pytest.mark.xfail(
        strict=True,
        reason="B224: plotDevicesOnImage scatters x=row.latitude, "
               "y=row.longitude, transposing the devices with respect to "
               "the image plotImage draws (x from left/right, i.e. "
               "longitudes; y from lower/upper, i.e. latitudes). The label "
               "offset compounds it, shifting a longitude by the latitude "
               "spread. See the consolidated findings issue.",
    )
    def test_the_markers_should_land_at_longitude_by_latitude(self, containedPresentation):
        _, ax = containedPresentation.plotDevicesOnImage(
            "Measurements", "T1", "sonic", "map1", scatterkwargs=dict())
        offsets = sorted(tuple(c.get_offsets().tolist()[0]) for c in ax.collections)
        assert offsets == [SONIC_1, SONIC_2]

    def test_the_markers_currently_land_transposed(self, containedPresentation):
        """Characterisation of B224."""
        _, ax = containedPresentation.plotDevicesOnImage(
            "Measurements", "T1", "sonic", "map1", scatterkwargs=dict())
        offsets = sorted(tuple(c.get_offsets().tolist()[0]) for c in ax.collections)
        assert offsets == [(SONIC_1[1], SONIC_1[0]), (SONIC_2[1], SONIC_2[0])]

    def test_the_label_offset_is_taken_from_the_latitude_spread(self, containedPresentation):
        """Characterisation of B224: the offset applied to the y (longitude)
        coordinate is a fraction of the latitude range."""
        _, ax = containedPresentation.plotDevicesOnImage(
            "Measurements", "T1", "sonic", "map1", scatterkwargs=dict(),
            textDeltaX=0, textDeltaY=1.0)
        latitudeSpread = SONIC_2[1] - SONIC_1[1]
        label = sorted(text.get_position() for text in ax.texts)[0]
        assert label[1] == pytest.approx(SONIC_1[0] + latitudeSpread)


@pytest.mark.unit
class TestPlotDevicesOnImageMutatesItsDefault:
    """B225: ``setdefault`` on a mutable default argument."""

    def test_the_scatter_keywords_default_to_a_shared_dictionary(self):
        """Characterisation of B225: the object that gets mutated."""
        parameter = inspect.signature(
            experimentPresentation.plotDevicesOnImage).parameters["scatterkwargs"]
        assert isinstance(parameter.default, dict)

    @pytest.mark.xfail(
        strict=True,
        reason="B225: plotDevicesOnImage calls setdefault on its "
               "scatterkwargs={} default, so a call that omits the argument "
               "permanently writes s=50 and c='r' into the shared default "
               "dict and no later caller can obtain an unstyled scatter by "
               "omitting it. See the consolidated findings issue.",
    )
    def test_omitting_the_scatter_keywords_should_leave_the_default_empty(
            self, containedPresentation):
        containedPresentation.plotDevicesOnImage("Measurements", "T1", "sonic", "map1")
        assert inspect.signature(
            experimentPresentation.plotDevicesOnImage).parameters["scatterkwargs"].default == {}

    def test_omitting_the_scatter_keywords_currently_pollutes_the_default(
            self, containedPresentation):
        """Characterisation of B225."""
        containedPresentation.plotDevicesOnImage("Measurements", "T1", "sonic", "map1")
        assert inspect.signature(
            experimentPresentation.plotDevicesOnImage
        ).parameters["scatterkwargs"].default == {"s": 50, "c": "r"}


# ---------------------------------------------------------------------------
# plotImage
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPlotImage:
    """B218: the extent is built from keys argos does not use."""

    def test_the_image_metadata_bounds_the_box_with_left_right_lower_upper(self, experiment):
        """Characterisation of B218: the dict the method indexes."""
        metadata = experiment.getImageMetadata("map1")
        assert sorted(k for k in metadata if k in
                      ("left", "right", "lower", "upper",
                       "xleft", "xright", "ybottom", "ytop")) == \
            ["left", "lower", "right", "upper"]

    def test_the_image_itself_loads_from_inside_the_experiment_zip(self, experiment):
        """So B218 is the only thing standing between plotImage and a plot."""
        image = experiment.getImage("map1")
        assert image.shape[:2] == (4, 4)

    @pytest.mark.xfail(
        strict=True,
        reason="B218: plotImage builds the imshow extent from "
               "metadata['xleft'/'xright'/'ybottom'/'ytop'], but argos's "
               "getImageMetadata returns the imageStandalone entry, whose "
               "bounding-box keys are left/right/lower/upper, so every call "
               "raises KeyError: 'xleft'. "
               "See the consolidated findings issue.",
    )
    def test_it_should_draw_the_image_over_its_geographic_extent(self, presentation):
        ax = presentation.plotImage("map1", withGrid=False)
        assert len(ax.images) == 1
        assert ax.images[0].get_extent() == (MAP_LEFT, MAP_RIGHT, MAP_LOWER, MAP_UPPER)

    def test_it_currently_raises_a_key_error(self, presentation):
        """Characterisation of B218."""
        with pytest.raises(KeyError, match="xleft"):
            presentation.plotImage("map1", withGrid=False)

    def test_a_caller_supplied_axes_reaches_the_same_failure(self, presentation):
        """Characterisation of B218: the ax branch is no better off."""
        _, ax = plt.subplots()
        with pytest.raises(KeyError, match="xleft"):
            presentation.plotImage("map1", ax=ax, withGrid=False)
        assert len(ax.images) == 0

    def test_an_unknown_image_name_raises_on_the_lookup_instead(self, presentation):
        with pytest.raises(KeyError, match="noSuchMap"):
            presentation.plotImage("noSuchMap")


# ---------------------------------------------------------------------------
# plotMap
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPlotMap:
    """B219: the trial sets are read off the wrong object."""

    def test_the_presentation_layer_has_no_trial_sets_of_its_own(self, presentation):
        """Characterisation of B219: the attribute the method reads."""
        assert not hasattr(presentation, "trialSet")
        assert hasattr(presentation.datalayer, "trialSet")

    def test_it_reads_names_that_are_not_its_parameters(self):
        """Characterisation of B219: ax, plot_kwargs, deviceType and
        toolkitDataSource are all read but never passed or defined."""
        parameters = list(inspect.signature(experimentPresentation.plotMap).parameters)
        assert parameters == ["self", "trialSetName", "trialName"]
        source = inspect.getsource(experimentPresentation.plotMap)
        for name in ("ax", "plot_kwargs", "deviceType", "toolkitDataSource"):
            assert name in source

    @pytest.mark.xfail(
        strict=True,
        reason="B219: plotMap reads self.trialSet, but trialSet lives on the "
               "datalayer, so the first statement raises AttributeError. "
               "Behind it: ax/plot_kwargs are read although neither is a "
               "parameter, deviceType and toolkitDataSource are undefined "
               "names, and the method returns nothing. "
               "See the consolidated findings issue.",
    )
    def test_it_should_plot_the_devices_of_the_trial_on_a_tile_map(self, presentation):
        result = presentation.plotMap("Measurements", "T1")
        assert result is not None

    def test_it_currently_raises_an_attribute_error(self, presentation):
        """Characterisation of B219."""
        with pytest.raises(AttributeError, match="trialSet"):
            presentation.plotMap("Measurements", "T1")


# ---------------------------------------------------------------------------
# _plotEntityLocationScatter / _plotEntityLocationNames
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPlotEntityLocationHelpers:
    """B220: both start from a datalayer attribute that does not exist."""

    def test_the_experiment_is_the_setup_and_has_no_experiment_setup_attribute(self, experiment):
        """Characterisation of B220: the attribute both methods read."""
        assert not hasattr(experiment, "experimentSetup")
        assert hasattr(experiment, "trialSet")

    def test_the_presentation_layer_defines_no_entity_markers(self):
        """Characterisation of B220: the next attribute the scatter needs."""
        assert not hasattr(experimentPresentation, "_entityMarkers")

    def test_the_module_defines_no_floor_constants(self):
        """Characterisation of B220: the names the label branch indexes."""
        import hera.measurements.experiment.presentation as presentationModule

        assert not hasattr(presentationModule, "FLOOR_PLATFORM")
        assert not hasattr(presentationModule, "FLOOR_CONCOURSE")

    @pytest.mark.xfail(
        strict=True,
        reason="B220: _plotEntityLocationScatter reads "
               "self.datalayer.experimentSetup.trialSet[...], but "
               "experimentSetupWithData has no experimentSetup attribute -- "
               "it exposes .trialSet directly (the same mistake as B165). "
               "Behind it: entitiesTable(status) calls an argos property, "
               "self._entityMarkers is never defined, and FLOOR_PLATFORM / "
               "FLOOR_CONCOURSE are undefined names. "
               "See the consolidated findings issue.",
    )
    def test_the_scatter_should_draw_the_entities_of_the_type(self, presentation):
        ax = presentation._plotEntityLocationScatter(
            "sonic", "Measurements", "T1", None, "map1")
        assert len(ax.collections) == 1

    def test_the_scatter_currently_raises_an_attribute_error(self, presentation):
        """Characterisation of B220."""
        with pytest.raises(AttributeError, match="experimentSetup"):
            presentation._plotEntityLocationScatter(
                "sonic", "Measurements", "T1", None, "map1")

    @pytest.mark.xfail(
        strict=True,
        reason="B220: _plotEntityLocationNames reads the same non-existent "
               "self.datalayer.experimentSetup, so it raises AttributeError "
               "before it can label anything. "
               "See the consolidated findings issue.",
    )
    def test_the_names_helper_should_label_the_entities(self, presentation):
        ax = presentation._plotEntityLocationNames(
            "sonic", "Measurements", "T1", None, "map1",
            plotNameMode=experimentPresentation.NAMES_LONG, text_kw=dict())
        assert len(ax.texts) == 2

    def test_the_names_helper_currently_raises_an_attribute_error(self, presentation):
        """Characterisation of B220."""
        with pytest.raises(AttributeError, match="experimentSetup"):
            presentation._plotEntityLocationNames(
                "sonic", "Measurements", "T1", None, "map1",
                plotNameMode=experimentPresentation.NAMES_LONG, text_kw=dict())

    def test_an_invalid_plot_name_mode_is_never_reached(self, presentation):
        """Characterisation of B220: the documented ValueError is dead too."""
        with pytest.raises(AttributeError, match="experimentSetup"):
            presentation._plotEntityLocationNames(
                "sonic", "Measurements", "T1", None, "map1",
                plotNameMode="neither short nor long", text_kw=dict())


# ---------------------------------------------------------------------------
# plotDevices
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPlotDevices:
    """B221: a misspelled local, then B219's mistake again."""

    def test_the_keyword_parameter_is_spelled_without_an_underscore(self):
        """Characterisation of B221: the name the body reads does not exist."""
        parameters = inspect.signature(experimentPresentation.plotDevices).parameters
        assert "plotkwargs" in parameters
        assert "plot_kwargs" not in parameters
        assert "plot_kwargs" in inspect.getsource(experimentPresentation.plotDevices)

    @pytest.mark.xfail(
        strict=True,
        reason="B221: plotDevices opens with plot_kwargs = plot_kwargs or {} "
               "while its parameter is spelled plotkwargs, so with the "
               "default ax=None it raises UnboundLocalError. Behind it: "
               "self.trialSet (B219's mistake), self.datalayer._process_row "
               "(dead per B160) and row.stationName, a column argos does not "
               "produce. See the consolidated findings issue.",
    )
    def test_it_should_plot_the_devices_of_the_trial(self, presentation):
        fig, ax = presentation.plotDevices("Measurements", "T1", "sonic", "map1")
        assert len(ax.collections) == 2

    def test_it_currently_raises_an_unbound_local_error(self, presentation):
        """Characterisation of B221."""
        with pytest.raises(UnboundLocalError, match="plot_kwargs"):
            presentation.plotDevices("Measurements", "T1", "sonic", "map1")

    def test_passing_an_axes_skips_that_line_and_fails_on_the_trial_sets(self, presentation):
        """Characterisation of B221: the second defect on the same path."""
        _, ax = plt.subplots()
        with pytest.raises(AttributeError, match="trialSet"):
            presentation.plotDevices("Measurements", "T1", "sonic", "map1", ax=ax)

    def test_the_plot_keywords_argument_cannot_be_used_either(self, presentation):
        """Characterisation of B221: naming the real parameter does not help,
        because the body still reads the misspelled local."""
        with pytest.raises(UnboundLocalError, match="plot_kwargs"):
            presentation.plotDevices("Measurements", "T1", "sonic", "map1",
                                     plotkwargs=dict(figsize=(3, 4)))


# ---------------------------------------------------------------------------
# plotDeviceTypeFunctionality
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPlotDeviceTypeFunctionality:
    """B226: dead because of the analysis defects B216 and B215."""

    def test_it_asks_the_analysis_layer_to_recompute_and_normalise(self):
        """Characterisation of B226: neither recalculate=False nor
        normalize=False is ever passed, so both dead branches are taken."""
        source = inspect.getsource(experimentPresentation.plotDeviceTypeFunctionality)
        assert "getDeviceTypeTransmissionFrequencyOfTrial" in source
        assert "normalize=True" in source
        assert "recalculate" not in source

    @pytest.mark.xfail(
        strict=True,
        reason="B226: plotDeviceTypeFunctionality calls "
               "getDeviceTypeTransmissionFrequencyOfTrial with normalize=True "
               "and without recalculate=False, so it always enters the "
               "recompute branch killed by B216 ('dict' object has no "
               "attribute 'query') and would then hit B215's missing "
               "getOptimalFrequencyHz. See the consolidated findings issue.",
    )
    def test_it_should_return_the_heatmap_axes_and_the_frequency_table(
            self, presentation, deviceData):
        ax, pvt = presentation.plotDeviceTypeFunctionality("sonic", "T1")
        assert sorted(pvt.columns) == ["sonic 1", "sonic 2"]
        assert ax.get_xlabel() == "Device name"

    def test_it_currently_raises_the_attribute_error_of_the_recompute_branch(
            self, presentation, deviceData):
        """Characterisation of B226."""
        with pytest.raises(AttributeError, match="'dict' object has no attribute 'query'"):
            presentation.plotDeviceTypeFunctionality("sonic", "T1")

    def test_without_raw_data_it_fails_even_earlier(self, presentation):
        """Characterisation of B226: the trial read raises first."""
        with pytest.raises(ValueError, match="There is no data for sonic"):
            presentation.plotDeviceTypeFunctionality("sonic", "T1")

    def test_no_default_trial_set_makes_the_lookup_fail_on_none(self, experimentFactory):
        """The trialSetName=None fallback is the datalayer's default, so with
        no default the trial-set lookup is given None."""
        presentation = experimentFactory(defaultTrialSetName=None).presentation
        with pytest.raises(KeyError):
            presentation.plotDeviceTypeFunctionality("sonic", "T1")
