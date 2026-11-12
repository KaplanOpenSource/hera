"""measurements/experiment/analysis.py: ``experimentAnalysis``'s
trial-facing methods.

``test_experiment_analysis_dataengine.py`` already covers
``getTurbulenceStatistics`` (B75) and the data-engine factory; this file
covers the six remaining methods: ``getDeviceLocations``, ``_splitName``,
``getDeviceTypeTransmissionFrequencyOfTrial``,
``getDeviceTypePlannedMessageCount``, ``addMetadata`` and
``addTrialProperties``.

Construction strategy
---------------------
Five of the six need a *loaded* experiment, so the fixtures reuse the
approach of ``test_experiment_experiment.py``: a genuine argos **v3.0.0**
setup dict is synthesised, written as
``<root>/runtimeExperimentData/Datasources_Configurations.json`` plus
``<root>/runtimeExperimentData/<name>.zip`` under ``tmp_path``, and loaded
through the real ``experimentSetupWithData`` constructor. Nothing is
stubbed -- the real ``argos`` package does the migration and the property
composition, and the analysis layer under test is the one the constructor
wires up (``experiment.analysis``). Device data is registered as ordinary
``type='Experiment_rawData'`` measurement documents holding pickled
frames.

``unit_project`` is requested by the factory fixture *before* the
experiment is built, because ``parquetDataEngineHera.__init__`` builds its
own ``datalayer.Project`` without forwarding a ``filesDirectory`` and
would otherwise create ``~/.hera/<projectName>``, which the unit layer's
``_no_home_writes`` guard fails.

Covered
-------
``_splitName`` (both branches, and its divergence from the identically
named helper in ``presentation.py``); ``getDeviceLocations`` (the
``trialSetName is None`` guard, and the two defects that kill the rest);
``getDeviceTypeTransmissionFrequencyOfTrial`` (the cache-hit path in wide
and long format, the default trial set, the cache key, and the two
defects that kill the compute and the normalisation branches);
``getDeviceTypePlannedMessageCount``; ``addMetadata`` (the missing-trial
path and the two defects); ``addTrialProperties`` (the working path, the
index-vs-column handling, the missing-``timestamp`` ``ValueError``, the
naive/aware mismatch, and the missing-``ReleaseStart`` defect).

Deliberately not covered
------------------------
``getTurbulenceStatistics`` -- already pinned as B75 in
``test_experiment_analysis_dataengine.py``.
The ``PANDASDB``/``DASKDB`` engine variants of the experiment: both open
raw ``pymongo`` connections to a real server, and their defects are
already pinned as B76/B77/B165.

Bugs pinned here
----------------
* B212 (fixed): ``getDeviceLocations`` and ``addMetadata`` both wrote
  ``...[trialName].entitiesTable()``, but ``argos.Trial.entitiesTable``
  is a ``property`` returning a ``DataFrame``. Every call therefore
  raised ``TypeError: 'DataFrame' object is not callable``, so both
  methods were dead. (``experiment.py`` gets this right -- it reads
  ``.entitiesTable`` without the parentheses.)
* B213 (fixed): ``getDeviceLocations`` filtered with
  ``query("entityType==@entityTypeName")``, but a *trial's* entities
  table has no ``entityType`` column -- argos names it
  ``deviceTypeName`` there, and ``entityType`` exists only on
  ``Experiment.entitiesTable`` / ``EntityType.entitiesTable``.
* B214 (fixed): ``addMetadata`` merged the trial metadata with
  ``right_on="entityName"``, a column that likewise does not exist on
  ``Trial.entitiesTable`` (argos spells it ``deviceItemName``). This is
  the same mismatch as B162 in ``experiment.py``, copy-pasted into the
  analysis layer.
* B215: ``getDeviceTypePlannedMessageCount`` calls
  ``self.getOptimalFrequencyHz(deviceType)``, a method that does not
  exist on ``experimentAnalysis`` or anywhere else in hera, so it always
  raises ``AttributeError``. Because
  ``getDeviceTypeTransmissionFrequencyOfTrial`` calls it whenever
  ``normalize=True`` -- its default -- the normalised frequency, which is
  what the docstring advertises, can never be produced.
* B216 (fixed): the recompute branch of
  ``getDeviceTypeTransmissionFrequencyOfTrial`` read
  ``experimentData.trialSet[trialSetName].trials.query("trialName == @trialName")``,
  but argos's ``TrialSet.trials`` is a **dict** of serialised trials, not
  a table (``TrialSet.trialsTable`` is the DataFrame, and it is indexed
  by trial name and has no ``trialName`` column either), so every call
  raised ``AttributeError: 'dict' object has no attribute 'query'``.
  Since ``recalculate`` defaults to ``True``, the method could only ever
  return anything at all on an explicit ``recalculate=False`` against an
  already-populated cache -- which nothing could populate. Two further
  defects sat on the same dead path: ``freq.set_index('timestamp')``
  discarded its result, and the ``groupby(...).resample("1min")``
  hard-coded a one-minute bin instead of using the ``samplingWindow``
  argument that the cache key and the normalisation both honour. The
  trial start/end now come from the trial's own properties.
* B217 (fixed): ``addTrialProperties`` reads
  ``properties.get("ReleaseStart", None)`` -- deliberately tolerating a
  trial with no release -- and then unconditionally computed
  ``tmp.timestamp - releaseTime``, which raised ``TypeError:
  unsupported operand type(s) for -: 'DatetimeArray' and 'NoneType'``.
  Every trial without a ``ReleaseStart`` property (the ordinary case for
  a measurement trial) therefore could not get its ``fromStart`` columns
  either, although those do not depend on the release at all. The
  release columns are now added only when there is a release.
"""
import inspect
import json
import os
import zipfile

import pandas
import pytest

from hera.datalayer import datatypes
from hera.measurements.experiment.experiment import experimentSetupWithData
from hera.tests.unit.conftest import UNIT_PROJECT_NAME

TRIAL_START = "2020-01-01 00:00:00"
TRIAL_END = "2020-01-01 00:09:00"
RELEASE_START = "2020-01-01 00:02:00"
TZ = "Asia/Jerusalem"


# ---------------------------------------------------------------------------
# A genuine argos v3.0.0 experiment, laid out on disk
# ---------------------------------------------------------------------------

def _deviceOnTrial(name, longitude, latitude):
    """One entry of a v3.0.0 ``devicesOnTrial`` list."""
    return {
        "deviceTypeName": "sonic",
        "deviceItemName": name,
        "location": {"name": "map1", "coordinates": [longitude, latitude]},
        "attributes": [{"name": "height", "value": "3"}],
    }


def _experimentSetup(withRelease=False):
    """A minimal but genuine argos v3.0.0 setup dict.

    withRelease
        also give the trial a ``ReleaseStart`` property, which is the only
        shape ``addTrialProperties`` can actually read (see B217).
    """
    trialAttributeTypes = [
        {"key": "TrialStart", "name": "TrialStart", "type": "datetime-local",
         "label": "TrialStart", "description": ""},
        {"key": "TrialEnd", "name": "TrialEnd", "type": "datetime-local",
         "label": "TrialEnd", "description": ""},
    ]
    trialProperties = [
        {"key": "TrialStart", "val": TRIAL_START},
        {"key": "TrialEnd", "val": TRIAL_END},
    ]
    if withRelease:
        trialAttributeTypes.append(
            {"key": "ReleaseStart", "name": "ReleaseStart", "type": "datetime-local",
             "label": "ReleaseStart", "description": ""})
        trialProperties.append({"key": "ReleaseStart", "val": RELEASE_START})

    return {
        "version": "3.0.0",
        "name": "UNITEXP",
        "description": "a synthetic experiment",
        "startDate": "2020-01-01T00:00:00.000Z",
        "endDate": "2020-01-02T00:00:00.000Z",
        "trialTypes": [{
            "name": "Measurements",
            "description": "the measurement trials",
            "attributeTypes": trialAttributeTypes,
            "trials": [{
                "name": "T1",
                "createdDate": "2020-01-01T00:00:00.000Z",
                "properties": trialProperties,
                "devicesOnTrial": [
                    _deviceOnTrial("sonic 1", 34.0, 32.0),
                    _deviceOnTrial("sonic 2", 34.5, 32.5),
                ],
            }],
        }],
        "deviceTypes": [{
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
        }],
        "imageStandalone": [{
            "name": "map1", "filename": "map1.png",
            "left": 34.0, "right": 35.0, "lower": 32.0, "upper": 33.0,
            "width": 100, "height": 100,
        }],
    }


@pytest.fixture()
def experimentFactory(tmp_path, unit_files_directory, unit_project):
    """Build a real ``experimentSetupWithData`` from a setup dict.

    ``unit_project`` is requested first on purpose: it persists the tmp
    files directory in the project config, so the ``Project`` that
    ``parquetDataEngineHera`` builds for itself (which forwards no
    ``filesDirectory``) does not fall back to ``~/.hera``.
    """
    counter = {"n": 0}

    def _build(setup=None, defaultTrialSetName="Measurements"):
        counter["n"] += 1
        root = str(tmp_path / f"exp{counter['n']}")
        runtime = os.path.join(root, "runtimeExperimentData")
        os.makedirs(runtime, exist_ok=True)
        with open(os.path.join(runtime, "Datasources_Configurations.json"), "w") as handle:
            json.dump({"experimentName": "UNITEXP"}, handle)
        with zipfile.ZipFile(os.path.join(runtime, "UNITEXP.zip"), "w") as archive:
            archive.writestr("data.json",
                             json.dumps(_experimentSetup() if setup is None else setup))
        return experimentSetupWithData(
            projectName=UNIT_PROJECT_NAME,
            pathToExperiment=root,
            filesDirectory=unit_files_directory,
            defaultTrialSetName=defaultTrialSetName,
        )

    return _build


@pytest.fixture()
def experiment(experimentFactory):
    return experimentFactory()


@pytest.fixture()
def analysis(experiment):
    """The analysis layer the experiment constructor wired up."""
    return experiment.analysis


@pytest.fixture()
def releaseExperiment(experimentFactory):
    """An experiment whose trial carries a ``ReleaseStart`` property."""
    return experimentFactory(_experimentSetup(withRelease=True))


def _deviceFrame():
    """Twenty minutes of alternating 'sonic 1'/'sonic 2' rows, tz-aware."""
    index = pandas.date_range("2020-01-01 00:00:00", periods=20, freq="1min", tz=TZ)
    return pandas.DataFrame(
        {
            "deviceName": ["sonic 1" if i % 2 == 0 else "sonic 2" for i in range(20)],
            "TC": [20.0 + i for i in range(20)],
        },
        index=index,
    ).rename_axis("timestamp")


@pytest.fixture()
def deviceData(unit_project, tmp_path):
    """The shared-file raw-data document for the sonic device type."""
    frame = _deviceFrame()
    path = tmp_path / "rawData_sonic.pkl"
    frame.to_pickle(str(path))
    unit_project.addMeasurementsDocument(
        resource=str(path),
        dataFormat="pickle",
        type="Experiment_rawData",
        desc=dict(experimentName="UNITEXP", deviceType="sonic"),
    )
    return frame


def _cacheKey(**overrides):
    """The description ``getDeviceTypeTransmissionFrequencyOfTrial`` looks up."""
    key = dict(deviceType="sonic", samplingWindow="1min", trialName="T1",
               trialSetName="Measurements", completeTimeSeries=True,
               completeDevices=True)
    key.update(overrides)
    return key


CACHED_FREQUENCY = pandas.DataFrame(
    {"sonic 1": [5.0, 5.0], "sonic 2": [4.0, 3.0]},
    index=pandas.DatetimeIndex(["2020-01-01 00:00:00", "2020-01-01 00:01:00"]),
)


@pytest.fixture()
def cachedFrequency(experiment):
    """Seed the frequency cache document the method reads on a cache hit."""

    def _seed(**overrides):
        experiment.addCacheDocument(
            type=experiment.analysis.TECHNICALDOC_FREQUENCY,
            dataFormat=datatypes.JSON_PANDAS,
            resource=CACHED_FREQUENCY.to_json(),
            desc=_cacheKey(**overrides),
        )

    return _seed


# ---------------------------------------------------------------------------
# _splitName
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestSplitName:
    def test_a_two_token_name_yields_the_second_token(self, analysis):
        assert analysis._splitName("sonic 12") == "12"

    def test_a_single_token_name_is_returned_whole(self, analysis):
        assert analysis._splitName("sonic") == "sonic"

    def test_only_the_second_token_is_taken_from_a_longer_name(self, analysis):
        assert analysis._splitName("sonic 3 north") == "3"

    def test_an_empty_name_yields_the_empty_string(self, analysis):
        assert analysis._splitName("") == ""

    def test_it_returns_a_string_where_the_presentation_helper_returns_an_int(self, analysis):
        """The two identically named helpers do not agree on the return type.

        ``presentation._splitName`` wraps the same token in ``int()``, so
        the analysis layer sorts device ids lexicographically while the
        presentation layer sorts them numerically.
        """
        from hera.measurements.experiment.presentation import experimentPresentation

        assert analysis._splitName("sonic 12") == "12"
        assert experimentPresentation._splitName(None, "sonic 12") == 12


# ---------------------------------------------------------------------------
# getDeviceLocations
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetDeviceLocationsTrialSetResolution:
    def test_an_omitted_trial_set_name_falls_back_to_the_default(self, analysis):
        """The default trial set is used, so the lookup succeeds."""
        locations = analysis.getDeviceLocations("sonic", "T1")
        assert sorted(locations.deviceItemName) == ["sonic 1", "sonic 2"]

    def test_no_trial_set_name_and_no_default_raises_value_error(self, experimentFactory):
        analysis = experimentFactory(defaultTrialSetName=None).analysis
        with pytest.raises(ValueError, match="trialSetName is None"):
            analysis.getDeviceLocations("sonic", "T1")

    def test_an_unknown_trial_set_name_raises_key_error(self, analysis):
        with pytest.raises(KeyError):
            analysis.getDeviceLocations("sonic", "T1", trialSetName="noSuchTrialSet")


@pytest.mark.unit
class TestGetDeviceLocationsReadsTheProperty:
    """B212: ``entitiesTable`` is a property, not a method."""

    def test_the_trial_entities_table_is_a_property_not_a_method(self, experiment):
        """Why the method reads it without parentheses."""
        trial = experiment.trialSet["Measurements"]["T1"]
        assert isinstance(inspect.getattr_static(trial, "entitiesTable"), property)
        assert isinstance(trial.entitiesTable, pandas.DataFrame)

    def test_it_returns_the_rows_of_the_requested_device_type(self, analysis):
        locations = analysis.getDeviceLocations("sonic", "T1", trialSetName="Measurements")
        assert sorted(locations.deviceItemName) == ["sonic 1", "sonic 2"]

    def test_an_unknown_device_type_yields_an_empty_table(self, analysis):
        assert len(analysis.getDeviceLocations("noSuchType", "T1",
                                               trialSetName="Measurements")) == 0


@pytest.mark.unit
class TestGetDeviceLocationsFiltersOnTheTrialColumn:
    """B213: the trial's table spells the device type ``deviceTypeName``."""

    def test_the_trial_entities_table_names_the_type_column_devicetypename(self, experiment):
        """Why the query filters on deviceTypeName."""
        table = experiment.trialSet["Measurements"]["T1"].entitiesTable
        assert "deviceTypeName" in table.columns
        assert "entityType" not in table.columns

    def test_the_query_the_method_builds_selects_the_device_type(self, experiment):
        table = experiment.trialSet["Measurements"]["T1"].entitiesTable
        entityTypeName = "sonic"
        assert len(table.query("deviceTypeName==@entityTypeName")) == 2


# ---------------------------------------------------------------------------
# addMetadata
# ---------------------------------------------------------------------------

def _dataToMerge():
    return pandas.DataFrame(
        {"deviceName": ["sonic 1", "sonic 2"], "TC": [20.0, 21.0]},
        index=pandas.DatetimeIndex(
            ["2020-01-01 00:00:00", "2020-01-01 00:01:00"], name="timestamp"),
    )


@pytest.mark.unit
class TestAddMetadataTrialSetResolution:
    def test_an_omitted_trial_set_name_falls_back_to_the_default(self, analysis):
        merged = analysis.addMetadata(_dataToMerge(), "T1")
        assert "height" in merged.columns

    def test_no_default_trial_set_makes_the_lookup_fail_on_none(self, experimentFactory):
        """Unlike getDeviceLocations, addMetadata has no None guard, so the
        None goes straight into the trial-set lookup."""
        analysis = experimentFactory(defaultTrialSetName=None).analysis
        with pytest.raises(KeyError):
            analysis.addMetadata(_dataToMerge(), "T1")

    def test_an_unknown_trial_name_raises_key_error(self, analysis):
        with pytest.raises(KeyError):
            analysis.addMetadata(_dataToMerge(), "noSuchTrial", trialSetName="Measurements")


@pytest.mark.unit
class TestAddMetadataReadsTheProperty:
    """B212 again, at its second call site."""

    def test_it_joins_the_trial_device_metadata_onto_the_data(self, analysis):
        merged = analysis.addMetadata(_dataToMerge(), "T1", trialSetName="Measurements")
        assert "height" in merged.columns
        assert len(merged) == 2
        assert merged.index.name == "timestamp"


@pytest.mark.unit
class TestAddMetadataMergesOnTheTrialColumn:
    """B214: the trial's table spells the device name ``deviceItemName``."""

    def test_the_trial_entities_table_names_the_device_column_deviceitemname(self, experiment):
        """Why the merge joins on deviceItemName."""
        table = experiment.trialSet["Measurements"]["T1"].entitiesTable
        assert "deviceItemName" in table.columns
        assert "entityName" not in table.columns

    def test_the_merge_the_method_builds_joins_on_the_device_name(self, experiment):
        table = experiment.trialSet["Measurements"]["T1"].entitiesTable
        merged = _dataToMerge().reset_index().merge(
            table, left_on="deviceName", right_on="deviceItemName")
        assert len(merged) == 2


# ---------------------------------------------------------------------------
# getDeviceTypePlannedMessageCount
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetDeviceTypePlannedMessageCount:
    """B215: the helper it depends on does not exist."""

    def test_the_analysis_layer_has_no_optimal_frequency_helper(self, analysis):
        """Characterisation of B215: the missing attribute."""
        from hera.measurements.experiment.analysis import experimentAnalysis

        assert not hasattr(experimentAnalysis, "getOptimalFrequencyHz")
        assert not hasattr(analysis, "getOptimalFrequencyHz")

    @pytest.mark.xfail(
        strict=True,
        reason="B215: getDeviceTypePlannedMessageCount calls "
               "self.getOptimalFrequencyHz(deviceType), a method that does "
               "not exist on experimentAnalysis or anywhere else in hera, so "
               "every call raises AttributeError -- and with it the "
               "normalize=True default of "
               "getDeviceTypeTransmissionFrequencyOfTrial. "
               "See the consolidated findings issue.",
    )
    def test_it_should_return_the_number_of_messages_expected_in_the_window(self, analysis):
        assert analysis.getDeviceTypePlannedMessageCount("sonic", samplingWindow="1min") > 0

    def test_it_currently_raises_an_attribute_error(self, analysis):
        """Characterisation of B215."""
        with pytest.raises(AttributeError, match="getOptimalFrequencyHz"):
            analysis.getDeviceTypePlannedMessageCount("sonic", samplingWindow="1min")

    def test_the_sampling_window_is_parsed_before_the_failure(self, analysis):
        """An unparseable window fails earlier, in pandas.to_timedelta."""
        with pytest.raises(ValueError):
            analysis.getDeviceTypePlannedMessageCount("sonic", samplingWindow="not a window")


# ---------------------------------------------------------------------------
# getDeviceTypeTransmissionFrequencyOfTrial
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestTransmissionFrequencyFromTheCache:
    """The one path that survives: an explicit ``recalculate=False`` hit."""

    def test_a_cache_hit_returns_the_stored_wide_table(self, analysis, cachedFrequency):
        cachedFrequency()
        pvt = analysis.getDeviceTypeTransmissionFrequencyOfTrial(
            "sonic", "T1", trialSetName="Measurements",
            normalize=False, recalculate=False)
        assert list(pvt.columns) == ["sonic 1", "sonic 2"]
        assert pvt.loc[pandas.Timestamp("2020-01-01 00:01:00"), "sonic 2"] == 3

    def test_the_long_format_melts_one_row_per_device_and_time(self, analysis, cachedFrequency):
        cachedFrequency()
        long = analysis.getDeviceTypeTransmissionFrequencyOfTrial(
            "sonic", "T1", trialSetName="Measurements",
            normalize=False, recalculate=False, wideFormat=False)
        assert list(long.columns) == ["timestamp", "deviceName", "Frequency"]
        assert len(long) == 4
        assert sorted(long.deviceName.unique()) == ["sonic 1", "sonic 2"]
        row = long.query("deviceName=='sonic 2'").iloc[-1]
        assert row.Frequency == 3

    def test_an_omitted_trial_set_name_looks_the_cache_up_under_the_default(
            self, analysis, cachedFrequency):
        cachedFrequency()
        pvt = analysis.getDeviceTypeTransmissionFrequencyOfTrial(
            "sonic", "T1", normalize=False, recalculate=False)
        assert len(pvt) == 2

    def test_a_cache_document_for_another_sampling_window_is_not_a_hit(
            self, analysis, cachedFrequency, deviceData):
        """The window is part of the cache key, so the lookup misses and the
        table is recomputed from the raw data instead."""
        cachedFrequency(samplingWindow="5min")
        pvt = analysis.getDeviceTypeTransmissionFrequencyOfTrial(
            "sonic", "T1", trialSetName="Measurements",
            normalize=False, recalculate=False)
        assert len(pvt) == 10


@pytest.mark.unit
class TestTransmissionFrequencyRecompute:
    """B216: the trial start/end come from the trial, not from a query."""

    def test_the_trial_set_exposes_its_trials_as_a_dict(self, experiment):
        """Why the method cannot .query() TrialSet.trials."""
        trialSet = experiment.trialSet["Measurements"]
        assert isinstance(trialSet.trials, dict)
        assert not hasattr(trialSet.trials, "query")

    def test_the_trials_table_has_no_trialname_column_either(self, experiment):
        """The obvious substitute is also wrong -- trialsTable is indexed
        by the trial name."""
        table = experiment.trialSet["Measurements"].trialsTable
        assert "trialName" not in table.columns
        assert list(table.index) == ["T1"]

    def test_it_computes_the_frequency_table_from_the_raw_data(self, analysis, deviceData):
        pvt = analysis.getDeviceTypeTransmissionFrequencyOfTrial(
            "sonic", "T1", trialSetName="Measurements", normalize=False)
        assert sorted(pvt.columns) == ["sonic 1", "sonic 2"]
        # one row per minute of the trial window, TrialStart..TrialEnd
        assert len(pvt) == 10

    def test_the_recompute_populates_the_cache(self, analysis, deviceData):
        """B216: nothing could write the cache before, so recalculate=False
        had nothing to read."""
        analysis.getDeviceTypeTransmissionFrequencyOfTrial(
            "sonic", "T1", trialSetName="Measurements", normalize=False)
        cached = analysis.datalayer.getCacheDocuments(
            type=analysis.TECHNICALDOC_FREQUENCY, **_cacheKey())
        assert len(cached) == 1
        assert sorted(cached[0].getData().columns) == ["sonic 1", "sonic 2"]

    def test_the_sampling_window_sets_the_bin_width(self, analysis, deviceData):
        """B216: the bin used to be hard-coded to one minute."""
        pvt = analysis.getDeviceTypeTransmissionFrequencyOfTrial(
            "sonic", "T1", trialSetName="Measurements",
            samplingWindow="5min", normalize=False)
        assert len(pvt) == 2

    def test_recalculate_is_true_by_default_so_a_cache_hit_is_refreshed(
            self, analysis, cachedFrequency, deviceData):
        """The default recomputes and overwrites the stored table."""
        cachedFrequency()
        pvt = analysis.getDeviceTypeTransmissionFrequencyOfTrial(
            "sonic", "T1", trialSetName="Measurements", normalize=False)
        assert len(pvt) == 10
        assert len(pvt) != len(CACHED_FREQUENCY)

    def test_the_recompute_branch_fails_earlier_when_there_is_no_raw_data(self, analysis):
        """Without documents the trial read raises before B216 is reached."""
        with pytest.raises(ValueError, match="There is no data for sonic"):
            analysis.getDeviceTypeTransmissionFrequencyOfTrial(
                "sonic", "T1", trialSetName="Measurements", normalize=False)


@pytest.mark.unit
class TestTransmissionFrequencyNormalisationIsDead:
    """B215 reached through the frequency method's own default."""

    @pytest.mark.xfail(
        strict=True,
        reason="B215: normalize=True divides by "
               "getDeviceTypePlannedMessageCount, which calls the "
               "non-existent self.getOptimalFrequencyHz, so the normalised "
               "frequency the docstring advertises can never be produced. "
               "See the consolidated findings issue.",
    )
    def test_it_should_normalise_the_cached_table_to_the_planned_rate(
            self, analysis, cachedFrequency):
        cachedFrequency()
        pvt = analysis.getDeviceTypeTransmissionFrequencyOfTrial(
            "sonic", "T1", trialSetName="Measurements", recalculate=False)
        assert pvt.max().max() <= 1

    def test_it_currently_raises_an_attribute_error(self, analysis, cachedFrequency):
        """Characterisation of B215 through the frequency method."""
        cachedFrequency()
        with pytest.raises(AttributeError, match="getOptimalFrequencyHz"):
            analysis.getDeviceTypeTransmissionFrequencyOfTrial(
                "sonic", "T1", trialSetName="Measurements", recalculate=False)

    def test_normalize_defaults_to_true(self):
        """So the failure above is what an ordinary caller gets."""
        from hera.measurements.experiment.analysis import experimentAnalysis

        parameters = inspect.signature(
            experimentAnalysis.getDeviceTypeTransmissionFrequencyOfTrial).parameters
        assert parameters["normalize"].default is True
        assert parameters["recalculate"].default is True


# ---------------------------------------------------------------------------
# addTrialProperties
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestAddTrialProperties:
    def test_it_adds_the_elapsed_columns_relative_to_start_and_release(self, releaseExperiment):
        data = _deviceFrame()
        result = releaseExperiment.analysis.addTrialProperties(data, "T1")
        assert set(["fromStart", "fromRelease", "fromStartSeconds",
                    "fromReleaseSeconds"]).issubset(result.columns)

    def test_the_seconds_are_measured_from_the_trial_start(self, releaseExperiment):
        result = releaseExperiment.analysis.addTrialProperties(_deviceFrame(), "T1")
        # the frame starts exactly at TrialStart and steps one minute
        assert result.fromStartSeconds.iloc[0] == 0
        assert result.fromStartSeconds.iloc[3] == 180

    def test_the_seconds_before_the_release_are_negative(self, releaseExperiment):
        result = releaseExperiment.analysis.addTrialProperties(_deviceFrame(), "T1")
        # ReleaseStart is two minutes after the first row
        assert result.fromReleaseSeconds.iloc[0] == -120
        assert result.fromReleaseSeconds.iloc[2] == 0

    def test_the_elapsed_columns_are_timedeltas(self, releaseExperiment):
        result = releaseExperiment.analysis.addTrialProperties(_deviceFrame(), "T1")
        assert result.fromStart.iloc[1] == pandas.Timedelta("1min")

    def test_a_timestamp_index_is_promoted_to_a_column(self, releaseExperiment):
        result = releaseExperiment.analysis.addTrialProperties(_deviceFrame(), "T1")
        assert "timestamp" in result.columns
        assert len(result) == 20

    def test_a_timestamp_column_is_used_as_it_stands(self, releaseExperiment):
        data = _deviceFrame().reset_index()
        result = releaseExperiment.analysis.addTrialProperties(data, "T1")
        assert result.fromStartSeconds.iloc[0] == 0

    def test_the_original_columns_are_preserved(self, releaseExperiment):
        result = releaseExperiment.analysis.addTrialProperties(_deviceFrame(), "T1")
        assert list(result.deviceName.head(2)) == ["sonic 1", "sonic 2"]
        assert result.TC.iloc[0] == 20.0

    def test_data_without_a_timestamp_anywhere_raises(self, releaseExperiment):
        data = _deviceFrame().reset_index(drop=True)
        with pytest.raises(ValueError, match="'timestamp' is not part of data"):
            releaseExperiment.analysis.addTrialProperties(data, "T1")

    def test_a_timezone_naive_frame_cannot_be_compared_to_the_aware_trial_times(
            self, releaseExperiment):
        naive = _deviceFrame().tz_localize(None, level=0)
        with pytest.raises(TypeError, match="tz-naive and tz-aware"):
            releaseExperiment.analysis.addTrialProperties(naive, "T1")

    def test_an_omitted_trial_set_name_falls_back_to_the_default(self, releaseExperiment):
        result = releaseExperiment.analysis.addTrialProperties(
            _deviceFrame(), "T1", trialSetName=None)
        assert result.fromStartSeconds.iloc[0] == 0

    def test_an_unknown_trial_name_raises_key_error(self, releaseExperiment):
        with pytest.raises(KeyError):
            releaseExperiment.analysis.addTrialProperties(_deviceFrame(), "noSuchTrial")


@pytest.mark.unit
class TestAddTrialPropertiesWithoutARelease:
    """B217: the optional release time is only subtracted when it exists."""

    def test_a_trial_without_a_release_has_no_releasestart_property(self, experiment):
        """The value the method fetches is None."""
        properties = experiment.trialSet["Measurements"]["T1"].properties
        assert "TrialStart" in properties
        assert properties.get("ReleaseStart", None) is None

    def test_it_still_adds_the_columns_measured_from_the_start(self, analysis):
        result = analysis.addTrialProperties(_deviceFrame(), "T1")
        assert result.fromStartSeconds.iloc[0] == 0
        assert result.fromStartSeconds.iloc[3] == 180

    def test_the_release_columns_are_omitted(self, analysis):
        """There is no release to measure from, so no release columns."""
        result = analysis.addTrialProperties(_deviceFrame(), "T1")
        assert "fromRelease" not in result.columns
        assert "fromReleaseSeconds" not in result.columns
