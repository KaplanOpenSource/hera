"""
Experiment Toolkit (Measurements)
=================================

This module implements the experiment toolkit layer used under
``hera.measurements``. It provides two main capabilities:

1) **experimentHome** (factory/home)
   A "home" object that lives under a project and acts as a factory for loading
   specific experiment toolkits dynamically from experiment folders registered
   as datasources in the project.

2) **experimentSetupWithData** (experiment toolkit)
   A toolkit that unifies an Argos experiment setup (metadata + trial structure)
   with a Hera data engine for efficient data access (parquet / pandas / dask).

Key Concepts
------------
- **Project**: The Hera project context (MongoDB + configuration + datasources).
- **Experiment datasource**: A datasource entry in the project that points to a
  folder on disk with an experiment structure (runtimeExperimentData, code, etc).
- **Dynamic loading**: Each experiment typically contains a Python class with the
  same name as the experiment, inside ``<experimentPath>/code``. We load it with
  ``pydoc.locate()`` after adding that folder to ``sys.path``.
- **Argos dependency**: Some classes inherit from Argos types. If Argos is not
  installed, runtime usage may be limited, and Sphinx should mock-import it.

Notes for Documentation (Sphinx)
--------------------------------
- If you build API docs with autodoc and Argos is not installed, configure:
  ``autodoc_mock_imports = ["argos"]`` in the relevant ``conf.py``.

"""

import os
import pydoc
import sys
import logging
import pandas as pd

from hera import toolkit
from .presentation import experimentPresentation
from .analysis import experimentAnalysis
from hera.measurements.GIS.utils import WSG84, ITM, convertCRS

try:
    from argos.experimentSetup import dataObjects as argosDataObjects
except ImportError:
    # Argos is optional; if it is not installed, experiment toolkits cannot be used.
    # For Sphinx docs, prefer mocking this import via autodoc_mock_imports.
    print("Must have argos installed and in the path. ")

from .dataEngine import dataEngineFactory, PARQUETHERA, PANDASDB, DASKDB
from hera.utils import loadJSON


# ---------------------------------------------------------------------------
# Argos trial time field names.
# These names must match the properties in the Argos web interface.
# Do not change unless the Argos schema changes.
# ---------------------------------------------------------------------------

TRIALSTART = 'TrialStart'
TRIALEND = 'TrialEnd'


class experimentHome(toolkit.abstractToolkit):
    """
    experimentHome
    --------------

    Factory / home object for experiment toolkits within a project.

    This class does **not** represent a specific experiment. Instead, it manages
    the *list of experiments registered as project datasources* and provides
    ``getExperiment()`` to dynamically load the appropriate toolkit class for a
    requested experiment.

    Typical usage
    -------------
    >>> from hera.measurements.experiment.experiment import experimentHome
    >>> eh = experimentHome(projectName="my_project")
    >>> exp = eh.getExperiment("MyExperiment")
    >>> df = exp.getDataFromDateRange(deviceType="IMS", startTime=..., endTime=...)

    Notes
    -----
    - Each experiment datasource is expected to expose a filesystem path.
      The toolkit class is expected at: ``<experimentPath>/code/<ExperimentName>.py``
      containing a class named ``<ExperimentName>``.
    - This home object uses the Hera datasource APIs inherited from
      ``abstractToolkit``: ``getDataSourceMap()``, ``getDataSourceTable()``,
      ``getDataSourceDocument()``.
    """

    DOCTYPE_ENTITIES = 'EntitiesData'
    CODE_DIRECTORY = 'code'

    def __init__(self, projectName, filesDirectory=None):
        """
        Initialize experimentHome.

        Parameters
        ----------
        projectName : str
            The project name to work with.
        filesDirectory : str, optional
            Directory to store cache / intermediate files for the toolkit.
            If None, uses the default behavior of ``abstractToolkit``.
        """
        super().__init__(projectName=projectName, toolkitName="experimentToolKit", filesDirectory=filesDirectory)
        self.logger = logging.getLogger()
        self.logger.info("Init experiment toolkit")

    @property
    def experimentMap(self):
        """
        Alias property for experiment map.

        Returns
        -------
        dict
            Mapping of experimentName -> datasource metadata.

        Notes
        -----
        The implementation in this file returns ``self.experimentMap()``, which is
        recursive. Prefer calling ``getExperimentsMap()`` directly.
        """
        return self.getExperimentsMap()

    def getExperimentsMap(self):
        """
        Get a dictionary mapping experiment names to datasource entries.

        Returns
        -------
        dict
            A dict where keys are experiment datasource names and values are
            datasource metadata dicts (as returned by ``getDataSourceMap()``).
        """
        M = dict()
        for experiment in self.getDataSourceMap():
            experimentName = experiment['datasourceName']
            M[experimentName] = experiment
        return M

    @property
    def experimentsTable(self):
        """
        Table view of experiment datasources.

        Returns
        -------
        pandas.DataFrame
            Datasource table for experiments (as returned by ``getDataSourceTable()``).
        """
        return self.getDataSourceTable()

    def getExperimentsTable(self):
        """
        Get the experiment datasource table.

        Returns
        -------
        pandas.DataFrame
            Datasource table for experiments.
        """
        return self.getDataSourceTable()

    def getExperiment(self, experimentName, filesDirectory=None):
        """
        Load and instantiate a specific experiment toolkit class by name.

        The loader:
        1) looks up the experiment datasource document in the project
        2) resolves its path on disk
        3) appends ``<experimentPath>/code`` to ``sys.path``
        4) loads the class ``<experimentName>.<experimentName>`` via ``pydoc.locate``
        5) instantiates and returns it

        Parameters
        ----------
        experimentName : str
            The name of the experiment datasource/toolkit to load.
        filesDirectory : str, optional
            Directory to store cache / intermediate files for the experiment toolkit.
            If None, experiment toolkit decides its default cache folder.

        Returns
        -------
        experimentSetupWithData
            An instantiated experiment toolkit (or a subclass), loaded dynamically.

        Raises
        ------
        ValueError
            If the datasource exists but the toolkit class cannot be located, or
            if the experiment is not registered in the project.
        """
        self.logger.info(f"Getting experiment {experimentName}")

        L = self.getDataSourceDocument(datasourceName=experimentName)
        if L:
            self.logger.info("Found experiment. Loading")
            experimentPath = L.getData()

            sys.path.append(os.path.join(experimentPath, self.CODE_DIRECTORY))
            self.logger.debug(f"Adding path {os.path.join(experimentPath, self.CODE_DIRECTORY)} to classpath")

            toolkitName = f"{experimentName}.{experimentName}"
            self.logger.debug(f"Loading toolkits: {toolkitName}")

            toolkitCls = pydoc.locate(toolkitName)
            if toolkitCls is None:
                err = f"Cannot find toolkit {toolkitName} in {os.path.join(experimentPath, self.CODE_DIRECTORY)}"
                self.logger.error(err)
                raise ValueError(err)

            return toolkitCls(
                projectName=self.projectName,
                pathToExperiment=experimentPath,
                filesDirectory=filesDirectory
            )

        err = (
            f"Experiment {experimentName} not found in Project {self.projectName}. "
            f"Please load the experiment to the project."
        )
        self.logger.error(err)
        raise ValueError(err)

    def keys(self):
        """
        List all experiment names registered in the project.

        Returns
        -------
        list[str]
            Experiment datasource names.
        """
        return [x for x in self.getExperimentsMap()]

    def __getitem__(self, item):
        """
        Dictionary-like access to ``getExperiment(item)``.

        Parameters
        ----------
        item : str
            Experiment name.

        Returns
        -------
        experimentSetupWithData
            Loaded experiment toolkit.
        """
        return self.getExperiment(item)

    def experimentDataType(self):
        """
        Return the configured experiment data type (if present).

        Returns
        -------
        Any
            The attribute ``_experimentDataType`` if defined by a subclass.

        Notes
        -----
        This base file does not define ``self._experimentDataType``.
        """
        return self._experimentDataType


class experimentSetupWithData(argosDataObjects.ExperimentZipFile, toolkit.abstractToolkit):
    """
    experimentSetupWithData
    -----------------------

    A toolkit that unifies:
    - Argos experiment setup (trial sets, entities, metadata in a zip)
    - Hera data engine (parquet/pandas/dask) for retrieving measurement data

    This class is typically the main entry point after loading a specific experiment
    (either directly or via ``experimentHome.getExperiment``).

    Attributes
    ----------
    configuration : dict
        Experiment runtime configuration loaded from
        ``runtimeExperimentData/Datasources_Configurations.json`` (and updated by overrides).
    analysis : experimentAnalysis
        Analysis layer attached to this experiment.
    presentation : experimentPresentation
        Presentation layer attached to this experiment.

    Notes
    -----
    - Requires Argos classes to be importable at runtime (unless you mock them for docs).
    - Builds a cache directory ``<filesDirectory>/experimentCache``.
    """

    _configuration = None
    entityType = None
    trialSet = None

    _analysis = None
    _presentation = None

    @property
    def analysis(self):
        """
        Analysis layer accessor.

        Returns
        -------
        experimentAnalysis
            Analysis layer instance.
        """
        return self._analysis

    @property
    def presentation(self):
        """
        Presentation layer accessor.

        Returns
        -------
        experimentPresentation
            Presentation layer instance.
        """
        return self._presentation

    @property
    def configuration(self):
        """
        Experiment configuration accessor.

        Returns
        -------
        dict
            Runtime experiment configuration.
        """
        return self._configuration

    @property
    def name(self):
        """
        Experiment name accessor.

        Returns
        -------
        str
            ``configuration["experimentName"]``.
        """
        return self.configuration['experimentName']

    def _initTrialSets(self):
        """
        Initialize trial set objects with data access.

        This populates ``self.trialSet`` with ``TrialSetWithData`` entries, keyed
        by the trial set name, based on the Argos setup structure.

        Notes
        -----
        Expects Argos setup structure to include:
        ``self.setup["trialSets"]``.
        """
        experimentSetup = self.setup
        for trialset in experimentSetup['trialSets']:
            self.trialSet[trialset['name']] = TrialSetWithData(
                experiment=self,
                TrialSetSetup=trialset,
                experimentData=self._experimentData
            )

    def _initEntitiesTypes(self):
        """
        Initialize entity type objects with data access.

        This populates ``self.entityType`` with ``EntityTypeWithData`` entries, keyed
        by the entity type name, based on the Argos setup structure.

        Notes
        -----
        Expects Argos setup structure to include:
        ``self.setup["entityTypes"]``.
        """
        experimentSetup = self.setup
        for entityType in experimentSetup['entityTypes']:
            self.entityType[entityType['name']] = EntityTypeWithData(
                experiment=self,
                metadata=entityType,
                experimentData=self._experimentData
            )

    def getExperimentData(self):
        """
        Return the experiment data engine.

        Accessing experiment data is typically done through this data engine using
        ``.getData(...)``.

        Returns
        -------
        object
            One of:
            - parquetDataEngineHera
            - pandasDataEngineDB
            - daskDataEngineDB
            (depending on the configured ``dataType``).

        Notes
        -----
        The exact concrete engine class is created by ``dataEngineFactory()``.
        """
        return self._experimentData

    def __init__(
        self,
        projectName,
        pathToExperiment,
        dataType=PARQUETHERA,
        dataSourceConfiguration=dict(),
        filesDirectory=None,
        defaultTrialSetName=None
    ):
        """
        Initialize the specific experiment toolkit.

        Parameters
        ----------
        projectName : str
            The project name to work with.
        pathToExperiment : str
            Absolute path to the experiment folder on disk.
            Must include ``runtimeExperimentData`` and experiment setup zip.
        dataType : str, optional
            Define how the data is retrieved:
            - PARQUETHERA (default)
            - PANDASDB
            - DASKDB
        dataSourceConfiguration : dict, optional
            Overrides/updates the experiment configuration loaded from
            ``Datasources_Configurations.json``.
        filesDirectory : str, optional
            Directory to store cache/intermediate files.
            If None, uses current working directory.
        defaultTrialSetName : str, optional
            Default trial set name to use when the user does not specify a trial set.

        Raises
        ------
        ValueError
            If the configuration file or setup zip cannot be found.
        """
        configurationFileName = os.path.join(
            pathToExperiment, 'runtimeExperimentData', "Datasources_Configurations.json"
        )
        if not os.path.isfile(configurationFileName):
            raise ValueError(f"The configuration file doesn't exist. Looking for {configurationFileName}")

        self._configuration = loadJSON(configurationFileName)

        dataSourceConfiguration = dict() if dataSourceConfiguration is None else dataSourceConfiguration
        self._configuration.update(dataSourceConfiguration)

        experimentName = self.configuration['experimentName']
        setupFile = os.path.join(pathToExperiment, 'runtimeExperimentData', f"{experimentName}.zip")
        if not os.path.isfile(setupFile):
            raise ValueError(f"The experiment setup file doesn't exist. Looking for {setupFile}")

        # Initialize the data engine.
        self._experimentData = dataEngineFactory().getDataEngine(
            projectName, self._configuration, experimentObj=self, dataType=dataType
        )

        self.entityType = dict()
        self.trialSet = dict()

        if filesDirectory is None:
            filesDirectory = os.getcwd()

        cacheDir = os.path.join(filesDirectory, "experimentCache")
        os.makedirs(cacheDir, exist_ok=True)

        # Init Argos setup + Hera toolkit base
        argosDataObjects.ExperimentZipFile.__init__(self, setupFile)
        toolkit.abstractToolkit.__init__(
            self,
            projectName=projectName,
            toolkitName=f"{experimentName}Toolkit",
            filesDirectory=cacheDir
        )

        self._defaultTrialSetName = defaultTrialSetName
        self._analysis = experimentAnalysis(self,)
        self._presentation = experimentPresentation(self, self.analysis)

    @property
    def defaultTrialSet(self):
        """
        Default trial set name.

        Returns
        -------
        str or None
            Default trial set name if configured.
        """
        return self._defaultTrialSetName

    @property
    def trialsOfDefaultTrialSet(self):
        """
        Convenience accessor for trials in the default trial set.

        Returns
        -------
        TrialSetWithData
            Trial set object corresponding to ``defaultTrialSet``.

        Raises
        ------
        KeyError
            If ``defaultTrialSet`` is None or not found in ``self.trialSet``.
        """
        return self.trialSet[self.defaultTrialSet]

    def _initAnalysisAndPresentation(self, analysisCLS, presentationCLS):
        """
        Initialize analysis and presentation layers using custom classes.

        Parameters
        ----------
        analysisCLS : type
            Analysis layer class. Recommended to inherit from ``experimentAnalysis``.
        presentationCLS : type
            Presentation layer class. Recommended to inherit from ``experimentPresentation``.

        Notes
        -----
        This method replaces the current analysis and presentation objects.
        """
        self._analysis = analysisCLS(self)
        self._presentation = presentationCLS(self, self._analysis)

    def getDataFromDateRange(self, deviceType, startTime, endTime, deviceName=None, withMetadata=True):
        """
        Retrieve device data within a date range, optionally enriched with metadata.

        Parameters
        ----------
        deviceType : str
            The device type (as defined by the experiment schema / Argos).
        startTime : datetime-like
            Range start time (inclusive). Passed to the data engine.
        endTime : datetime-like
            Range end time (inclusive/exclusive depending on engine implementation).
        deviceName : str, optional
            Optional device name filter.
        withMetadata : bool, optional
            If True, merges device metadata (entities table) into the returned data.

        Returns
        -------
        pandas.DataFrame
            Data frame indexed by timestamp (after merge) containing measurement data,
            possibly with metadata columns added.

        Raises
        ------
        ValueError
            If no data is found in the specified date range.
        """
        data = self._experimentData.getData(
            deviceType=deviceType,
            deviceName=deviceName,
            startTime=startTime,
            endTime=endTime
        )

        if len(data) == 0:
            raise ValueError(f"There is no data for {deviceType} between the dates {startTime} and {endTime}")

        if withMetadata:
            devicemetadata = self.entitiesTable()
            if len(devicemetadata) > 0:
                data = (
                    data.reset_index()
                    .merge(devicemetadata, left_on="deviceName", right_on="entityName")
                    .set_index("timestamp")
                )

        return data

    def _process_row(self, row):
        """
        Convert a single row's WGS84 coordinates into ITM coordinates.

        Parameters
        ----------
        row : pandas.Series
            A row that contains ``row.Longitude`` and ``row.Latitude``.

        Returns
        -------
        pandas.Series
            Two values: [ITM_x, ITM_y] (as returned by ``convertCRS``).
        """
        pp = convertCRS([[row.Longitude, row.Latitude]], inputCRS=WSG84, outputCRS=ITM)
        return pd.Series([pp.x[0], pp.y[0]])

    def get_devices_image_coordinates(self, trialSetName, trialName, deviceType, outputCRS=ITM):
        """
        Compute a bounding box around device coordinates for a trial.

        This is useful for requesting map tiles or images that include all devices.

        Parameters
        ----------
        trialSetName : str
            Trial set name.
        trialName : str
            Trial name inside the trial set.
        deviceType : str
            Device type name used for filtering entities metadata.
        outputCRS : int, optional
            Desired CRS for output coordinates.
            - If ITM: converts lat/lon to ITM and returns bbox in meters.
            - Else: returns bbox in original WGS84 values.

        Returns
        -------
        tuple
            (min_latitude, min_longitude, max_latitude, max_longitude)
            Note: "latitude/longitude" names remain for backward compatibility,
            even when outputCRS is ITM (then values are ITM coordinates).

        Notes
        -----
        - Uses ``entitiesTable`` for the given trial.
        - For ITM conversion, creates columns: ``ITM_Latitude``, ``ITM_Longitude``.
        """
        devices_df = self.trialSet[trialSetName][trialName].entitiesTable.query("deviceTypeName==@deviceType")

        if outputCRS == ITM:
            devices_df[['ITM_Latitude', 'ITM_Longitude']] = devices_df.apply(self._process_row, axis=1)
            latitudes = devices_df['ITM_Latitude']
            longitudes = devices_df['ITM_Longitude']
        else:
            latitudes = devices_df['Latitude']
            longitudes = devices_df['Longitude']

        min_latitude, max_latitude = min(latitudes), max(latitudes)
        min_longitude, max_longitude = min(longitudes), max(longitudes)
        return min_latitude, min_longitude, max_latitude, max_longitude


class TrialSetWithData(argosDataObjects.TrialSet):
    """
    TrialSetWithData
    ----------------

    Argos TrialSet extended with a Hera data engine reference.

    Each TrialSet holds multiple trials. This wrapper ensures each created trial
    is a ``TrialWithdata`` that can fetch its measurement data using the engine.
    """

    def _initTrials(self):
        """
        Initialize trials within this trial set as ``TrialWithdata`` objects.

        Notes
        -----
        Expects Argos metadata structure to include:
        ``self._metadata["trials"]``.
        """
        for trial in self._metadata['trials']:
            self[trial['name']] = TrialWithdata(
                trialSet=self,
                metadata=trial,
                experimentData=self._experimentData
            )

    def __init__(self, experiment: experimentSetupWithData, TrialSetSetup: dict, experimentData: dataEngineFactory):
        """
        Initialize the trial set with access to the experiment's data engine.

        Parameters
        ----------
        experiment : experimentSetupWithData
            Parent experiment toolkit instance.
        TrialSetSetup : dict
            Trial set setup metadata.
        experimentData : object
            Data engine instance responsible for retrieving data.
        """
        self._experimentData = experimentData
        super().__init__(experiment, TrialSetSetup)


class TrialWithdata(argosDataObjects.Trial):
    """
    TrialWithdata
    -------------

    Argos Trial extended with a Hera data engine reference.

    Adds ``getData(...)`` with sensible defaults for start/end time based on the
    trial properties (TRIALSTART/TRIALEND).
    """

    def getData(self, deviceType, deviceName=None, startTime=None, endTime=None, withMetadata=False):
        """
        Retrieve data for this trial.

        Parameters
        ----------
        deviceType : str
            Device type to query.
        deviceName : str, optional
            Specific device to query. If None, returns all devices of this type.
        startTime : datetime-like, optional
            If None, defaults to the trial's ``TRIALSTART`` property.
        endTime : datetime-like, optional
            If None, defaults to the trial's ``TRIALEND`` property.
        withMetadata : bool, optional
            If True, merges entities metadata into the returned data.

        Returns
        -------
        pandas.DataFrame
            Data frame indexed by timestamp (after merge) containing trial data.

        Raises
        ------
        ValueError
            If no data exists for this trial in the requested interval.
        """
        startTime = self.properties[TRIALSTART] if startTime is None else startTime
        endTime = self.properties[TRIALEND] if endTime is None else endTime

        data = self._experimentData.getData(
            deviceType=deviceType,
            deviceName=deviceName,
            startTime=startTime,
            endTime=endTime
        )

        if len(data) == 0:
            raise ValueError(f"There is no data for {deviceType} between the dates {startTime} and {endTime}")

        if withMetadata:
            devicemetadata = self.entitiesTable()
            if len(devicemetadata) > 0:
                data = (
                    data.reset_index()
                    .merge(devicemetadata, left_on="deviceName", right_on="entityName")
                    .set_index("timestamp")
                )

        return data

    def __init__(self, trialSet: TrialSetWithData, metadata: dict, experimentData: dataEngineFactory):
        """
        Initialize the trial wrapper.

        Parameters
        ----------
        trialSet : TrialSetWithData
            Parent trial set object.
        metadata : dict
            Trial metadata (Argos schema).
        experimentData : object
            Data engine instance used to fetch trial data.
        """
        self._experimentData = experimentData
        super().__init__(trialSet, metadata)


class EntityTypeWithData(argosDataObjects.EntityType):
    """
    EntityTypeWithData
    ------------------

    Argos EntityType extended with a Hera data engine reference.

    Each entity type contains multiple entities (devices). This wrapper ensures
    each created entity is an ``EntityWithData`` that can query the data engine.
    """

    def _initEntities(self):
        """
        Initialize entities of this entity type as ``EntityWithData`` objects.

        Notes
        -----
        Expects Argos metadata structure to include:
        ``self._metadata["entities"]``.
        """
        for entity in self._metadata['entities']:
            self[entity['name']] = EntityWithData(
                entityType=self,
                metadata=entity,
                experimentData=self._experimentData
            )

    def __init__(self, experiment: experimentSetupWithData, metadata: dict, experimentData: dataEngineFactory):
        """
        Initialize entity type wrapper with access to the data engine.

        Parameters
        ----------
        experiment : experimentSetupWithData
            Parent experiment toolkit instance.
        metadata : dict
            Entity type metadata (Argos schema).
        experimentData : object
            Data engine instance used to fetch device data.
        """
        self._experimentData = experimentData
        super().__init__(experiment, metadata)

    def getData(self, startTime=None, endTime=None):
        """
        Retrieve data for this entity type (device type) over a time range.

        Parameters
        ----------
        startTime : datetime-like, optional
            Start time for query.
        endTime : datetime-like, optional
            End time for query.

        Returns
        -------
        pandas.DataFrame
            Data for this entity type from the data engine.
        """
        return self._experimentData.getData(self.name, startTime=startTime, endTime=endTime)

    def getDataTrial(self, trialSetName, trialName):
        """
        Retrieve data for a specific entity within a specific trial.

        This method:
        1) locates the trial object
        2) uses trial start/end times (TRIALSTART/TRIALEND)
        3) queries the data engine, optionally per-device depending on metadata

        Parameters
        ----------
        trialSetName : str
            Trial set name.
        trialName : str
            Trial name.

        Returns
        -------
        pandas.DataFrame
            Retrieved measurement data for this entity within the trial time window.

        Notes
        -----
        Uses ``self.properties["StoreDataPerDevice"]`` to choose per-device storage mode.
        """
        trial = self.experiment.trialSet[trialSetName][trialName]
        startTime = trial.properties[TRIALSTART]
        endTime = trial.properties[TRIALEND]

        StoreDataPerDevice = self.properties['StoreDataPerDevice']
        data = self._experimentData.getData(
            deviceType=self.entityType,
            deviceName=self.name,
            startTime=startTime,
            endTime=endTime,
            perDevice=StoreDataPerDevice
        )
        return data


class EntityWithData(argosDataObjects.Entity):
    """
    EntityWithData
    --------------

    Argos Entity extended with a Hera data engine reference.

    Provides ``getData`` to retrieve this specific device/entity measurements.
    """

    def __init__(self, entityType: EntityTypeWithData, metadata: dict, experimentData):
        """
        Initialize entity wrapper.

        Parameters
        ----------
        entityType : EntityTypeWithData
            Parent entity type object.
        metadata : dict
            Entity metadata (Argos schema).
        experimentData : object
            Data engine instance used to fetch device data.
        """
        self._experimentData = experimentData
        super().__init__(entityType, metadata)

    def getData(self, startTime=None, endTime=None):
        """
        Retrieve data for this entity (device) over a time range.

        Parameters
        ----------
        startTime : datetime-like, optional
            Query start time.
        endTime : datetime-like, optional
            Query end time.

        Returns
        -------
        pandas.DataFrame
            Retrieved data for this entity from the data engine.

        Notes
        -----
        Uses ``self.properties["StoreDataPerDevice"]`` to choose per-device storage mode.
        """
        StoreDataPerDevice = self.properties['StoreDataPerDevice']
        return self._experimentData.getData(
            deviceType=self.entityType,
            deviceName=self.name,
            startTime=startTime,
            endTime=endTime,
            perDevice=StoreDataPerDevice
        )
