import os
import sys
import logging
import pydoc
import pandas as pd

from hera import toolkit
from .presentation import experimentPresentation
from .analysis import experimentAnalysis
from hera.measurements.GIS.utils import WSG84, ITM, convertCRS

try:
    from argos.experimentSetup import dataObjects as argosDataObjects
except ImportError:
    # Argos is optional; if it is not installed, experiment toolkits cannot be used.
    print("Must have argos installed and in the path.")

from .dataEngine import dataEngineFactory, PARQUETHERA, PANDASDB, DASKDB
from hera.utils import loadJSON

# The name of the properties. These must match the Argos web interface.
# Do not change unless the Argos schema changes.
TRIALSTART = "TrialStart"
TRIALEND = "TrialEnd"


class experimentHome(toolkit.abstractToolkit):
    """
    This object functions as a factory/home for the other experiments.
    It is responsible for getting the right toolkit for the requested experiment.
    """

    DOCTYPE_ENTITIES = "EntitiesData"
    CODE_DIRECTORY = "code"

    def __init__(self, projectName, filesDirectory=None):
        super().__init__(
            projectName=projectName,
            toolkitName="experimentToolKit",
            filesDirectory=filesDirectory,
        )
        self.logger = logging.getLogger()
        self.logger.info("Init experiment toolkit")

    @property
    def experimentMap(self):
        """
        Backward-compatibility alias.
        Historically this was a property, so we keep the interface,
        even though today the real logic is in getExperimentsMap().
        """
        return self.getExperimentsMap()

    def getExperimentsMap(self):
        """
        Get a dictionary mapping experiment name -> datasource document.

        Returns
        -------
        dict
            Keys are experiment names (datasourceName),
            values are the matching datasource documents.
        """
        experiments_map = {}
        for experiment in self.getDataSourceMap():
            experimentName = experiment["datasourceName"]
            experiments_map[experimentName] = experiment
        return experiments_map

    @property
    def experimentsTable(self):
        """
        Return a tabular view (DataFrame-like) of experiment datasources.
        """
        return self.getDataSourceTable()

    def getExperimentsTable(self):
        """
        Backward-compatible alias for experimentsTable.
        """
        return self.getDataSourceTable()

    def getExperiment(self, experimentName, filesDirectory=None):
        """
        Get the specific experiment class.

        Parameters
        ----------
        experimentName : str
            The name of the experiment.
        filesDirectory : str, optional
            The directory to save cache/intermediate files.
            If None, uses the current working directory / 'experimentCache'.

        Returns
        -------
        experimentSetupWithData
        """
        self.logger.info(f"Getting experiment {experimentName}")
        ds_doc = self.getDataSourceDocument(datasourceName=experimentName)

        if ds_doc:
            self.logger.info("Found experiment. Loading")
            experimentPath = ds_doc.getData()

            # Add experiment's 'code' directory to sys.path so we can import its toolkit.
            sys.path.append(os.path.join(experimentPath, self.CODE_DIRECTORY))
            self.logger.debug(
                f"Adding path {os.path.join(experimentPath, self.CODE_DIRECTORY)} to sys.path"
            )

            toolkitName = f"{experimentName}.{experimentName}"
            self.logger.debug(f"Loading toolkit: {toolkitName}")

            toolkitCls = pydoc.locate(toolkitName)
            if toolkitCls is None:
                err = (
                    f"Cannot find toolkit {toolkitName} in "
                    f"{os.path.join(experimentPath, self.CODE_DIRECTORY)}"
                )
                self.logger.error(err)
                raise ValueError(err)

            return toolkitCls(
                projectName=self.projectName,
                pathToExperiment=experimentPath,
                filesDirectory=filesDirectory,
            )

        err = (
            f"Experiment {experimentName} not found in Project {self.projectName}. "
            f"Please load the experiment to the project."
        )
        self.logger.error(err)
        raise ValueError(err)

    def keys(self):
        """
        Get the experiment names of the project.

        Returns
        -------
        list of str
        """
        return [name for name in self.getExperimentsMap()]

    def __getitem__(self, item):
        """
        Allow experimentHome['expName'] syntax to return a specific experiment.
        """
        return self.getExperiment(item)

    def experimentDataType(self):
        """
        Backward-compatibility hook for experiment data type.
        """
        return getattr(self, "_experimentDataType", None)


class experimentSetupWithData(argosDataObjects.ExperimentZipFile, toolkit.abstractToolkit):
    """
    A class that unifies the argos.experiment setup with the data.
    """

    _configuration = None
    entityType = None
    trialSet = None

    _analysis = None
    _presentation = None

    @property
    def analysis(self):
        return self._analysis

    @property
    def presentation(self):
        return self._presentation

    @property
    def configuration(self):
        return self._configuration

    @property
    def name(self):
        return self.configuration["experimentName"]

    def _initTrialSets(self):
        """
        Initialize trial sets from the experiment setup metadata.
        """
        experimentSetup = self.setup
        for trialset in experimentSetup["trialSets"]:
            self.trialSet[trialset["name"]] = TrialSetWithData(
                experiment=self,
                TrialSetSetup=trialset,
                experimentData=self._experimentData,
            )

    def _initEntitiesTypes(self):
        """
        Initialize entity types from the experiment setup metadata.
        """
        experimentSetup = self.setup
        for entityType in experimentSetup["entityTypes"]:
            self.entityType[entityType["name"]] = EntityTypeWithData(
                experiment=self,
                metadata=entityType,
                experimentData=self._experimentData,
            )

    def getExperimentData(self):
        """
        Get the Data Engine of the experiment.
        Accessing experiment data is done through this object (using .getData()).

        Returns
        -------
        parquetDataEngineHera or pandasDataEngineDB or daskDataEngineDB
        """
        return self._experimentData

    def __init__(
        self,
        projectName,
        pathToExperiment,
        dataType=PARQUETHERA,
        dataSourceConfiguration=dict(),
        filesDirectory=None,
        defaultTrialSetName=None,
    ):
        """
        Initialize the specific experiment toolkit.

        Parameters
        ----------
        projectName : str
            The project name to work with.
        pathToExperiment : str
            The path to the experiment data.
        dataType : str
            How data is retrieved: dask/pandas from MongoDB, or parquet.
        dataSourceConfiguration : dict
            Override datasources configuration of the experiment.
        filesDirectory : str, optional
            Directory to save cache/intermediate files.
            If None, uses [current directory]/experimentCache.
        defaultTrialSetName : str, optional
            Default trial set to use if not supplied.
        """
        # Locate configuration file
        configurationFileName = os.path.join(
            pathToExperiment, "runtimeExperimentData", "Datasources_Configurations.json"
        )

        if not os.path.isfile(configurationFileName):
            raise ValueError(
                f"The configuration file does not exist. Looking for {configurationFileName}"
            )
        self._configuration = loadJSON(configurationFileName)

        dataSourceConfiguration = dict() if dataSourceConfiguration is None else dataSourceConfiguration
        self._configuration.update(dataSourceConfiguration)

        experimentName = self.configuration["experimentName"]
        setupFile = os.path.join(
            pathToExperiment, "runtimeExperimentData", f"{experimentName}.zip"
        )
        if not os.path.isfile(setupFile):
            raise ValueError(
                f"The experiment setup file does not exist. Looking for {setupFile}"
            )

        # Initialize the data engine
        self._experimentData = dataEngineFactory().getDataEngine(
            projectName,
            self._configuration,
            experimentObj=self,
            dataType=dataType,
        )
        self.entityType = dict()
        self.trialSet = dict()

        if filesDirectory is None:
            filesDirectory = os.getcwd()

        cacheDir = os.path.join(filesDirectory, "experimentCache")
        os.makedirs(cacheDir, exist_ok=True)

        argosDataObjects.ExperimentZipFile.__init__(self, setupFile)
        toolkit.abstractToolkit.__init__(
            self,
            projectName=projectName,
            toolkitName=f"{experimentName}Toolkit",
            filesDirectory=cacheDir,
        )

        self._defaultTrialSetName = defaultTrialSetName
        self._analysis = experimentAnalysis(self)
        self._presentation = experimentPresentation(self, self.analysis)

    @property
    def defaultTrialSet(self):
        return self._defaultTrialSetName

    @property
    def trialsOfDefaultTrialSet(self):
        return self.trialSet[self.defaultTrialSet]

    def _initAnalysisAndPresentation(self, analysisCLS, presentationCLS):
        """
        Initialize the analysis and presentation classes and set the data layer.

        Parameters
        ----------
        analysisCLS : class
            The analysis class, recommended to inherit from .analysis.experimentAnalysis.
        presentationCLS : class
            The presentation class, recommended to inherit from .presentation.experimentPresentation.
        """
        self._analysis = analysisCLS(self)
        self._presentation = presentationCLS(self, self._analysis)

    def getDataFromDateRange(
        self,
        deviceType,
        startTime,
        endTime,
        deviceName=None,
        withMetadata=True,
    ):
        """
        Retrieve data for a given device type and time range.

        Parameters
        ----------
        deviceType : str
        startTime : datetime-like
        endTime : datetime-like
        deviceName : str, optional
        withMetadata : bool

        Returns
        -------
        DataFrame
        """
        data = self._experimentData.getData(
            deviceType=deviceType,
            deviceName=deviceName,
            startTime=startTime,
            endTime=endTime,
        )

        if len(data) == 0:
            raise ValueError(
                f"There is no data for {deviceType} between the dates {startTime} and {endTime}"
            )

        if withMetadata:
            devicemetadata = self.entitiesTable()
            if len(devicemetadata) > 0:
                data = (
                    data.reset_index()
                    .merge(
                        devicemetadata,
                        left_on="deviceName",
                        right_on="entityName",
                    )
                    .set_index("timestamp")
                )

        return data

    def _process_row(self, row):
        """
        Helper for coordinate conversion of a single row (Longitude, Latitude).
        """
        pp = convertCRS([[row.Longitude, row.Latitude]], inputCRS=WSG84, outputCRS=ITM)
        return pd.Series([pp.x[0], pp.y[0]])

    def get_devices_image_coordinates(
        self,
        trialSetName,
        trialName,
        deviceType,
        outputCRS=ITM,
    ):
        """
        Compute bounding box of devices in image coordinates (ITM or original).

        Returns
        -------
        (min_latitude, min_longitude, max_latitude, max_longitude)
        """
        devices_df = self.trialSet[trialSetName][trialName].entitiesTable.query(
            "deviceTypeName==@deviceType"
        )

        if outputCRS == ITM:
            devices_df[["ITM_Latitude", "ITM_Longitude"]] = devices_df.apply(
                self._process_row, axis=1
            )
            latitudes = devices_df["ITM_Latitude"]
            longitudes = devices_df["ITM_Longitude"]
        else:
            latitudes = devices_df["Latitude"]
            longitudes = devices_df["Longitude"]

        min_latitude, max_latitude = min(latitudes), max(latitudes)
        min_longitude, max_longitude = min(longitudes), max(longitudes)
        return min_latitude, min_longitude, max_latitude, max_longitude


class TrialSetWithData(argosDataObjects.TrialSet):
    """
    TrialSet that is aware of the experiment data engine.
    """

    def _initTrials(self):
        for trial in self._metadata["trials"]:
            self[trial["name"]] = TrialWithdata(
                trialSet=self,
                metadata=trial,
                experimentData=self._experimentData,
            )

    def __init__(self, experiment: experimentSetupWithData, TrialSetSetup: dict, experimentData: dataEngineFactory):
        """
        Initialize the TrialSet with a link to the shared experiment data engine.
        """
        self._experimentData = experimentData
        super().__init__(experiment, TrialSetSetup)


class TrialWithdata(argosDataObjects.Trial):
    """
    Trial object that knows how to pull its data from the shared experiment data engine.
    """

    def getData(
        self,
        deviceType,
        deviceName=None,
        startTime=None,
        endTime=None,
        withMetadata=False,
    ):
        """
        Retrieve trial data for a given device type and time range.
        If startTime/endTime are not provided, the trial properties TRIALSTART/TRIALEND are used.
        """
        startTime = self.properties[TRIALSTART] if startTime is None else startTime
        endTime = self.properties[TRIALEND] if endTime is None else endTime

        data = self._experimentData.getData(
            deviceType=deviceType,
            deviceName=deviceName,
            startTime=startTime,
            endTime=endTime,
        )

        if len(data) == 0:
            raise ValueError(
                f"There is no data for {deviceType} between the dates {startTime} and {endTime}"
            )

        if withMetadata:
            devicemetadata = self.entitiesTable()
            if len(devicemetadata) > 0:
                data = (
                    data.reset_index()
                    .merge(
                        devicemetadata,
                        left_on="deviceName",
                        right_on="entityName",
                    )
                    .set_index("timestamp")
                )

        return data

    def __init__(
        self,
        trialSet: TrialSetWithData,
        metadata: dict,
        experimentData: dataEngineFactory,
    ):
        self._experimentData = experimentData
        super().__init__(trialSet, metadata)


class EntityTypeWithData(argosDataObjects.EntityType):
    """
    EntityType that knows how to pull its data from the shared experiment data engine.
    """

    def _initEntities(self):
        for entity in self._metadata["entities"]:
            self[entity["name"]] = EntityWithData(
                entityType=self,
                metadata=entity,
                experimentData=self._experimentData,
            )

    def __init__(
        self,
        experiment: experimentSetupWithData,
        metadata: dict,
        experimentData: dataEngineFactory,
    ):
        """
        Initialize the EntityType with a link to the shared experiment data engine.
        """
        self._experimentData = experimentData
        super().__init__(experiment, metadata)

    def getData(self, startTime=None, endTime=None):
        """
        Retrieve all data for this entity type (optionally in a time range).
        """
        return self._experimentData.getData(
            self.name,
            startTime=startTime,
            endTime=endTime,
        )

    def getDataTrial(self, trialSetName, trialName):
        """
        Return the device data for this entity type in a specific trial.

        Parameters
        ----------
        trialSetName : str
        trialName : str

        Returns
        -------
        DataFrame
        """
        trial = self.experiment.trialSet[trialSetName][trialName]
        startTime = trial.properties[TRIALSTART]
        endTime = trial.properties[TRIALEND]

        StoreDataPerDevice = self.properties["StoreDataPerDevice"]
        data = self._experimentData.getData(
            deviceType=self.entityType,
            deviceName=self.name,
            startTime=startTime,
            endTime=endTime,
            perDevice=StoreDataPerDevice,
        )
        return data


class EntityWithData(argosDataObjects.Entity):
    """
    Entity that knows how to pull its data from the shared experiment data engine.
    """

    def __init__(
        self,
        entityType: EntityTypeWithData,
        metadata: dict,
        experimentData,
    ):
        self._experimentData = experimentData
        super().__init__(entityType, metadata)

    def getData(self, startTime=None, endTime=None):
        """
        Retrieve data for this specific entity (device).
        """
        StoreDataPerDevice = self.properties["StoreDataPerDevice"]

        return self._experimentData.getData(
            deviceType=self.entityType,
            deviceName=self.name,
            startTime=startTime,
            endTime=endTime,
            perDevice=StoreDataPerDevice,
        )
