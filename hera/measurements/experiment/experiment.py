import json
import os
import pydoc
import sys
from ... import toolkit,toolkitHome,datalayer
from .presentation import experimentPresentation
from .analysis import experimentAnalysis
from hera.measurements.GIS.utils import WSG84,ITM,convertCRS
import pandas as pd

try:
    from argos.experimentSetup import dataObjects as argosDataObjects
except ImportError:
    print("Must have argos installed and in the path. ")

from .dataEngine import dataEngineFactory, PARQUETHERA, PANDASDB,DASKDB
from ...utils  import loadJSON
import logging

class experimentHome(toolkit.abstractToolkit):
    """
        This is the object that function as a factory/home to the other experiments.
        It is responsible for getting the right toolkit for the requested experiment.

    """

    DOCTYPE_ENTITIES = 'EntitiesData'
    CODE_DIRECTORY = 'code'

    def __init__(self, projectName, filesDirectory=None):
        super().__init__(projectName=projectName, toolkitName="experimentToolKit", filesDirectory=filesDirectory)
        self.logger = logging.getLogger()
        self.logger.info("Init experiment toolkit")


    @property
    def experimentMap(self):
        return self.experimentMap()

    #
    def getExperimentsMap(self):
        """
        Get dictionary of experiments map of project.

        Returns
        -------
            dict
        """
        M=dict()
        for experiment in self.getDataSourceMap():
            experimentName=experiment['datasourceName']
            M[experimentName]=experiment

        return M

    @property
    def experimentsTable(self):
        return self.getDataSourceTable()

    def getExperimentsTable(self):
        return self.getDataSourceTable()

    def getExperiment(self,experimentName,filesDirectory=None):
        """
        Get the specific experiment class.

        Parameters
        ----------
        experimentName : str
                The name of the experimen
        filesDirectory: str
                The directory to save the cache/intermediate files.
                If None, use the [current directory]/experimentCache.

        Returns
        -------
            experimentSetupWithData
        """

        self.logger.info(f"Getting experiment {experimentName}")
        L = self.getDataSourceDocument(datasourceName=experimentName)
        if L:
            self.logger.info(f"Found experiment. Loading")
            experimentPath=L.getData()
            sys.path.append(os.path.join(experimentPath,self.CODE_DIRECTORY))
            self.logger.debug(f"Adding path {os.path.join(experimentPath,self.CODE_DIRECTORY)} to classpath")
            toolkitName = f"{experimentName}.{experimentName}"
            self.logger.debug(f"Loading toolkits: {toolkitName}")

            toolkitCls = pydoc.locate(toolkitName)
            if toolkitCls is None:
                err = f"Cannot find toolkit {toolkitName} in {os.path.join(experimentPath,self.CODE_DIRECTORY)}"
                self.logger.error(err)
                raise ValueError(err)

            return toolkitCls(projectName=self.projectName,
                              pathToExperiment=experimentPath,filesDirectory=filesDirectory)
        else:
            err = f"Experiment {experimentName} not found in Project {self.projectName}. Please load the experiment to the project. "
            self.logger.error(err)
            raise ValueError(err)


    def keys(self):
        """
        Get the experiments names of project.

        Returns
        -------
            list
        """
        return [x for x in self.getExperimentsMap()]

    def __getitem__(self, item):
        return self.getExperiment(item)

    def experimentDataType(self):
        return self._experimentDataType

class experimentSetupWithData(argosDataObjects.ExperimentZipFile,toolkit.abstractToolkit):
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
        return self.configuration['experimentName']

    def _initTrialSets(self):
        experimentSetup = self.setup
        for trialset in experimentSetup['trialSets']:
            self.trialSet[trialset['name']] = TrialSetWithData(experiment = self, TrialSetSetup=trialset,experimentData= self._experimentData)

    def _initEntitiesTypes(self):
        experimentSetup = self.setup
        for entityType in experimentSetup['entityTypes']:
            self.entityType[entityType['name']] = EntityTypeWithData(experiment=self, metadata = entityType, experimentData= self._experimentData)

    def getExperimentData(self):
        """
        Get the parquet Data Engine of experiment. Acessing data of experiment is through this class (using .getData()).

        Returns
        -------
            parquetDataEngineHera , pandasDataEngineDB or daskDataEngineDB.
        """
        return self._experimentData

    def __init__(self, projectName, pathToExperiment, dataType=PARQUETHERA, dataSourceConfiguration=dict(), filesDirectory=None,defaultTrialSetName=None):
        """
            Initializes the specific experiment toolkit.

        Parameters
        ----------
        projectName: str
                The project name to work with.

        pathToExperiment:
                The path to the experiment data.

        dataType: str
                Define how the data is retrieved: dask or pandas directly from the mongoDB, or through the
                                                  parquet.

        dataSourceConfiguration : dict
                overwrite the datasources configuration of the experiment.
                See ... for structure.

        filesDirectory: str
                The directory to save the cache/intermediate files.
                If None, use the [current directory]/experimentCache.

        defaultTrialSet: str
                A default trialset to use if not supplied.

        """
        # setup the configuration file name
        configurationFileName = os.path.join(pathToExperiment, 'runtimeExperimentData', "Datasources_Configurations.json")

        if not os.path.isfile(configurationFileName):
            raise ValueError(f" The configuration file doesn't exist. Looking for {configurationFileName}")
        self._configuration = loadJSON(configurationFileName)


        dataSourceConfiguration = dict() if dataSourceConfiguration is None else dataSourceConfiguration
        self._configuration.update(dataSourceConfiguration)

        experimentName = self.configuration['experimentName']
        setupFile = os.path.join(pathToExperiment, 'runtimeExperimentData', f"{experimentName}.zip" )
        if not os.path.isfile(setupFile):
            raise ValueError(f"The experiment setup file doesn't exist. Looking for {setupFile}  ")

        # Now initialize the data engine.
        self._experimentData = dataEngineFactory().getDataEngine(projectName,self._configuration,experimentObj=self, dataType = dataType)
        self.entityType = dict()
        self.trialSet = dict()

        if filesDirectory is None:
            filesDirectory = os.getcwd()

        cacheDir = os.path.join(filesDirectory, "experimentCache")
        os.makedirs(cacheDir,exist_ok=True)

        argosDataObjects.ExperimentZipFile.__init__(self,setupFile)
        toolkit.abstractToolkit.__init__(self,projectName=projectName,toolkitName=f"{experimentName}Toolkit",filesDirectory=cacheDir)

        self._defaultTrialSetName = defaultTrialSetName
        self._analysis = experimentAnalysis(self,)
        self._presentation = experimentPresentation(self,self.analysis)

    @property
    def defaultTrialSet(self):
        return self._defaultTrialSetName

    @property
    def trialsOfDefaultTrialSet(self):
        return self.trialSet[self.defaultTrialSet]

    def _initAnalysisAndPresentation(self,analysisCLS,presentationCLS):
        """
            Initializes the analysis and the presentation classes
            and sets the datalayer.

        Parameters
        ----------
        analysisCLS :  class
                The analysis class. It is recommended that it will inherit from
                .analysis.experimentAnalysis

        presentationCLS : class
                The presentation class. It is recommended that it will inherit from
                .presentation.experimentPresentation


        Returns
        -------

        """
        self._analysis = analysisCLS(self)
        self._presentation = presentationCLS(self,self._analysis)

    def getDataFromTrial(self,deviceType,deviceName = None,startTime = None, endTime = None,withMetadata = True):

        startTime = self.properties['TrialStart'] if startTime is None else startTime
        endTime = self.properties['TrialEnd'] if endTime is None else endTime

        data = self._experimentData.getData(deviceType=deviceType,deviceName=deviceName,startTime=startTime,endTime=endTime)

        if len(data) == 0:
            raise ValueError(f"There is no data for {deviceType} between the dates {startTime} and {endTime}")

        if withMetadata:
            devicemetadata = self.entitiesTable()
            if len(devicemetadata) > 0:
                data = data.reset_index().merge(devicemetadata, left_on="deviceName", right_on="entityName").set_index(
                    "timestamp")

        return data

    def _process_row(self,row):
        pp = convertCRS([[row.Longitude, row.Latitude]], inputCRS=WSG84, outputCRS=ITM)
        return pd.Series([pp.x[0], pp.y[0]])

    def get_devices_image_coordinates(self,trialSetName,trialName,deviceType,outputCRS=ITM):
        devices_df = self.trialSet[trialSetName][trialName].entitiesTable.query("deviceTypeName==@deviceType")

        if outputCRS==ITM:
            devices_df[['ITM_Latitude', 'ITM_Longitude']] = devices_df.apply(self._process_row, axis=1)
            latitudes = devices_df['ITM_Latitude']
            longitudes = devices_df['ITM_Longitude']
        else:
            latitudes = devices_df['Latitude']
            longitudes = devices_df['Longitude']
        min_latitude, max_latitude = min(latitudes), max(latitudes)
        min_longitude, max_longitude = min(longitudes), max(longitudes)
        return min_latitude,min_longitude,max_latitude,max_longitude

class TrialSetWithData(argosDataObjects.TrialSet):

    def _initTrials(self):
        for trial in self._metadata['trials']:
            self[trial['name']] = TrialWithdata(trialSet=self,metadata=trial, experimentData =self._experimentData )

    def __init__(self, experiment:experimentSetupWithData, TrialSetSetup: dict, experimentData: dataEngineFactory):
        self._experimentData = experimentData
        super().__init__(experiment, TrialSetSetup)


class TrialWithdata(argosDataObjects.Trial):

    def getData(self,deviceType,deviceName = None,startTime = None, endTime = None,withMetadata = False):

        startTime = self.properties['TrialStart'] if startTime is None else startTime
        endTime = self.properties['TrialEnd'] if endTime is None else endTime

        data = self._experimentData.getData(deviceType=deviceType,deviceName=deviceName,startTime=startTime,endTime=endTime)

        if len(data) == 0:
            raise ValueError(f"There is no data for {deviceType} between the dates {startTime} and {endTime}")

        if withMetadata:
            devicemetadata = self.entitiesTable()
            if len(devicemetadata) > 0:
                data = data.reset_index().merge(devicemetadata, left_on="deviceName", right_on="entityName").set_index("timestamp")

        return data


    def __init__(self, trialSet: TrialSetWithData, metadata: dict, experimentData: dataEngineFactory):
        self._experimentData = experimentData
        super().__init__(trialSet, metadata)


class EntityTypeWithData(argosDataObjects.EntityType):

    def _initEntities(self):
        for entity in self._metadata['entities']:
            self[entity['name']] = EntityWithData(entityType=self, metadata=entity,experimentData =self._experimentData)

    def __init__(self, experiment:experimentSetupWithData, metadata: dict, experimentData: dataEngineFactory):
        self._experimentData = experimentData
        super().__init__(experiment, metadata)

    def getData(self, startTime=None, endTime=None):
        return self._experimentData.getData(self.name,startTime = startTime,endTime = endTime )

class EntityWithData(argosDataObjects.Entity):

    def __init__(self, entityType: EntityTypeWithData, metadata: dict, experimentData):
        self._experimentData = experimentData
        super().__init__(entityType, metadata)


    def getData(self,startTime=None, endTime=None):
        return self._experimentData.getData(self.entityType.name,self.name,startTime,endTime)