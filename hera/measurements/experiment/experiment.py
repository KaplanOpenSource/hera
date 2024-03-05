import json
import os
import pydoc
import sys
from ... import toolkit,toolkitHome,datalayer
try:
    from argos.experimentSetup import dataObjects as argosDataObjects
    from argos import DESIGN,DEPLOY
except ImportError:
    print("Must have argos installed and in the path. ")

from .dataEngine import dataEngineFactory, PARQUETHERA, PANDASDB,DASKDB
from hera.utils.jsonutils import loadJSON


class experimentHome(toolkit.abstractToolkit):
    """
        This is the object that function as a factory/home to the other experiments.
        It is responsible for getting the right toolkit for the requested experiment.

    """

    DOCTYPE_ENTITIES = 'EntitiesData'
    CODE_DIRECTORY = 'code'

    def __init__(self, projectName, filesDirectory=None):

        super().__init__(projectName=projectName, toolkitName="experimentToolKit", filesDirectory=filesDirectory)
        self.logger.info("Init experiment toolkit")

    #
    def getExperimentsMap(self):
        """

            Returns a map of the names

            experiment name -> experiment data dict.


        :return:
        """
        M=dict()
        for experiment in self.getDataSourceMap():
            experimentName=experiment['datasourceName']
            M[experimentName]=experiment

        return M

    def getExperimentsTable(self):
        return self.getDataSourceTable()

    def getExperiment(self,
                      experimentName,
                      filesDirectory=None):
        """
        get the experiment data source.
        get the experiemt path from data source
        add it to the python path (sys.path.append)
        get the handler class name
        get the handler with pydoc.

        Parameters
        ----------
        experimentName : str
                The name of the experiment

        experimentDataType : str
                Can be either EXPERIMENTDATALAYER_HERA or EXPERIMENTDATALAYER_DB

        dataSourceConfiguration: dict
                overwrite the dataSourceConfiguration of the experiment.
                see ... for details on its structure.

        filesDirectory: str
                The directory to save the cache/intermediate files.
                If None, use the [current directory]/experimentCache.

        defaultTrialSetName:str
                The default trialSetName to use. (in procedures that require trialSetName).

        Returns
        -------
            Instance of the hera.measurements.old.experiment.experiment.experimentSetupWithData
            that was derived for the specific experiment.
        """
        self.logger.info(f"Getting experiment {experimentName}")
        L = self.getDatasourceDocument(datasourceName=experimentName)
        if L:
            self.logger.info(f"Found experiment. Loading")
            experimentPath=L.desc['experimentPath']
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
        return [x for x in self.getExperimentsMap()]

    # def parserList(self):
    #     """
    #         Return the list of parsers.
    #
    #     Returns
    #     -------
    #         list of str
    #     """
    #     className = ".".join(__class__.__module__.split(".")[:-1])
    #     parserPath = f"{className}.parsers"
    #     mod = pydoc.locate(parserPath)
    #     return [x.split("_")[1] for x in dir(mod) if x.startswith("Parser_")]

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
        self._experimentData = dataEngineFactory().getDataEngine(projectName,self._configuration, dataType = dataType)
        self.entityType = dict()
        self.trialSet = dict()

        if filesDirectory is None:
            filesDirectory = os.getcwd()

        cacheDir = os.path.join(filesDirectory, "experimentCache")
        os.makedirs(cacheDir,exist_ok=True)

        argosDataObjects.ExperimentZipFile.__init__(self,setupFile)
        toolkit.abstractToolkit.__init__(self,projectName=projectName,toolkitName=f"{experimentName}Toolkit",filesDirectory=cacheDir)

        self._sonicHighFreqToolkit = toolkitHome.getToolkit(toolkitName=toolkitHome.METEOROLOGY_HIGHFREQ,
                                                            projectName=projectName,
                                                            filesDirectory=cacheDir)

        self._defaultTrialSetName = defaultTrialSetName

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

    def getDataFromTrial(self,deviceType,deviceName = None,trialState = DEPLOY,startTime = None, endTime = None,withMetadata = True):

        startTime = self.properties['TrialStart'] if startTime is None else startTime
        endTime = self.properties['TrialEnd'] if endTime is None else endTime

        data = self._experimentData.getData(deviceType=deviceType,deviceName=deviceName,startTime=startTime,endTime=endTime)

        if len(data) == 0:
            raise ValueError(f"There is no data for {deviceType} between the dates {startTime} and {endTime}")

        if withMetadata:
            devicemetadata = self.entitiesTable(trialState)
            if len(devicemetadata) > 0:
                data = data.reset_index().merge(devicemetadata, left_on="deviceName", right_on="entityName").set_index(
                    "timestamp")

        return data

    def updateDataWithVersion(self,olddata,newdata):
        """
            This procedure gets an old data with version <xx> in field 'version'
            and creates a copy with version+1, then updates all the  data from the new data into it (including
            data that did not exist).

            For example, lets assume that the old data is

            deviceName, 

        Parameters
        ----------
        olddata
        newdata

        Returns
        -------

        """

class TrialSetWithData(argosDataObjects.TrialSet):

    def _initTrials(self):
        for trial in self._metadata['trials']:
            self[trial['name']] = TrialWithdata(trialSet=self,metadata=trial, experimentData =self._experimentData )

    def __init__(self, experiment:experimentSetupWithData, TrialSetSetup: dict, experimentData: dataEngineFactory):
        self._experimentData = experimentData
        super().__init__(experiment, TrialSetSetup)


class TrialWithdata(argosDataObjects.Trial):

    def getData(self,deviceType,deviceName = None,trialState = DEPLOY,startTime = None, endTime = None,withMetadata = False):

        startTime = self.properties['TrialStart'] if startTime is None else startTime
        endTime = self.properties['TrialEnd'] if endTime is None else endTime

        data = self._experimentData.getData(deviceType=deviceType,deviceName=deviceName,startTime=startTime,endTime=endTime)

        if len(data) == 0:
            raise ValueError(f"There is no data for {deviceType} between the dates {startTime} and {endTime}")

        if withMetadata:
            devicemetadata = self.entitiesTable(trialState)
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




#def loadData(self, fileNameOrData, parser, saveMode=None, **kwargs):
    #
    #     """
    #
    #     Parameters:
    #     -----------
    #     fileNameOrData: str
    #         path to parse.
    #
    #     parser: str
    #         The name of the parser to use.
    #
    #     saveMode: str
    #         Ingored in this method.
    #     kwargs: additional parameters
    #             overWrite - if true, overWrite the database records, else throws exception if exists.
    #                         default : False.
    #
    #     :return:
    #     """
    #
    #     if parser not in self.parserList():
    #         raise ValueError(f"{parser} is not a valid parser. Must be {','.join(self.parserList())}")
    #
    #     if isinstance(fileNameOrData,str):
    #         pathToData = os.path.abspath(fileNameOrData)
    #         className = ".".join(__class__.__module__.split(".")[:-1])
    #         parserPath = f"{className}.parsers.Parser_{parser}"
    #         parserCls = pydoc.locate(parserPath)
    #         parser = parserCls()
    #
    #         experimentMetadata = parser.parse(pathToData=pathToData)
    #
    #     # elif isinstance(fileNameOrData,pandas.DataFrame) or isinstance(fileNameOrData,dask.dataframe):
    #     #     data = fileNameOrData
    #     else:
    #         raise ValueError("fileNameOrData must be a path to a metadata file (or a file), all other posibilities not impolemented yet")
    #
    #     overWrite = kwargs.get("overWrite", False)
    #
    #     # experimentMetadata = parser.parse(pathToData=pathToData)
    #
    #
    #     # gets the meta data and adds to the DB.
    #     # A dictionary with the following:
    #     #
    #     #   Stations: { the metadata for the stations } (JSON, or pandas)
    #     #   Devices: [ A list of the metadata of each device, remember to add the folder name of the parquet  (the resource). ].
    #
    #     # adds to the db.  (with self.addMeasurementDocument).
    #
    #     for experiment, experimentData in experimentMetadata.items():
    #
    #         D = experimentData['devices']
    #
    #         # Create the data source of the experiment itself.
    #         L=self.getDatasourceDocument(datasourceName=experiment)
    #
    #         if not (L) or (L and overWrite):
    #             delList=self.deleteDataSourceDocuments(datasourceName=experiment)
    #
    #             for delDoc in delList:
    #
    #                 self.logger.info("deleted experiment document:")
    #                 self.logger.info(json.dumps(delDoc, indent=4, sort_keys=True))
    #
    #             self.addDataSource(dataSourceName=experiment,
    #                                resource=experimentData.get('properties',{}),
    #                                dataFormat=datalayer.datatypes.DICT,
    #                                handlerPath=fileNameOrData,
    #                                handlerClass=experimentData['toolkitHandler'])
    #
    #         # create the assets
    #         L = self.getMeasurementsDocuments(type=self.DOCTYPE_ENTITIES,
    #                                           experiment=experiment)
    #
    #         if not (L) or (L and overWrite):
    #
    #             delList = self.deleteMeasurementsDocuments(type=self.DOCTYPE_ENTITIES,
    #                                                        experiment=experiment)
    #             for deldoc in delList:
    #                 self.logger.info("deleted entities document:")
    #                 self.logger.info(json.dumps(deldoc, indent=4, sort_keys=True))
    #
    #             desc = dict()
    #
    #             desc['experiment'] = experiment
    #             self.addMeasurementsDocument(resource=experimentMetadata[experiment]['entities'],
    #                                          type=self.DOCTYPE_ENTITIES,
    #                                          dataFormat=datalayer.datatypes.DICT,
    #                                          desc=desc)
    #
    #         # create the devices (with the data from the stations).
    #
    #         for entityDict in D:
    #             entity=entityDict["deviceName"]
    #             path = entityDict['deviceDataPath']
    #
    #             if not (os.path.exists(path)):
    #                 raise FileNotFoundError("Wrong folder path")
    #             else:
    #
    #                 print(f'{path} path to {entity} data  is O.K')
    #
    #             L = self.getMeasurementsDocuments(type=self.DOCTYPE_DEVICES,
    #                                               deviceName=entity,
    #                                               experiment=experiment)
    #
    #             if not (L) or (L and overWrite):
    #
    #                 delList = self.deleteMeasurementsDocuments(type=self.DOCTYPE_DEVICES,
    #                                                            deviceName=entity,
    #                                                            experiment=experiment)
    #                 for deldoc in delList:
    #                     print("delete devices document")
    #                     print(json.dumps(deldoc, indent=4, sort_keys=True))
    #
    #                 print(f'Loading {entity} data to DB')
    #                 entityDict['experiment'] = experiment
    #
    #                 self.addMeasurementsDocument(type=self.DOCTYPE_DEVICES,
    #                                              dataFormat=datalayer.datatypes.PARQUET,
    #                                              resource=path,
    #                                              desc=entityDict)
    #             else:
    #                 print(f'{entityDict["deviceName"]} DB list is not empty, but overWrite flag is {overWrite}')
    #                 print(f'{entityDict["deviceName"]} not stored to DB ')
