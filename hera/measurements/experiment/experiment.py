from hera import toolkit
import json
from hera import datalayer
import os
import pydoc
import pandas
import dask
import sys
import argos.experimentSetup
from argos.experimentSetup import dataObjects as argosDataObjects
from .datalayer import dataEngine, EXPERIMENTDATALAYER_HERA, EXPERIMENTDATALAYER_DB
import logging
from hera.utils.jsonutils import loadJSON

EXPERIMENT_SETUPFILE = '/runtimeExperimentData/experiment.json'
CONFIGURATION = '/runtimeExperimentData/Datasources_Configurations.json'
TOOLKIT_FILE = 'toolkit'

class experimentToolKit(toolkit.abstractToolkit):

    DOCTYPE_ENTITIES = 'EntitiesData'

    _experimentDataType = None
    _experimentSetupSource = None
    _configuration = None

    @property
    def trialSet(self):
        return self._trialSet

    @trialSet.setter
    def trialSet(self, value):
        self._trialSet = value

    def __init__(self, projectName, filesDirectory=None):

        super().__init__(projectName=projectName, toolkitName="experimentToolKit", filesDirectory=filesDirectory)
        self.logger.info("Init experiment toolkit")

    def loadData(self, fileNameOrData, parser, saveMode=None, **kwargs):

        """

        Parameters:
        -----------
        fileNameOrData: str
            path to parse.

        parser: str
            The name of the parser to use.

        saveMode: str
            Ingored in this method.
        kwargs: additional parameters
                overWrite - if true, overWrite the database records, else throws exception if exists.
                            default : False.

        :return:
        """

        if parser not in self.parserList():
            raise ValueError(f"{parser} is not a valid parser. Must be {','.join(self.parserList())}")

        if isinstance(fileNameOrData,str):
            pathToData = os.path.abspath(fileNameOrData)
            className = ".".join(__class__.__module__.split(".")[:-1])
            parserPath = f"{className}.parsers.Parser_{parser}"
            parserCls = pydoc.locate(parserPath)
            parser = parserCls()

            experimentMetadata = parser.parse(pathToData=pathToData)

        # elif isinstance(fileNameOrData,pandas.DataFrame) or isinstance(fileNameOrData,dask.dataframe):
        #     data = fileNameOrData
        else:
            raise ValueError("fileNameOrData must be a path to a metadata file (or a file), all other posibilities not impolemented yet")

        overWrite = kwargs.get("overWrite", False)

        # experimentMetadata = parser.parse(pathToData=pathToData)


        # gets the meta data and adds to the DB.
        # A dictionary with the following:
        #
        #   Stations: { the metadata for the stations } (JSON, or pandas)
        #   Devices: [ A list of the metadata of each device, remember to add the folder name of the parquet  (the resource). ].

        # adds to the db.  (with self.addMeasurementDocument).

        for experiment, experimentData in experimentMetadata.items():

            D = experimentData['devices']

            # Create the data source of the experiment itself.
            L=self.getDatasourceDocument(datasourceName=experiment)

            if not (L) or (L and overWrite):
                delList=self.deleteDataSourceDocuments(datasourceName=experiment)

                for delDoc in delList:

                    self.logger.info("deleted experiment document:")
                    self.logger.info(json.dumps(delDoc, indent=4, sort_keys=True))

                self.addDataSource(dataSourceName=experiment,
                                   resource=experimentData.get('properties',{}),
                                   dataFormat=datalayer.datatypes.DICT,
                                   handlerPath=fileNameOrData,
                                   handlerClass=experimentData['toolkitHandler'])

            # create the assets
            L = self.getMeasurementsDocuments(type=self.DOCTYPE_ENTITIES,
                                              experiment=experiment)

            if not (L) or (L and overWrite):

                delList = self.deleteMeasurementsDocuments(type=self.DOCTYPE_ENTITIES,
                                                           experiment=experiment)
                for deldoc in delList:
                    self.logger.info("deleted entities document:")
                    self.logger.info(json.dumps(deldoc, indent=4, sort_keys=True))

                desc = dict()

                desc['experiment'] = experiment
                self.addMeasurementsDocument(resource=experimentMetadata[experiment]['entities'],
                                             type=self.DOCTYPE_ENTITIES,
                                             dataFormat=datalayer.datatypes.DICT,
                                             desc=desc)

            # create the devices (with the data from the stations).

            for entityDict in D:
                entity=entityDict["deviceName"]
                path = entityDict['deviceDataPath']

                if not (os.path.exists(path)):
                    raise FileNotFoundError("Wrong folder path")
                else:

                    print(f'{path} path to {entity} data  is O.K')

                L = self.getMeasurementsDocuments(type=self.DOCTYPE_DEVICES,
                                                  deviceName=entity,
                                                  experiment=experiment)

                if not (L) or (L and overWrite):

                    delList = self.deleteMeasurementsDocuments(type=self.DOCTYPE_DEVICES,
                                                               deviceName=entity,
                                                               experiment=experiment)
                    for deldoc in delList:
                        print("delete devices document")
                        print(json.dumps(deldoc, indent=4, sort_keys=True))

                    print(f'Loading {entity} data to DB')
                    entityDict['experiment'] = experiment

                    self.addMeasurementsDocument(type=self.DOCTYPE_DEVICES,
                                                 dataFormat=datalayer.datatypes.PARQUET,
                                                 resource=path,
                                                 desc=entityDict)
                else:
                    print(
                        f' {entityDict["deviceName"]} DB list is not empty, but overWrite flag is {overWrite}')
                    print(f'{entityDict["deviceName"]} not stored to DB ')

            # for asset,assetDict in experimentData['assets'].items():
            #     for trial in assetDict['trials']:
            #         for entity in assetDict['trials'][trial]['entities']:
            #
            #             entityDict = list(filter(lambda x: x['deviceName'] == entity, D))[0]
            #
            #             path = entityDict['deviceDataPath']
            #
            #             if not (os.path.exists(path)):
            #                 raise FileNotFoundError("Wrong folder path")
            #             else:
            #                 print(f'{path} path to {entity} data from {assetDict["name"]} is O.K')
            #
            #             L = self.getMeasurementsDocuments(type=self.DOCTYPE_DEVICES, deviceName=entity)
            #
            #             if not (L) or (L and overWrite):
            #
            #                 delList = self.deleteMeasurementsDocuments(type=self.DOCTYPE_DEVICES,
            #                                                            deviceName=entityDict['deviceName'])
            #                 for deldoc in delList:
            #                     print("delete document")
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
            #                 print(
            #                     f' {entityDict["deviceName"]} DB list is not empty, but overWrite flag is {overWrite}')
            #                 print(f'{entityDict["deviceName"]} not stored to DB ')

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


    def getExperiment(self,experimentName, experimentSetupSource = None, experimentDataType = EXPERIMENTDATALAYER_HERA, configuration = None):
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

        experimentSetupSource :

        experimentDataType : str
                Can be either EXPERIMENTDATALAYER_HERA or EXPERIMENTDATALAYER_DB

        Returns
        -------
            Instance of the hera.measurements.experiment.experiment.experimentSetupWithData
            that was derived for the specific experiment.
        """

        if experimentDataType not in [EXPERIMENTDATALAYER_HERA,EXPERIMENTDATALAYER_DB]:
            raise ValueError(f"experimentData type must be EXPERIMENTDATALAYER_HERA ({EXPERIMENTDATALAYER_HERA}) or EXPERIMENTDATALAYER_DB ({EXPERIMENTDATALAYER_DB})")

        self._experimentDataType = experimentDataType

        L = self.getDatasourceDocument(datasourceName=experimentName)
        if L:
            path=L.desc['handlerPath']
            sys.path.append(path)

            if  os.path.isfile(path+'/runtimeExperimentData/'+experimentName+ '/experiment.json'):
                self._experimentSetupSource = path+'/runtimeExperimentData/'+experimentName+ '/experiment.json'
            else:
                raise ValueError(f" The experiment setup file doesn't exist in the target folder")

            if  os.path.isfile( path + CONFIGURATION):
                self._configuration = path + CONFIGURATION
            else:
                raise ValueError(f" The configuration file doesn't exist in the target folder")

            toolkit=TOOLKIT_FILE +'.'+ L.desc['className'] #L.desc['handlerClass']
            toolkitCls=pydoc.locate(toolkit)

            return toolkitCls(self)
        else:
            raise ValueError(f"Please load first the experiment datasource to hera")


    def keys(self):
        return [x for x in self.getExperimentsMap()]

    def parserList(self):
        """
            Return the list of parsers_old.

        Returns
        -------
            list of str
        """
        className = ".".join(__class__.__module__.split(".")[:-1])
        parserPath = f"{className}.parsers"
        mod = pydoc.locate(parserPath)
        return [x.split("_")[1] for x in dir(mod) if x.startswith("Parser_")]

    def __getitem__(self, item):
        return self.getExperiment(item)

    def experimentDataType(self):
        return self._experimentDataType

    def getConfiguration(self):
        return loadJSON(self._configuration)

    def getExperimentSetup(self):
        return loadJSON(self._experimentSetupSource)



class experimentSetupWithData(argosDataObjects.Experiment):
    """
        A class that unifies the argos.experiment setup with the data.
    """

    def _initTrialSets(self):
        experimentSetup = self._toolkit.experimentSetup()
        for trialset in experimentSetup['trialSets']:
            self._trialSetsDict[trialset['name']] = TrialSetWithData(experiment = self, TrialSetSetup=trialset,experimentData= self._experimentData)

    def _initEntitiesTypes(self):
        experimentSetup = self._toolkit.experimentSetup()
        for entityType in experimentSetup['entitiesTypes']:
            self._entitiesTypesDict[entityType['name']] = EntityTypeWithData(experiment=self, metadata = entityType, experimentData= self._experimentData)

    def getExperimentData(self):
        return self._experimentData

    def __init__(self, toolkit): # The setup here is the dict already (choose the source in getExperiment)
        self._toolkit = toolkit
        self._experimentData = None#dataEngine(toolkit,toolkit.experimentDataType())
        super().__init__(toolkit.experimentSetup())
        #self.logger.info("Init experiment data")


    def getData(self):

        pass

class TrialSetWithData(argosDataObjects.TrialSet):

    def _initTrials(self):
        for trial in self._metadata['trials']:
            self[trial['name']] = TrialWithdata(trialSet=self,metadata=trial, experimentData =self._experimentData )

    def __init__(self,experiment:experimentSetupWithData, TrialSetSetup: dict, experimentData: dataEngine):
        self._experimentData = experimentData
        super().__init__(experiment, TrialSetSetup)



class TrialWithdata(argosDataObjects.Trial):

    def getData(self,deviceType = None,deviceName = None,trialState = 'deploy',startTime = None, endTime = None,withMetadata = True):
        return self._experimentData.getDataFromTrial(deviceType = deviceType,trialName = self.name(),trialSet = self.trialSet.name(),deviceName=deviceName,
                                                     trialState=trialState,startTime=startTime,endTime=endTime, withMetadata=withMetadata)


    def __init__(self, trialSet: TrialSetWithData, metadata: dict, experimentData: dataEngine):
        self._experimentData = experimentData
        super().__init__(trialSet, metadata)


class EntityTypeWithData(argosDataObjects.EntityType):

    def _initEntities(self):
        for entity in self._metadata['entities']:
            self[entity['name']] = EntityWithData(entityType=self, metadata=entity,experimentData =self._experimentData)

    def __init__(self, experiment:experimentSetupWithData, metadata: dict, experimentData: dataEngine):
        self._experimentData = experimentData
        super().__init__(experiment, metadata)

    def getData(self, startTime=None, endTime=None):
        return self._experimentData.getData(self.name(),startTime = startTime,endTime = endTime )

class EntityWithData(argosDataObjects.Entity):

    def __init__(self, entityType: EntityTypeWithData, metadata: dict, experimentData: dataEngine):
        self._experimentData = experimentData
        super().__init__(entityType, metadata)


    def getData(self,startTime=None, endTime=None):
        return self._experimentData.getData(self.entityType().name(),self.name(),startTime,endTime)


if __name__ =="__main__":

    tool = experimentSetupWithData('ds','testDS', '/home/hadas/Projects/2022/DS2022/experiment.json',EXPERIMENTDATALAYER_DB)
