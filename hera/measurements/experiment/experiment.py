from ... import toolkit
import json
from hera import datalayer
import os
import pydoc
import pandas
import dask
import sys



class experimentToolKit(toolkit.abstractToolkit):

    DOCTYPE_ASSETS = 'AssetsData'
    DOCTYPE_DEVICES = 'DevicesData'

    def __init__(self,projectName,FilesDirectory=None):

        super().__init__(projectName=projectName, toolkitName="experimentData",FilesDirectory=FilesDirectory)
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
            L = self.getMeasurementsDocuments(type=self.DOCTYPE_ASSETS,
                                              experiment=experiment)

            if not (L) or (L and overWrite):

                delList = self.deleteMeasurementsDocuments(type=self.DOCTYPE_ASSETS,
                                                           experiment=experiment)
                for deldoc in delList:
                    self.logger.info("deleted assets document:")
                    self.logger.info(json.dumps(deldoc, indent=4, sort_keys=True))

                desc = dict()

                desc['experiment'] = experiment
                self.addMeasurementsDocument(resource=experimentMetadata[experiment]['assets'],
                                             type=self.DOCTYPE_ASSETS,
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

    def getExperiment(self,experimentName):

        try:
            L = self.getDatasourceDocument(datasourceName=experimentName)
            path=L.desc['handlerPath']
            sys.path.append(path)

            toolkit=L.desc['handlerClass']
            toolkitCls=pydoc.locate(toolkit)
            return toolkitCls(self.projectName)

        except ImportError as e:
            print(f"Cannot load the specific experiment handler, falling to the default one. The exception is {e}")
            return experimentDataLayer(self.projectName, experimentName=experimentName)


        # get the experiment data source.
        # get the experiemt path from data source
        # add it to the python path (sys.path.append)
        # get the handler class name
        # get the handler with pydoc.

        # if anything fails:
        #       return the experimentDataLayer

        # return experimentDataLayer(self.projectName, experimentName=experimentName)

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


class experimentDataLayer(datalayer.Project):

    DOCTYPE_ASSETS = 'AssetsData'
    DOCTYPE_DEVICES = 'DevicesData'


    def __init__(self, projectName, experimentName):

        super().__init__(projectName=projectName)
        self.logger.info("Init experiment data")
        self.experimentName = experimentName


    def getAssetsMap(self):
        L = self.getMeasurementsDocuments(type=self.DOCTYPE_ASSETS,
                                          experiment=self.experimentName)
        return L[0].getData()

    def getAssetsTable(self):
        Table=[]
        for Asset,assetdata in self.getAssetsMap().items():
            table=pandas.json_normalize(assetdata)
            Table.append(table)

        return pandas.concat((Table))


    def getDevices(self):

        devicesDict=dict()
        L = self.getMeasurementsDocuments(type=self.DOCTYPE_DEVICES,
                                          experiment=self.experimentName)
        for device in L:
            d=dict(device.desc)
            key=d['deviceName']
            del d['deviceName']
            devicesDict[key]=d

        return devicesDict

        pass

    def getDevicesTable(self):

        Table = []
        for Device, devicedata in self.getDevices().items():
            table = pandas.json_normalize(devicedata)
            Table.append(table)

        return pandas.concat((Table))


    def getDevicesDocument(self,deviceName,**kwargs):
        """

        deviceName:
        :param kwargs:
        :return:
        """
        self.getMeasurementsDocuments()
        pass




class analysis():

    pass


class presentation():

    pass
