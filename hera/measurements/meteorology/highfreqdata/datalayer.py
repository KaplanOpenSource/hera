import pandas
import dask
import pydoc
import os
import json
from .analysis.analysislayer import RawdataAnalysis
from hera import datalayer

from ....datalayer import Project
from .... import toolkit

class HighFreqToolKit(toolkit.abstractToolkit):
    """
        Manages the loading and storing of high frequency sonic data

        The data can be in the formats:

        - CampbellBinary
        - TOA5 (not implemented yet)


        TODO:
            Complete the other parsers from the older versions.

    """

    DOCTYPE_STATIONS = 'StationsData'
    DOCTYPE_MEASUREMENTS = 'MeasurementsData'

    def __init__(self, projectName,campaignName=None,FilesDirectory=None):
        """
            Initializes a datalayer for the highfreqdata data.


        Parameters
        ----------

        projectName: str
                The project name
        """
        super().__init__(projectName=projectName,toolkitName="highFreqMeteorology",FilesDirectory=FilesDirectory)
        self.logger.info("Init High frequency data")
        self.campaignName = campaignName
        self._analysis = RawdataAnalysis(self)
       # self._presentation = presenation(self,self.analysis)


        self._np_size = "100MB"


    @property
    def docType(self):
        return f"{self.toolkitName}_HighFreqData"


    def loadData(self,fileNameOrData,parser,saveMode=None,**kwargs):

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


        if self.campaignName is None:
            print('go to experimet loadData')
            pass

        pathToData = os.path.abspath(fileNameOrData)
        className = ".".join(__class__.__module__.split(".")[:-1])
        parserPath = f"{className}.parsers.Parser_{parser}"
        parserCls = pydoc.locate(parserPath)
        parser = parserCls()

        overWrite = kwargs.get("overWrite",False)

        ExperimentMetadata = parser.parse(pathToData=pathToData)
        import pdb
        pdb.set_trace()
        # gets the meta data and adds to the DB.
        # A dictionary with the following:
        #
        #   Stations: { the metadata for the stations } (JSON, or pandas)
        #   Devices: [ A list of the metadata of each device, remember to add the folder name of the parquet  (the resource). ].

        # adds to the db.  (with self.addMeasurementDocument).

        for campaign,campaignData in ExperimentMetadata.items():
            for assetDict in campaignData['Assets']:
                for trial in assetDict['trials']:
                    for entity in assetDict['trials'][trial]['entities']:

                        D=campaignData['Devices']
                        entityDict = list(filter(lambda x: x['deviceName'] == entity, D))[0]

                        path = entityDict['deviceDataPath']

                        if not (os.path.exists(path)):
                            raise FileNotFoundError("Wrong folder path")
                        else:
                            print(f'{path} path to {entity} data from {assetDict["name"]} is O.K')

                        L = self.getMeasurementsDocuments(type=self.DOCTYPE_MEASUREMENTS,deviceName=entity)

                        if not (L) or (L and overWrite):

                            delList = self.deleteMeasurementsDocuments(type=self.DOCTYPE_MEASUREMENTS,
                                                                       deviceName=entityDict['deviceName'])
                            for deldoc in delList:
                                print("delete document")
                                print(json.dumps(deldoc, indent=4, sort_keys=True))

                            print(f'Laoding {entity} data to DB')
                            entityDict['Campaign'] = campaign

                            self.addMeasurementsDocument(type=self.DOCTYPE_MEASUREMENTS,
                                                         dataFormat=datalayer.datatypes.PARQUET,
                                                         resource=path,
                                                         desc=entityDict)
                        else:
                            print(f' {entityDict["deviceName"]} DB list is not empty, but overWrite flag is {overWrite}')
                            print(f'{entityDict["deviceName"]} not stored to DB ')

            # for deviceDict in campaignData['Devices']:

            #     location = deviceDict['trials']['deploy']['deviceLocation']
            #     instType = deviceDict['trials'][0]['deploy']['deviceType']
            #     hgt = deviceDict['trials'][0]['deploy']['deviceHeight']
            #
            #     path = os.path.join(self._dataBasePath, location, instType, hgt)
            #
            #     if not (os.path.exists(path)):
            #         raise FileNotFoundError("Wrong folder path")
            #     else:
            #         print(f'{deviceDict["deviceName"]} path is O.K')
            #
            #     L=self.getMeasurementsDocuments(type = self.DOCTYPE_MEASUREMENTS,
            #                                     deviceName = deviceDict['deviceName'])
            #
            #     if not(L) or (L and overWrite):
            #
            #         delList=self.deleteMeasurementsDocuments(type = self.DOCTYPE_MEASUREMENTS,
            #                                                  deviceName = deviceDict['deviceName'])
            #         for deldoc in delList:
            #             print("delete document")
            #             print(json.dumps(deldoc, indent=4, sort_keys=True))
            #
            #         print(f'Laoding {deviceDict["deviceName"]} to DB')
            #         deviceDict['Campaign']=self.CAMPAIGN
            #
            #         self.addMeasurementsDocument(type        = self.DOCTYPE_MEASUREMENTS,
            #                                      dataFormat  = datalayer.datatypes.PARQUET,
            #                                      resource    = path,
            #                                      desc        = deviceDict)
            #     else:
            #         print(f' {deviceDict["deviceName"]} DB list is not empty, but overWrite flag is {args.overWrite}')
            #         print(f'{deviceDict["deviceName"]} not stored to DB ')

    def getDevicesDesc(self):
        L = self.getMeasurementsDocuments(type=self.DOCTYPE_MEASUREMENTS)

        STATIONS=[]

        for entity in L:
            df= pandas.json_normalize(entity.desc)
            df[['lat','lon']]=df['trials.deploy.deviceCoords'].to_list()
            df=df[['deviceName','trials.name','trials.deploy.deviceLocation','lat','lon','trials.deploy.deviceHeight','Campaign']]
            df=df.rename(columns={"trials.name": "trialName",
                                  "deploy.name":"entity",
                                  "trials.deploy.deviceLocation":"asset",
                                  "trials.deploy.deviceHeight":'deviceHeight'})
            STATIONS.append(df)

        stationsDF=pandas.concat(STATIONS)
        return stationsDF



    def _getRawData(self,stationName,height=None,start=None,end=None,inmemory=False,**kwargs):
        """
            Returns raw data of the requested type.

        Parameters
        ----------

        stationName: str
                The name of the station

        dataType: str
                The type of data to retrieve.
                Either 'rawsonic' or 'temperature'.

        height: float
                The height of the sonic from the mast base.

        start:  pandas.Timestamp or str
                Start time for date range.

        end:  pandas.Timestamp or str
                End time for date range.

        Returns
        -------

        dask.dataframe or pandas.Dataframe.
            The raw data.
        """
        if type(start) is str:
            start = pandas.Timestamp(start)

        if type(end) is str:
            end = pandas.Timestamp(end)

        # should add a document type='SonicData' or something.

        ## THIS SHOULD BE CHANGED ACCORDING TO THE NEW METADATA FORMAT (type, deviceType).
        qry = dict(type=self.dataType,
                   stationName=stationName,
                   height=height)

        kwargs.update(qry)


        docList = self.getMeasurementsDocuments(**kwargs)
        rawData = dask.dataframe.concat([doc.getData() for doc in docList])[start:end]

        return rawData.compute() if inmemory else rawData


    def getSonicData(self,stationName,height=None,start=None,end=None,inmemory=False):
        """
            Implement using the getRawData.
        :param stationName:
        :param height:
        :param start:
        :param end:
        :param inmemory:
        :return:
        """
        RAW=self._getRawData(stationName=stationName,height=height,start=start,end=end,inmemory=inmemory)

        return RAW

        pass



class RawData(Project):

    _dataType = None

    @property
    def dataType(self):
        return self._dataType


    def __init__(self, projectName,dataType,databaseNameList=None,useAll=True): # databaseNameList=None, useAll=False

        super().__init__(projectName=projectName,
                         publicProjectName=f"{dataType}_highfreq",
                         databaseNameList=databaseNameList,
                         useAll=useAll)


    def _getRawData(self,stationName,height=None,start=None,end=None,inmemory=False,**kwargs):
        """
            Returns raw data of the requested type.

        Parameters
        ----------

        stationName: str
                The name of the station

        dataType: str
                The type of data to retrieve.
                Either 'rawsonic' or 'temperature'.

        height: float
                The height of the sonic from the mast base.

        start:  pandas.Timestamp or str
                Start time for date range.

        end:  pandas.Timestamp or str
                End time for date range.

        Returns
        -------

        dask.dataframe or pandas.Dataframe.
            The raw data.
        """
        if type(start) is str:
            start = pandas.Timestamp(start)

        if type(end) is str:
            end = pandas.Timestamp(end)

        # shoudl add a document type='SonicData' or something.
        qry = dict(type=self.dataType,
                   stationName=stationName,
                   height=height)

        kwargs.update(qry)


        docList = self.getMeasurementsDocuments(**kwargs)
        rawData = dask.dataframe.concat([doc.getData() for doc in docList])[start:end]

        return rawData.compute() if inmemory else rawData


class RawSonic(RawData):


    _analysis = None

    @property
    def analysis(self):
        return self._analysis


    def __init__(self, projectName, databaseNameList=None, useAll=False):

        super().__init__(projectName=projectName,
                         publicProjectName="RawSonic",
                         databaseNameList= databaseNameList,
                         useAll = useAll)

        self._analysis = RawdataAnalysis(self)



    def getData(self,stationName,height=None,start=None,end=None,inmemory=False,**kwargs):
        """
            Return the sonic raw data either in dask or loaded to memeory (pandas).

        Parameters
        ----------

        stationName: str
                The name of the station
        height: float
                The height of the sonic from the mast base.

        start:  pandas.Timestamp or str
                Start time for date range.

        end:  pandas.Timestamp or str
                End time for date range.

        kwargs : dict,
            Additional query

        Returns
        -------

        dask.dataframe or pandas.Dataframe.
            The raw data.

        """

        self._getRawData(stationName=stationName,
                         dataType='RawSonic',
                         height=height,
                         start=start,
                         end=end,
                         inmemory=inmemory, **kwargs)


class Temperature(RawData):



    def __init__(self, projectName, databaseNameList=None, useAll=False):

        super().__init__(projectName=projectName,
                         publicProjectName="RawSonic",
                         databaseNameList= databaseNameList,
                         useAll = useAll)


    def getData(self,stationName,height=None,start=None,end=None,inmemory=False,**kwargs):
            """
                Return the temperature raw data either in dask or loaded to memeory (pandas).


                Parameters
                ----------

                stationName: str
                        The name of the station
                height: float
                        The height of the sonic from the mast base.

                start:  pandas.Timestamp or str
                        Start time for date range.

                end:  pandas.Timestamp or str
                        End time for date range.

                kwargs : dict,
                    Additional query

            """
            self._getRawData(stationName=stationName,
                             dataType='Temperature',
                             height=height,
                             start=start,
                             end=end,
                             inmemory=inmemory, **kwargs)


