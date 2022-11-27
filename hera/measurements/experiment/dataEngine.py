from hera import toolkit
import json
from hera import datalayer
import pymongo
import pandas
import dask_mongo
from argos.manager import experimentSetup, DEPLOY, DESIGN

PARQUETHERA = 'parquetDataEngingHera'
PANDASDB = 'pandasDataEngineDB'
DASKDB = 'daskDataEngineDB'

class dataEngineFactory:

    def __init__(self):
        pass

    def getDataEngine(self,projectName, datasourceConfiguration, dataType = PARQUETHERA):

        if dataType == PANDASDB:
            return pandasDataEngineDB(projectName,datasourceConfiguration)
        elif dataType == PARQUETHERA:
            return parquetDataEngineHera(projectName,datasourceConfiguration)
        elif dataType == DASKDB:
            return daskDataEngineDB(projectName,datasourceConfiguration)
        else:
            raise NotImplementedError(f"Hera datalayer {dataType} is not implemeneted yet. ")

class pandasDataEngineDB:

    _mongo_client = None
    _dbConfiguration = None

    @property
    def DBconfiguration(self):
        self._dbConfiguration

    def dbConnect(self):
        # print("mongodb://" + self._DB['login']['username'] + ":" + self._DB['login']['password'] + "@" +
        #      self._DB['login']['ip'] + ":" + self._DB['login']['port'] + "/")
        try:
            self._mongo_client = pymongo.MongoClient(f"mongodb://{self._dbConfiguration['login']['username']}:{self._dbConfiguration['login']['password']}@{self._dbConfiguration['login']['ip']}:{self._dbConfiguration['login']['port']}/")
            try:
                self._mongo_client.server_info()
            except OperationFailure as e:
                return e
        except Exception as e:
            return e

    def __init__(self, projectName,datasourceConfiguration):

        if 'DB' not in datasourceConfiguration:
            raise ValueError(f"The configuration file does not have 'DB' definitions\n Got: {json.dumps(datasourceConfiguration,indent=4)}")

        self._dbConfiguration = datasourceConfiguration['DB']
        self.dbConnect()

    def getDataFromTrial(self, deviceType, trialName, trialSet=None, deviceName=None, withMetadata=True,
                         trialState=DEPLOY, startTime=None, endTime=None):
        """
            Return the device data from the set. Use the default trial set if it is None.

        deviceType: str
                The device type KAIJO,PT1000,GC,NDIR

        trialName: str
                The name of the trial

        trialSet: str [optional]
                The trial set. Use default if None

        deviceName: str [optional]
                Filter according to one specific Device.

        startTime: str [optional]
                Use as start time if exists, otherwise use the trial start Time.

        endTime: str [optional]
                Use as end time if exists, otherwise use the trial end Time.

        withMetadata: bool
                If true adds the following Fields:
                * TimeFromReleaseStart:
                        The time elasped from the release start
                * TimeFromReleaseEnd:
                        The time elasped from the release end
                * longitude
                        The position X of the device.
                            Use the trial state to determine the location
                * latitute
                        The position Y of the device
                            Use the trial state to determine the location

        trialState : DEPLOY or DESIGN

        :return: pandas
                Pandas with the data

        """
        trialSet = self.experimentObj.trialSet if trialSet is None else trialSet

        trial = self._experimentObj.experimentSetup.trialSet[trialSet][trialName]
        startTime = trial.properties['TrialStart'] if startTime is None else startTime
        endTime = trial.properties['TrialEnd'] if endTime is None else endTime

        data = self.getData(deviceType=deviceType,
                            deviceName=deviceName,
                            startTime=startTime,
                            endTime=endTime)

        if len(data) == 0:
            raise ValueError(f"There is no data for {deviceType} between the dates {startTime} and {endTime}")

        if withMetadata:
            devicemetadata = self.experimentObj.experimentSetup.trialSet[trialSet][trialName].entitiesTable(trialState)
            if len(devicemetadata) > 0:
                data = data.reset_index().merge(devicemetadata, left_on="deviceName", right_on="entityName").set_index(
                    "timestamp")

        return data

    def getData(self, deviceType, deviceName=None, startTime=None, endTime=None):

        collectionList = [x['name'] for x in self._mongo_client[self._dbConfiguration['db_name']].list_collections()]
        if deviceType not in collectionList:
            raise ValueError(f"device type {deviceType} not found. Should be one of {','.join(collectionList)}")

        collection = self._mongo_client[self._dbConfiguration['db_name']][deviceType]

        full_qry = {"deviceName": deviceName} if deviceName is not None else {}

        timeFilter = {}
        if startTime is not None:
            if not isinstance(startTime, float):
                try:
                    startTime = pandas.to_datetime(startTime).timestamp() * 1000
                except Exception:
                    raise ValueError(f"{startTime} is not a valid date")

            timeFilter['$gte'] = startTime

        if endTime is not None:
            if not isinstance(endTime, float):
                try:
                    endTime = pandas.to_datetime(endTime).timestamp() * 1000
                except Exception:
                    raise ValueError(f"{endTime} is not a valid date")

            timeFilter['$lte'] = endTime

        if len(timeFilter) > 0:
            full_qry["timestamp"] = timeFilter

        ret = pandas.DataFrame(list(collection.find(full_qry)))
        if not ret.empty:
            # When the ts is saved in UTC time.  use this line:
            # ret.timestamp.apply(lambda x: pandas.to_datetime(x, unit="ms",utc=True).tz_convert("israel"))
            ret.timestamp = ret.timestamp.apply(lambda x: pandas.to_datetime(x, unit="ms",utc=True).tz_convert("israel"))
            ret = ret.set_index("timestamp")

        return ret

    def getDeviceList(self, device):

        devicesList = self._experimentObj.experimentSetup.entityType
        return devicesList[device]

    def getDeviceTable(self, device):

        Table = self._experimentObj.experimentSetup.entityTypeTable
        return Table[device]

class daskDataEngineDB:

    _mongo_client = None
    _dbConfiguration = None

    @property
    def DBconfiguration(self):
        self._dbConfiguration

    @property
    def connectionString(self):
        return f"mongodb://{self._dbConfiguration['login']['username']}:{self._dbConfiguration['login']['password']}@{self._dbConfiguration['login']['ip']}:{self._dbConfiguration['login']['port']}/"

    def __init__(self, projectName, datasourceConfiguration):

        if 'DB' not in datasourceConfiguration:
            raise ValueError(f"The configuration file does not have 'DB' definitions\n Got: {json.dumps(datasourceConfiguration,indent=4)}")

        self._dbConfiguration = datasourceConfiguration['DB']

    def getDataFromTrial(self, deviceType, trialName, trialSet=None, deviceName=None, withMetadata=True,
                         trialState=DEPLOY, startTime=None, endTime=None):
        """
            Return the device data from the set. Use the default trial set if it is None.

        deviceType: str
                The device type KAIJO,PT1000,GC,NDIR

        trialName: str
                The name of the trial

        trialSet: str [optional]
                The trial set. Use default if None

        deviceName: str [optional]
                Filter according to one specific Device.

        startTime: str [optional]
                Use as start time if exists, otherwise use the trial start Time.

        endTime: str [optional]
                Use as end time if exists, otherwise use the trial end Time.

        withMetadata: bool
                If true adds the following Fields:
                * TimeFromReleaseStart:
                        The time elasped from the release start
                * TimeFromReleaseEnd:
                        The time elasped from the release end
                * longitude
                        The position X of the device.
                            Use the trial state to determine the location
                * latitute
                        The position Y of the device
                            Use the trial state to determine the location

        trialState : DEPLOY or DESIGN

        :return: pandas
                Pandas with the data

        """
        trialSet = self.experimentObj.trialSet if trialSet is None else trialSet

        trial = self._experimentObj.experimentSetup.trialSet[trialSet][trialName]
        startTime = trial.properties['TrialStart'] if startTime is None else startTime
        endTime = trial.properties['TrialEnd'] if endTime is None else endTime

        data = self.getData(deviceType=deviceType,
                            deviceName=deviceName,
                            startTime=startTime,
                            endTime=endTime)

        if len(data) == 0:
            raise ValueError(f"There is no data for {deviceType} between the dates {startTime} and {endTime}")

        if withMetadata:
            devicemetadata = self.experimentObj.experimentSetup.trialSet[trialSet][trialName].entitiesTable(trialState)
            if len(devicemetadata) > 0:
                data = data.reset_index().merge(devicemetadata, left_on="deviceName", right_on="entityName").set_index(
                    "timestamp")

        return data

    def getData(self, deviceType, deviceName=None, startTime=None, endTime=None):

        # collectionList = [x['name'] for x in self._mongo_client[self._dbConfiguration['db_name']].list_collections()]
        # if deviceType not in collectionList:
        #     raise ValueError(f"device type {deviceType} not found. Should be one of {','.join(collectionList)}")

        full_qry = {"deviceName": deviceName} if deviceName is not None else {}

        timeFilter = {}
        if startTime is not None:
            if not isinstance(startTime, float):
                try:
                    startTime = pandas.to_datetime(startTime).timestamp() * 1000
                except Exception:
                    raise ValueError(f"{startTime} is not a valid date")

            timeFilter['$gte'] = startTime

        if endTime is not None:
            if not isinstance(endTime, float):
                try:
                    endTime = pandas.to_datetime(endTime).timestamp() * 1000
                except Exception:
                    raise ValueError(f"{endTime} is not a valid date")

            timeFilter['$lte'] = endTime

        if len(timeFilter) > 0:
            full_qry["timestamp"] = timeFilter

        try:
            ret = dask_mongo.read_mongo(database="expradmin",
                                        collection=deviceType,
                                        connection_kwargs=dict(host=self.connectionString),
                                        chunksize=10,
                                        match=full_qry
                                        )
        except StopIteration:
            ret = None

        return ret

    def getDeviceList(self, device):

        devicesList = self._experimentObj.experimentSetup.entityType
        return devicesList[device]

    def getDeviceTable(self, device):

        Table = self._experimentObj.experimentSetup.entityTypeTable
        return Table[device]

class parquetDataEngineHera(datalayer.Project):

    def __init__(self, projectName, datasourceConfiguration):
        super().__init__(projectName=projectName)

    def getDataFromTrial(self, deviceType, trialName, trialSet=None, deviceName=None, withMetadata=True,
                         trialState=DEPLOY, startTime=None, endTime=None):
        """
            Return the device data from the set. Use the default trial set if it is None.

        deviceType: str
                The device type KAIJO,PT1000,GC,NDIR

        trialName: str
                The name of the trial

        trialSet: str [optional]
                The trial set. Use default if None

        deviceName: str [optional]
                Filter according to one specific Device.

        startTime: str [optional]
                Use as start time if exists, otherwise use the trial start Time.

        endTime: str [optional]
                Use as end time if exists, otherwise use the trial end Time.

        withMetadata: bool
                If true adds the following Fields:
                * TimeFromReleaseStart:
                        The time elasped from the release start
                * TimeFromReleaseEnd:
                        The time elasped from the release end
                * longitude
                        The position X of the device.
                            Use the trial state to determine the location
                * latitute
                        The position Y of the device
                            Use the trial state to determine the location

        trialState : DEPLOY or DESIGN

        :return: pandas
                Pandas with the data

        """
        trialSet = self.experimentObj.trialSet if trialSet is None else trialSet

        trial = self._experimentObj.experimentSetup.trialSet[trialSet][trialName]
        startTime = trial.properties['TrialStart'] if startTime is None else startTime
        endTime = trial.properties['TrialEnd'] if endTime is None else endTime

        data = self.getData(deviceType=deviceType,
                            deviceName=deviceName,
                            startTime=startTime,
                            endTime=endTime)

        if len(data) == 0:
            raise ValueError(f"There is no data for {deviceType} between the dates {startTime} and {endTime}")

        if withMetadata:
            devicemetadata = self.experimentObj.experimentSetup.trialSet[trialSet][trialName].entitiesTable(
                trialState)
            if len(devicemetadata) > 0:
                data = data.reset_index().merge(devicemetadata, left_on="deviceName",
                                                right_on="entityName").set_index(
                    "timestamp")

        return data

    def getData(self, deviceType, deviceName=None, startTime=None, endTime=None):

        collection = self.getMeasurementsDocuments(type = 'rawData',deviceType = deviceType)
        if len(collection) == 0:
            return pandas.DataFrame()

        data = collection[0].getData()

        if deviceName is not None:
            data = data.query(f"deviceName == '{deviceName}'")

        # ### This is just because the current pandas are not TZ-aware.
        # ### Once they are, we will remove this!!!
        # ### This removes the tz!.
        # if startTime is not None:
        #     if isinstance(startTime,str):
        #         startTime = pandas.to_datetime(startTime)
        #     startTime = startTime.strftime("%Y-%m-%d %H:%M:%S")
        #
        #
        # if endTime is not None:
        #     if isinstance(endTime,str):
        #         endTime = pandas.to_datetime(endTime)
        #     endTime = endTime.strftime("%Y-%m-%d %H:%M:%S")
        # ## ------------------------------------------------------

        data = data.loc[slice(startTime,endTime)]
        return data.compute()

