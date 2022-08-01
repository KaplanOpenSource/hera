from hera import toolkit
import json
from hera import datalayer
import pymongo
import os
import pydoc
import pandas
import dask
import sys
import argos.experimentSetup
import logging
from argos.manager import experimentSetup, DEPLOY, DESIGN

EXPERIMENTDATALAYER_HERA = 'parquet'
EXPERIMENTDATALAYER_DB = 'MongoDB'


class dataEngine():

    def __init__(self, toolkit, dataType = EXPERIMENTDATALAYER_HERA):

        if dataType not in [EXPERIMENTDATALAYER_DB, EXPERIMENTDATALAYER_HERA]:
            raise ValueError(
                f"datalayer type {dataType} must be {EXPERIMENTDATALAYER_DB} or {EXPERIMENTDATALAYER_HERA}")

        if dataType == EXPERIMENTDATALAYER_DB:
            return experimentDataLayerDB(toolkit,toolkit.configuration())
        elif dataType == EXPERIMENTDATALAYER_HERA:
            return experimentDataLayerParqet(toolkit)

        else:
            raise NotImplementedError("Hera datalayer is not implemeneted yet. ")

class experimentDataLayerDB():

    _mongo_client = None
    _DB = None
    _toolkit = None

    @property
    def FilesDirectory(self):
        return self.toolkit.FilesDirectory

    @property
    def projectName(self):
        return self.toolkit.projectName

    @property
    def configurationDB(self):
        return self.toolkit.configuration()["DB"]

    @property
    def toolkit(self):
        return self._toolkit


    def dbConnect(self):

        self._DB = self.configurationDB
        print("Connecting to mongoDB")
        # print("mongodb://" + self._DB['login']['username'] + ":" + self._DB['login']['password'] + "@" +
        #      self._DB['login']['ip'] + ":" + self._DB['login']['port'] + "/")
        try:
            self._mongo_client = pymongo.MongoClient(
                f"mongodb://{self._DB['login']['username']}:{self._DB['login']['password']}@{self._DB['login']['ip']}:{self._DB['login']['port']}/")
            try:
                self._mongo_client.server_info()
                return "Connection Successfull"
            except OperationFailure as e:
                return e
        except Exception as e:
            return e

        print(dbConnect())


    def __init__(self, toolkit, configuration):

        self._toolkit = toolkit
        # self._analysis = experimentDataAnalysis(self)
        # self._setupPresentation = setupPresentationLayer(self, self._analysis)
        # self._dataPresentationLayer = dataPresetationLayer(self, self._technicalAnalysis)

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
        trialSet = self.toolkit.trialSet if trialSet is None else trialSet

        trial = self._toolkit.experimentSetup.trialSet[trialSet][trialName]
        startTime = trial.properties['TrialStart'] if startTime is None else startTime
        endTime = trial.properties['TrialEnd'] if endTime is None else endTime

        data = self.getData(deviceType=deviceType,
                            deviceName=deviceName,
                            startTime=startTime,
                            endTime=endTime)

        if len(data) == 0:
            raise ValueError(f"There is no data for {deviceType} between the dates {startTime} and {endTime}")

        if withMetadata:
            devicemetadata = self.toolkit.experimentSetup.trialSet[trialSet][trialName].entitiesTable(trialState)
            if len(devicemetadata) > 0:
                data = data.reset_index().merge(devicemetadata, left_on="deviceName", right_on="entityName").set_index(
                    "timestamp")

        return data

    def getData(self, deviceType, deviceName=None, startTime=None, endTime=None):

        collectionList = [x['name'] for x in self._mongo_client[self._DB['db_name']].list_collections()]
        if deviceType not in collectionList:
            raise ValueError(f"device type {deviceType} not found. Should be one of {','.join(collectionList)}")

        collection = self._mongo_client[self._DB['db_name']][deviceType]

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
            ret.timestamp = ret.timestamp.apply(lambda x: pandas.to_datetime(x, unit="ms"))
            ret = ret.set_index("timestamp")

        return ret

    def getDeviceList(self, device):

        devicesList = self._toolkit.experimentSetup.entityType
        return devicesList[device]

    def getDeviceTable(self, device):

        Table = self._toolkit.experimentSetup.entityTypeTable
        return Table[device]

class experimentDataLayerParqet(datalayer.Project):
    pass

class analysis():

    pass


class presentation():

    pass
