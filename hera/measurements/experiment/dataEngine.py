import pandas
from ... import datalayer
from ...utils.logging import get_classMethod_logger


PARQUETHERA = 'parquetDataEngingHera'
PANDASDB = 'pandasDataEngineDB'
DASKDB = 'daskDataEngineDB'

class dataEngineFactory:
    """Factory for creating experiment data engine instances."""

    def __init__(self):
        """Initialize the data engine factory."""
        pass

    def getDataEngine(self,projectName, datasourceConfiguration,experimentObj, dataType = PARQUETHERA):
        """Create and return a data engine of the requested type.

        Parameters
        ----------
        projectName : str
            Name of the project.
        datasourceConfiguration : dict
            Data source configuration dictionary.
        experimentObj : object
            The parent experiment object.
        dataType : str, optional
            Engine type identifier (PARQUETHERA, PANDASDB, or DASKDB).

        Returns
        -------
        pandasDataEngineDB or parquetDataEngineHera or daskDataEngineDB
        """

        if dataType == PARQUETHERA:
            return parquetDataEngineHera(projectName,datasourceConfiguration,experimentObj)
        elif dataType in (PANDASDB, DASKDB):
            raise NotImplementedError(f"Hera datalayer {dataType} is not implemeneted yet. ")
        else:
            raise NotImplementedError(f"Hera datalayer {dataType} is not implemeneted yet. ")

class parquetDataEngineHera(datalayer.Project):
    """Data engine that retrieves experiment data from Hera parquet storage."""

    experimentName = None
    experimentObj = None


    def __init__(self, projectName, datasourceConfiguration,experimentObj):
        """Initialize the parquet data engine.

        Parameters
        ----------
        projectName : str
            Name of the project.
        datasourceConfiguration : dict
            Configuration dictionary with experiment name and settings.
        experimentObj : object
            The parent experiment object.
        """
        logger = get_classMethod_logger(self,"parquetDataEngineHera")

        super().__init__(projectName=projectName)

        self.experimentName = datasourceConfiguration['experimentName']
        self.experimentObj = experimentObj

    def getDataFromTrial(self, deviceType, trialName, trialSet=None, deviceName=None, withMetadata=True, startTime=None, endTime=None,**query):
        """
        Return the device data from the set. Use the default trial set if it is None.

        Parameters
        ----------
        deviceType: str
            The device type.

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
                    - TimeFromReleaseStart: The time elasped from the release start.
                    - TimeFromReleaseEnd: The time elasped from the release end.
                    - longitude: The position X of the device. Use the trial state to determine the location
                    - latitute: The position Y of the device. Use the trial state to determine the location


        Returns
        -------
            pd.DataFrame

        """
        #3logger = get_classMethod_logger(self, "parquetDataEngineHera")
        trialSet = self.experimentObj.trialSet if trialSet is None else trialSet

        trial = self.experimentObj.experimentSetup.trialSet[trialSet][trialName]
        startTime = trial.properties['TrialStart'] if startTime is None else startTime
        endTime = trial.properties['TrialEnd'] if endTime is None else endTime

        data = self.getData(deviceType=deviceType,
                            deviceName=deviceName,
                            startTime=startTime,
                            endTime=endTime)

        if len(data) == 0:
            raise ValueError(f"There is no data for {deviceType} between the dates {startTime} and {endTime}")

        if withMetadata:
            devicemetadata = self.experimentObj.experimentSetup.trialSet[trialSet][trialName].entitiesTable()
            if len(devicemetadata) > 0:
                data = data.reset_index().merge(devicemetadata, left_on="deviceName",
                                                right_on="entityName").set_index(
                    "timestamp")

        return data

    def getData(self, deviceType, deviceName=None, startTime=None, endTime=None,autoCompute=False,perDevice=False,**query):
        """
        Returns the data from of the device type. Queries on device if it exists.

        Parameters
        ----------
        deviceType : str
            The device type.

        deviceName : str, default=None
            The device specific name.

        startTime : datetime, default=None
            The requested starting time to bring the data.

        endTime : datetime, default=None
            The requested ending time to bring the data.

        autoCompute  : bool, default=False
            If true, compute and return the pandas. Else return dask.

        perDevice: bool, default=False
            Is the data organized per device. If false, assumes all data from same type is in the same file. If true, device name should also be defined.


        Returns
        -------
            dask.DataFrame, dask.Pandas. 
        """
        logger = get_classMethod_logger(self, "getData")
        logger.execution("------- Start --------")
        logger.debug(f"Getting {deviceType} with device name {deviceName} from {startTime} to {endTime}. Autocompute? {autoCompute}")

        if perDevice:
            assert deviceName, "If perDeivce=True then deviceName should be defined!"

            collection = self.getMeasurementsDocuments(type='Experiment_rawData', experimentName=self.experimentName,
                                                        deviceType=deviceType, deviceName=deviceName, **query)
        else:
            collection = self.getMeasurementsDocuments(type = 'Experiment_rawData',experimentName=self.experimentName,deviceType = deviceType,**query)

        if len(collection) == 0:
                return pandas.DataFrame()

        data = collection[0].getData()

        if deviceName is not None and not perDevice:
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

        if startTime is not None or endTime is not None:
            data = data.loc[slice(startTime,endTime)]

        if autoCompute:
            data =data.compute()

        return data

