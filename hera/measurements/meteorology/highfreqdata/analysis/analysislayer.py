import pandas
import dask
from .turbulencestatistics import singlePointTurbulenceStatistics
from .meandatacalculator import AveragingCalculator, MeanDataCalculator


class RawdataAnalysis:
    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, datalayer):
        self._datalayer = datalayer

    def singlePointTurbulenceStatistics(self,
                                        deviceNameOrData,
                                        samplingWindow,
                                        start,
                                        end,
                                        height,
                                        buildingHeight,
                                        averagedHeight,
                                        inmemory=False,
                                        isMissingData=False,
                                        **kwargs):
        """
        This method loads the raw data that corresponds to the requirements (projectName, station, instrument.. ) and
        creates a turbulence calculator with the desirable sampling window.


        Parameters
        ----------
        deviceNameOrData : str / pandas.DataFrame / dask.Dataframe / None
            the data to process.

            if str, queries the database with the deviceName as input.

             if NOne, query the database only with kwargs.

        samplingWindow : str
            The desirable sampling window.

        start : str/pandas.Timestamp
            Datetime of the begin.

        end : str/pandas.Timestamp
            Datetime of the end.

        inmemory : bool, positional, default False
            A flag of whether or not to use pandas.

        isMissingData : bool, positional, default False
            A flag if there is a missing data to compute accordingly.

        kwargs :
            Other query arguments for the database.

        Returns
        -------
        singlePointTurbulenceStatistics
            A turbulence calculator of the loaded raw data.
        """

        identifier = {'projectName': self.datalayer.projectName,
                      'samplingWindow': samplingWindow,
                      'height': height,
                      'buildingHeight': buildingHeight,
                      'averagedHeight': averagedHeight,
                      'start': start,
                      'end': end,
                      "isMissingData": isMissingData,
                      "filters": None,
                      "dataSource1": None,
                      "dataSource2": None
                      }
        identifier.update(kwargs)

        if isinstance(deviceNameOrData, pandas.DataFrame) or isinstance(deviceNameOrData, dask.dataframe.DataFrame):

            rawData = deviceNameOrData

        else:
            raise ValueError("deviceNameOrData must be a dask/pandas dataframe")

        return singlePointTurbulenceStatistics(rawData=rawData, metadata=identifier)

    def AveragingCalculator(self,
                            deviceNameOrData,
                            samplingWindow,
                            start,
                            end,
                            height,
                            buildingHeight,
                            averagedHeight,
                            inmemory=False,
                            isMissingData=False,
                            **kwargs):
        """
        This method loads the raw data that corresponds to the requirements (projectName, station, instrument.. ) and
        creates a TRH calculator with the desirable sampling window. It then uses the calculator to return a mean temperature
        pandas dataframe.


        Parameters
        ----------
        deviceNameOrData : str / pandas.DataFrame / dask.Dataframe / None
            the data to process.

            if str, queries the database with the deviceName as input.

             if NOne, query the database only with kwargs.

        samplingWindow : str
            The desirable sampling window.

        start : str/pandas.Timestamp
            Datetime of the begin.

        end : str/pandas.Timestamp
            Datetime of the end.

        inmemory : bool, positional, default False
            A flag of whether or not to use pandas.

        isMissingData : bool, positional, default False
            A flag if there is a missing data to compute accordingly.

        kwargs :
            Other query arguments for the database.

        Returns
        -------
        singlePointTurbulenceStatistics
            A turbulence calculator of the loaded raw data.
        """

        identifier = {'projectName': self.datalayer.projectName,
                      'samplingWindow': samplingWindow,
                      'height': height,
                      'buildingHeight': buildingHeight,
                      'averagedHeight': averagedHeight,
                      'start': start,
                      'end': end,
                      "filters": None,
                      "dataSource1": None,
                      "dataSource2": None
                      }
        identifier.update(kwargs)

        if isinstance(deviceNameOrData, pandas.DataFrame) or isinstance(deviceNameOrData, dask.dataframe.DataFrame):
            rawData = deviceNameOrData
        else:
            raise ValueError("deviceNameOrData must be a dask/pandas dataframe")

        calculator = AveragingCalculator(rawData=rawData, metadata=identifier)

        return calculator

    def MeanDataCalculator(self, TurbCalcOrData=None, compute_mode_turb='not_from_db_and_not_save',
                           AverageCalcOrData=None, compute_mode_AverageCalc=None, **metadata):

        return MeanDataCalculator(TurbCalcOrData, compute_mode_turb, AverageCalcOrData, compute_mode_AverageCalc,
                                  **metadata)
