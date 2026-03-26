import pandas
import dask
from .turbulencestatistics import singlePointTurbulenceStatistics
from .meandatacalculator import AveragingCalculator, MeanDataCalculator


class RawdataAnalysis:
    """Analysis layer for the high-frequency meteorology toolkit.

    Provides factory methods for creating turbulence calculators,
    averaging calculators, and mean-data calculators from raw sonic
    anemometer data.
    """

    _datalayer = None

    @property
    def datalayer(self):
        """HighFreqToolKit : The parent toolkit instance."""
        return self._datalayer

    def __init__(self, datalayer):
        """Initialise with a reference to the parent toolkit.

        Parameters
        ----------
        datalayer : HighFreqToolKit
            The high-frequency meteorology toolkit.
        """
        self._datalayer = datalayer

    def singlePointTurbulenceStatistics(self,
                                        sonicData,
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
        sonicData : str / pandas.DataFrame / dask.Dataframe / None
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

        if isinstance(sonicData, pandas.DataFrame) or isinstance(sonicData, dask.dataframe.DataFrame):
         rawData = sonicData
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
        """Create a MeanDataCalculator for deriving mean-field turbulence statistics.

        Parameters
        ----------
        TurbCalcOrData : singlePointTurbulenceStatistics or pandas.DataFrame or dask.DataFrame, optional
            Turbulence calculator or pre-computed second-moment data.
        compute_mode_turb : str
            Compute mode for the turbulence calculator (default ``'not_from_db_and_not_save'``).
        AverageCalcOrData : AveragingCalculator or pandas.DataFrame or dask.DataFrame, optional
            Averaging calculator or pre-computed mean data.
        compute_mode_AverageCalc : str, optional
            Compute mode for the averaging calculator. Defaults to *compute_mode_turb*.
        **metadata
            Additional metadata forwarded to the calculator (``start``, ``end``, etc.).

        Returns
        -------
        MeanDataCalculator
            A calculator pre-populated with second moments and optional mean fields.
        """
        return MeanDataCalculator(TurbCalcOrData, compute_mode_turb, AverageCalcOrData, compute_mode_AverageCalc,
                                  **metadata)
