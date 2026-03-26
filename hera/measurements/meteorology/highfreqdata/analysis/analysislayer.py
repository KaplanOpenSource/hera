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

    def _resolveData(self, dataOrName):
        """Resolve a data source name or DataFrame to a DataFrame.

        Parameters
        ----------
        dataOrName : str or pandas.DataFrame or dask.dataframe.DataFrame
            If a string, loads the named data source from the toolkit.
            If a DataFrame, returns it directly.

        Returns
        -------
        pandas.DataFrame or dask.dataframe.DataFrame

        Raises
        ------
        ValueError
            If the data source name is not found or the type is unsupported.
        """
        if isinstance(dataOrName, (pandas.DataFrame, dask.dataframe.DataFrame)):
            return dataOrName
        if isinstance(dataOrName, str):
            data = self.datalayer.getDataSourceData(dataOrName)
            if data is None:
                raise ValueError(
                    f"Data source '{dataOrName}' not found in project "
                    f"'{self.datalayer.projectName}'."
                )
            return data
        raise ValueError(
            f"Expected a data source name (str) or DataFrame, got {type(dataOrName)}."
        )

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
        """Create a single-point turbulence statistics calculator.

        Parameters
        ----------
        sonicData : str or pandas.DataFrame or dask.dataframe.DataFrame
            Sonic anemometer data. If a string, the named data source is
            loaded from the project automatically.
        samplingWindow : str
            Resampling window (e.g. ``'30min'``, ``'10S'``).
        start : str or pandas.Timestamp
            Start of the analysis period.
        end : str or pandas.Timestamp
            End of the analysis period.
        height : float
            Instrument height above ground (m).
        buildingHeight : float
            Height of the building the instrument is mounted on (m).
        averagedHeight : float
            Area-averaged building height in the surroundings (m).
        inmemory : bool, optional
            Whether to use pandas (in-memory) mode. Default ``False``.
        isMissingData : bool, optional
            Set ``True`` if the data has gaps. Default ``False``.
        **kwargs
            Additional metadata stored in the calculator's identifier.

        Returns
        -------
        singlePointTurbulenceStatistics
            A turbulence calculator ready for method-chained analysis.
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

        rawData = self._resolveData(sonicData)

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
        """Create an averaging calculator for time-windowed mean fields.

        Parameters
        ----------
        deviceNameOrData : str or pandas.DataFrame or dask.dataframe.DataFrame
            Raw data or data source name. If a string, the named data
            source is loaded from the project automatically.
        samplingWindow : str
            Resampling window (e.g. ``'30min'``, ``'5min'``).
        start : str or pandas.Timestamp
            Start of the analysis period.
        end : str or pandas.Timestamp
            End of the analysis period.
        height : float
            Instrument height above ground (m).
        buildingHeight : float
            Height of the building the instrument is mounted on (m).
        averagedHeight : float
            Area-averaged building height in the surroundings (m).
        inmemory : bool, optional
            Whether to use pandas (in-memory) mode. Default ``False``.
        isMissingData : bool, optional
            Set ``True`` if the data has gaps. Default ``False``.
        **kwargs
            Additional metadata stored in the calculator's identifier.

        Returns
        -------
        AveragingCalculator
            A calculator with pre-computed windowed means.
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

        rawData = self._resolveData(deviceNameOrData)

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
