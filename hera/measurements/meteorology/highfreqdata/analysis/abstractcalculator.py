import dask
import pandas
from copy import deepcopy
from hera.datalayer import collection as datalayer


class AbstractCalculator(object):
    """Base class for high-frequency meteorological data calculators.

    Provides the common infrastructure for computing turbulence statistics
    from raw sonic anemometer data. Manages raw data, temporary computed
    columns, in-memory averaged references, and optional database
    persistence of results.

    Subclasses (``singlePointTurbulenceStatistics``, ``AveragingCalculator``)
    add domain-specific calculation methods that populate ``_TemporaryData``
    and ``_CalculatedParams``.
    """

    _RawData = None
    _metadata = None
    _DataType = None
    _TemporaryData = None
    _CalculatedParams = None
    _AllCalculatedParams = None
    _InMemoryAvgRef = None
    _Karman = 0.4
    _saveProperties = {'dataFormat': None}

    def __init__(self, rawData, metadata):
        """Initialise with raw data and metadata.

        Parameters
        ----------
        rawData : pandas.DataFrame or dask.dataframe.DataFrame
            High-frequency time-series data (u, v, w, T columns expected).
        metadata : dict
            Experiment metadata including ``samplingWindow``, ``projectName``,
            ``height``, ``buildingHeight``, ``averagedHeight``, ``start``,
            ``end``, etc.

        Raises
        ------
        ValueError
            If *rawData* is neither a pandas nor a dask DataFrame.
        """
        if isinstance(rawData,pandas.DataFrame):
            self._DataType = 'pandas'
        elif isinstance(rawData,dask.dataframe.DataFrame):
            self._DataType = 'dask'
        else:
            raise ValueError("'rawData' type must be 'pandas.DataFrame' or 'dask.dataframe.core.DataFrame'.\nGot '%s'." % type(rawData))

        self._RawData = rawData
        self._metadata = metadata
        self._TemporaryData = pandas.DataFrame()
        self._CalculatedParams = []
        self._AllCalculatedParams = []
        self._joinmethod = "left"

    @property
    def JoinMethod(self):
        """str : Join method used when merging computed columns (default ``'left'``)."""
        return self._joinmethod

    @property
    def RawData(self):
        """pandas.DataFrame or dask.dataframe.DataFrame : The original raw data."""
        return self._RawData


    @property
    def TemporaryData(self):
        """pandas.DataFrame : Intermediate computed columns before final aggregation."""
        return self._TemporaryData

    @property
    def metaData(self):
        """dict : Experiment metadata dictionary."""
        return self._metadata

    @property
    def SamplingWindow(self):
        """str : Resampling window from metadata (e.g. ``'30min'``, ``'10min'``)."""
        return self._metadata['samplingWindow']

    @property
    def Karman(self):
        """float : Von Kármán constant (0.4)."""
        return self._Karman

    def set_saveProperties(self, dataFormat, **kwargs):
        """
        Setting the parameters for handling the saving part.

        :param dataFormat: The format to save the data to.
        :param kwargs: Other arguments required for saving the data to the specific dataFormat.
        :return:
        """
        self._saveProperties['dataFormat'] = dataFormat
        self._saveProperties.update(kwargs)

    def _updateInMemoryAvgRef(self, df):
        from .turbulencestatistics import InMemoryAvgData
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = InMemoryAvgData(df, turbulenceCalculator=self)
        else:
            self._InMemoryAvgRef = InMemoryAvgData(pandas.concat([df, self._InMemoryAvgRef], axis=1),
                                                   turbulenceCalculator=self)

        self._CalculatedParams = []

    def _compute(self):
        self._joinmethod = "left"

        if self._DataType == 'dask':
            try:
                df = self._TemporaryData[[x[0] for x in self._CalculatedParams]].compute()
            except ValueError as valueError:
                errorMessage = """A value error occurred while computing the data.
                                                This is probably because one of the problems bellow:
                                                1.The time index of the data is not divisible in the sampling window.
                                                  Please ensure that the sampling window and the time index are matching and try again.
                                                2.Data is missing. Try using the keyword argument \'isMissingData\' (isMissingData=True).

                                                The error that was raised: %s""" % valueError

                raise ValueError(errorMessage)
        else:
            df = self._TemporaryData[[x[0] for x in self._CalculatedParams]]

        self._updateInMemoryAvgRef(df)

    def compute(self, mode='not_from_db_and_not_save'):
        """Execute all pending calculations and return the aggregated result.

        Parameters
        ----------
        mode : str
            Execution mode controlling database interaction:

            - ``'not_from_db_and_not_save'`` — compute locally, do not persist.
            - ``'from_db_and_save'`` — check DB first; compute and save if missing.
            - ``'from_db_and_not_save'`` — check DB first; compute if missing, don't save.
            - ``'not_from_db_and_save'`` — compute locally and save to DB.

        Returns
        -------
        InMemoryAvgData
            Aggregated result containing all computed parameters.

        Raises
        ------
        ValueError
            If no parameters have been calculated yet.
        """
        if self._TemporaryData.columns.empty:
            raise ValueError("Parameters have not been calculated yet.")

        self._AllCalculatedParams.extend(self._CalculatedParams)

        getattr(self, '_compute_%s' % mode)()

        return self._InMemoryAvgRef

    def _compute_from_db_and_save(self):
        query = deepcopy(self.metaData)
        query.pop("projectName")
        query.pop("start")
        query.pop("end")
        docExist = list(
            datalayer.Cache.getDocuments(params__all=self._AllCalculatedParams, start__lte=self.metaData['end'],
                                         end__gte=self.metaData['start'], **query))


        if docExist:
            df = docExist[-1].getDocFromDB(usePandas=True)[[x[0] for x in self._CalculatedParams]]
            self._updateInMemoryAvgRef(df)
        else:
            self._compute()
            self._save_to_db(self._AllCalculatedParams)

    def _compute_from_db_and_not_save(self):
        query = deepcopy(self.metaData)
        query.pop("start")
        query.pop("end")
        query.pop("projectName")

        docExist = list(
            datalayer.Cache.getDocuments(params__all=self._AllCalculatedParams, start__lte=self.metaData['end'],
                                         end__gte=self.metaData['start'], **query))

        if docExist:
            df = docExist[-1].getDocFromDB(usePandas=True)[[x[0] for x in self._CalculatedParams]]
            self._updateInMemoryAvgRef(df)
        else:
            self._compute()

    def _compute_not_from_db_and_save(self):
        self._compute()
        self._save_to_db(self._AllCalculatedParams)

    def _compute_not_from_db_and_not_save(self):
        self._compute()

    def _save_to_db(self, params):
        if self._saveProperties['dataFormat'] is None:
            raise AttributeError('No save properties are set. Please use set_saveProperties function')
        else:
            doc = {}
            desc = deepcopy(self.metaData)
            doc['projectName'] = desc.pop('projectName')
            doc['dataFormat'] = self._saveProperties['dataFormat']
            doc['type'] = 'meteorology'
            doc['desc'] = desc
            doc['desc']['params'] = params
            doc['resource'] = getSaveData(data=self._InMemoryAvgRef, **self._saveProperties)



            datalayer.Cache.addDocument(**doc)


def getSaveData(dataFormat, **kwargs):
    """Dispatch data saving to the appropriate ``SaveDataHandler`` method.

    Parameters
    ----------
    dataFormat : str
        Target format name (``'HDF'``, ``'JSON_pandas'``, ``'parquet'``).
    **kwargs
        Arguments forwarded to the handler (``data``, ``path``, etc.).

    Returns
    -------
    str or dict
        Path or metadata dict returned by the handler.
    """
    return getattr(globals()['SaveDataHandler'], 'getSaveData_%s' % dataFormat)(**kwargs)


class SaveDataHandler(object):
    """Static methods for persisting computed data in various formats."""

    @staticmethod
    def getSaveData_HDF(data, path, key):
        """Save data to HDF5 format.

        Parameters
        ----------
        data : pandas.DataFrame
            Data to save.
        path : str
            File path.
        key : str
            HDF5 group key.

        Returns
        -------
        dict
            ``{'path': path, 'key': key}``.
        """
        data.to_HDF(path, key)
        return dict(path=path,
                    key=key
                    )

    @staticmethod
    def getSaveData_JSON_pandas(data, path=None):
        """Save data as JSON.

        Parameters
        ----------
        data : pandas.DataFrame
            Data to save.
        path : str, optional
            File path. If ``None``, returns JSON string directly.

        Returns
        -------
        str
            JSON string (if *path* is ``None``) or file path.
        """
        if path is None:
            return data.to_json()
        else:
            data.to_json(path)
            return path

    @staticmethod
    def getSaveData_parquet(data, path):
        """Save data as Parquet files.

        Parameters
        ----------
        data : dask.dataframe.DataFrame
            Data to save (repartitioned to 100 MB chunks).
        path : str
            Output directory path.

        Returns
        -------
        str
            The output path.
        """
        data = data.repartition(partition_size='100MB')
        data.to_parquet(path)
        return path

