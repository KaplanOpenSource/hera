import dask.dataframe.core
import pandas

from hera.measurements.meteorological.highfreqdata import InMemoryAvgData
from .....  import datalayer


class AbstractCalculator(object):
    _RawData = None
    _Metadata = None
    _DataType = None
    _TemporaryData = None
    _Identifier = None
    _CalculatedParams = None
    _AllCalculatedParams = None
    _InMemoryAvgRef = None
    _Karman = 0.4
    _saveProperties = {'dataFormat': None}

    def __init__(self, rawData, metadata, identifier):
        if type(rawData) == pandas.DataFrame:
            self._DataType = 'pandas'
        elif type(rawData) == dask.dataframe.core.DataFrame:
            self._DataType = 'dask'
        else:
            raise ValueError("'rawData' type must be 'pandas.DataFrame' or 'dask.dataframe.core.DataFrame'.\nGot '%s'." % type(rawData))

        self._RawData = rawData
        self._Metadata = metadata
        self._TemporaryData = pandas.DataFrame()
        self._Identifier = identifier
        self._CalculatedParams = []
        self._AllCalculatedParams = []
        self._joinmethod = "left"

    @property
    def JoinMethod(self):
        return self._joinmethod

    @property
    def RawData(self):
        return self._RawData

    @property
    def Metadata(self):
        return self._Metadata

    @property
    def TemporaryData(self):
        return self._TemporaryData

    @property
    def Identifier(self):
        return self._Identifier

    @property
    def SamplingWindow(self):
        return self._Identifier['samplingWindow']

    @property
    def Karman(self):
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
        if self._TemporaryData.columns.empty:
            raise ValueError("Parameters have not been calculated yet.")

        self._AllCalculatedParams.extend(self._CalculatedParams)

        getattr(self, '_compute_%s' % mode)()

        return self._InMemoryAvgRef

    def _compute_from_db_and_save(self):
        query = self._query()
        docExist = list(
            datalayer.Cache.getDocuments(params__all=self._AllCalculatedParams, start__lt=self.Identifier['end'],
                                         end__gt=self.Identifier['start'], **query))


        if docExist:
            df = docExist[-1].getDocFromDB(usePandas=True)[[x[0] for x in self._CalculatedParams]]
            self._updateInMemoryAvgRef(df)
        else:
            self._compute()
            self._save_to_db(self._AllCalculatedParams, query)

    def _compute_from_db_and_not_save(self):
        query = self._query()

        docExist = list(
            datalayer.Cache.getDocuments(params__all=self._AllCalculatedParams, start__lt=self.Identifier['end'],
                                         end__gt=self.Identifier['start'], **query))

        if docExist:
            df = docExist[-1].getDocFromDB(usePandas=True)[[x[0] for x in self._CalculatedParams]]
            self._updateInMemoryAvgRef(df)
        else:
            self._compute()

    def _compute_not_from_db_and_save(self):
        query = self._query()
        self._compute()
        self._save_to_db(self._AllCalculatedParams, query)

    def _compute_not_from_db_and_not_save(self):
        self._compute()

    def _query(self):
        query = dict(projectName=self.Identifier['projectName'],
                     start=self.Identifier['start'],
                     end=self.Identifier['end'],
                     samplingWindow=self.SamplingWindow,
                     station=self.Identifier['station'],
                     instrument=self.Identifier['instrument'],
                     height=self.Identifier['height']
                     )
        return query

    def _save_to_db(self, params, query):
        if self._saveProperties['dataFormat'] is None:
            raise AttributeError('No save properties are set. Please use set_saveProperties function')
        else:
            doc = {}
            doc['projectName'] = query.pop('projectName')
            doc['dataFormat'] = self._saveProperties['dataFormat']
            doc['type'] = 'meteorological'
            doc['desc'] = query
            doc['desc']['start'] = self.Identifier['start']
            doc['desc']['end'] = self.Identifier['end']
            doc['desc']['samplingWindow'] = self.SamplingWindow
            doc['desc']['params'] = params
            doc['resource'] = getSaveData(data=self._InMemoryAvgRef, **self._saveProperties)
            datalayer.Cache.addDocument(**doc)


def getSaveData(dataFormat, **kwargs):
    return getattr(globals()['SaveDataHandler'], 'getSaveData_%s' % dataFormat)(**kwargs)


class SaveDataHandler(object):

    @staticmethod
    def getSaveData_HDF(data, path, key):
        data.to_HDF(path, key)
        return dict(path=path,
                    key=key
                    )

    @staticmethod
    def getSaveData_JSON_pandas(data, path=None):
        if path is None:
            return data.to_json()
        else:
            data.to_json(path)
            return path

    @staticmethod
    def getSaveData_parquet(data, path):
        data = data.repartition(partition_size='100MB')
        data.to_parquet(path)
        return path

