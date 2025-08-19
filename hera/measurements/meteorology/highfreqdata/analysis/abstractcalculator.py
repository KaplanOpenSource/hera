import dask
import pandas
from copy import deepcopy
from hera.datalayer import collection as datalayer


class AbstractCalculator(object):
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
        return self._joinmethod

    @property
    def RawData(self):
        return self._RawData


    @property
    def TemporaryData(self):
        return self._TemporaryData

    @property
    def metaData(self):
        return self._metadata

    @property
    def SamplingWindow(self):
        return self._metadata['samplingWindow']

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

