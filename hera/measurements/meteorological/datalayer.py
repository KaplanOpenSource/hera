import os
import json
import pandas
from .parserClasses import *
from ...datalayer import Measurements
from ...datalayer.document.metadataDocument import nonDBMetadata

class DataLayer(object):
    _DataSource = None

    def __init__(self, DataSource):
        self._DataSource = DataSource

    def getDocFromDB(self, projectName, resource=None, dataFormat=None, **desc):
        desc['DataSource'] = self._DataSource
        return Measurements.getDocuments(projectName=projectName,
                                         resource=resource,
                                         dataFormat=dataFormat,
                                         type='meteorological',
                                         **desc
                                         )

    def getDocFromFile(self):
        pass

    def extractTrnasformLoad(self):
        pass

class DataLayer_IMS(DataLayer):

    def __init__(self):
        super().__init__(DataSource='IMS')

    def getDocFromDB(self, projectName, resource=None, dataFormat=None, StationName=None, **kwargs):

        """
        This function returns a list of 'doc' objects from the database that matches the requested query

        parameters
        ----------
        projectName : String
            The project to which the data is associated
        resource : String/dict/JSON
            The resource of the data
        dataFormat: String
            The data format
        StationName : String
            The name of the requested station. default None
        kwargs : dict
            Other properties for query

        returns
        -------
        docList : List

        """

        desc = kwargs.copy()
        desc['StationName'] = StationName

        docList = super().getDocFromDB(projectName=projectName,
                                       resource=resource,
                                       dataFormat=dataFormat,
                                       desc=desc
                                       )
        return docList

    def getDocFromFile(self,path, time_coloumn='time_obs', **kwargs):

        """
        Reads data from file and returns a 'metadata like' object

        parameters
        ----------

        path : The path to the data file
        time_coloumn : The name of the Time column for indexing. default ‘time_obs’
        kwargs :

        returns
        -------
        nonDBMetadata : list

        """

        # dl = DataLayer()

        loaded_dask, _ = self._getFromDir(path, time_coloumn)
        return [nonDBMetadata(loaded_dask, **kwargs)]


class InMemoryRawData(pandas.DataFrame):
    _Attrs = None

    def __init__(self, data=None, index=None, columns=None, dtype=None, copy=False):
        super(InMemoryRawData, self).__init__(data=data, index=index, columns=columns, dtype=dtype, copy=copy)
        self._Attrs = {}

    def append(self, other, ignore_index=False, verify_integrity=False):
        ret = super(InMemoryRawData, self).append(other, ignore_index=ignore_index, verify_integrity=verify_integrity)
        ret = InMemoryRawData(ret)
        ret._Attrs = other._Attrs
        ret._Attrs.update(self._Attrs)

        return ret

    @classmethod
    def read_hdf(cls, path_or_buf, key=None, **kwargs):
        ret = InMemoryRawData(pandas.read_hdf(path_or_buf, key, **kwargs))
        path_or_buf = '%s%s' % (path_or_buf.rpartition('.')[0], '.json')

        if os.path.isfile(path_or_buf):
            with open(path_or_buf, 'r') as jsonFile:
                ret._Attrs = json.load(jsonFile)

        return ret

    def to_hdf(self, path_or_buf, key, **kwargs):
        pandasCopy = self.copy()
        path_or_buf = '%s%s' % (path_or_buf.rpartition('.')[0], '.hdf')
        pandasCopy.to_hdf(path_or_buf, key, **kwargs)
        path_or_buf = '%s%s' % (path_or_buf.rpartition('.')[0], '.json')
        attrsToSave = self._Attrs

        if len(self._Attrs) > 0:
            if os.path.isfile(path_or_buf):
                with open(path_or_buf, 'r') as jsonFile:
                    attrsFile = json.load(jsonFile)
                    attrsFile.update(attrsToSave)
                    attrsToSave = attrsFile

            with open(path_or_buf, 'w') as jsonFile:
                json.dump(attrsToSave, jsonFile, indent=4, sort_keys=True)


class InMemoryAvgData(InMemoryRawData):
    _TurbulenceCalculator = None

    def __init__(self, data = None, index = None, columns = None, dtype = None, copy = False, turbulenceCalculator = None):
        super(InMemoryAvgData, self).__init__(data = data, index = index, columns = columns, dtype = dtype, copy = copy)
        self._TurbulenceCalculator = turbulenceCalculator
        self._Attrs['samplingWindow'] = turbulenceCalculator.SamplingWindow

    def __getattr__(self, item):
        if self._TurbulenceCalculator is None:
            raise AttributeError("The attribute '_TurbulenceCalculator' is None.")
        elif not item in dir(self._TurbulenceCalculator):
            raise NotImplementedError("The attribute '%s' is not implemented." % item)
        elif item == 'compute':
            ret = getattr(self._TurbulenceCalculator, item)
        else:
            ret = lambda *args, **kwargs: getattr(self._TurbulenceCalculator, item)(inMemory = self, *args, **kwargs)

        return ret
