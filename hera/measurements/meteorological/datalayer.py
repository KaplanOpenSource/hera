import os
import glob
import json
import pandas
import dask.dataframe as dd
from .parserClasses import *
from ... import datalayer
from ...datalayer.document.metadataDocument import nonDBMetadata
import warnings

class DataLayer(object):
    _DataSource = None
    _parser = None
    _docType = 'meteorological'
    _np_size = "100Mb"

    def __init__(self, DataSource):
        self._DataSource = DataSource
        self._parser = globals()['Parser_%s' % self._DataSource]

    def getDocFromDB(self, projectName, resource=None, dataFormat=None, **desc):
        desc['DataSource'] = self._DataSource
        docList = datalayer.Measurements.getDocuments(projectName=projectName,
                                                      resource=resource,
                                                      dataFormat=dataFormat,
                                                      type=self._docType,
                                                      **desc
                                                      )
        return docList

    def getDocFromFile(self):
        pass

    def _parse(self, **kwargs):
        return self._parser().parse(**kwargs)

    def load(self):
        pass


class DataLayer_IMS(DataLayer):
    _HebRenameDict = None
    _hebStnRename = None

    def __init__(self):
        super().__init__(DataSource='IMS')

        self._HebRenameDict = {"שם תחנה": 'Station_name',
                               "תאריך": "Date",
                               "שעה- LST": "Time_(LST)",
                               "טמפרטורה(C°)": "Temperature_(°C)",
                               "טמפרטורת מקסימום(C°)": "Maximum_Temperature_(°C)",
                               "טמפרטורת מינימום(C°)": "Minimum_Temperature_(°C)",
                               "טמפרטורה ליד הקרקע(C°)": "Ground_Temperature_(°C)",
                               "לחות יחסית(%)": "Relative_humidity_(%)",
                               "לחץ בגובה התחנה(hPa)": "Pressure_at_station_height_(hPa)",
                               "קרינה גלובלית(W/m²)": "Global_radiation_(W/m²)",
                               "קרינה ישירה(W/m²)": "Direct Radiation_(W/m²)",
                               "קרינה מפוזרת(W/m²)": "scattered radiation_(W/m²)",
                               'כמות גשם(מ"מ)': "Rain_(mm)",
                               "מהירות הרוח(m/s)": "wind_speed_(m/s)",
                               "כיוון הרוח(מעלות)": "wind_direction_(deg)",
                               "סטיית התקן של כיוון הרוח(מעלות)": "wind_direction_std_(deg)",
                               "מהירות המשב העליון(m/s)": "upper_gust_(m/s)",
                               "כיוון המשב העליון(מעלות)": "upper_gust_direction_(deg)",
                               'מהירות רוח דקתית מקסימלית(m/s)': "maximum_wind_1minute(m/s)",
                               "מהירות רוח 10 דקתית מקסימלית(m/s)": "maximum_wind_10minute(m/s)",
                               "זמן סיום 10 הדקות המקסימליות()": "maximum_wind_10minute_time"

                               }
        self._hebStnRenameDict = {'בית דגן                                           ': "Bet_Dagan"

                                  }


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

    def getDocFromFile(self, path, station_column='stn_name', time_coloumn='time_obs', **kwargs):

        """
        Reads data from file/directory and returns a 'metadata like' object

        parameters
        ----------

        path : The path to the data file
        time_coloumn : The name of the Time column for indexing. default ‘time_obs’
        kwargs :

        returns
        -------
        nonDBMetadata : list

        """

        loaded_dask, _ = self._parse(path=path, station_column=station_column, time_coloumn=time_coloumn)
        return [nonDBMetadata(loaded_dask, **kwargs)]

    def LoadData(self, newdata_path, outputpath, projectname, metadatafile=None, station_column='stn_name', time_coloumn='time_obs', **metadata):

        """
            This function load data from file to database:


        Parameters
        ----------

        newdata_path : string
            the path to the new data. in future might also be a web address.
        outputpath : string
            Destination folder path for saving files
        projectname : string
            The project to which the data is associated. Will be saved in Matadata
        metadatafile : string
            The path to a metadata file, if exist
        station_column : string
            The name of the 'Station Name' column, for the groupby method.  default 'stn_name'
        time_column : string
            The name of the Time column for indexing. default 'time_obs'
        metadata : dict, optional
            These parameters will be added into the metadata desc.

        """

        metadata.update(dict(DataSource=self._DataSource))

        # 1- load the data #

        loaded_dask, metadata_dict = self._parse(path=newdata_path,
                                                 station_column=station_column,
                                                 time_coloumn=time_coloumn,
                                                 metadatafile=metadatafile,
                                                 **metadata
                                                 )

        groupby_data = loaded_dask.groupby(station_column)

        for stnname in metadata_dict:
            stn_dask = groupby_data.get_group(stnname)

            filtered_stnname = metadata_dict['StationName']
            print('updating %s data' % filtered_stnname)

            # 2- check if station exist in DataBase #

            docList = datalayer.Measurements.getDocuments(projectName=projectname,
                                                          type=type,
                                                          DataSource=self._DataSource,
                                                          StationName=filtered_stnname
                                                          )

            dir_path = os.path.join(outputpath, filtered_stnname).replace(' ', '_')

            if docList:
                if len(docList)>1:
                    raise ValueError("the list should be at max length of 1. Check your query.")
                else:

                    # 3- get current data from database
                    doc = docList[0]
                    stn_db = doc.getData()
                    data = [stn_db.reset_index(), stn_dask]
                    new_Data = dd.concat(data, interleave_partitions=True)\
                                 .set_index(time_coloumn)\
                                 .drop_duplicates()\
                                 .repartition(partition_size=self._np_size)

                    new_Data.to_parquet(doc.resource, engine='pyarrow')

                    if doc.resource != dir_path:
                        warnings.warn('The outputpath argument does not match the resource of the matching data '
                                      'in the database.\nThe new data is saved in the resource of the matching '
                                      'old data: %s' % doc.resource,
                                      ResourceWarning)

            else:
                os.makedirs(dir_path, exist_ok=True)

                # 4- create meta data
                desc = metadata_dict[stnname]

                new_Data = stn_dask.repartition(partition_size=self._np_size)
                new_Data.to_parquet(dir_path, engine='pyarrow')

                datalayer.Measurements.addDocument(projectName=projectname,
                                                   resource=dir_path,
                                                   dataFormat='parquet',
                                                   type=type,
                                                   desc=desc
                                                   )


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
