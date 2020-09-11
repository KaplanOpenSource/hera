import os
import dask.dataframe
from ..datalayer import DataLayer as meteorological_datalayer
from ....datalayer import nonDBMetadata
import warnings


class IMSDatalayer(meteorological_datalayer):
    """
        Loads the IMS data format.
    """
    _HebRenameDict = None
    _hebStnRename = None

    def __init__(self, projectName, databaseNameList=None, useAll=False):
        """
            Initializes a datalayer for the IMS data.

            Also looks up the 'IMSData' in the public database.

        Parameters
        ----------

        projectName: str
                The project name
        databaseNameList: list
            The list of database naems to look for data.
            if None, uses the database with the name of the user.

        useAll: bool
            Whether or not to return the query results after found in one DB.
        """
        super().__init__(DataSource='IMS',
                         projectName=projectName,
                         publicProjectName="IMSData",
                         databaseNameList=databaseNameList,
                         useAll=useAll)

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

    def getDocFromDB(self, resource=None, dataFormat=None, StationName=None, **kwargs):

        """
        This function returns a list of 'doc' objects from the database that matches the requested query

        parameters
        ----------

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

        if StationName is not None:
            kwargs['StationName'] = StationName

        docList = super().getDocFromDB(resource=resource,
                                       dataFormat=dataFormat,
                                       desc=kwargs)
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

        loaded_dask, _ = self.parse(path=path, station_column=station_column, time_coloumn=time_coloumn)
        return [nonDBMetadata(loaded_dask, **kwargs)]

    def loadData(self, newdata_path, outputpath, projectName, metadatafile=None, station_column='stn_name',
                 time_column='time_obs', **metadata):

        """
            This function load data from file to database:


        Parameters
        ----------

        newdata_path : string
            the path to the new data. in future might also be a web address.
        outputpath : string
            Destination folder path for saving files
        projectName : string
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

        metadata['DataSource'] = self._DataSource

        # 1- load the data #

        loaded_dask, metadata_dict = self.parse(path=newdata_path,
                                                station_column=station_column,
                                                time_column=time_column,
                                                metadatafile=metadatafile,
                                                **metadata
                                                )

        groupby_data = loaded_dask.groupby(station_column)

        for stnname in metadata_dict:
            stn_dask = groupby_data.get_group(stnname)

            filtered_stnname = metadata_dict[stnname]['StationName']
            print('updating %s data' % filtered_stnname)

            # 2- check if station exist in DataBase #

            docList = self.get(projectName=projectName,
                               type=self._docType,
                               DataSource=self._DataSource,
                               StationName=filtered_stnname
                               )

            dir_path = os.path.join(outputpath, filtered_stnname).replace(' ', '_')

            if docList:
                if len(docList) > 1:
                    raise ValueError("the list should be at max length of 1. Check your query.")
                else:

                    # 3- get current data from database
                    doc = docList[0]
                    stn_db = doc.getDocFromDB()
                    data = [stn_db.reset_index(), stn_dask]
                    new_Data = dask.dataframe.concat(data, interleave_partitions=True) \
                        .set_index(time_column) \
                        .drop_duplicates() \
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

                self.addMeasurementsDocument(resource=dir_path,
                                             dataFormat='parquet',
                                             type=self._docType,
                                             desc=desc)


class CampbellBinary(meteorological_datalayer):
    """
        A datalayer that loads the Campbell binary data.
    """

    def __init__(self, projectName, databaseNameList=None, useAll=False):
        super().__init__(DataSource='IMS',
                         projectName=projectName,
                         publicProjectName="IMSData",
                         databaseNameList=databaseNameList,
                         useAll=useAll)

    def getDocFromDB(self, resource=None, dataFormat=None, station=None, instrument=None, height=None, **kwargs):

        """
        This function returns a list of 'doc' objects from the database that matches the requested query

        parameters
        ----------

        resource: String/dict/JSON
            The resource of the data
        dataFormat: String
            The data format
        station: String
            The name of the requested station. default None
        instrument: String
            The name of the requested instrument. default None
        height: String/int
            The requested height. default None
        kwargs : dict
            Other properties for query

        returns
        -------
        docList : List

        """
        if station is not None:
            kwargs['station'] = station
        if instrument is not None:
            kwargs['instrument'] = instrument
        if height is not None:
            kwargs['height'] = int(height)

        docList = super().getDocFromDB(resource=resource,
                                       dataFormat=dataFormat,
                                       desc=kwargs
                                       )
        return docList

    def getDocFromFile(self, path, **kwargs):

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

        loaded_dask, _ = self.parse(path=path)
        return [nonDBMetadata(loaded_dask, **kwargs)]

    def loadData(self, newdata_path, outputpath, projectName, **metadata):

        """
            This function load data from file to database:


        Parameters
        ----------

        newdata_path : string
            the path to the new data.
        outputpath : string
            Destination folder path for saving files
        projectName : string
            The project to which the data is associated. Will be saved in Matadata
        metadata : dict, optional
            These parameters will be added into the metadata desc.

        """

        metadata['DataSource'] = self._DataSource

        # 1- load the data #

        loaded_dask, metadata_dict = self.parse(path=newdata_path, **metadata)

        groupby_data = loaded_dask.groupby(['station', 'instrument', 'height'])

        for station in metadata_dict:
            for instrument in metadata_dict[station]:
                for height in metadata_dict[station][instrument]:
                    dir_path = os.path.join(outputpath, station, instrument, str(height))
                    docList = self.getMeasurementsDocuments(projectName=projectName,
                                                            type=self._docType,
                                                            DataSource=self._DataSource,
                                                            station=station,
                                                            instrument=instrument,
                                                            height=height
                                                            )

                    new_dask = groupby_data.get_group((station, instrument, height)) \
                        .drop(columns=['station', 'instrument', 'height'])

                    for col in new_dask.columns:
                        if new_dask[col].isnull().all().compute():
                            new_dask = new_dask.drop(col, axis=1)

                    if docList:
                        if len(docList) > 1:
                            raise ValueError("the list should be at max length of 1. Check your query.")
                        else:

                            # 3- get current data from database
                            doc = docList[0]
                            db_dask = doc.getDocFromDB()
                            data = [db_dask, new_dask]
                            new_Data = dask.dataframe.concat(data, interleave_partitions=True) \
                                .reset_index() \
                                .drop_duplicates() \
                                .set_index('index') \
                                .repartition(partition_size=self._np_size)

                            new_Data.to_parquet(doc.resource, engine='pyarrow')

                            if doc.resource != dir_path:
                                warnings.warn(
                                    'The outputpath argument does not match the resource of the matching data '
                                    'in the database.\nThe new data is saved in the resource of the matching '
                                    'old data: %s' % doc.resource,
                                    ResourceWarning)

                    else:
                        os.makedirs(dir_path, exist_ok=True)

                        # 4- create meta data
                        desc = metadata_dict[station][instrument][height]

                        new_Data = new_dask.repartition(partition_size=self._np_size)
                        new_Data.to_parquet(dir_path, engine='pyarrow')

                        self.addMeasurementsDocument(projectName=projectName,
                                                     resource=dir_path,
                                                     dataFormat='parquet',
                                                     type=self._docType,
                                                     desc=desc
                                                     )
