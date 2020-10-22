import os
import dask.dataframe
import warnings
import pydoc

from .parsers import JSONIMS

from ....datalayer.document import nonDBMetadataFrame
from ....datalayer import ProjectMultiDBPublic,datatypes


class datalayer(ProjectMultiDBPublic):
    """
        Loads the lowfreqdata data format.
    """
    _HebRenameDict = None
    _hebStnRename = None
    _np_size = None

    _doc_type = None

    _outputPath = None

    def __init__(self, projectName):
        """
            Initializes a datalayer for the lowfreqdata data.

            Also looks up the 'IMSData' in the public database.

        Parameters
        ----------

        projectName: str
                The project name
        """
        super().__init__(projectName=projectName,
                         publicProjectName='lowfreqdata',
                         databaseNameList=None,
                         useAll=True)

        self.logger.info("Init Low frequency data")

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
        self._hebStnRenameDict = {'בית דגן                                           ': "Bet_Dagan"}
        self._np_size = "100MB"

        self._doc_type = "IMSData"

    @property
    def docType(self):
        return self._doc_type


    @property
    def outputPath(self):
        return self._outputPath

    @outputPath.setter
    def outputPath(self, value):
        if value is not None:
            finalPath = os.path.join(value,"parquet")
            os.makedirs(finalPath, exist_ok=True)
        self._outputPath = value



    def getStationDataFromDB(self,StationName=None, **kwargs):
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

        docList = self.getMeasurementsDocuments(dataFormat=datatypes.PARQUET,type="IMSData",**kwargs)
        return docList

    def getStationDataFromFile(self,
                               fileName,
                               fileFormat=JSONIMS,
                               storeDB = True,
                               **kwargs):
        """
        Reads data from file/directory and returns a 'metadata like' object

        parameters
        ----------

        fileName : str
            The path to the data file
        fileName : str

        kwargs :

        returns
        -------
        nonDBMetadataFrame : list

        """

        classPath = ".".join(__class__.split(".")[:-1])
        parserPath = classPath + ".parsers.Parser_" + fileFormat

        parserCls = pydoc.locate(parserPath)
        parser = parserCls()

        loaded_dask, metadata_dict  = parser.parse(fileName=fileName)

        if storeDB:
            groupby_data = loaded_dask.groupby(parser.station_column)
            for stnname in metadata_dict:
                stn_dask = groupby_data.get_group(stnname)

                filtered_stnname = metadata_dict[stnname]['StationName']
                self.logger.execution('updating the data of %s' % filtered_stnname)

                # 2- check if station exist in DataBase #
                docList = self.getStationDataFromDB(StationName=filtered_stnname)

                if len(docList) == 0:

                    if self.outputPath is None:
                        self.outputPath = os.path.basedir(fileName)

                    dir_path = os.path.join(self.outputpath, filtered_stnname).replace(' ', '_')


                    # 4- create meta data
                    desc = metadata_dict[stnname]

                    new_Data = stn_dask.repartition(partition_size=self._np_size)
                    new_Data.to_parquet(dir_path, engine='pyarrow')

                    doc = self.addMeasurementsDocument(resource=dir_path,
                                                       dataFormat=datatypes.PARQUET,
                                                       type=self.docType,
                                                       desc=desc)

                    ret = [doc]

                else:

                    # 3- get current data from database
                    doc = docList[0]

                    data = [doc.getData().reset_index(), stn_dask]
                    new_Data = dask.dataframe.concat(data, interleave_partitions=True) \
                        .set_index(parser.time_column) \
                        .drop_duplicates() \
                        .repartition(partition_size=self._np_size)

                    new_Data.to_parquet(doc.resource, engine='pyarrow')

                    ret = docList[0]

            else:
                ret = [nonDBMetadataFrame(loaded_dask, **kwargs)]

        return ret



class CampbellBinary(meteorological_datalayer):
    """
        A datalayer that loads the Campbell binary data.
    """

    def __init__(self, projectName, databaseNameList=None, useAll=False):
        super().__init__(DataSource='lowfreqdata',
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
        nonDBMetadataFrame : list

        """

        loaded_dask, _ = self.parse(path=path)
        return [nonDBMetadataFrame(loaded_dask, **kwargs)]

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
