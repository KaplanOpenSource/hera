import os
import dask.dataframe
import warnings
import pydoc
import pandas


from ....datalayer.document import nonDBMetadataFrame
from ....datalayer import ProjectMultiDBPublic,datatypes


class DataLayer(ProjectMultiDBPublic):
    """
        Loads the lowfreqdata data format.
    """
    _HebRenameDict = None
    _hebStnRename = None
    _np_size = None

    _doc_type = None

    _outputPath = None

    JSONIMS = "JSONIMS"

    STATIONNAME = 'StationName'

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
        database entry or None

        """
        if StationName is not None:
            kwargs['StationName'] = StationName

        docList = self.getMeasurementsDocuments(dataFormat=datatypes.PARQUET,type="IMSData",**kwargs)
        return None if len(docList) == 0 else docList[0]

    def listStations(self):
        """
            Return the list of stations that are currently loaded in the DB (for the project that was initialized)
        Returns
        -------

        list
            A list of station names.
        """
        docList = self.getMeasurementsDocuments(type="IMSData")
        return [x.desc[self.STATIONNAME] for x in docList]

    def listStationsMetaData(self):
        """
            Return all the metadata of the stations.

        Returns
        --------
        pandas.DataFrame
            The metadata of the padas dataframe
        """
        docList = self.getMeasurementsDocuments(type="IMSData")
        return pandas.DataFrame([x.desc for x in docList])

    def getStationDataFromFile(self,
                               pathToData,
                               fileFormat=None,
                               storeDB = True,
                               storeTo = None,
                               **kwargs):
        """
                Reads data from file/directory and returns a 'metadata like' object

        parameters
        ----------

        pathToData: str
                The file name or path that contains the json files.
                if pathToData is a directory, then load all the .json in the directory.

        fileFormat: str
                The format of the file.
                if fileFormat is None, use JSONIMS

        storeDB: bool
                If true, store the data in the database. If the

        storeTo: str
                The name of the database to store the loaded data.
                When None use the default database (the current username).

        kwargs: general parameters.
                Not used now.

        Returns
        -------
            Database Document of station


        """
        fileFormat = self.JSONIMS if fileFormat is None else fileFormat

        self.logger.info(f"getStationDataFromFile: {pathToData} with format {fileFormat}. {'storing in DB' if storeDB else 'without storing'} (storeDB={storeDB})")
        className = ".".join(__class__.__module__.split(".")[:-1])
        parserPath = f"{className}.parsers.Parser_{fileFormat}"

        parserCls = pydoc.locate(parserPath)
        parser = parserCls()

        loaded_dask, metadata_dict  = parser.parse(pathToData=pathToData)
        self.logger.debug(f"Loaded file: {loaded_dask.head()}")

        if storeDB:
            groupby_data = loaded_dask.groupby(parser.station_column)
            for stnname in metadata_dict:
                stn_dask = groupby_data.get_group(stnname)

                filtered_stnname = metadata_dict[stnname]['StationName']
                self.logger.execution('updating the data of %s' % filtered_stnname)

                # 2- check if station exist in DataBase #
                doc = self.getStationDataFromDB(StationName=filtered_stnname)

                if doc is None:

                    if self.outputPath is None:
                        self.outputPath = os.path.dirname(os.path.abspath(pathToData))

                    dir_path = os.path.join(self.outputPath, f"{filtered_stnname}.parquet").replace(' ', '_')


                    # 4- create meta data
                    desc = metadata_dict[stnname]

                    new_Data = stn_dask.repartition(partition_size=self._np_size)
                    new_Data.to_parquet(dir_path, engine='fastparquet')

                    doc = self.addMeasurementsDocument(resource=dir_path,
                                                       dataFormat=datatypes.PARQUET,
                                                       type=self.docType,
                                                       desc=desc,
                                                       users=storeTo)

                    ret = doc

                    #addMeasurementsDocument(self, resource="", dataFormat="string", type="", desc={}, users=None)

                else:
                    data = [doc.getData().reset_index(), stn_dask.reset_index()]
                    new_Data = dask.dataframe.concat(data, interleave_partitions=True) \
                        .set_index(parser.time_column) \
                        .drop_duplicates() \
                        .repartition(partition_size=self._np_size)

                    new_Data.to_parquet(doc.resource, engine='fastparquet')

                    ret = doc

        else:
            self.logger.debug("Creating nonDBMetadataFrame with loaded_dask")

            ret = nonDBMetadataFrame(loaded_dask)

        return ret


