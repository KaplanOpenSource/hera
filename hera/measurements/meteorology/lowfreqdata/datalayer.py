import os
import dask.dataframe
import warnings
import pydoc
import pandas


from ....datalayer import datatypes
from ....datalayer.document import nonDBMetadataFrame
from .... import toolkit

from .analysis import analysis
from .presentationLayer import presenation

class lowFreqToolKit(toolkit.abstractToolkit):
    """
        Manages the loading and storing of low frequency meteorological data

        The data can be in the formats:

        - TOAA
        - IMS_JSON
        -

        TODO:
            Complete the other parsers_old from the older versions.

    """

    _np_size = None
    _doc_type = None


    STATIONNAME = 'StationName'

    def __init__(self, projectName, filesDirectory=None):
        """
            Initializes a datalayer for the lowfreqdata data.

            Also looks up the 'IMSData' in the public database.

        Parameters
        ----------

        projectName: str
                The project name
        """
        super().__init__(projectName=projectName, toolkitName="lowFreqMeteorology", filesDirectory=filesDirectory)
        self.logger.info("Init Low frequency data")

        self._analysis = analysis(self)
        self._presentation = presenation(self,self.analysis)


        self._np_size = "100MB"


    @property
    def docType(self):
        return f"{self.toolkitName}_LowFreqData"



    def loadData(self, fileNameOrData, parser, saveMode=toolkit.TOOLKIT_SAVEMODE_NOSAVE, additional_data=dict()):
        """
            Loads a low frequency meteorological data to a dask/pandas dataframe.

            If saveMode is not NO_SAVE the data is saved to a parquet file and possibly
            stores it to the DB.


        Parameters
        ----------

        fileNameOrData: str, dask.dataframe, pandas.dataframe
            if  str, the file name or the path to the files.
            if  dask.dataframe or pandas.dataframe it is the data.

        saveMode: str
                Can be either:

                    - TOOLKIT_SAVEMODE_NOSAVE   : Just load the data from file and return the datafile

                    - TOOLKIT_SAVEMODE_ONLYFILE : Loads the data from file and save to a file.
                                                  raise exception if file exists.

                    - TOOLKIT_SAVEMODE_ONLYFILE_REPLACE: Loads the data from file and save to a file.
                                                  Replace the file if it exists.

                    - TOOLKIT_SAVEMODE_FILEANDDB : Loads the data from file and save to a file and store to the DB as a source.
                                                    Raise exception if the entry exists.

                    - TOOLKIT_SAVEMODE_FILEANDDB_REPLACE: Loads the data from file and save to a file and store to the DB as a source.
                                                    Replace the entry in the DB if it exists.

        parser: str
            The name of the parser to use

            current parsers_old available:
                - IMS

        additional_data: dictdata to store in the DB metadata.
            additional
        Returns
        -------
            Returns a list of documents, one for each station.

        """
        if parser not in self.parserList():
            raise ValueError(f"{parser} is not a valid parser. Must be {','.join(self.parserList())}")

        if isinstance(fileNameOrData,str):
            pathToData = os.path.abspath(fileNameOrData)
            className = ".".join(__class__.__module__.split(".")[:-1])
            parserPath = f"{className}.parsers_old.Parser_{parser}"
            parserCls = pydoc.locate(parserPath)
            parser = parserCls()

            data = parser.parse(pathToData=pathToData)

        elif isinstance(fileNameOrData,pandas.DataFrame) or isinstance(fileNameOrData,dask.dataframe):
            data = fileNameOrData
        else:
            raise ValueError("fileNameOrData must be a path to a data (or a file), or a dask/pandas dataframe")


        if saveMode in [toolkit.TOOLKIT_SAVEMODE_ONLYFILE,
                        toolkit.TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,
                        toolkit.TOOLKIT_SAVEMODE_FILEANDDB,
                        toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE]:

            groupby_data = data.groupby(parser.station_column)

            ret = []
            for stnname,stndata in groupby_data:

                filteredName = parser.filterStationName(stnname)

                doc = self._storeStation(datasourceName=filteredName,
                                           stationData=stndata,
                                           parser=parser,
                                           metadata=dict())

                ret.append(doc)



        else:
            self.logger.debug("Creating nonDBMetadataFrame with loaded_dask")
            ret = [nonDBMetadataFrame(data)]

        return ret


    def _storeStation(self,datasourceName,stationData,parser,metadata=dict()):

        doc = self.getDatasourceData(datasourceName=datasourceName)

        if doc is None:

            outputFile = os.path.join(self.FilesDirectory, "parquet",f"{datasourceName}.parquet".replace(' ', '_'))


            new_Data = stationData.repartition(partition_size=self._np_size)
            new_Data.to_parquet(outputFile, engine='fastparquet')

            doc = self.addDataSource(dataSourceName=datasourceName,
                                     resource=outputFile,
                                     dataFormat=datatypes.PARQUET,
                                     **metadata)

            ret = doc

        else:
            data = [doc.getData().reset_index(), stationData.reset_index()]
            new_Data = dask.dataframe.concat(data, interleave_partitions=True) \
                .set_index(parser.time_column) \
                .drop_duplicates() \
                .repartition(partition_size=self._np_size)

            new_Data.to_parquet(doc.resource, engine='fastparquet')
            ret = doc
        return ret

    def parserList(self):
        """
            Return the list of parsers_old.

        Returns
        -------
            list of str
        """
        className = ".".join(__class__.__module__.split(".")[:-1])
        parserPath = f"{className}.parsers_old"
        mod = pydoc.locate(parserPath)
        return [x.split("_")[1] for x in dir(mod) if x.startswith("Parser_")]



