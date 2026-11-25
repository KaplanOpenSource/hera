import os
import dask.dataframe
import warnings
import pydoc
import pandas
import logging

from ....datalayer import datatypes
from ....datalayer.document import nonDBMetadataFrame
from .... import toolkit

from .analysis import analysis
from .presentationLayer import presenation

class lowFreqToolKit(toolkit.abstractToolkit):
    """
        Low frequency meteorological data.

        The toolkit has no ingestion method of its own. Station data is
        registered as ordinary datasources (typically parquet, through a
        repository JSON or addDataSource) and read with getDataSourceData.
    """

    _np_size = None
    _doc_type = None


    STATIONNAME = 'StationName'

    def __init__(self, projectName, filesDirectory=None,connectionName=None):
        """
            Initializes a datalayer for the lowfreqdata data.

            Also looks up the 'IMSData' in the public database.

        Parameters
        ----------

        projectName: str
                The project name
        """
        super().__init__(projectName=projectName, toolkitName="lowFreqMeteorology", filesDirectory=filesDirectory,connectionName=connectionName)
        self.logger = logging.getLogger()
        self.logger.info("Init Low frequency data")

        self._analysis = analysis(self)
        self._presentation = presenation(self,self.analysis)


        self._np_size = "100MB"


    @property
    def docType(self):
        return f"{self.toolkitName}_LowFreqData"
