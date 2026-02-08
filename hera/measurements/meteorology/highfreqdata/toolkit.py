from .analysis.analysislayer import RawdataAnalysis
from hera import toolkit
from hera.utils.logging import get_classMethod_logger
from .parsers.CampbellBinary import Parser
from .parsers.TOA5 import ASCIIParser

class HighFreqToolKit(toolkit.abstractToolkit):
    """
    Toolkit for managing high-frequency meteorological sonic anemometer data.

    The HighFreqToolKit handles loading, parsing, and storing high-frequency
    meteorological measurements from sonic anemometers. It supports multiple data
    formats and provides analysis capabilities for turbulence and wind measurements.

    Key Features
    ------------
    - Support for Campbell Scientific binary format
    - TOA5 ASCII format support (planned)
    - Conversion to Parquet format for efficient storage
    - Time-based data filtering
    - Analysis layer for turbulence statistics

    Data Sources
    ------------
    Data sources are station measurement documents stored in the project database.
    Each station can have multiple measurement time series.

    Supported Formats
    ----------------
    - CampbellBinary: Binary format from Campbell Scientific dataloggers
    - TOA5: ASCII format (not yet implemented)

    Examples
    --------
    >>> from hera import toolkitHome
    >>> meteo_tk = toolkitHome.getToolkit("MeteoHighFreq", projectName="my_project")
    >>> 
    >>> # Convert Campbell binary file to Parquet
    >>> df = meteo_tk.campbelToParquet(
    ...     binaryFile="/path/to/data.dat",
    ...     fromTime="2024-01-01",
    ...     toTime="2024-01-31"
    ... )
    """

    DOCTYPE_STATIONS = 'StationsData'
    DOCTYPE_MEASUREMENTS = 'MeasurementsData'

    def __init__(self, projectName, filesDirectory=None):
        """
            Initializes a datalayer for the highfreqdata data.


        Parameters
        ----------

        projectName: str
                The project name
        """
        super().__init__(projectName=projectName, toolkitName="highFreqMeteorology", filesDirectory=filesDirectory)
        logger = get_classMethod_logger(self,"init")
        logger.info("Init High frequency data")
        self._analysis = RawdataAnalysis(self)
       # self._presentation = presenation(self,self.datalayer)

    @property
    def docType(self):
        return f"{self.toolkitName}_HighFreqData"


    def campbelToParquet(self, binaryFile, fromTime=None, toTime=None, chunkSize=10000):
        """
        Read Campbell Scientific binary file and convert to Dask DataFrame.

        This method parses Campbell Scientific binary data files and returns a
        Dask DataFrame for efficient handling of large datasets. The data can be
        filtered by time range.

        Parameters
        ----------
        binaryFile : str
            Path to the Campbell Scientific binary data file (.dat file).

        fromTime : str or datetime, optional
            Start time for data extraction. If None, reads from the beginning
            of the file. Can be ISO format string or datetime object.
            Default is None.

        toTime : str or datetime, optional
            End time for data extraction. If None, reads to the end of the file.
            Can be ISO format string or datetime object. Default is None.

        chunkSize : int, optional
            Number of rows to process in each chunk. Larger values use more
            memory but may be faster. Default is 10000.

        Returns
        -------
        dask.dataframe.DataFrame
            Dask DataFrame containing the parsed meteorological data with
            columns for wind components, temperature, and other measured variables.

        Examples
        --------
        >>> # Read entire file
        >>> df = meteo_tk.campbelToParquet("/path/to/data.dat")
        >>> 
        >>> # Read specific time range
        >>> df = meteo_tk.campbelToParquet(
        ...     "/path/to/data.dat",
        ...     fromTime="2024-01-01 00:00:00",
        ...     toTime="2024-01-31 23:59:59"
        ... )
        """
        campelParser = Parser(chunkSize = chunkSize)
        return campelParser.parse(path=binaryFile,fromTime=fromTime,toTime=toTime)

    def asciiToParquet(self, path, fromTime=None, toTime=None):
        """
        Read TOA5 ASCII file and convert to list of DataFrames.

        Parses TOA5 format ASCII files containing meteorological data. Returns
        a list of DataFrames, one for each device/station in the file.

        Parameters
        ----------
        path : str
            Path to the TOA5 ASCII data file.

        fromTime : str or datetime, optional
            Start time for data extraction. If None, reads from the beginning.
            Default is None.

        toTime : str or datetime, optional
            End time for data extraction. If None, reads to the end.
            Default is None.

        Returns
        -------
        list of pandas.DataFrame
            List of DataFrames, one per device/station in the file. Each DataFrame
            contains time series data for that device.

        Examples
        --------
        >>> # Read entire file
        >>> dfs = meteo_tk.asciiToParquet("/path/to/data.txt")
        >>> 
        >>> # Read specific time range
        >>> dfs = meteo_tk.asciiToParquet(
        ...     "/path/to/data.txt",
        ...     fromTime="2024-01-01",
        ...     toTime="2024-01-31"
        ... )
        """
        asciiParser = ASCIIParser()
        return asciiParser.parse(path=path, fromTime=fromTime, toTime=toTime)
