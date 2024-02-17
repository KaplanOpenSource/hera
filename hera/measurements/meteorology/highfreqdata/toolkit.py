from .analysis.analysislayer import RawdataAnalysis
from .... import toolkit
from .... import get_classMethod_logger
from .parsers.CampbellBinary import Parser

class HighFreqToolKit(toolkit.abstractToolkit):
    """
        Manages the loading and storing of high frequency sonic data

        The data can be in the formats:

        - CampbellBinary
        - TOA5 (not implemented yet)


        TODO:
            Complete the other parsers from the older versions.

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


    def campbelToParquet(self,binaryFile,fromTime=None, toTime=None,chunkSize = 10000):
        """
            Reads the cambell binary file from fromTime to toTime and return a dask dataframe.


        Parameters
        ----------
        binaryFile
        fromTime
        toTime
        existintData

        Returns
        -------

        """
        campelParser = Parser(chunkSize = chunkSize)
        return campelParser.parse(path=binaryFile,fromTime=fromTime,toTime=toTime)