from ...toolkit import abstractToolkit
from .Sigma import BriggsRural
from .gasCloud import abstractGasCloud
from .Meteorology import MeteorologyFactory
from ...utils import *

class gaussianToolkit(abstractToolkit):
    """
    Toolkit for Gaussian plume dispersion modeling.

    The gaussianToolkit implements Gaussian dispersion models for atmospheric
    releases. It supports various sigma (dispersion coefficient) formulations
    and can generate gas cloud objects for further analysis.

    Key Features
    ------------
    - Multiple sigma type formulations (Briggs rural/urban)
    - Meteorology factory for creating meteorological conditions
    - Gas cloud generation for instantaneous and continuous releases
    - Integration with unum units for dimensional consistency

    Supported Sigma Types
    --------------------
    - briggsRural: Briggs rural dispersion coefficients

    Examples
    --------
    >>> from hera import toolkitHome
    >>> from unum.units import *
    >>> 
    >>> gauss_tk = toolkitHome.getToolkit("GaussianDispersion", projectName="my_project")
    >>> 
    >>> # Get meteorology
    >>> meteo = gauss_tk.getMeteorologyFromU10(
    ...     u10=5*m/s,
    ...     inversion=100*m,
    ...     stability="D"
    ... )
    >>> 
    >>> # Create gas cloud
    >>> cloud = gauss_tk.getGasCloud(
    ...     sourceQ=100*kg/s,
    ...     sourceHeight=10*m,
    ...     initialCloudSize=(1*m, 1*m, 1*m),
    ...     sigmaTypeName="briggsRural"
    ... )
    """
    _sigmaDict = None

    def __init__(self, projectName: str, filesDirectory: str = None):
        """
        Initialize the Gaussian dispersion toolkit.

        Parameters
        ----------
        projectName : str
            Name of the Hera project to work with.
        filesDirectory : str, optional
            Directory for storing output files. Default is None.
        """
        super().__init__(projectName=projectName, toolkitName="gaussianToolkit", filesDirectory=filesDirectory)
        self._sigmaDict = dict(briggsRural=BriggsRural)

    def getSigmaType(self,sigmaName):
        """

        Parameters
        ----------
        sigmaName

        Returns
        -------

        """
        try:
            sigmaCls = self._sigmaDict[sigmaName]
        except KeyError:
            err = f"The type {sigmaName} is not found. Must be one of {','.join(self.listSigmaTypes())}"
            raise ValueError(err)
        return sigmaCls()

    def listSigmaTypes(self):
        """
            Print the list of sigma types
        Returns
        -------


        """
        return [x for x in self._sigmaDict.keys()]


    def getMeteorologyFromU10(self, u10, inversion, verticalProfileType="log", temperature=20*celsius, stability="D",
                              z0=0.1*m, ustar=0.3*m/s, skinSurfaceTemperature=35*celsius):
        return MeteorologyFactory().getMeteorologyFromU10(u10=u10, inversion=inversion, verticalProfileType=verticalProfileType,
                    temperature=temperature, stability=stability, z0=z0, ustar=ustar, skinSurfaceTemperature=skinSurfaceTemperature)


    def getMeteorologyFromURefHeight(self, u, refHeight, inversion, verticalProfileType="log", temperature=20*celsius, stability="D",
                              z0=0.1*m, ustar=0.3*m/s, skinSurfaceTemperature=35*celsius):
        return MeteorologyFactory().getMeteorologyFromURefHeight(u=u, refHeight=refHeight,  inversion=inversion,
                    verticalProfileType=verticalProfileType, temperature=temperature, stability=stability, z0=z0,
                    ustar=ustar,skinSurfaceTemperature=skinSurfaceTemperature)


    def getGasCloud(self, sourceQ, sourceHeight, initialCloudSize, sigmaTypeName="briggsRural"):
        """

        Parameters
        ----------
        sourceQ : unum, method
            If unum:
                The unit determine the release time.
                [mass] - Instantaneous
                [mass/time] - Continuous
            else
                Continuous (not implementaed yet.)

        sourceHeight : unum
        initialCloudSize : 3-touple unum, the sigmas in each axis.
        sigmaTypeName : Name of the sigma type, for example from Briggs, rural/urban.

        Returns
        -------
        An instance of the class gadCloud

        """
        sigmaType = self.getSigmaType(sigmaTypeName)
        gascloud = abstractGasCloud.createGasCloud(sourceQ=sourceQ,sourceHeight=sourceHeight,initialCloudSize=initialCloudSize,sigmaType=sigmaType)
        return gascloud























