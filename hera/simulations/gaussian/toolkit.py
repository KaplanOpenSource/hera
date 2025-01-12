from ...toolkit import abstractToolkit
from .Sigma import BriggsRural
from .gasCloud import abstractGasCloud
from .Meteorology import MeteorologyFactory
from ...utils import *

class gaussianToolkit(abstractToolkit):

    _sigmaDict = None

    def __init__(self, projectName: str, filesDirectory: str = None):
        """
            Initializes the toolkit
        Parameters
        ----------
        projectName
        filesDirectory
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























