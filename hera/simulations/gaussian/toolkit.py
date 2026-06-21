from hera.toolkit import abstractToolkit
from hera.simulations.gaussian.Sigma import BriggsRural
from hera.simulations.gaussian.gasCloud import abstractGasCloud
from hera.simulations.gaussian.Meteorology import MeteorologyFactory
from hera.utils.unitHandler import ureg, unumToPint, celsius, K
import matplotlib.pyplot as plt



class gaussianToolkit(abstractToolkit):

    _sigmaDict = None
    _presentation = None

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

        self._presentation = presentationLayer()

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






    def getMeteorologyFromU10(self, u10, inversion, verticalProfileType="log", temperature=20*ureg.degC, stability="D",
                              z0=0.1*ureg.m, ustar=0.3*ureg.m/ureg.s, skinSurfaceTemperature=35*ureg.degC):
        return MeteorologyFactory().getMeteorologyFromU10(u10=u10, inversion=inversion, verticalProfileType=verticalProfileType,
                    temperature=temperature, stability=stability, z0=z0, ustar=ustar, skinSurfaceTemperature=skinSurfaceTemperature)


    def getMeteorologyFromURefHeight(self, u, refHeight, inversion, verticalProfileType="log", temperature=20*ureg.degC, stability="D",
                              z0=0.1*ureg.m, ustar=0.3*ureg.m/ureg.s, skinSurfaceTemperature=35*ureg.degC):
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


    @property
    def presentation(self):

        return self._presentation







#--------------------------------------------------PRESENTATION_LAYER--------------------------------------------------

class presentationLayer:


    #--------------------------------Concentration--------------------------------

    def plotMaxConcentrationOverTime(self,C, y, z, x_min=None,x_max=None):
        """
        Given a line along the X-axis, this function plots the maximum concentration (over time).
        :param C: xarray of the concentrations in units of [mass/volume]
        :param y: Y-coordinate of the line
        :param z: Z-coordinate of the line
        :return:
        """
        y = unumToPint(y).m_as(ureg.m)
        z = unumToPint(z).m_as(ureg.m)

        conc_x_inst = C.sel(y=y, z=z, method="nearest").max(dim="time")
        x_array = conc_x_inst.squeeze().x
        plt.plot(x_array, conc_x_inst)
        plt.xlabel("Distance from source $[m]$")

        unitsTemp = str(C.attrs['Q'])
        units = unitsTemp.split(' ')[1]

        plt.ylabel(f"Concentration {units}")
        plt.title(f"Maximum concentration over time. y={y}[m], z={z}[m]")
        plt.grid()
        if x_min is not None and x_max is not None:
            x_min = unumToPint(x_min).m_as(ureg.m)
            x_max = unumToPint(x_max).m_as(ureg.m)
            plt.xlim(x_min, x_max)
        plt.show()


    def plotFixedPointConcentrationOverTime(self, C, x, y, z, t_min=None,t_max=None):
        """
        Given a point in space, this function plots the concentration as a function of time.
        :param C: xarray of the concentrations in units of [mass/volume]
        :param x: X-coordinate of the point
        :param y: Y-coordinate of the point
        :param z: Z-coordinate of the point
        :return:
        """

        x = unumToPint(x).m_as(ureg.m)
        y = unumToPint(y).m_as(ureg.m)
        z = unumToPint(z).m_as(ureg.m)

        conc_inst_t = C.sel(x=x, y=y, z=z, method="nearest")
        time_array = conc_inst_t.squeeze().time

        unitsTemp = str(C.attrs['Q'])
        units = unitsTemp.split(' ')[1]

        plt.plot(time_array, conc_inst_t)
        plt.xlabel("Time from release $[min]$")
        plt.ylabel(f"Concentration {units}")
        plt.title(f"Receptor at x={x}[m], y={y}[m], z={z}[m].")
        plt.grid()
        if t_min is not None and t_max is not None:
            t_min = unumToPint(t_min).m_as(ureg.min)
            t_max = unumToPint(t_max).m_as(ureg.min)
            plt.xlim(t_min, t_max)
        plt.show()


#--------------------------------Dosage--------------------------------

    def plotDosagePerDistance(self,D, y, z, time):
        """
        Given a line along the X-axis, this function plots the dosage over distance.
        :param C: xarray of the concentrations in units of [mass*time/volume]
        :param y: Y-coordinate of the line
        :param z: Z-coordinate of the line
        :return:
        """
        y = unumToPint(y).m_as(ureg.m)
        z = unumToPint(z).m_as(ureg.m)
        time = unumToPint(time).m_as(ureg.min)

        dos_x_inst = D.sel(y=y, z=z, time=time, method="nearest")
        x_array = dos_x_inst .squeeze().x
        plt.plot(x_array, dos_x_inst )
        plt.xlabel("Distance from source $[m]$")

        unitsTemp = str(D.attrs['Q'])
        units = unitsTemp.split(' ')[1]

        plt.ylabel(f"Dosage  {units}")
        plt.title(f"Dosage per distance. y={y}[m], z={z}[m], time={time}[min]")
        plt.grid()
        plt.show()


    def plotFixedPointDosageOverTime(self, D, x, y, z):
        """
        Given a point in space, this function plots the dosage as a function of time.
        :param C: xarray of the dosage in units of [mass*time/volume]
        :param x: X-coordinate of the point
        :param y: Y-coordinate of the point
        :param z: Z-coordinate of the point
        :return:
        """

        x = unumToPint(x).m_as(ureg.m)
        y = unumToPint(y).m_as(ureg.m)
        z = unumToPint(z).m_as(ureg.m)

        dos_inst_t = D.sel(x=x, y=y, z=z, method="nearest")
        time_array = dos_inst_t .squeeze().time

        unitsTemp = str(D.attrs['Q'])
        units = unitsTemp.split(' ')[1]

        plt.plot(time_array, dos_inst_t )
        plt.xlabel("Time from release $[min]$")
        plt.ylabel(f"Dosage over time {units}")
        plt.title(f"Receptor at x={x}[m], y={y}[m], z={z}[m].")
        plt.grid()
        plt.show()



#--------------------------------For Shlomi--------------------------------


    def plotMaxConcentrationOverTime_noQ(self,C, y, z, x_min=None,x_max=None):
        """
        Given a line along the X-axis, this function plots the maximum concentration (over time).
        :param C: xarray of the concentrations in units of [1/volume]
        :param y: Y-coordinate of the line
        :param z: Z-coordinate of the line
        :return:
        """
        y = unumToPint(y).m_as(ureg.m)
        z = unumToPint(z).m_as(ureg.m)

        conc_x_inst = C.sel(y=y, z=z, method="nearest").max(dim="time")
        x_array = conc_x_inst.squeeze().x
        plt.plot(x_array, conc_x_inst)
        plt.xlabel("Distance from source $[m]$")

        plt.ylabel(r"Concentration $\left[\frac{1}{m^3}\right]$")
        plt.title(f"Maximum concentration over time. y={y}[m], z={z}[m]")
        plt.grid()
        if x_min is not None and x_max is not None:
            x_min = unumToPint(x_min).m_as(ureg.m)
            x_max = unumToPint(x_max).m_as(ureg.m)
            plt.xlim(x_min, x_max)
        plt.show()


    def plotFixedPointConcentrationOverTime_noQ(self, C, x, y, z, t_min=None,t_max=None):
        """
        Given a point in space, this function plots the concentration as a function of time.
        :param C: xarray of the concentrations in units of [1/volume]
        :param x: X-coordinate of the point
        :param y: Y-coordinate of the point
        :param z: Z-coordinate of the point
        :return:
        """

        x = unumToPint(x).m_as(ureg.m)
        y = unumToPint(y).m_as(ureg.m)
        z = unumToPint(z).m_as(ureg.m)

        conc_inst_t = C.sel(x=x, y=y, z=z, method="nearest")
        time_array = conc_inst_t.squeeze().time

        plt.plot(time_array, conc_inst_t)
        plt.xlabel("Time from release $[min]$")
        plt.ylabel(r"Concentration $\left[\frac{1}{m^3}\right]$")
        plt.title(f"Receptor at x={x}[m], y={y}[m], z={z}[m].")
        plt.grid()
        if t_min is not None and t_max is not None:
            t_min = unumToPint(t_min).m_as(ureg.min)
            t_max = unumToPint(t_max).m_as(ureg.min)
            plt.xlim(t_min, t_max)
        plt.show()



    def plotTIACPerDistance_noQ(self,TIAC, y, z, time):
        """
        Given a line along the X-axis, this function plots the dosage over distance.
        :param TIAC: xarray of the TIAC in units of [1*time/volume]
        :param y: Y-coordinate of the line
        :param z: Z-coordinate of the line
        :return:
        """
        y = unumToPint(y).m_as(ureg.m)
        z = unumToPint(z).m_as(ureg.m)
        time = unumToPint(time).m_as(ureg.min)

        dos_x_inst = TIAC.sel(y=y, z=z, time=time, method="nearest")
        x_array = dos_x_inst .squeeze().x
        plt.plot(x_array, dos_x_inst )
        plt.xlabel("Distance from source $[m]$")

        plt.ylabel(r"Dosage $\left[\frac{1}{m^3} \cdot min\right]$")
        plt.title(f"Dosage per distance. y={y}[m], z={z}[m], time={time}[min]")
        plt.grid()
        plt.show()


    def plotFixedPointTIACOverTime_noQ(self, TIAC, x, y, z):
        """
        Given a point in space, this function plots the dosage as a function of time.
        :param TIAC: xarray of the TIAC in units of [1*time/volume]
        :param x: X-coordinate of the point
        :param y: Y-coordinate of the point
        :param z: Z-coordinate of the point
        :return:
        """

        x = unumToPint(x).m_as(ureg.m)
        y = unumToPint(y).m_as(ureg.m)
        z = unumToPint(z).m_as(ureg.m)

        dos_inst_t = TIAC.sel(x=x, y=y, z=z, method="nearest")
        time_array = dos_inst_t .squeeze().time

        plt.plot(time_array, dos_inst_t )
        plt.xlabel("Time from release $[min]$")
        plt.ylabel(r"Dosage over time $\left[\frac{1}{m^3} \cdot min\right]$")
        plt.title(f"Receptor at x={x}[m], y={y}[m], z={z}[m].")
        plt.grid()
        plt.show()




