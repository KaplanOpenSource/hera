from hera.toolkit import abstractToolkit
from hera.simulations.gaussian.Sigma import BriggsRural
from hera.simulations.gaussian.gasCloud import abstractGasCloud
from hera.simulations.gaussian.Meteorology import MeteorologyFactory
from hera.utils.unitHandler import ureg, unumToPint, celsius, K
import matplotlib.pyplot as plt
import scipy
import numpy



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




    def getMeteorologyFromU10(self, u10, inversion, verticalProfileType="powerLaw", temperature=ureg.Quantity(20, ureg.degC), stability="D",
                              z0=0.1*ureg.m, ustar=0.3*ureg.m/ureg.s, skinSurfaceTemperature=ureg.Quantity(35, ureg.degC)):
        return MeteorologyFactory().getMeteorologyFromU10(u10=u10, inversion=inversion, verticalProfileType=verticalProfileType,
                    temperature=temperature, stability=stability, z0=z0, ustar=ustar, skinSurfaceTemperature=skinSurfaceTemperature)


    def getMeteorologyFromURefHeight(self, u, refHeight, inversion, verticalProfileType="powerLaw", temperature=ureg.Quantity(20, ureg.degC), stability="D",
                              z0=0.1*ureg.m, ustar=0.3*ureg.m/ureg.s, skinSurfaceTemperature=ureg.Quantity(35, ureg.degC)):
        return MeteorologyFactory().getMeteorologyFromURefHeight(u=u, refHeight=refHeight,  inversion=inversion,
                    verticalProfileType=verticalProfileType, temperature=temperature, stability=stability, z0=z0,
                    ustar=ustar,skinSurfaceTemperature=skinSurfaceTemperature)





    def getSpaceTime(self, meteorology, sourceHeight, wind_profile_type, maxx, maxz, dt, dz, dxdy_multiplier, minimal_maxy, initialCloudSize):

        # setting the wind speed at the source height according to the Hot Spot's manual, or Hera's built-in functions.
        if wind_profile_type == 'HotSpot':
            windSpeed = meteorology.getWindVelocity_hotSpot(height=sourceHeight)
        elif wind_profile_type == 'default':
            windSpeed = meteorology.getWindVelocity(height=sourceHeight)
        else:
            raise ValueError("wind_profile_type must be either 'default' or 'HotSpot'")

        # it is preferable to take dxdy as a multiple of dt*windSpeed, so that the cloud center is directly above a grid-point.
        # the dxdy_multiplier makes the cloud above every dxdy_multiplier-th grid-point (on the datetime axis).
        dxdy = (dt * windSpeed * dxdy_multiplier).to(ureg.m)

        # the Y-span should be the smallest multiple of dxdy that is greater/equal to the initially given maxy.
        # the purpose is to limit the grid along the y-axis, to reduce runtime.
        maxy = numpy.ceil(minimal_maxy.m_as(ureg.m) / dxdy.m_as(ureg.m)) * dxdy

        # the time span should be such that the cloud has passed far enough beyond the maximum X (maxx),
        # yet not too far as to calculate for times that don't affect the concentraion up to maxx.
        # therefore we take a time span for which maxx is 3*sigmas away from the cloud center:
        def find_x_for_timeSpan(maxx):
            maxx = maxx.m_as(ureg.m)

            def get_function(x, maxx):
                sigmaX = self.getSigmaType(sigmaName='briggsRural').getSigma(x=x, stability=meteorology.stability,
                                                                           sigma0=initialCloudSize,
                                                                           units=False)['sigmaX'][0]
                return x - 3 * sigmaX - maxx

            root = scipy.optimize.root_scalar(get_function, bracket=[0, 10 ** 6], args=maxx).root
            return root

        x_timeSpan = find_x_for_timeSpan(maxx) * ureg.m
        timeSpan = x_timeSpan / (windSpeed * 60 * (ureg.s / ureg.min))
        timeSpan = numpy.ceil(timeSpan.m_as(ureg.min)) * ureg.min

        spaceTime = {
            'minx': 0 * ureg.m, 'maxx': maxx,
            'miny': -maxy, 'maxy': maxy + 1 * ureg.m,
            'minz': 0 * ureg.m, 'maxz': maxz,
            'dxdy': dxdy, 'dz': dz,
            'timeSpan': timeSpan, 'dt': dt}

        return spaceTime




    def getGasCloud(self, sourceQ, sourceHeight, initialCloudSize, meteorology, wind_profile_type,
                    spaceTime, deposition_velocity, sigmaTypeName="briggsRural"):
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
        wind_profile_type : str, gets either "HotSpot" or "default".
        sigmaTypeName : Name of the sigma type, for example from Briggs, rural/urban.

        Returns
        -------
        An instance of the class gadCloud

        """
        sigmaType = self.getSigmaType(sigmaTypeName)
        gascloud = abstractGasCloud.createGasCloud(sourceQ=sourceQ,sourceHeight=sourceHeight,
                                                   initialCloudSize=initialCloudSize,meteorology=meteorology,
                                                   wind_profile_type=wind_profile_type, spaceTime=spaceTime,
                                                   sigmaType=sigmaType ,deposition_velocity=deposition_velocity)
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
        plt.title(f"Detector at x={x}[m], y={y}[m], z={z}[m].")
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
        plt.ylabel(f"TIAC over time {units}")
        plt.title(f"Detector at x={x}[m], y={y}[m], z={z}[m].")
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

        plt.ylabel(rf"Concentration $\left[{C.attrs['Q'].units:~L}\right]$")
        plt.title(f"Maximum concentration over time. y={y}[m], z={z}[m]")
        plt.grid()

        x_min = (C.x[0].values)*ureg.meter if x_min is None else x_min
        x_max = (C.x[-1].values)*ureg.meter if x_max is None else x_max
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
        plt.ylabel(rf"Concentration $\left[{C.attrs['Q'].units:~L}\right]$")
        plt.title(f"Detector at x={x}[m], y={y}[m], z={z}[m].")
        plt.grid()
        t_min = (C.time[0].values)*ureg.min if t_min is None else t_min
        t_max = (C.time[-1].values)*ureg.min if t_max is None else t_max
        t_min = unumToPint(t_min).m_as(ureg.min)
        t_max = unumToPint(t_max).m_as(ureg.min)
        plt.xlim(t_min, t_max)
        plt.show()



    def plotTIACPerDistance_noQ(self,TIAC, y, z, time, x_min=None,x_max=None):
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

        plt.ylabel(rf"TIAC $\left[{TIAC.attrs['Q'].units:~L}\right]$")
        plt.title(f"TIAC per distance. y={y}[m], z={z}[m], time={time}[min]")
        plt.grid()
        x_min = (TIAC.x[0].values)*ureg.meter if x_min is None else x_min
        x_max = (TIAC.x[-1].values)*ureg.meter if x_max is None else x_max
        x_min = unumToPint(x_min).m_as(ureg.m)
        x_max = unumToPint(x_max).m_as(ureg.m)
        plt.xlim(x_min, x_max)
        plt.show()


    def plotFixedPointTIACOverTime_noQ(self, TIAC, x, y, z, t_min=None,t_max=None):
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
        plt.ylabel(rf"TIAC $\left[{TIAC.attrs['Q'].units:~L}\right]$")
        plt.title(f"Detector at x={x}[m], y={y}[m], z={z}[m].")
        plt.grid()

        t_min = (TIAC.time[0].values)*ureg.min if t_min is None else t_min
        t_max = (TIAC.time[-1].values)*ureg.min if t_max is None else t_max
        t_min = unumToPint(t_min).m_as(ureg.min)
        t_max = unumToPint(t_max).m_as(ureg.min)
        plt.xlim(t_min, t_max)
        plt.show()




