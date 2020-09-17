import pandas
import numpy
from ..utils import *
import pydoc
from ..gaussian import Meteorology as metmodule

class StandardMeteorolgyConstant(object):
    """
        Implements a standard meteorology.

    """
    _temperature = None
    _stability   = None  # stability, for the wind profile.
    _z0          = None  # roughness, for the wind profile.
    _wind_p      = None  # the calculated roughness.
    _wind_mathematical_angle = None
    _ustar       = None

    _refHeight = None
    _u_refHeight = None

    @property
    def ustar(self):
        return self._ustar

    @ustar.setter
    def ustar(self,value):
        self._ustar = toUnum(value,m/s)


    @property
    def refHeight(self):
        return self._refHeight

    @refHeight.setter
    def refHeight(self,value):
        self._refHeight = toUnum(value,m)


    @property
    def wind_p(self):
        return self._wind_p

    @property
    def u10(self):
        return self._u_refHeight

    @u10.setter
    def u10(self,value):
        self._u_refHeight = toUnum(value,m/s)
        self._refHeight = 10

    @property
    def u_refHeight(self):
        return self._u_refHeight

    @u_refHeight.setter
    def u_refHeight(self,value):
        self._u_refHeight = toUnum(value,m/s)

    @property
    def stability(self):
        return self._stability

    @stability.setter
    def stability(self,value):
        try:
            value = value.upper()
            if value not in ["A","B","C","D","E","F"]:
                raise ValueError("Stability %s is not valid. should be a letter from A to F")
        except:
                raise ValueError("Stability %s is not valid. shoudl be a letter from A to F")
        self._stability = value
        self._setPvalues()

    @property
    def z0(self):
        return self._z0

    @z0.setter
    def z0(self,value):
        self._z0 = toUnum(value,m)
        self._setPvalues()

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self,value):
        self._temperature = toUnum(value,celsius)

    @property
    def skinSurfaceTemperature(self):
        return self._skinSurfaceTemperature

    @skinSurfaceTemperature.setter
    def skinSurfaceTemperature(self,value):
        self._skinSurfaceTemperature = toUnum(value, celsius)



    #================================  Wind profile calcluation
    ## Calculate the wind profile based on stability and z0.
    #  These are the const. See getWindP for details.
    _pvalues = pandas.DataFrame({
                        "A": [0, 0.05, 0.08, 0.17, 0.27],
                        "B": [0, 0.06, 0.09, 0.17, 0.28],
                        "C": [0, 0.06, 0.11, 0.2 , 0.31],
                        "D": [0, 0.12, 0.16, 0.27, 0.37],
                        "E": [0, 0.34, 0.32, 0.38, 0.47],
                        "F": [0, 0.53, 0.54, 0.61, 0.69]}, index=[0,0.01,0.1,1,3]) # roughness, in [m]


    def __init__(self,**kwargs):
        """
            Define the base parameters of the meteorology

        :param kwargs:
            temperature - The temperature on the ground. [C]
                          default is 19C

            z0          - units are [m] unless specified otherwise.
                        The default is 10cm.
            stability   - The default is stability D.

            u/refHeight - The wind velocity at refHeight.
                          if not specified, check if u10 is present.

            u10         - The wind velocity at 10m.
                          default units are m/s.
                          default value is 4m/s

            ustar       - The u* of the meteorological conditions.
                          the default value is 30cm/s

            skinSurfaceTemperature - The temperature of the surface skin.
                                     The default value is 30C
        """
        self.temperature = kwargs.get("temperature",20)
        self.stability = kwargs.get("stability", "D")
        self.z0        = kwargs.get("z0", 0.1) # 10cm

        if "u" in kwargs:
            if "refHeight" not in kwargs:
                raise ValueError("Can not set u and not send refHeight")
            else:
                self.refHeight = kwargs['refHeight']
                self.u_refHeight = kwargs['u']
        else:
            self.u10       = kwargs.get("u10", 4) # 10m/s

        self.ustar     = kwargs.get("ustar",0.3) #
        self.skinSurfaceTemperature = kwargs.get("skinSurfaceTemperature",35) # day

    def getAirPressure(self,height):
        """
            Return the air pressure at requested height.

        :param height:
            The height.
            Default unit [m].

        :return:
            The air pressure at mmHg units.
        """
        height = toUnum(height,m)

        return 760.*numpy.exp(-1.186e-4*height.asNumber(m))*mmHg

    def getTKE(self,height):
        """
            A simplistic model for TKE in the atmosphere.

            Currently implemented only neutral conditions.

        :param height:
        :return:
            the tke [m^2/s^2].
        """
        return (3.25*self.ustar)**2.

    def getAirTemperature(self,height):
        """
            Return the air temperature.

            The air temperature drops at 6.5C/km

        :param height:
                default unit m
        :return:
            air temperature at C.
        """
        return self._temperature - 6.5e-3*toNumber(height,m)*celsius


    def getAirDensity(self,height):
        """
            Calculate the air density

            \begin{equation}
                \rho_{air} =  \frac{1.701316e-6*P}{(1+0.00367*T)}*[g\cdot cm^{-3}]
            \end{equation}

            Where P [mmHg] is the air pressure and T [celsius] is the temperature.

            Then convert to the mks system.

        :param height:
               default m
        :return:
               air density in kg/m**3
        """
        P = self.getAirPressure(height).asNumber(mmHg)
        T = self.getAirTemperature(height).asNumber(celsius)

        density = 1.701316e-6*P/(1+0.00367*T)*g/cm**3
        return density.asUnit(kg/m**3)


    def getAirDynamicViscosity(self,height):
        """
            Calculate the dynamic viscosity

            \begin{equation}
                \nu = 1e-6*(170.27 + 0.911409*T - 0.00786742*T**2)*[dyne*s/cm**2]
            \end{equation}

        :param height:
                The height in [m]
        :return:
                The air viscosity in dyne/sec/m**2.
        """
        T = self.getAirTemperature(height).asNumber(celsius)
        return (1e-6*(170.27 + 0.911409*T - 0.00786742*T**2)*dyne*s/cm**2).asUnit(dyne*s/m**2)


    def _setPvalues(self):
        """
            Return the p-values (exponent coefficient of the wind profile).
            Taken from Irwin JS "A theoretical variation of the wind profile power-law exponent as a function of surface roughness and stability" 1984.

        :return:
            The coefficient (dimensionless).
        """
        if (self.z0 is None or self.stability is None):
            return
        pstab        = self._pvalues[self.stability]
        self._wind_p = numpy.interp(self.z0.asNumber(m),pstab.index,pstab)

    def getWindVeclocity(self,height):
        """
            Return the wind velocity defined as:
            \begin{equation}
                    u(z) = u_{x [m]}\cdot \left(\frac{height [m]}{x [m]}\right)^{pconst}
            \end{equation}

            where pconst is calculated with getPvalues.

        :param height:
                default units [m]
        :return:
            The wind velocity at the requested height.
        """
        height = toNumber(height,m)
        refHeight = toNumber(self.refHeight,m)
        height = numpy.min([numpy.max([height,0]),300])

        return self.u_refHeight*(height/refHeight)**self.wind_p



class MeteorologyProfile(StandardMeteorolgyConstant):
    """
        Gets a profile of the wind velocity and the wind direction.

    """
    pass





#################################################################################################
#                                Factory                                                        #
#################################################################################################

class MeteorolgyFactory(object):

        def getMeteorology(self,name,**kwargs):
            """
               Creating a meteorology object.

            :param kwargs:
                    name: The meteorology object name.

                    Other kwparams are passed to the meteorology object.

            :return:
                The meteorology object.
            """
            return getattr(metmodule,name)(**kwargs)



meteorologyFactory = MeteorolgyFactory()

