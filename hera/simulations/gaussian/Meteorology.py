import pandas
import numpy
from hera.utils import *
from hera.utils.unitHandler import ureg, unumToPint
# from hera.utils import tounit, tonumber


class StandardMeteorolgyConstant_powerLaw:
    """
        Implements a standard meteorology.

    """
    _temperature = None
    _stability = None  # stability, for the wind profile.
    _z0 = None  # roughness, for the wind profile.
    _wind_p = None  # the calculated roughness.
    _wind_mathematical_angle = None
    _ustar = None

    _refHeight = None
    _u_refHeight = None

    # inversion = None
    # @property
    # def inversion(self):
    #     return self.inversion

    @property
    def ustar(self):
        return self._ustar

    @ustar.setter
    def ustar(self, value):
        self._ustar = tounit(value, ureg.m / ureg.s)

    @property
    def refHeight(self):
        return self._refHeight

    @refHeight.setter
    def refHeight(self, value):
        self._refHeight = tounit(value, ureg.m)

    @property
    def wind_p(self):
        return self._wind_p

    @property
    def u10(self):
        return self._u_refHeight

    @u10.setter
    def u10(self, value):
        self._u_refHeight = tounit(value, ureg.m / ureg.s)
        self._refHeight = 10

    @property
    def u_refHeight(self):
        return self._u_refHeight

    @u_refHeight.setter
    def u_refHeight(self, value):
        self._u_refHeight = tounit(value, ureg.m / ureg.s)

    @property
    def stability(self):
        return self._stability

    @stability.setter
    def stability(self, value):
        try:
            value = value.upper()
            if value not in ["A", "B", "C", "D", "E", "F"]:
                raise ValueError("Stability %s is not valid. should be a letter from A to F")
        except:
            raise ValueError("Stability %s is not valid. shoudl be a letter from A to F")
        self._stability = value
        self._setPvalues()

    @property
    def z0(self):
        return self._z0

    @z0.setter
    def z0(self, value):
        self._z0 = tounit(value, ureg.m)
        self._setPvalues()

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self, value):
        self._temperature = tounit(value, ureg.K)

    @property
    def skinSurfaceTemperature(self):
        return self._skinSurfaceTemperature

    @skinSurfaceTemperature.setter
    def skinSurfaceTemperature(self, value):
        self._skinSurfaceTemperature = tounit(value, ureg.K)

    # ================================  Wind profile calcluation
    ## Calculate the wind profile based on stability and z0.
    #  These are the const. See getWindP for details.
    _pvalues = pandas.DataFrame({
        "A": [0, 0.05, 0.08, 0.17, 0.27],
        "B": [0, 0.06, 0.09, 0.17, 0.28],
        "C": [0, 0.06, 0.11, 0.2, 0.31],
        "D": [0, 0.12, 0.16, 0.27, 0.37],
        "E": [0, 0.34, 0.32, 0.38, 0.47],
        "F": [0, 0.53, 0.54, 0.61, 0.69]}, index=[0, 0.01, 0.1, 1, 3])  # roughness, in [m]

    def __init__(self, inversion, temperature=20, stability="D", z0=0.1, ustar=0.3, skinSurfaceTemperature=35, **kwargs):
        """
            Define the base parameters of the meteorology

            temperature - The temperature on the ground. [C]
                          default is 19C

            z0          - units are [m] unless specified otherwise.
                        The default is 10cm.
            stability   - The default is stability D.


            ustar       - The u* of the meteorology conditions.
                          the default value is 30cm/s

            skinSurfaceTemperature - The temperature of the surface skin.
                                     The default value is 30C

            kwargs:
                u/refHeight - The wind velocity at refHeight.
                              if not specified, check if u10 is present.

                u10         - The wind velocity at 10m.
                              default units are m/s.
                              default value is 4m/s

        """
        self.temperature = temperature
        self.stability = stability
        self.z0 = z0
        self.inversion = inversion

        if "u" in kwargs:
            if "refHeight" not in kwargs:
                raise ValueError("Can not set u and not send refHeight")
            else:
                self.refHeight = kwargs['refHeight']
                self.u_refHeight = kwargs['u']
        else:
            self.u10 = kwargs.get("u10", 4)  # 10m/s

        self.ustar = ustar
        self.skinSurfaceTemperature = skinSurfaceTemperature  # day

    def getAirPressure(self, height):
        """
            Return the air pressure at requested height.

        :param height:
            The height.
            Default unit [m].

        :return:
            The air pressure at mmHg units.
        """
        height = unumToPint(height).to(ureg.m)

        return 760. * numpy.exp(-1.186e-4 * height.m_as(ureg.m)) * ureg.mmHg

    def getTKE(self, height):
        """
            A simplistic model for TKE in the atmosphere.

            Currently implemented only neutral conditions.

        :param height:
        :return:
            the tke [m^2/s^2].
        """
        return (3.25 * self.ustar) ** 2.

    def getAirTemperature(self, height):
        """
            Return the air temperature.

            The air temperature drops at 6.5C/km

        :param height:
                default unit m
        :return:
            air temperature at C.
        r"""
        return self._temperature - 6.5e-3 * tonumber(height, ureg.m) * ureg.K

    def getAirDensity(self, height):
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
        r"""
        P = unumToPint(self.getAirPressure(height)).m_as(ureg.mmHg)
        T = unumToPint(self.getAirTemperature(height)).m_as(ureg.degC)

        density = 1.701316e-6 * P / (1 + 0.00367 * T) * ureg.g / ureg.cm ** 3
        return density.to(ureg.kg / ureg.m ** 3)

    def getAirDynamicViscosity(self, height):
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
        T = unumToPint(self.getAirTemperature(height)).m_as(ureg.degC)
        return (1e-6 * (170.27 + 0.911409 * T - 0.00786742 * T ** 2) * ureg.dyne * ureg.s / ureg.cm ** 2).to(ureg.dyne * ureg.s / ureg.m ** 2)

    def _setPvalues(self):
        """
            Return the p-values (exponent coefficient of the wind profile).
            Taken from Irwin JS "A theoretical variation of the wind profile power-law exponent as a function of surface roughness and stability" 1984.

        :return:
            The coefficient (dimensionless).
        r"""
        if (self.z0 is None or self.stability is None):
            return
        pstab = self._pvalues[self.stability]
        self._wind_p = numpy.interp(unumToPint(self.z0).m_as(ureg.m), pstab.index, pstab)

    def getWindVelocity(self, height):
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
        r"""
        height = tonumber(height, ureg.m)
        refHeight = tonumber(self.refHeight, ureg.m)
        height = numpy.min([numpy.max([height, 0]), 300])

        return self.u_refHeight * (height / refHeight) ** self.wind_p


    def getWindVelocity_hotSpot(self, height):
        facfor_dict = dict(A=0.07, B=0.07, C=0.1, D=0.15, E=0.35, F=0.55)
        u_H = self.u10*(height/(10*ureg.m))**facfor_dict[self.stability]  # The value 10*m in the denominator is the reference height for u10.
        # u_H_float = round(u_H.m_as(ureg.m/ureg.s), 2)
        # ret = u_H_float*ureg.m/ureg.s
        return u_H

class StandardMeteorolgyConstant_log(StandardMeteorolgyConstant_powerLaw):

    def getWindVelocity(self, height):
        """
            Return the wind velocity defined as:
            \begin{equation}
                    u(z) = u_{x [m]}\cdot log(\frac{height [m]}{z_0 [m]}\right)
            \end{equation}

            where pconst is calculated with getPvalues.

        :param height:
                default units [m]
        :return:
            The wind velocity at the requested height.
        """

        z0 = tonumber(self.z0, ureg.m)
        height = tonumber(height, ureg.m)
        height = numpy.min([numpy.max([height, 0]), 300])
        refHeight = tonumber(self.refHeight, ureg.m)
        ustar_over_kappa = self.u_refHeight / numpy.log(refHeight / z0)

        if (height <= z0):
            u = 0 * ureg.m / ureg.s
        else:
            u = ustar_over_kappa * numpy.log(height / z0)

        return u


class StandardMeteorolgyConstant_uniformWind(StandardMeteorolgyConstant_powerLaw):

    def getWindVelocity(self, height):
        """
            Constant wind
        :param height:
                default units [m]
        :return:
            The wind velocity at the requested height.
        """
        return self.u_refHeight


class MeteorologyProfile(StandardMeteorolgyConstant_powerLaw):
    """
        Gets a profile of the wind velocity and the wind direction.

    """
    pass


#################################################################################################
#                                Factory                                                        #
#################################################################################################

class MeteorologyFactory:

    def __init__(self):
        self.meteorology = dict(powerLaw=StandardMeteorolgyConstant_powerLaw,
                                log=StandardMeteorolgyConstant_log,
                                uniformWind=StandardMeteorolgyConstant_uniformWind)

    def getMeteorologyFromU10(self, u10, inversion, verticalProfileType="log", temperature=ureg.Quantity(20, ureg.degC), stability="D", z0=0.1*ureg.m, ustar=0.3*ureg.m/ureg.s, skinSurfaceTemperature=ureg.Quantity(35, ureg.degC)):
        """
           Creating a meteorology object.

        :param kwargs:
                name: The meteorology object name.

                Other kwparams are passed to the meteorology object.

        :return:
            The meteorology object.
        """
        return self.meteorology[verticalProfileType](u10=u10, inversion=inversion, temperature=temperature, stability=stability,
                                                     z0=z0, ustar=ustar, skinSurfaceTemperature=skinSurfaceTemperature)

    def getMeteorologyFromURefHeight(self, u, inversion, refHeight, verticalProfileType="log", temperature=ureg.Quantity(20, ureg.degC), stability="D",
                                     z0=0.1*ureg.m, ustar=0.3*ureg.m/ureg.s, skinSurfaceTemperature=ureg.Quantity(35, ureg.degC)):
        """
           Creating a meteorology object.

        :param kwargs:
                name: The meteorology object name.

                Other kwparams are passed to the meteorology object.

        :return:
            The meteorology object.
        """
        return self.meteorology[verticalProfileType](u=u, inversion=inversion, refHeight=refHeight, temperature=temperature,
                                                     stability=stability, z0=z0, ustar=ustar, skinSurfaceTemperature=skinSurfaceTemperature)


