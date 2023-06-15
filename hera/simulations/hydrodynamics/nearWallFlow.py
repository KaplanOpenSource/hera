"""
    Algorithm for calculating the ustar in a channel flow.

    Taken from Boundary-Layer Theory 9th edition (schlichting) pg 537.

"""
from ...utils.unum import tounit
from scipy.optimize import fsolve
from unum.units import *
import numpy

class functionG:
    """
        The G(Lambda,D) is defined implicitly (Eqn. 17.60, page 537):

        \begin{equation}
            \frac{\Lambda}{G} + 2\ln\left(\frac{\Lambda}{G}\right) -D = \Lambda
        \end{equation}

    """
    def __init__(self):
        pass


    def _implicitG(self,x,Lambda,D):

        return Lambda/x + 2*numpy.log(Lambda/x) -D -x

    def solve(self,Lambda,D):
        """
            Solves the G function implicitely.

        :param Lambda:
        :param D:
        :return:
        """
        return fsolve(self._implicitG,1,args=(Lambda,D))[0]

class nearWallFlow:

    _Ra = None # The arithmetic mean deviations.

    _nu = None # kinematic viscosity of the fluid.

    @property
    def kinematicViscosity(self):
        return self._nu

    @property
    def technical_roughness(self):
        """
            Estimate the technical roughness from Ra (page 534,  Schlichting).

        :return:
        """
        return 3.5*self._Ra

    @property
    def height(self):
        return self._channelHeight.asNumber(m)

    @property
    def hydraulicHeight(self):
        return 2*self.height

    def technical_roughness_plus(self,ustar):
        """
            The dimensionless roughness .
        :param ustar:
        :return:
        """
        return (ustar*self.technical_roughness/self.kinematicViscosity).asNumber()


    def Cplus(self,ustar):
        """
            Cplus (C^+) is the correction to the velocity profile due to roughness. .

            The equation for C^+ in rough channel is given in Eqn. 17.40


        :return: float
            The C^+.
        """
        kappa = 0.41 # von karman constant.
        kplus = (ustar*self.technical_roughness/self.kinematicViscosity).asNumber()

        if kplus <= 5:
            ret = 5
        elif (kplus >5) and (kplus < 70):
            ret = 8-1/kappa*numpy.log(3.4+kplus)
        else:
            ret = 8

        return ret

    def ReynoldsUm(self,ustar,Um):
        """
            Return the reynolds based on the Um and hydraulic diameter.
        :return:
        """
        return ustar*self.hydraulicHeight/self.kinematicViscosity


    def __init__(self,Ra,nu):
        """

        Parameters
        ----------

        Ra: float/unit
            Arithmetic mean deviations of surface asperities. default unit [m]

        nu: float/ unit
            The viscosity of the fluid. default unit [m^2/s]
        """
        self._nu = tounit(nu,m**2/s)
        self._Ra = tounit(Ra,m)



class channelFlow(nearWallFlow):
    """
        Defines a list of function fr the anaytical management of rough cannel.
        See chapter 16 and 17. 

        Units are mks.

        This is the model 'flow' for the atmosphere.

    """
    def __init__(self,Ra,nu,channelHeight):
        """

        Parameters
        ----------

        Ra: float/unit
            Arithmetic mean deviations. default unit [m]

        nu: float/ unit
            The viscosity of the fluid. default unit [m^2/s]

        channelHeight: float/unit
            The height of the channel. default unit [m]

            Note that in schlichting H is channelHeight/2.

        """
        super().__init__(Ra,nu)
        self._channelHeight = tounit(channelHeight,m)
        self._functionG = functionG()
        self.C_bar_plus_C_barbar = -1.7 # Cbar + Cbar bar (equation 17.91).
        self.C_bar = 0.94  # Cbar + Cbar bar (equation 17.91).



    def skin_friction(self,ustar,Um):
        """
            Skin friction (c_f) is given by eqn 17.92.


        ustar: float/unum.
            The guess for the current Ustar. units [m**2/s]

        Um: float/unum
            The flow velocity at the center of the channel. (u_m).

        :return: float
            The skin friction.
        """
        kappa = 0.41  # von karman constant.
        Re_dh = self.ReynoldsUm(ustar, Um)
        Lambda = 2*numpy.log()
        Cplus  = self.Cplus(ustar)
        D      = 0.82*Cplus - 4.56
        G      = self._functionG.solve(Lambda,D)

        ret = 2*(kappa/numpy.log(Re_dh)*G)**2
        return ret


    def Re_tau(self,ustar):
        """
            Returns the Reynolds number with ustar and half channel height as measurements.old.

        :param ustar:
        :return:
        """
        ustar = tounit(ustar, m / s)
        return ((ustar*self._channelHeight/2)/self.kinematicViscosity).asNumber()


    def get_Umean_from_Ustar(self,ustar):
        """
            Calculate the U mean that is needed in order to obtain the requested Ustar.

        Parameters
        -----------

        ustar: float, unum
            The ustar, default units [m/s]

        Returns
        -------
        the velocity in the center of the channel (Um) [m/s]
        """
        kappa = 0.41  # von karman constant.

        ustar = tounit(ustar,m/s)

        Re_tau = self.Re_tau(ustar)

        Um_plus = 1/kappa*numpy.log(Re_tau) + self.Cplus(ustar) + self.C_bar_plus_C_barbar

        Um = Um_plus*ustar ## see Eqn. 17.86

        return Um

    def get_Ucenter_from_Ustar(self,ustar):
        """
            Calculate the U at the center of the channel that is needed in order to obtain the requested Ustar.

            Eqn 17.53.

        Parameters
        -----------

        ustar: float, unum
            The ustar, default units [m/s]

        Returns
        -------
        the velocity in the center of the channel (Um) [m/s]
        """
        kappa = 0.41  # von karman constant.

        ustar = tounit(ustar,m/s)

        Re_tau = self.Re_tau(ustar)

        Um_plus = 1/kappa*numpy.log(Re_tau) + self.Cplus(ustar) + self.C_bar

        Um = Um_plus*ustar ## see Eqn. 17.86

        return Um


class couetteFlow(nearWallFlow):
    """
        Defines a list of function fr the anaytical management of rough cannel.
        See chapter 16 and 17.

        Units are mks.

        This is the model 'flow' for indoors.

    """
    def __init__(self, Ra, nu, channelHeight):
        """

        Parameters
        ----------

        Ra: float/unit
            Arithmetic mean deviations. default unit [m]

        nu: float/ unit
            The viscosity of the fluid. default unit [m^2/s]

        channelHeight: float/unit
            The height of the channel. default unit [m]

            Note that in schlichting H is channelHeight/2.

        """
        super().__init__(Ra,nu)
        self._channelHeight = tounit(channelHeight, m)
        self._functionG = functionG()


    def skin_friction(self, ustar, Um):
        """
            Skin friction (c_f) is given by eqn 17.92.


        ustar: float/unum.
            The guess for the current Ustar. units [m**2/s]

        Um: float/unum
            The flow velocity at the center of the channel. (u_m).

        :return: float
            The skin friction.
        """
        kappa = 0.41  # von karman constant.
        Re_dh = self.ReynoldsUm(ustar, Um)
        Lambda = 2 * numpy.log()
        Cplus = self.Cplus(ustar)
        D = 0.82 * Cplus - 4.56
        G = self._functionG.solve(Lambda, D)

        ret = 2 * (kappa / numpy.log(Re_dh) * G) ** 2
        return ret

    def Re_tau(self, ustar):
        """
            Returns the Reynolds number with ustar and half channel height as measurements.old.

        :param ustar:
        :return:
        """
        ustar = tounit(ustar, m / s)
        return ((ustar * self._channelHeight / 2) / self.kinematicViscosity).asNumber()

    def get_Um_from_Ustar(self, ustar):
        """
            Calculate the Um that is needed in order to obtain the requested Ustar.

            17.21

        Parameters
        -----------

        ustar: float, unum
            The ustar, default units [m/s]

        Returns
        -------
        the velocity in the center of the channel (Um) [m/s]
        """
        kappa = 0.41  # von karman constant.

        ustar = tounit(ustar, m / s)

        Re_tau = self.Re_tau(ustar)

        Um_plus = 1 / kappa * numpy.log(Re_tau) + self.Cplus(ustar)

        Um = Um_plus * ustar  ## see Eqn. 17.86

        return Um
