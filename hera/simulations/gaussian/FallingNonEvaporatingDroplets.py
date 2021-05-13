import numpy
import pandas
try:
	from scipy.integrate import solve_ivp
except ImportError:
	print("wrong scipy, can't solve the system")

from unum.units import *
from ..utils import toUnum,toNumber
from scipy.optimize import root
from scipy.constants import g as gconst
from .Meteorology import meteorologyFactory
from .Sigma import briggsRural

from scipy.constants import g as gconst

class FallingNonEvaporatingDroplets(object):
    """
        Gaussian model for falling non Evaporating Droplets.

        Here we assume that the particles are in the terminal velocity.
        We neglect the exit velocity.

    """
    _particleDiameter = None

    _position = None
    _meteorology = None

    _rho_p  = None

    _particleMass   = None

    _cloudQ = None # The total mass in the cloud.
    _cloudSigma = None


    @property
    def beta(self):
        return 5.

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self,value):
        try:
            if len(value) != 3:
                raise ValueError("Must be a 3 tuple")
        except:
                raise ValueError("Must be a 3 tuple. ")

        self._position=[toUnum(x,m) for x in value]

    @property
    def cloudSigma(self):
        return self._cloudSigma

    @cloudSigma.setter
    def cloudSigma(self, value):
        try:
            if len(value) != 3:
                raise ValueError("Must be a 3 tuple")
        except:
            raise ValueError("Must be a 3 tuple. ")

        self._cloudSigma = [toUnum(x, m) for x in value]

    @property
    def cloudSigmaX(self):
        return self._cloudSigma[0]

    @property
    def cloudSigmaY(self):
        return self._cloudSigma[1]

    @property
    def cloudSigmaZ(self):
        return self._cloudSigma[2]

    @property
    def Q(self):
        return self._cloudQ
    @Q.setter
    def Q(self,value):
        self._cloudQ = toUnum(value,kg)

    @property
    def particleMass(self):
        return self._particleMass

    @property
    def meteorology(self):
        return self._meteorology

    @property
    def x(self):
        return self._position[0]

    @property
    def y(self):
        return self._position[1]

    @property
    def z(self):
        return self._position[2]

    @property
    def particleDiameter(self):
        return self._particleDiameter

    @particleDiameter.setter
    def particleDiameter(self,value):
        self._particleDiameter = toUnum(value,um)
        self.__setParticleMass()

    @property
    def rho_p(self):
        return self._rho_p

    @rho_p.setter
    def rho_p(self, value):
        self._rho_p = toUnum(value, g / cm ** 3)
        self.__setParticleMass()

    @property 
    def N(self): 
        """
           The total number of particles.  
        """
        return (self.Q/self.particleMass).asNumber()

    @property
    def SpreadFactor(self): 
        """
          The spread factor of the droplets. Assuming oil and taken from the VX parameters. 
        """        
        return 4.5

    @property
    def g(self):
        return gconst*m/s**2

    @property
    def AreaOnSurface(self): 
        """
          The total surface area that the droplet occupies after its spread. 
            Q                                    3 * Q * SpreadFactor 
          ------ * pi*(d/2)**2*SpreadFactor =   --------------------- 
           Mass                                      2 * rho_l*d
          
        """
        return (self.Q*self.SpreadFactor*numpy.pi*(self.particleDiameter/2.)**2)/self.particleMass

    def __setParticleMass(self):
        if (self.rho_p is None or self.particleDiameter is None):
            return
        self._particleMass = (1 / 6. * self.rho_p * numpy.pi * self.particleDiameter ** 3).asUnit(kg)

    def __init__(self,particleDiameter,Q=1*kg,position=(0*m,0*m,0*m),meteorologyName="StandardMeteorolgyConstant",**kwargs):
        """
            Initializes the particle cloud.
            Currently initialized to point source.

	    See Aroesty - Atmospheric diffusion of droplet clouds for details. 

	    Plume (Aroesty pg 19): Plume diffusion refers to circumstances where material release and smaplling times are long compared with travel time from the source. 
	    Puff  		     : material release and samplling times are short compared with travel time. 

		Since the fall even from 500m is O(min) ~ O(release time) than we are always in the plume regime. 


        :param particleDiameter:
                particle diameter.
                default unit [um].
        :param Q:
                The total mass of the cloud.
                default units [kg]
        :param cloudDimensions:
                The initial cloud dimensions.

        :param position:
                A 3 tuple of unum length  that indicates the initial position of the cloud in space.
                default unit [m].
        :param meteorology:
                The name meteorology class to use

	:param dispersionType: 
		The dispersion is a plume or a puff.
		
		The differance is taken from Aerosty - Atmospheric diffusion of droplet clouds. 
		
		
		
		According to Aerosty if it is a plume (release time 

        :param kwargs:
                * all Passed to the meteorology factory.

                * dragCoeffFunc: The name of the drag coefficient function.
                                 Default Kelbaliyev15.

                                 Note that IK is not a good approximation for Re > 1000.

                * correctionCloudFunc: The correction to the cloud sigma.
                                  Either Plume or Puff.

                * rho_l : The density of the liquid.
        """
        cloudSigma = (0 * m, 0 * m, 0 * m)
        self._meteorology = meteorologyFactory.getMeteorology(meteorologyName,**kwargs)

        dragCoeffFunc   = kwargs.get("dragCoeffFunc","Haugen")
        self._dragfunc = getattr(self,"_DragCoefficient_%s" % dragCoeffFunc.title())

        cloudCorrectionFunc = kwargs.get("correctionCloudFunc","Plume")
        self._correctionfunc = getattr(self,"correctionCloud_%s" % cloudCorrectionFunc.title())

        self.particleDiameter = particleDiameter
        self.rho_p = kwargs.get("rho_l", 0.9 * g / cm ** 3) # oil
        self.Q = Q
        self.cloudSigma = cloudSigma
        self.position   = position

    def getTerminalVelocity(self,height=None):
        Res = root(self._VTFunc,1,args=(height,))
        return Res.x[0]*m/s


    def _DragCoefficient_Ik(self,Re):
        """
            Taken from IK. (ask yehudaa for details.)

        :return:
            Drag coefficient.
        """

        if Re < 2500:
            if Re < 0.037:
                f = 1
            elif Re < 1.5:
                f = 1 + 0.000294018 - 0.140373135*Re + 0.021993847*Re**2
            elif (Re < 370):
                f = -0.051616279 + 3.38774356/Re**0.5-7.1500312/Re + 11.207882/Re**1.5-7.30089/Re**2
            else:
                # 370 < f < 2500:
                f = -0.0643804358 + 6.303699848/(Re**0.5) - 170.035459/Re + 3317.9984639/Re**1.5 - 22352.02823/Re**2

            CD = 24.0/(f*Re)
        else:
            if Re < 2905:
                CD = 0.21667354 + 1.44e-4*Re
            else:
                CD = 0.11174    + 1.80124e-4*Re

        return CD

    def _DragCoefficient_Kelbaliyev15(self, Re):
        """
            Return the Cd to the Re.

            Ref 15.  Kelbaliyev Drag Coefficients of Variously Shaped Solid Particles,Drops, and Bubbles
        """
        Cd = 0
        if Re < 0.01:
            Cd = 24. / Re
        elif Re <= 20:
            Cd = (24. / Re) * (1 + 0.1315 * Re ** (0.82 - 0.05 * numpy.log10(Re)))
        elif Re <= 260:
            Cd = (24. / Re) * (1 + 0.1935 * Re ** 0.6305)
        elif Re <= 1500:
            F = 1.6425 - 1.12421 * numpy.log10(Re) + 0.1558 * (numpy.log10(Re) ** 2)
            Cd = 10 ** F
        elif Re <= 1.2e4:
            F = -2.4571 + 2.55581 * numpy.log10(Re) - 0.92951 * (numpy.log10(Re) ** 2) + 0.10491 * (
                        numpy.log10(Re) ** 3)
            Cd = 10 ** F
        else:
            Cd = 0.44

        return Cd

    def _DragCoefficient_Fan(self,Re):
        """
            The drag coefficient.

            Taken from principles of gas-solid flows. Fan and Zhu (1998).

        :param Re:
        :return:
            Drag coefficient
        """
        if Re < 2:
            Cd = 24 / Re
        elif Re < 500:
            Cd = 18.5 / (Re ** 0.6)
        else:
            Cd = 0.44
        return Cd

    def _DragCoefficient_Haugen(self,Re):
        """
            Haugen 2010 - Particle impaction on a cylinder in a crossflow as function of Stokes and Reynolds numbers
            j fluid mesh 2010.

        :param Re:
        :return:
        """
        if Re > 1000:
            return 0.44
        else:
            return 24/Re *(1+ 0.15*Re**0.687)


    def _VTFunc(self,Vt,height=None):

        height = toUnum(height if height is not None else self.z,m)


        rho_air = self.meteorology.getAirDensity(height)
        nu_air  = self.meteorology.getAirDynamicViscosity(height)

        uVt = Vt[0]*m/s
        Re  = (uVt*self.particleDiameter/nu_air).asNumber()
        CD  = self._dragfunc(Re)

        GuessVt = 4 * gconst * (m/s**2) * (self.rho_p - rho_air) * self.particleDiameter ** 2 / (3 * nu_air * CD * Re)

        return (GuessVt - uVt).asNumber(m/s)

    def correctionCloud_None(self):
        return 1

    def correctionCloud_Plume(self):
        Ubar = numpy.mean([self.meteorology.getWindVeclocity(z) for z in numpy.arange(0, self.z.asNumber(m))])
        return ((1 + (self.beta * self.getTerminalVelocity() / Ubar) ** 2) ** -0.25).asNumber()  # taking beta=5

    def correctionCloud_Puff(self):
        Ubar = numpy.mean([self.meteorology.getWindVeclocity(z) for z in numpy.arange(0, self.z.asNumber(m))])
        return ((1 + (self.beta * self.getTerminalVelocity() / Ubar) ** 2) ** -0.5).asNumber()  # taking beta=5


    def solveToTime(self, T):
        """
            Solves the equations to the requested time.

        :param T:
                The time to solve to.
        :return:
                A pandas with the times and the columns:
                    - position
                    - Q
                    - sigma
                    - particle velocity
        """
        # y0: x,y,z,dist,u,w
        y0 = numpy.array([self.x.asNumber(m),
                          self.y.asNumber(m),
                          self.z.asNumber(m),
                          0,
                          0,
                          0])

        sol = solve_ivp(self._fallingParticle, [0, T] , y0=y0,events=hit_ground())
        ret = pandas.DataFrame({'t': sol.t, 'x': sol.y[0, :],'y' : sol.y[1,:], 'z': sol.y[2, :],'distance' : sol.y[3,:], 'u_p' : sol.y[4,:], 'w_p' : sol.y[5,:]})

        ret = ret.merge(briggsRural.getSigma(ret['distance'], self.meteorology.stability),
                        left_on='distance', right_on='distance')

        ## Followig aroesty: csanady plume correction or ... puff correction.
        correction = self._correctionfunc()
        ret = ret.assign(sigmaXCorrected=ret['sigmaX']*correction)\
           .assign(sigmaZCorrected=ret['sigmaZ']*correction)\
           .assign(Q=self.Q.asUnit(kg))\
           .assign(diameter=self.particleDiameter.asUnit(um))
        return ret


    def _fallingParticle(self, t, y):
        """
            solving
                dx/dt = u_p
                dz/dt = w_p
                ds/dt = (u_p-u)**2+w_p**2)**0.2
                dU/dt = Cd(Re)*A*rho_air/2* |U-U_p|(u-u_p)
                dw/dt = Cd(Re)*A*rho_air/2* |U-U_p|(w_p) + g

                dist - the distance of the flight.
                U = (u,w)
                U_p = (u_p,w_p)
        """
        x, yc, z, dist, u_p, w_p = y

        z = numpy.max(z,0)
        z = toUnum(z,m)

        u_p = toUnum(u_p,m/s)
        w_p = toUnum(w_p, m / s)

        try:
            U = self.meteorology.getWindVeclocity(z)
        except RuntimeWarning:
            print("warning",z)
        rho_air = self.meteorology.getAirDensity(z)
        nu_air = self.meteorology.getAirDynamicViscosity(z)

        # 1. Calculate |U-U_p|
        Uabs = ((u_p - U)**2 + w_p**2)**0.5

        # 2. Calculate Re
        Re = (rho_air*Uabs * self.particleDiameter / nu_air).asNumber()

        # 3. Calculate Cd
        Cd = self._dragfunc(Re)

        # 4. Calculate Coeff = rho*A/2 / mp
        #                    = rho/2 (pi * d**2/4) / (rho_p*pi*d**3/6) = rho * (2/3)/d = 2rho/(3d*rho_p).
        Coeff = 2*rho_air/(3*self.particleDiameter*self.rho_p)

        new_u_p = Cd*Coeff*Uabs*(U-u_p)
        new_w_p = -Cd*Coeff*Uabs*(w_p)- self.g

        newy = numpy.zeros(y.shape)
        newy[0] = u_p.asNumber(m/s)                   # x
        newy[1] = 0                                   # y
        newy[2] = w_p.asNumber(m/s)                   # z
        newy[3] = Uabs.asNumber(m/s)                  # dist - distance,
        newy[4] = new_u_p.asNumber()    # u_p
        newy[5] = new_w_p.asNumber()    # w_p
        return newy

class hit_ground(object):
    def __call__(self, t, y):
        return y[2]-1e-3
    terminal= True

    #particleDiameter, Q = 1 * kg, cloudSigma = (1 * m, 1 * m, 1 * m), position = (0 * m, 0 * m,0 * m), meteorologyName = "StandardMeteorolgyConstant", ** kwargs):

    #particle = FallingNonEvaporatingDroplets(particleDiameter=3*mm,position = (0 * m, 0 * m,300 * m))


