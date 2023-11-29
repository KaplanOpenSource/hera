from unum.units import *
import pandas
import numpy
from scipy.stats import lognorm
from ...utils import tounit,tonumber
from .FallingNonEvaporatingDroplets import FallingNonEvaporatingDroplets

class FixedPositionDropletsCloud(object):
    """
        Holds a list of FallingNonEvaporatingDroplets
        that were created using the lognormal ditribution.

    """

    _dropletList = None

    @property
    def dropletList(self):
        return self._dropletList

    def __init__(self, mmd, geometricstd, position, Q, clouds=30, meteorologyname="logNormal",met_kwargs={}, **kwargs):
        """
            Creates a list of clouds (discretization according to the number of clouds).

        :param mmd:
                The mmd of the droplets.
        :param geometricstd:
                The particle distribution geometric std.
        :param position:
                The initial position of the cloud.
        :param Q:
                The Q of the cloud.
                default units [kg]
        :param clouds:
            The number of clouds to generate.
        :param meteorologyname:
            The name of the meteorology
        :param kwargs:
            Parameters to pass to the droplets.
        """
        self._dropletList = []
        self._initDropletPosition(mmd, geometricstd, position, Q, clouds=clouds, meteorologyname=meteorologyname,met_kwargs=met_kwargs, **kwargs)

    def getGround(self,T):
        """
     		Returns the ground concentration using the pancake model (Aroesty) 
                Also return the total number of particles for each class (N) 
		and calculate their relative surface area. 

		
        """
        retList = []
        #print("Total of %s dropletClouds" % len(self._dropletList))
        total = len(self._dropletList)
        for i,droplet in enumerate(self._dropletList):
            print(f"solving droplets {i}/{total}" ,end="\r")
            res = droplet.solveToTime(T)
            res['N'] = droplet.N 
            res['DropletArea'] = droplet.AreaOnSurface.asNumber(um**2)
            retList.append(res.iloc[-1].to_frame().T)

        ret = pandas.concat(retList, ignore_index=True)
        return ret


    def _initDropletPosition(self, mmd, geometricstd, position, Q, clouds=30, meteorologyname="logNormal",met_kwargs={}, **kwargs):

        rv = lognorm(numpy.log(geometricstd), scale=tonumber(mmd, m))
        lower = rv.ppf(1e-4)
        upper = rv.ppf(1-1e-4)

        met_params = dict(z0=10*cm)
        met_params.update(met_kwargs)

        interval = numpy.logspace(numpy.log(lower),numpy.log(upper),clouds,base=numpy.e)
        dh = numpy.diff(numpy.log(interval))[0]

        massFractionVector = numpy.diff(rv.cdf(interval))*tonumber(Q,kg)
        diameterVector     = numpy.exp(numpy.log(interval[:-1])+dh/2.)  # [m]

        for dropletDiam,dropletQ in zip(diameterVector,massFractionVector):
            droplets = FallingNonEvaporatingDroplets(particleDiameter=dropletDiam*m,
                                                     Q=dropletQ*kg,
                                                     position=position,
                                                     meteorologyName=meteorologyname,
                                                     met_kwargs=met_params,
                                                     **kwargs)
            self._dropletList.append(droplets)


class LinePositionDropletsCloud(FixedPositionDropletsCloud):
    """
        A line across the wind.
    """
    def __init__(self, mmd, geometricstd, position, Q, linelength, clouds=30,linepositions=100, meteorologyname="StandardMeteorolgyConstant", **kwargs):
        """
            Creates a list of clouds (discretization according to the number of clouds).

        :param mmd:
                The mmd of the droplets.
        :param geometricstd:
                The particle distribution geometric std.
        :param position:
                The initial position of the cloud.
        :param Q:
                The Q of the cloud.
                default units [kg]
        :param linelength:
                The length of the line long with the initial clouds are dispersed.
        :param clouds:
            The number of clouds to generate.
        :param linepositions:
            The number of clouds along the line to generate.

        :param meteorologyname:
            The name of the meteorology
        :param kwargs:
            Parameters to pass to the droplets.
        """
        self._dropletList = []

        qCloud = Q/linepositions
        for Ypos in numpy.linspace(0,tonumber(linelength,m),linepositions):
            curpos = (position[0],tounit(Ypos,m),position[2])
            self._initDropletPosition(mmd, geometricstd, curpos, qCloud, clouds=clouds,meteorologyname=meteorologyname, **kwargs)



class FixedPointClippedDropletCloud(FixedPositionDropletsCloud):

    _clippedDiameter = None

    @property
    def clippedDiameter(self):
        return self._clippedDiameter

    @clippedDiameter.setter
    def clippedDiameter(self, value):

        if isinstance(value,str):
            self._clippedDiameter = eval(value)
        else:
            self._clippedDiameter = tounit(value,mm)


    def __init__(self,mmd,geometricstd,position,Q,clippedDiameter,clouds=30, meteorologyname="StandardMeteorolgyConstant",**kwargs):

        super().__init__(mmd=mmd,
                         geometricstd=geometricstd,
                         position=position,
                         Q=Q,
                         clouds=clouds,
                         meteorologyname=meteorologyname
                         )

        self.clippedDiameter = clippedDiameter


    def _initDropletPosition(self, mmd, geometricstd, position, Q, clouds=30, meteorologyname="StandardMeteorolgyConstant", **kwargs):

        clippedDiameter = tonumber(self.clippedDiameter,m)

        rv = lognorm(numpy.log(geometricstd), scale=tonumber(mmd, m))
        lower = rv.ppf(1e-4)
        upper = clippedDiameter

        interval = numpy.logspace(numpy.log(lower),numpy.log(upper),clouds,base=numpy.e)
        dh = numpy.diff(numpy.log(interval))[0]

        maxMass = rv.cdf(interval[-1])

        massFractionVector = numpy.diff(rv.cdf(interval))*tonumber(Q,kg)
        diameterVector     = numpy.exp(numpy.log(interval[:-1])+dh/2.)  # [m]

        massFractionVector /= maxMass

        for dropletDiam,dropletQ in zip(diameterVector,massFractionVector):
            droplets = FallingNonEvaporatingDroplets(particleDiameter=dropletDiam*m, Q=dropletQ*kg, position=position, meteorologyName=meteorologyname, **kwargs)
            self._dropletList.append(droplets)



class CirclePositionClippedDropletsCloud(FixedPointClippedDropletCloud):
    """
        A line across the wind.
    """
    def __init__(self, mmd, geometricstd, position, Q, clippedDiameter, clouds=30,radius=10*m,circlepositions=4, meteorologyname="StandardMeteorolgyConstant", **kwargs):
        """
            Creates a list of clouds (discretization according to the number of clouds).

        :param mmd:
                The mmd of the droplets.
        :param geometricstd:
                The particle distribution geometric std.
        :param position:
                The initial position of the cloud.
        :param Q:
                The Q of the cloud.
                default units [kg]
        :param linelength:
                The length of the line long with the initial clouds are dispersed.
        :param clouds:
            The number of clouds to generate.
        :param linepositions:
            The number of clouds along the line to generate.

        :param meteorologyname:
            The name of the meteorology
        :param kwargs:
            Parameters to pass to the droplets.
        """
        self._dropletList = []

        self.clippedDiameter = clippedDiameter
        radius_m = tounit(radius,m)

        qCloud = Q/circlepositions
        for deg in numpy.linspace(0,2*numpy.pi,circlepositions):

            if radius is None:
                curpos = position
            else:
                curpos = (position[0]+ radius_m*numpy.cos(deg),position[0]+ radius_m*numpy.sin(deg),position[2])

            self._initDropletPosition(mmd, geometricstd, curpos, qCloud, clouds=clouds,meteorologyname=meteorologyname, **kwargs)

