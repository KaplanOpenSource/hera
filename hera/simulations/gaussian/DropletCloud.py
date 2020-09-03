from unum.units import *
import pandas
import numpy
from scipy.stats import lognorm
from ..utils import toUnum,toNumber
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

    def __init__(self, mmd, geometricstd, position, Q, clouds=30, meteorologyname="StandardMeteorolgyConstant", **kwargs):
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
        self._initDropletPosition(mmd, geometricstd, position, Q, clouds=clouds, meteorologyname=meteorologyname, **kwargs)

    def getGround(self,T):
        """
     		Returns the ground concentration using the pancake model (Aroesty) 
                Also return the total number of particles for each class (N) 
		and calculate their relative surface area. 

		
        """
        retList = []
        #print("Total of %s dropletClouds" % len(self._dropletList))
        for i,droplet in enumerate(self._dropletList):
            #print("solving droplets %s" %i )
            res = droplet.solveToTime(T)
            res['N'] = droplet.N 
            res['DropletArea'] = droplet.AreaOnSurface
            retList.append(res.iloc[-1].to_frame().T)

        ret = pandas.concat(retList, ignore_index=True)
        return ret


    def _initDropletPosition(self, mmd, geometricstd, position, Q, clouds=30, meteorologyname="StandardMeteorolgyConstant", **kwargs):

        rv = lognorm(numpy.log(geometricstd), scale=toNumber(mmd, m))
        lower = rv.ppf(1e-4)
        upper = rv.ppf(1-1e-4)

        interval = numpy.logspace(numpy.log(lower),numpy.log(upper),clouds,base=numpy.e)
        dh = numpy.diff(numpy.log(interval))[0]

        massFractionVector = numpy.diff(rv.cdf(interval))*toNumber(Q,kg)
        diameterVector     = numpy.exp(numpy.log(interval[:-1])+dh/2.)  # [m]

        for dropletDiam,dropletQ in zip(diameterVector,massFractionVector):
            droplets = FallingNonEvaporatingDroplets(particleDiameter=dropletDiam*m, Q=dropletQ*kg, position=position, meteorologyName=meteorologyname, **kwargs)
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
        for Ypos in numpy.linspace(0,toNumber(linelength,m),linepositions):
            curpos = (position[0],toUnum(Ypos,m),position[2])
            self._initDropletPosition(mmd, geometricstd, curpos, qCloud, clouds=clouds,meteorologyname=meteorologyname, **kwargs)
