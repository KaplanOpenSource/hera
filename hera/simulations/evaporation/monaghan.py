import pandas
import numpy
from unum.units import *
from ..gaussian.Meteorology import StandardMeteorolgyConstant
from ..utils import toUnum,toNumber
from pyriskassessment.agents.Agents import Agent
from unum import Unum

class MonaghanConstantConditions(object):
    """
        The Monaghan-McPhreson-Nusse model for fixed conditions.

        Take from

        Monaghan J & McPhreson WR: A mathematical model for prediction vapour disages on and downsid of contrmainated areas of grassland.
        Ref 281 (Yellow system).

        P0, Mass ratio, U2P1 and wind height are the constant from the monaghan report.

    """
    _agent = None # The agent that is being evaporated.
    _standardMeteorolgyConstant = None

    @property
    def agent(self):
        return self._agent

    @agent.setter
    def agent(self,value):

        if not isinstance(value,Agent):
            raise ValueError("Agent must be a pyriskassessment.agents.Agents.Agent")
        self._agent = value

    @property
    def P0(self):
        return 5.0/60*cm/s

    @property
    def U2P1(self):
        return 0.01

    @property
    def MASSRATIO(self):
        return 10.

    @property
    def WINDHEIGHT(self):
        return 200*cm


    @property
    def meteorology(self):
        return self._standardMeteorolgyConstant

    @meteorology.setter
    def meteorology(self,value):
        if not isinstance(value,StandardMeteorolgyConstant):
            raise ValueError("Meteorolgoy must be a pynumericalmodels.gaussian.Meteorology.StandardMeteorolgyConstant")
        self._standardMeteorolgyConstant = value



    def __init__(self,agent,meteorology):
        self.agent = agent
        self.meteorology = meteorology

    def Q(self,timelist,dropletRadius):
        """
            Evaluates the evaporation rate in constant wind from droplets on grass.
            Returns the cumulative evaporation.
            Note that it will not sum up to 1 due to irreversible absorption to the ground.  

        :param timelist:
                The time steps in which we need to estimate the results.

        :param dropletRadius:
                The radius of the droplets.

        :return:
                The (cumulative) fraction that was evaporated.

        """

        evaporation = [0]*9    # There are 8 constants. But since the index in NUSSE begins with 1 we define 9 vars and the first one is not used.
        evaporation_nc = [0]*4 # There are 3 constants. But since the index in NUSSE begins with 1 we define 9 vars and the first one is not used.

        sorption_p      = self.agent.physicalproperties.getSorptionCoefficient()
        sf              = self.agent.physicalproperties.getSpreadFactor()
        agentvolatility = self.agent.physicalproperties.getVolatility(self.meteorology.skinSurfaceTemperature)
        agentdensity    = self.agent.physicalproperties.getDensity(self.meteorology.skinSurfaceTemperature)

        p1 = self.meteorology.getWindVeclocity(self.WINDHEIGHT) * self.U2P1

        evaporation[1] = (sorption_p + p1) / (p1 - self.P0)
        evaporation[4] = (sorption_p + self.P0) / (p1-self.P0)
        evaporation[5] = self.MASSRATIO - evaporation[4]*(1-self.MASSRATIO)
        evaporation[2] = 2* sorption_p * evaporation[5]/(sorption_p + p1)
        evaporation[3] = (self.P0**2-sorption_p**2)/(sorption_p + p1)**2
        evaporation[6] = (1-self.MASSRATIO)* evaporation[1]
        evaporation[7] = (p1+self.P0)/(p1+sorption_p)
        evaporation[8] = 1./(1+self.MASSRATIO)

        evaporation_nc[1] = 4*self.MASSRATIO / (3* (self.MASSRATIO+1)) * agentdensity / (sf**2*agentvolatility*(self.P0+sorption_p))
        evaporation_nc[2] = (sorption_p+self.P0) / ((sorption_p+p1)*self.MASSRATIO) # equivalent to variable nc2 in the code of Hezi, but not multiplied by the droplet radius.
        evaporation_nc[3] = evaporation_nc[2]*evaporation_nc[1] # equivalent to variable nc3 in the code of Hezi

        evaporation = [Unum.coerceToUnum(x).asNumber() for x in evaporation]
        evaporation_nc = [Unum.coerceToUnum(x) for x in evaporation_nc]

        ret = numpy.zeros(len(timelist))

        sstime  = (2*toUnum(dropletRadius,um)*evaporation_nc[3]).asNumber(s)
        tau     = (2 * toUnum(dropletRadius, um) * evaporation_nc[1]).asNumber(s)

        sstime_timedelta = pandas.to_timedelta(sstime,"s")
        tau_timedelta = pandas.to_timedelta(tau, "s")

        Iss = timelist < timelist[0] + sstime_timedelta
        Ibetween_sstime_and_tau = (timelist > timelist[0] + sstime_timedelta) & (timelist < timelist[0] + tau_timedelta)
        Iabove_tau = (timelist > timelist[0] + tau_timedelta)

        Tlocal = (timelist[Iss]-timelist[0]).total_seconds()/sstime
        ret[Iss] = evaporation[8]*evaporation[7]*Tlocal

        Tlocal = (timelist[Ibetween_sstime_and_tau]-timelist[0]).total_seconds() / (sstime)
        ret[Ibetween_sstime_and_tau] = evaporation[8]*(evaporation[7] + evaporation[1]*(evaporation[5]**2 * (1/(evaporation[6]-Tlocal)-1/(evaporation[6]-1))\
                                                                                        - evaporation[2]*numpy.log((evaporation[6]-Tlocal)/(evaporation[6]-1)) + evaporation[3]*(1-Tlocal)))

        Tlocal = numpy.ones(len(timelist[Iabove_tau])) * (tau/sstime)
        ret[Iabove_tau] = evaporation[8]*(evaporation[7] + evaporation[1]*(evaporation[5]**2 * (1/(evaporation[6]-Tlocal)-1/(evaporation[6]-1))\
                                                                                        - evaporation[2]*numpy.log((evaporation[6]-Tlocal)/(evaporation[6]-1)) + evaporation[3]*(1-Tlocal)))


        return ret
