from ...datalayer import project
from ...riskassessment import AgentHome
from ..utils import toNumber
from unum.units import *
import numpy

class evaporationModels(object):

    _agent = None
    _Mair = None
    _Vair = None
    _Vagent = None
    _dinamicViscocityModel = None
    _molecularDiffusionModel = None
    _evaporationModel = None

    @property
    def agent(self):
        return self._agent

    @property
    def dinamicViscocityModel(self):
        return self._dinamicViscocityModel

    @dinamicViscocityModel.setter
    def dinamicViscocityModel(self,newModel):
        self._dinamicViscocityModel=newModel

    @property
    def molecularDiffusionModel(self):
        return self._molecularDiffusionModel

    @molecularDiffusionModel.setter
    def molecularDiffusionModel(self,newModel):
        self._molecularDiffusionModel=newModel

    @property
    def evaporationModel(self):
        return self._evaporationModel

    @evaporationModel.setter
    def evaporationModel(self,newModel):
        self._evaporationModel=newModel

    @property
    def agent(self):
        return self._agent

    @agent.setter
    def agent(self,newAgent):
        self._agent = AgentHome.getAgent(newAgent)

    @property
    def Mair(self):
        return self._Mair

    @property
    def Vair(self):
        return self._Vair

    @property
    def Magent(self):
        return self._agent.physicalproperties.molecularWeight.asNumber(g/mol)

    @property
    def Vagent(self):
        return self._Vagent

    def __init__(self,agent, evaporationModel="US",dinamicViscocityModel="powerLaw",molecularDiffusionModel="FSG"):

        self._agent=AgentHome.getAgent(agent)
        self._Mair = 28.967
        self._Vair = 20.1
        self._molecularDiffusionModel = molecularDiffusionModel
        self._dinamicViscocityModel = dinamicViscocityModel
        self._evaporationModel = evaporationModel
        try:
            self._Vagent = self._agent.physicalproperties.molecularVolume.asNumber(cm ** 3 / mol)  # agent.physicalproperties.molecularVolume.asNumber(ml/mol)
        except:
            print("Note that the molecular volume is not given, therefore, the molecular diffusion is set to EPA.")
            self._molecularDiffusionModel = "EPA"
            self._Vagent = None

    def molecularDiffusion(self,temperature):
        temperature = toNumber(temperature, K)
        return getattr(self,f"molecularDiffusion_{self._molecularDiffusionModel}")(temperature)

    def molecularDiffusion_FSG(self,temperature):
        return 0.0000001*(temperature**1.75)*numpy.sqrt(1/self.Mair+1/self.Magent)/((numpy.cbrt(self.Vair)+numpy.cbrt(self.Vagent))**2)

    def molecularDiffusion_EPA(self,temperature):
        return 0.00000000409*(temperature**1.9)*numpy.sqrt(1/self.Mair+1/self.Magent)/numpy.cbrt(self.Magent)

    def dynamicViscocityAir(self,temperature):
        temperature = toNumber(temperature, K)
        return getattr(self,f"dynamicViscocityAir_{self._dinamicViscocityModel}")(temperature)

    def dynamicViscocityAir_powerLaw(self,temperature):
        return 1.8205*10**(-2)*numpy.sqrt(temperature/293)

    def Reynolds(self,diameter,velocity,temperature):
        density = self.Mair / (temperature * 8.205 * 10 ** (-5))
        return diameter*velocity*density/self.dynamicViscocityAir(temperature=temperature)

    def Schmidt(self,temperature):
        density = self.Mair / (temperature *  8.205 * 10 ** (-5)) # (g/(m^3))
        return self.dynamicViscocityAir(temperature=temperature)/(density*self.molecularDiffusion(temperature))

    def flux(self,diameter,velocity,temperature):
        return getattr(self, f"flux_{self._evaporationModel}")(diameter,velocity,temperature)

    def flux_US(self, diameter,velocity,temperature):
        temperature = toNumber(temperature, K)
        diameter = toNumber(diameter, m)
        velocity = toNumber(velocity,m/s)
        Re = self.Reynolds(diameter=diameter,velocity=velocity,temperature=temperature)
        Sc = self.Schmidt(temperature=temperature)
        Km = 0.664*(Re**(-0.5))*(Sc**(-2/3))*velocity if Re<=20000 else 0.0366*(Re**(-0.2))*(Sc**(-2/3))*velocity
        R = 8.205 * 10 ** (-5) # gas constant
        return Km*self.Magent*self.agent.physicalproperties.vaporPressure(temperature)/(temperature * R)
