from ...datalayer import project
#from ...riskassessment import Agent
from ..utils import toNumber
from unum.units import *
import numpy

class evaporationModels(object):

    _agent = None
    _Mair = None
    _Vair = None
    _Magent = None
    _Vagent = None
    _dinamicViscocityModel = None
    _molecularDiffusionModel = None
    _vaporPressure = None
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
        self._agent = newAgent

    @property
    def Mair(self):
        return self._Mair

    @property
    def Vair(self):
        return self._Vair

    @property
    def Magent(self):
        return self._Magent

    @property
    def Vagent(self):
        return self._Vagent

    @property
    def vaporPressureCoeffs(self):
        return self._vaporPressure

    def __init__(self,agent, evaporationModel="US",dinamicViscocityModel="powerLaw",molecularDiffusionModel="FSG"):

        #self._agent=Agent(agent)
        self._agent = agent
        self._Mair = 28.967
        self._Vair = 20.1
        self._Magent = agent["physicalproperties"]["molecularWeight"].asNumber(g/mol)#agent.physicalproperties.molecularWeight.asNumber(g/mol)
        self._molecularDiffusionModel = molecularDiffusionModel
        self._dinamicViscocityModel = dinamicViscocityModel
        self._vaporPressure = agent["physicalproperties"]["vaporPressure"]
        self._evaporationModel = evaporationModel
        try:
            self._Vagent = agent["physicalproperties"]["molecularVolume"].asNumber(cm ** 3 / mol)  # agent.physicalproperties.molecularVolume.asNumber(ml/mol)
        except:
            print("Note that the molecular volume is not given, therefore, the molecular diffusion is set to EPA.")
            self._molecularDiffusionModel = "EPA"
            self._Vagent = None

    def molecularDiffusion(self,temperature):
        temperature = toNumber(temperature, K)
        return getattr(self,f"molecularDiffusion_{self._molecularDiffusionModel}")(temperature)

    def molecularDiffusion_FSG(self,temperature):
        return 0.001*(temperature**1.75)*numpy.sqrt(1/self.Mair+1/self.Magent)/((numpy.cbrt(self.Vair)+numpy.cbrt(self.Vagent))**2)

    def molecularDiffusion_EPA(self,temperature):
        return 0.0000409*(temperature**1.9)*numpy.sqrt(1/self.Mair+1/self.Magent)/numpy.cbrt(self.Magent)

    def dynamicViscocityAir(self,temperature):
        temperature = toNumber(temperature, K)
        return getattr(self,f"dynamicViscocityAir_{self._dinamicViscocityModel}")(temperature)

    def dynamicViscocityAir_powerLaw(self,temperature):
        return 1.8205*10**(-5)*numpy.sqrt(temperature/293)

    def Reynolds(self,diameter,velocity,temperature):
        density = self.Mair / (temperature * 8.205 * 10 ** (-2))
        return diameter*velocity*density/self.dynamicViscocityAir(temperature=temperature)

    def Schmidt(self,temperature):
        density = self.Mair / (temperature * 8.205 * 10 ** (-2))
        return self.dynamicViscocityAir(temperature=temperature)/(density*self.molecularDiffusion(temperature))

    def vaporPressure(self,temperature):
        return getattr(self,f"vaporPressure_{self.vaporPressureCoeffs['model']}")(temperature)

    def vaporPressure_Antoine(self,temperature):
        temperature = toNumber(temperature, K)
        A = self.vaporPressureCoeffs["A"]
        B = self.vaporPressureCoeffs["B"]
        C = self.vaporPressureCoeffs["C"]
        units = self.vaporPressureCoeffs["units"]
        return toNumber(10**(A-B/(temperature-C))*units,bar)

    def flux(self,diameter,velocity,temperature):
        return getattr(self, f"flux_{self._evaporationModel}")(diameter,velocity,temperature)

    def flux_US(self, diameter,velocity,temperature):

        Re = self.Reynolds(diameter=diameter,velocity=velocity,temperature=temperature)
        Sc = self.Schmidt(temperature=temperature)
        Km = 0.664*(Re**(-0.5))*(Sc**(-2/3))*velocity if Re<=20000 else 0.0366*(Re**(-0.2))*(Sc**(-2/3))*velocity
        return Km*self.Magent*self.vaporPressure(temperature)/(temperature * 8.205 * 10 ** (-2))
