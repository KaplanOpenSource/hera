import os 
import json 
from  pyriskassessment.agents import effects
from  hera.simulations.utils import toNumber,toUnum
from unum.units import *
from ..heraDatalayer.heraDatalayer import heraDatalayer

class Agent(object): 

	_effects = None 

	_effectsParameters = None 

	def effects(self):
		return [x for x in self._effects.keys()]

	def  __getitem__(self,name): 
		return self._effects[name]

	@property
	def physicalproperties(self):
		return self._physicalproperties

	@property
	def fullDescription(self):
		return self._agentconfig

	@property 
	def effectproperties(self): 
		return self._effectParameters

	@property
	def tenbergeCoefficient(self):
		return self._effectParameters.get("tenbergeCoefficient",1)

	@tenbergeCoefficient.setter
	def tenbergeCoefficient(self,value):
		self._effectParameters["tenbergeCoefficient"] = float(value)
		for effectname,effectconfig in self._agentconfig["effects"].items():
			self._effects[effectname] = effects.injuryfactory.getInjury(effectname,effectconfig,**self._effectParameters)

	def __init__(self,name, projectName="AgentsCollection"):
		"""
			Reads the agent json. 
			The name of the agent is the configuration file. 

			{
				"effectParameters" : { 
					TenBergeCoefficient and ect. 
				},
				"effects": { 
					"effect name" : { effect data (+ injury levels) } 
					

				}
			}
		"""
		# conffile = os.path.join(os.path.expanduser("~"),".pyriskassessment","%s.json" % name)
		# with open(conffile,'r') as config:
		# 	self._agentconfig = json.load(config)
		self._agentconfig = heraDatalayer().getAgent(projectName=projectName, Agent=name)
		
		self._effectParameters = self._agentconfig.get("effectParameters",{})

		self._effects = {}
		for effectname,effectconfig in self._agentconfig["effects"].items():
			self._effects[effectname] = effects.injuryfactory.getInjury(effectname,effectconfig,**self._effectParameters)

		self.__dict__.update(self._effects)

		self._physicalproperties = PhysicalPropeties(self._agentconfig)


class PhysicalPropeties(object):
	"""
		Implements the approximation of the physical properties of an agent.

	"""
	_params = None
	_molecularWeight = None
	_sorptionCoefficient = None
	_spreadFactor= None

	@property
	def _volatilityConst(self):
		return self._params["volatilityConstants"]

	@property
	def _densityConst(self):
		return self._params["densityConstants"]


	@property
	def molecularWeight(self):
		return self._molecularWeight

	@molecularWeight.setter
	def molecularWeight(self,value):
		self._molecularWeight = toUnum(eval(value),g/mol)


	@property
	def sorptionCoefficient(self):
		return self._sorptionCoefficient

	@sorptionCoefficient.setter
	def sorptionCoefficient(self,value):
		self._sorptionCoefficient = toUnum(eval(value),cm/s)

	@property
	def spreadFactor(self):
		return self._spreadFactor

	@spreadFactor.setter
	def spreadFactor(self,value):
		self._spreadFactor = float(value)

	def getMolecularWeight(self):
		return self._molecularWeight

	def getSpreadFactor(self):
		return self._spreadFactor

	def getSorptionCoefficient(self):
		return self._sorptionCoefficient


	def getVolatility(self,temperature):
		"""
			Return the vapor saturation concentration [g/cm**3]

		:param temperature:
						The temperature in [C]
		:return:
			The vapor saturation as Unum.
		"""
		temperature = toNumber(temperature,celsius)
		MW = self.getMolecularWeight().asNumber(g/mol)

		a,b,c,d = self._volatilityConst

		V1 = 10**(a-b/(temperature+c))

		return 1.585287951807229e-5*MW*V1/(temperature+273.16)*g/cm**3


	def getDensity(self, temperature):
		"""
			Return the density of the object (g/cm**3).

			:param temperature:
				The temperature in [C]

			:return:
				The density as Unum
		"""
		temperature = toNumber(temperature,celsius)
		a,b,c = self._densityConst

		return (a-b*(temperature-c))*g/cm**3


	def __init__(self,configJSON):

		if "physicalProperties" in configJSON:
			self._params 		 = configJSON["physicalProperties"]
			self.molecularWeight = self._params["molecularWeight"]
			self.sorptionCoefficient = self._params["sorptionCoefficient"]
			self.spreadFactor    = self._params["spreadFactor"]
	


