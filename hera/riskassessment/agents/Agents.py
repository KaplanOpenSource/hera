from unum.units import *
from .effects import  injuryfactory
from ...utils import tonumber,tounit
import numpy
import json



class Agent:

	_effects = None 

	_effectsParameters = None 

	@property
	def effectNames(self):
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
			self._effects[effectname] = injuryfactory.getInjury(effectname,effectconfig,**self._effectParameters)


	@property
	def name(self):
	    return self._agentconfig['name']


	def __init__(self,descriptor):
		"""
			Constructor of the Agent project.

			Initializes an agent.


		Parameters
		-----------
		descriptor: JSON
			A JSON object (dict) that holds all the information on the agent.

			{
				"name" : [the name of the agent],
				"effectParameters" : {
					TenBergeCoefficient and ect.
				},
				"effects": {
					"effect name" : { effect data (+ injury levels) }


				}
			}

		"""
		self._agentconfig = descriptor
		self._effectParameters = self._agentconfig.get("effectParameters",{})

		self._effects = {}
		for effectname,effectconfig in self._agentconfig["effects"].items():
			self._effects[effectname] = injuryfactory.getInjury(effectname,effectconfig,**self._effectParameters)

		self.__dict__.update(self._effects)

		self._physicalproperties = PhysicalPropeties(self._agentconfig)


	def toJSON(self):
		ret = dict(name=self.name,
				   physicalProperties=self.physicalproperties.toJSON(),
				   effect={})

		for effect in self.effectNames:
			ret['effect'][effect] = self[effect].toJSON()

		return ret


	def __str__(self):
		return json.dumps(self.toJSON(),indent=4)




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
	def _vaporConst(self):
		return self._params["vaporPressure"]

	@property
	def molecularWeight(self):
		return self._molecularWeight

	@property
	def molecularVolume(self):
		return self._params["molecularVolume"]

	@molecularWeight.setter
	def molecularWeight(self,value):
		self._molecularWeight = tounit(eval(value),g/mol)

	@property
	def sorptionCoefficient(self):
		return self._sorptionCoefficient

	@sorptionCoefficient.setter
	def sorptionCoefficient(self,value):
		self._sorptionCoefficient = tounit(eval(value),cm/s)

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
		temperature = tonumber(temperature,celsius)
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
		temperature = tonumber(temperature,celsius)
		a,b,c = self._densityConst

		return (a-b*(temperature-c))*g/cm**3

	def vaporPressure(self, temperature):
		temperature = tonumber(temperature, K)
		A = self._vaporConst["A"]
		B = self._vaporConst["B"]
		C = self._vaporConst["C"]
		D = self._vaporConst["D"] if "D" in self._vaporConst.keys() else 0
		E = self._vaporConst["E"] if "E" in self._vaporConst.keys() else 0
		F = self._vaporConst["F"] if "F" in self._vaporConst.keys() else 0
		units = self._vaporConst["units"]
		return tonumber((10 ** (A - B / (temperature - C) + D * numpy.log10(
			temperature) + E * temperature + F * temperature ** 2)) * units, bar)


	def __init__(self,configJSON):

		if "physicalProperties" in configJSON:
			self._params 		 = configJSON["physicalProperties"]
			self.molecularWeight = self._params.get("molecularWeight","1")
			self.sorptionCoefficient = self._params.get("sorptionCoefficient","1")
			self.spreadFactor    = self._params.get("spreadFactor",1)
	

	def toJSON(self):
		return self._params

