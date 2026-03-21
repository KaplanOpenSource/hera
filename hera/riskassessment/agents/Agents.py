from hera.utils.unitHandler import ureg, Unum, Quantity, celsius, K, tonumber
from hera.riskassessment.agents.effects import  injuryfactory
import numpy
import json



class Agent:
	"""
	Represents a hazardous agent with its associated injury effect models.

	An agent is initialized from a JSON descriptor that defines its name,
	effect parameters (e.g. Ten Berge coefficient), and one or more injury
	effects (e.g. inhalation, thermal). Each effect is accessible by name
	via dictionary-style access (``agent["inhalation"]``).
	"""

	_effects = None

	_effectsParameters = None

	@property
	def effectNames(self):
		"""
		List of effect names defined for this agent.

		Returns
		-------
		list of str
		"""
		return [x for x in self._effects.keys()]

	def  __getitem__(self,name):
		"""
		Access an effect by name.

		Parameters
		----------
		name : str
			The effect name.

		Returns
		-------
		Injury
			The injury effect object.
		"""
		return self._effects[name]

	@property
	def physicalproperties(self):
		"""
		Physical properties of the agent (molecular weight, density, vapor pressure, etc.).

		Returns
		-------
		PhysicalPropeties
		"""
		return self._physicalproperties

	@property
	def fullDescription(self):
		"""
		The full JSON descriptor used to initialize this agent.

		Returns
		-------
		dict
		"""
		return self._agentconfig

	@property
	def effectproperties(self):
		"""
		The effect parameters dictionary (e.g. tenbergeCoefficient).

		Returns
		-------
		dict
		"""
		return self._effectParameters

	@property
	def tenbergeCoefficient(self):
		"""
		The Ten Berge coefficient for dose-response calculations.

		Returns
		-------
		float
		"""
		return self._effectParameters.get("tenbergeCoefficient",1)

	@tenbergeCoefficient.setter
	def tenbergeCoefficient(self,value):
		self._effectParameters["tenbergeCoefficient"] = float(value)
		for effectname,effectconfig in self._agentconfig["effects"].items():
			self._effects[effectname] = injuryfactory.getInjury(effectname,effectconfig,**self._effectParameters)


	@property
	def name(self):
		"""
		The name of the agent.

		Returns
		-------
		str
		"""
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
		"""
		Serialize the agent to a JSON-compatible dictionary.

		Returns
		-------
		dict
		"""
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
		"""Molecular weight in g/mol."""
		return self._molecularWeight

	@property
	def molecularVolume(self):
		"""Molecular volume from the physical properties descriptor."""
		return self._params["molecularVolume"]

	@molecularWeight.setter
	def molecularWeight(self,value):
		self._molecularWeight = ureg(value).to("g/mol")

	@property
	def sorptionCoefficient(self):
		"""Sorption coefficient in cm/s."""
		return self._sorptionCoefficient

	@sorptionCoefficient.setter
	def sorptionCoefficient(self,value):
		self._sorptionCoefficient = ureg(value).to("cm/s")

	@property
	def spreadFactor(self):
		"""Spread factor (dimensionless)."""
		return self._spreadFactor

	@spreadFactor.setter
	def spreadFactor(self,value):
		self._spreadFactor = float(value)

	def getMolecularWeight(self):
		"""Return the molecular weight."""
		return self._molecularWeight

	def getSpreadFactor(self):
		"""Return the spread factor."""
		return self._spreadFactor

	def getSorptionCoefficient(self):
		"""Return the sorption coefficient."""
		return self._sorptionCoefficient


	def getVolatility(self,temperature):
		"""
			Return the vapor saturation concentration [g/cm**3]

		:param temperature:
						The temperature in [C]
		:return:
			The vapor saturation as Unum.
		"""
		if isinstance(temperature, Unum):
			temperature = temperature.asNumber(celsius)
		elif isinstance(temperature, Quantity):
			temperature = temperature.m_as(ureg.celsius)
		MW = self.getMolecularWeight().to(ureg.g/ureg.mol)

		a,b,c,d = self._volatilityConst

		V1 = 10**(a-b/(temperature+c))

		return 1.585287951807229e-5*MW*V1/(temperature+273.16)*ureg.g/ureg.cm**3


	def getDensity(self, temperature):
		"""
			Return the density of the object (g/cm**3).

			:param temperature:
				The temperature in [C]

			:return:
				The density as Unum
		"""
		if isinstance(temperature, Unum):
			temperature = temperature.asNumber(celsius)
		elif isinstance(temperature, Quantity):
			temperature = temperature.m_as(ureg.celsius)
		a,b,c = self._densityConst

		return (a-b*(temperature-c))*ureg.g/ureg.cm**3

	def vaporPressure(self, temperature):
		"""
		Calculate the vapor pressure at a given temperature.

		Parameters
		----------
		temperature : float or Unum or Quantity
			The temperature (in Kelvin or with units).

		Returns
		-------
		float
			Vapor pressure in bar.
		"""
		if isinstance(temperature, Unum):
			temperature = temperature.asNumber(K)
		elif isinstance(temperature, Quantity):
			temperature = temperature.m_as(ureg.kelvin)
		A = self._vaporConst["A"]
		B = self._vaporConst["B"]
		C = self._vaporConst["C"]
		D = self._vaporConst["D"] if "D" in self._vaporConst.keys() else 0
		E = self._vaporConst["E"] if "E" in self._vaporConst.keys() else 0
		F = self._vaporConst["F"] if "F" in self._vaporConst.keys() else 0
		units = ureg(self._vaporConst["units"])
		return ((10 ** (A - B / (temperature - C) + D * numpy.log10(
			temperature) + E * temperature + F * temperature ** 2)) * units).m_as(ureg.bar)


	def __init__(self,configJSON):
		"""
		Initialize physical properties from an agent config JSON.

		Parameters
		----------
		configJSON : dict
			The full agent descriptor. If it contains a ``physicalProperties``
			key, those values are used to set molecular weight, sorption
			coefficient, and spread factor.
		"""
		if "physicalProperties" in configJSON:
			self._params 		 = configJSON["physicalProperties"]
			self.molecularWeight = self._params.get("molecularWeight","1*g/mol")
			self.sorptionCoefficient = self._params.get("sorptionCoefficient","1*cm/s")
			self.spreadFactor    = self._params.get("spreadFactor",1)
	

	def toJSON(self):
		"""
		Serialize physical properties to a JSON-compatible dictionary.

		Returns
		-------
		dict
		"""
		return self._params

