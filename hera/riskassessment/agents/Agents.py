from unum.units import *
from .effects import  injuryfactory
from ...simulations.utils import toNumber,toUnum
import numpy

from ...datalayer import ProjectMultiDBPublic


class AgentHome(ProjectMultiDBPublic):
	"""
		A class to load and get agents.

		Supports retrieval from the DB and initializing from a descriptor.

	"""

	def __init__(self, projectName="AgentsCollection"):
		super().__init__(projectName=projectName, publicProjectName="AgentsCollection")



	def getAgent(self, nameOrDesc):
		"""
			Initialize the agents.

		:param nameOrDesc: str or JSON.
			Can be either the name of the agent (str) or
			the descriptor

			{
				"name" : [the name of the agent],
				"effectParameters" : {
					TenBergeCoefficient and ect.
				},
				"effects": {
					"effect name" : { effect data (+ injury levels) }


				}
			}


		:param projectName: str
				The name of the project in the local DB that will be searched for the agent.
		:return:
		"""
		if isinstance(nameOrDesc,str):
			configList = self.getMeasurementsDocuments(type='Agent',Agent=nameOrDesc)
			if len(configList)==0:
				raise ValueError(f"Agent {nameOrDesc} is not found. Load it with hera-risk-agent load")
			descriptor = configList[0].desc
		elif isinstance(nameOrDesc,dict):
			descriptor = nameOrDesc
		else:
			raise ValueError("nameOrDesc must be the agent name (str) or its JSON description (dict) ")

		return Agent(descriptor)


	def listAgents(self):
		"""
			Lists the agents that are currently loaded in the DB (both local and public).

		:return: list
			A list of agent names.

		"""
		configList = self.getMeasurementsDocuments(type='Agent')
		return [x.desc['Agent'] for x in configList]


	def loadAgent(self,name,agentDescription,public=True):
		"""
			Adds the agent to the DB. Either to the public or to the local DB.

		:param name: str
				Agent name
		:param agentDescription: dict
				The agent description

		:return:
				None
		"""

		agentDescription['name'] = name

		destDB = "public" if public else None

		self.addMeasurementsDocument(resource="",
									 type="Agent",
									 dataFormat="string",
									 users=destDB ,
									 desc=agentDescription)

class Agent:

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
			self._effects[effectname] = injuryfactory.getInjury(effectname,effectconfig,**self._effectParameters)


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

	def vaporPressure(self, temperature):
		temperature = toNumber(temperature, K)
		A = self._vaporConst["A"]
		B = self._vaporConst["B"]
		C = self._vaporConst["C"]
		D = self._vaporConst["D"] if "D" in self._vaporConst.keys() else 0
		E = self._vaporConst["E"] if "E" in self._vaporConst.keys() else 0
		F = self._vaporConst["F"] if "F" in self._vaporConst.keys() else 0
		units = self._vaporConst["units"]
		return toNumber((10 ** (A - B / (temperature - C) + D * numpy.log10(
			temperature) + E * temperature + F * temperature ** 2)) * units, bar)


	def __init__(self,configJSON):

		if "physicalProperties" in configJSON:
			self._params 		 = configJSON["physicalProperties"]
			self.molecularWeight = self._params["molecularWeight"]
			self.sorptionCoefficient = self._params["sorptionCoefficient"] if "sorptionCoefficient" in self._params.keys() else "1"
			self.spreadFactor    = self._params["spreadFactor"] if "spreadFactor" in self._params.keys() else "1"
	


