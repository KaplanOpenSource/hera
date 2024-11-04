import json
import os
from ..toolkit import abstractToolkit,TOOLKIT_SAVEMODE_FILEANDDB,TOOLKIT_SAVEMODE_FILEANDDB_REPLACE,TOOLKIT_SAVEMODE_NOSAVE
from ..datalayer import datatypes, nonDBMetadataFrame
from .agents.Agents import Agent
from .presentation.casualtiesFigs import casualtiesPlot
from .protectionpolicy.ProtectionPolicy import ProtectionPolicy
from ..simulations.LSM.toolkit import LSMToolkit
import geopandas
from unum.units import *

class RiskToolkit(abstractToolkit):
    """
        A class to load and get agents.

        Supports retrieval from the DB and initializing from a descriptor.

    """
    _presentation = None
    _protectionPolicy = None
    _analysis = None

    @property
    def analysis(self):
        return self._analysis

    @property
    def ProtectionPolicy(self):
        return self._protectionPolicy

    @property
    def presentation(self):
        return self._presentation

    def __init__(self, projectName, filesDirectory=None):
        super().__init__(projectName=projectName, filesDirectory=filesDirectory, toolkitName="RiskAssessment")
        self._presentation = casualtiesPlot()
        self._protectionPolicy = ProtectionPolicy
        self._analysis = analysis(self)

    def getAgent(self, nameOrDesc, version=None):
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
        if isinstance(nameOrDesc, str):
            descriptor = self.getDataSourceData(nameOrDesc, version=version)
            if descriptor is None:
                raise ValueError(f"Agent {nameOrDesc} is not found. Load it with hera-risk-agent load")

        elif isinstance(nameOrDesc, dict):
            descriptor = nameOrDesc
        else:
            raise ValueError("nameOrDesc must be the agent name (str) or its JSON description (dict) ")

        return Agent(descriptor)


    def listAgentsNames(self):
        """
            Lists the agents that are currently loaded in the DB (both local and public).

        :return: list
            A list of agent names.

        """
        return [x.desc["datasourceName"] for x in self.getDataSourceDocumentsList()]

    def loadAgent(self, name, agentDescription, version,saveMode=TOOLKIT_SAVEMODE_FILEANDDB):
        """
			Adds the agent to the DB. Either to the public or to the local DB.
			Equivalent to loadData

        :param name: str
                Agent name
        :param agentDescription: dict
                The agent description

        :return:
                None
        """
        agentDescription['name'] = name
        agentDescription['version'] = version
        return self.loadData(agentDescription,saveMode=saveMode)

    def loadData(self, fileNameOrData, saveMode=TOOLKIT_SAVEMODE_FILEANDDB,**kwargs):
        """
            Abstract loading a data from file. Manages the parsing of the
            datafile.

			Equivalent to loadData

        Parameters
        ----------
        fileNameOrData: str
                If str , the datafile to load
                If other objects - convert the
        parser: str
                The name of the parser to use

        :param saveMode: str
                Can be either:

                    - TOOLKIT_SAVEMODE_NOSAVE   : Just load the data from file and return the datafile


                    - TOOLKIT_SAVEMODE_FILEANDDB : Loads the data from file and save to a file and store to the DB as a source.
                                                    Raise exception if the entry exists.

                    - TOOLKIT_SAVEMODE_FILEANDDB_REPLACE: Loads the data from file and save to a file and store to the DB as a source.
                                                    Replace the entry in the DB if it exists.

        """
        if isinstance(fileNameOrData, str):
            if os.path.isfile(fileNameOrData):
                with open(fileNameOrData,'r') as readFile:
                    agentDescription = json.load(readFile)
            else:
                agentDescription = json.loads(fileNameOrData)
        elif isinstance(fileNameOrData,dict):
            agentDescription = fileNameOrData
        else:
            raise ValueError("fileNameOrData must be a file, JSON str or JSON object (dict)")

        name = agentDescription['name']
        version = agentDescription.get('version',None)

        agentDoc = self.getDataSourceDocument(datasourceName=name, version=version)

        if agentDoc is None:
            self.addDataSource(name, resource=json.dumps(agentDescription), dataFormat=datatypes.JSON_DICT, **agentDescription)

        elif saveMode == TOOLKIT_SAVEMODE_FILEANDDB:
            raise ValueError(f"Agent {name} version {agentDoc.desc.get('version',None)} in the database.")
        else:
            agentDoc.resource = json.dumps(agentDescription)
            agentDoc.desc['version']  = version
            agentDoc.save()
        return nonDBMetadataFrame(agentDescription) if agentDoc is None else agentDoc

class analysis():

    _datalayer = None
    _LSM = None

    @property
    def LSM(self):
        return self._LSM

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, dataLayer):
        self._datalayer = dataLayer
        self._LSM = LSMToolkit(projectName=self.datalayer.projectName)

    def getRiskAreas(self, tenbergeCoefficient, levels,Q=1*kg,protectionPolicy=None,LSMfile=None, **LSMparams):
        """
        returns the bounds and polygons of risk areas from dispersion of agents with different Ten Berge
        coefficients and hazardous levels.
        params:
            Q = the total mass of dispersed particles; default is 1 kg
            tenbergeCoefficients = a list of Ten Berge coefficients, or a single value.
            levels = a list of levels for which the risk areas are calculated, should be equal in length to tenbergeCoefficients
        """
        levels = levels if isinstance(levels, list) else [levels]
        lsmSim = self.LSM.getSimulations(**LSMparams)[0] if LSMfile is None else self.LSM.singleSimulation(LSMfile)
        Concentration = lsmSim.getConcentration(Q=Q).isel(z=0)
        if protectionPolicy is not None:
            Concentration = self.datalayer.ProtectionPolicy(protectionPolicy).compute(Concentration)
        description = {"effectParameters": {"tenbergeCoefficient": tenbergeCoefficient},
                       "effects": {"RegularPopulation": {"type": "Lognormal10",
                       "calculator": {"TenBerge": {"breathingRate": 10}},
                       "parameters": {"type": "Lognormal10DoseResponse",
                        "levels": ["Severe"],"parameters": {"Severe": {"TL_50": 10,"sigma": 0.5}}}}}}
        dumbAgent = self.datalayer.getAgent(description)
        boundsList = []
        polygonList = []
        for level in levels:
            ToxicLoad = dumbAgent.RegularPopulation.calculateRaw(Concentration,"C",isel={"datetime":-1})
            ToxicLoad = ToxicLoad.to_dataframe().reset_index()
            ToxicLoad = ToxicLoad.loc[ToxicLoad.C>=level]
            data = geopandas.geodataframe.GeoDataFrame(ToxicLoad, geometry=geopandas.points_from_xy(ToxicLoad.x, ToxicLoad.y))
            boundsList.append(data.unary_union.bounds)
            polygonList.append(data.unary_union.convex_hull)

        return geopandas.geodataframe.GeoDataFrame({"tenBergCoefficient":tenbergeCoefficient, "level":levels, "bounds":boundsList, "geometry":polygonList})

