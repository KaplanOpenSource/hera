import json
import os

from hera import get_classMethod_logger
from hera.utils import ureg
from ..toolkit import abstractToolkit,TOOLKIT_SAVEMODE_FILEANDDB,TOOLKIT_SAVEMODE_FILEANDDB_REPLACE,TOOLKIT_SAVEMODE_NOSAVE
from ..datalayer import datatypes, nonDBMetadataFrame
from .agents.Agents import Agent
from .presentation.casualtiesFigs import casualtiesPlot
from .protectionpolicy.ProtectionPolicy import ProtectionPolicy
from ..simulations.LSM.toolkit import LSMToolkit
import geopandas

class RiskToolkit(abstractToolkit):
    """
    Toolkit for agent-based risk assessment.

    Manages hazardous agents (chemical, thermal, blast), their injury
    effect models, protection policies, and casualty estimation.
    Agents are stored as versioned data sources and can be loaded
    from JSON descriptors or from the database.
    """
    _presentation = None
    _protectionPolicy = None
    _analysis = None

    @property
    def analysis(self):
        """
        Access the risk analysis layer.

        Returns
        -------
        analysis
            Object providing risk area calculation and LSM integration.
        """
        return self._analysis

    @property
    def ProtectionPolicy(self):
        """
        Access the ProtectionPolicy class for building protection action pipelines.

        Returns
        -------
        type
            The ProtectionPolicy class (not an instance).
        """
        return self._protectionPolicy

    @property
    def presentation(self):
        """
        Access the presentation layer for casualty visualizations.

        Returns
        -------
        casualtiesPlot
        """
        return self._presentation

    def __init__(self, projectName, filesDirectory=None, connectionName=None):
        """
        Initialize the RiskToolkit.

        Parameters
        ----------
        projectName : str
            The project name.
        filesDirectory : str, optional
            Directory for file outputs.
        connectionName : str, optional
            The DB connection name.
        """
        super().__init__(projectName=projectName, filesDirectory=filesDirectory, toolkitName="RiskAssessment", connectionName=connectionName)
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
        logger = get_classMethod_logger(self, "getAgent")
        if isinstance(nameOrDesc, str):
            descriptor = self.getDataSourceData(nameOrDesc, version=version)
            if descriptor is None:
                raise ValueError(f"Agent {nameOrDesc} is not found. Load it with hera-risk-agent load")

            agentDataSource = self.getDataSourceDocument(nameOrDesc, version=version)
            if 'name' not in descriptor:
                logger.warning("name not specified in the agent descriptor, using the datasource name as a default")
                descriptor['name'] = agentDataSource.desc['datasourceName']
            if 'version' not in descriptor:
                logger.warning("version not specified in the agent descriptor, using the datasource version as a default")
                descriptor['version'] = agentDataSource.desc['version']

        elif isinstance(nameOrDesc, dict):
            descriptor = nameOrDesc
        else:
            raise ValueError("nameOrDesc must be the agent name (str) or its JSON description (dict) ")

        return Agent(descriptor)


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
    """
    Risk analysis layer providing risk area calculations and LSM integration.
    """

    _datalayer = None
    _LSM = None

    @property
    def LSM(self):
        """
        Access the LSM toolkit for Lagrangian dispersion data.

        Returns
        -------
        LSMToolkit
        """
        return self._LSM

    @property
    def datalayer(self):
        """
        Access the parent RiskToolkit (data layer).

        Returns
        -------
        RiskToolkit
        """
        return self._datalayer

    def __init__(self, dataLayer):
        """
        Initialize the analysis layer.

        Parameters
        ----------
        dataLayer : RiskToolkit
            The parent toolkit providing data access.
        """
        self._datalayer = dataLayer
        self._LSM = LSMToolkit(projectName=self.datalayer.projectName)

    def getRiskAreas(self, tenbergeCoefficient, levels,Q=1*ureg.kg,protectionPolicy=None,LSMfile=None, **LSMparams):
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

