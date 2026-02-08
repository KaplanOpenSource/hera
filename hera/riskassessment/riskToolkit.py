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
    Toolkit for risk assessment and agent-based hazard analysis.

    The RiskToolkit manages toxic agents, their effects, and risk calculations
    based on dispersion simulations. It integrates with LSM simulations to
    calculate exposure and health impacts.

    Key Features
    ------------
    - Agent management (toxic substances, effect parameters)
    - Risk area calculation from dispersion results
    - Protection policy application
    - Casualty estimation and visualization
    - Integration with LSM dispersion simulations

    Data Sources
    ------------
    Data sources are agent definitions stored as JSON documents. Each agent
    defines:
    - Effect parameters (TenBerge coefficients, etc.)
    - Dose-response relationships
    - Injury levels and thresholds

    Examples
    --------
    >>> from hera import toolkitHome
    >>> risk_tk = toolkitHome.getToolkit("RiskAssessment", projectName="my_project")
    >>> 
    >>> # Get an agent
    >>> agent = risk_tk.getAgent("chlorine")
    >>> 
    >>> # List available agents
    >>> agents = risk_tk.listAgentsNames()
    >>> 
    >>> # Calculate risk areas
    >>> risk_areas = risk_tk.analysis.getRiskAreas(
    ...     tenbergeCoefficient=1.0,
    ...     levels=[0.1, 0.5, 1.0],
    ...     Q=100*kg
    ... )
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
        Retrieve an agent by name or create from descriptor dictionary.

        This method loads an agent definition from the database or creates an agent
        object from a provided descriptor dictionary. Agents define toxic substances
        with their effect parameters and dose-response relationships.

        Parameters
        ----------
        nameOrDesc : str or dict
            If str: Name of the agent to load from the database.
            If dict: Agent descriptor dictionary with the following structure:
            
            {
                "name": "agent_name",
                "effectParameters": {
                    "tenbergeCoefficient": 1.0,
                    # ... other effect parameters
                },
                "effects": {
                    "effect_name": {
                        "type": "Lognormal10",
                        "calculator": {
                            "TenBerge": {
                                "breathingRate": 10
                            }
                        },
                        "parameters": {
                            "type": "Lognormal10DoseResponse",
                            "levels": ["Severe"],
                            "parameters": {
                                "Severe": {
                                    "TL_50": 10,
                                    "sigma": 0.5
                                }
                            }
                        }
                    }
                }
            }

        version : tuple, optional
            Specific version of the agent to retrieve when nameOrDesc is a string.
            If None, uses the default version. Default is None.

        Returns
        -------
        Agent
            Agent object with effect calculation capabilities.

        Raises
        ------
        ValueError
            If agent name is not found in database, or if nameOrDesc has invalid type.

        Examples
        --------
        >>> # Get agent by name
        >>> agent = risk_tk.getAgent("chlorine")
        >>> 
        >>> # Get specific version
        >>> agent = risk_tk.getAgent("chlorine", version=(1, 0, 0))
        >>> 
        >>> # Create agent from descriptor
        >>> agent_desc = {
        ...     "name": "custom_agent",
        ...     "effectParameters": {"tenbergeCoefficient": 1.0},
        ...     "effects": {
        ...         "RegularPopulation": {
        ...             "type": "Lognormal10",
        ...             "calculator": {"TenBerge": {"breathingRate": 10}},
        ...             "parameters": {
        ...                 "type": "Lognormal10DoseResponse",
        ...                 "levels": ["Severe"],
        ...                 "parameters": {
        ...                     "Severe": {"TL_50": 10, "sigma": 0.5}
        ...                 }
        ...             }
        ...         }
        ...     }
        ... }
        >>> agent = risk_tk.getAgent(agent_desc)
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
        List all agent names available in the database.

        Returns all agent datasource names registered in both the local project
        database and the public database.

        Returns
        -------
        list of str
            List of agent names (datasource names) available in the database.

        Examples
        --------
        >>> # List all available agents
        >>> agents = risk_tk.listAgentsNames()
        >>> print(f"Available agents: {agents}")
        >>> # Output: ['chlorine', 'ammonia', 'sulfur_dioxide', ...]
        """
        return [x.desc["datasourceName"] for x in self.getDataSourceDocumentsList()]

    def loadAgent(self, name, agentDescription, version, saveMode=TOOLKIT_SAVEMODE_FILEANDDB):
        """
        Register an agent in the database.

        Adds a new agent definition to the database or updates an existing one.
        This is equivalent to calling loadData() with agent-specific parameters.

        Parameters
        ----------
        name : str
            Unique name identifier for the agent. This will be used as the
            datasource name in the database.

        agentDescription : dict
            Complete agent descriptor dictionary. Should include all effect
            parameters and dose-response relationships. The 'name' and 'version'
            fields will be automatically added/overwritten.

        version : tuple
            Version tuple (major, minor, patch) for the agent. Example: (1, 0, 0).

        saveMode : str, optional
            Save mode for the operation. Options:
            - TOOLKIT_SAVEMODE_FILEANDDB: Save to file and database, raise error if exists
            - TOOLKIT_SAVEMODE_FILEANDDB_REPLACE: Save to file and database, overwrite if exists
            - TOOLKIT_SAVEMODE_NOSAVE: Don't save, just return document object
            Default is TOOLKIT_SAVEMODE_FILEANDDB.

        Returns
        -------
        MeasurementDocument or nonDBMetadataFrame
            The created or updated document, or a non-DB metadata frame if
            saveMode is NOSAVE.

        Raises
        ------
        ValueError
            If agent already exists and saveMode doesn't allow overwrite.

        Examples
        --------
        >>> # Load a new agent
        >>> agent_desc = {
        ...     "effectParameters": {"tenbergeCoefficient": 1.0},
        ...     "effects": {...}
        ... }
        >>> doc = risk_tk.loadAgent(
        ...     name="new_agent",
        ...     agentDescription=agent_desc,
        ...     version=(1, 0, 0),
        ...     saveMode=risk_tk.TOOLKIT_SAVEMODE_FILEANDDB
        ... )
        """
        agentDescription['name'] = name
        agentDescription['version'] = version
        return self.loadData(agentDescription,saveMode=saveMode)

    def loadData(self, fileNameOrData, saveMode=TOOLKIT_SAVEMODE_FILEANDDB, **kwargs):
        """
        Load agent data from file, JSON string, or dictionary.

        This method handles loading agent definitions from various sources and
        storing them in the database. It supports JSON files, JSON strings, and
        Python dictionaries.

        Parameters
        ----------
        fileNameOrData : str or dict
            Agent definition source. Can be:
            - File path (str): Path to a JSON file containing agent definition
            - JSON string (str): Valid JSON string with agent definition
            - Dictionary (dict): Python dictionary with agent definition
            The agent definition must include 'name' and optionally 'version' fields.

        saveMode : str, optional
            Save mode for the operation. Options:
            - TOOLKIT_SAVEMODE_NOSAVE: Load and return without saving
            - TOOLKIT_SAVEMODE_FILEANDDB: Save to file and database, raise error if exists
            - TOOLKIT_SAVEMODE_FILEANDDB_REPLACE: Save to file and database, overwrite if exists
            Default is TOOLKIT_SAVEMODE_FILEANDDB.

        **kwargs : dict
            Additional keyword arguments passed to addDataSource() if creating
            a new document.

        Returns
        -------
        MeasurementDocument or nonDBMetadataFrame
            The created or updated document, or a non-DB metadata frame if
            saveMode is NOSAVE or document doesn't exist.

        Raises
        ------
        ValueError
            If fileNameOrData has invalid type, or if agent exists and saveMode
            doesn't allow overwrite.
        FileNotFoundError
            If fileNameOrData is a file path that doesn't exist.

        Examples
        --------
        >>> # Load from JSON file
        >>> doc = risk_tk.loadData(
        ...     fileNameOrData="/path/to/agent.json",
        ...     saveMode=risk_tk.TOOLKIT_SAVEMODE_FILEANDDB
        ... )
        >>> 
        >>> # Load from dictionary
        >>> agent_dict = {
        ...     "name": "agent_name",
        ...     "version": (1, 0, 0),
        ...     "effectParameters": {...},
        ...     "effects": {...}
        ... }
        >>> doc = risk_tk.loadData(agent_dict)
        >>> 
        >>> # Load from JSON string
        >>> json_str = '{"name": "agent", "effectParameters": {...}}'
        >>> doc = risk_tk.loadData(json_str)
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
    Analysis layer for risk assessment calculations.

    Provides methods for calculating risk areas, toxic loads, and exposure
    assessments based on dispersion simulation results and agent properties.

    Parameters
    ----------
    dataLayer : RiskToolkit
        The RiskToolkit instance to provide analysis for.

    Attributes
    ----------
    datalayer : RiskToolkit
        Reference to the toolkit instance.
    LSM : LSMToolkit
        LSM toolkit instance for accessing dispersion simulations.
    """

    _datalayer = None
    _LSM = None

    @property
    def LSM(self):
        """LSM toolkit instance for dispersion simulation access."""
        return self._LSM

    @property
    def datalayer(self):
        """RiskToolkit instance."""
        return self._datalayer

    def __init__(self, dataLayer):
        """
        Initialize the risk analysis layer.

        Parameters
        ----------
        dataLayer : RiskToolkit
            The RiskToolkit instance to provide analysis for.
        """
        self._datalayer = dataLayer
        self._LSM = LSMToolkit(projectName=self.datalayer.projectName)

    def getRiskAreas(self, tenbergeCoefficient, levels, Q=1*kg, protectionPolicy=None, LSMfile=None, **LSMparams):
        """
        Calculate risk area polygons for different toxic load levels.

        This method computes risk areas by calculating toxic loads from dispersion
        concentration data and identifying regions where toxic load exceeds specified
        threshold levels. Returns bounding boxes and convex hull polygons for each level.

        Parameters
        ----------
        tenbergeCoefficient : float
            TenBerge coefficient for toxic load calculation. This parameter relates
            concentration and exposure time to toxic load.

        levels : list of float or float
            Toxic load threshold levels for which risk areas are calculated.
            If a single float is provided, it is converted to a list. Each level
            represents a different severity threshold.

        Q : unum.units quantity, optional
            Total mass of dispersed particles. Must be a unum quantity with mass units.
            Default is 1*kg.

        protectionPolicy : dict, optional
            Protection policy dictionary to apply to concentration data before
            calculating toxic load. If None, raw concentration is used. Default is None.

        LSMfile : str, optional
            Path to LSM simulation file. If provided, this file is used instead of
            querying the database. If None, LSMparams are used to query simulations.
            Default is None.

        **LSMparams : dict
            Query parameters for finding LSM simulations in the database. These are
            passed to LSM.getSimulations(). Common parameters include:
            - releaseHeight: Release height in meters
            - releaseLocation: (x, y) tuple of release coordinates
            - simulationName: Name of specific simulation

        Returns
        -------
        geopandas.GeoDataFrame
            GeoDataFrame with columns:
            - 'tenBergCoefficient': TenBerge coefficient used
            - 'level': Toxic load level threshold
            - 'bounds': Tuple of (minx, miny, maxx, maxy) bounding box
            - 'geometry': Convex hull polygon of the risk area

        Notes
        -----
        - Creates a temporary agent with the specified TenBerge coefficient
        - Uses the last time step from the concentration data
        - Risk areas are calculated as convex hulls of points exceeding the threshold

        Examples
        --------
        >>> from unum.units import *
        >>> 
        >>> # Calculate risk areas for multiple levels
        >>> risk_areas = risk_tk.analysis.getRiskAreas(
        ...     tenbergeCoefficient=1.0,
        ...     levels=[0.1, 0.5, 1.0, 2.0],
        ...     Q=100*kg,
        ...     releaseHeight=10*m,
        ...     releaseLocation=(1000, 2000)
        ... )
        >>> 
        >>> # Use specific LSM file
        >>> risk_areas = risk_tk.analysis.getRiskAreas(
        ...     tenbergeCoefficient=1.5,
        ...     levels=[0.5],
        ...     Q=50*kg,
        ...     LSMfile="/path/to/simulation.nc"
        ... )
        >>> 
        >>> # With protection policy
        >>> policy = {"type": "evacuation", "threshold": 0.1}
        >>> risk_areas = risk_tk.analysis.getRiskAreas(
        ...     tenbergeCoefficient=1.0,
        ...     levels=[0.1, 0.5],
        ...     protectionPolicy=policy,
        ...     releaseHeight=15*m
        ... )
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

