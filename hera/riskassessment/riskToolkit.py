import json
import os
from ..toolkit import abstractToolkit,TOOLKIT_SAVEMODE_FILEANDDB,TOOLKIT_SAVEMODE_FILEANDDB_REPLACE
from ..datalayer import datatypes, nonDBMetadataFrame
from .agents.Agents import Agent
from .presentation.casualtiesFigs import casualtiesPlot
from .protectionpolicy.ProtectionPolicy import ProtectionPolicy


class RiskToolkit(abstractToolkit):
    """
        A class to load and get agents.

        Supports retrieval from the DB and initializing from a descriptor.

    """
    _presentation = None
    _protectionPolicy = None

    @property
    def ProtectionPolicy(self):
        return self._protectionPolicy

    @property
    def presentation(self):
        return self._presentation

    def __init__(self, projectName):
        super().__init__(projectName=projectName, toolkitName="RiskAssessment")
        self._presentation = casualtiesPlot()
        self._protectionPolicy = ProtectionPolicy

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
            descriptor = self.getDatasourceData(nameOrDesc,version=version)
            if descriptor is None:
                raise ValueError(f"Agent {nameOrDesc} is not found. Load it with hera-risk-agent load")

        elif isinstance(nameOrDesc, dict):
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
        return self.getDatasourceDocumentsList()

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

    def loadData(self, fileNameOrData, saveMode=TOOLKIT_SAVEMODE_FILEANDDB):
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

        agentDoc = self.getDatasourceDocument(datasourceName=name,version=version)

        if agentDoc is None:
            self.addDataSource(name, resource=json.dumps(agentDescription), dataFormat=datatypes.JSON_DICT, **agentDescription)

        elif saveMode == TOOLKIT_SAVEMODE_FILEANDDB:
            raise ValueError(f"Agent {name} version {agentDoc.desc.get('version',None)} in the database.")

        else:
            agentDoc.resource = agentDescription
            agentDoc.desc['version']  = version
            agentDoc.save()
        return nonDBMetadataFrame(agentDescription) if agentDoc is None else agentDoc
