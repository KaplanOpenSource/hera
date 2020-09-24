from hera.datalayer.collection import Measurements_Collection

class riskDatalayer():

    _userDB = None
    _publicDB = None

    def __init__(self):

        self._userDB = Measurements_Collection()
        self._publicDB = Measurements_Collection(user="public")

    def addAgent(self, Agent, description,  projectName="AgentsCollection", public=False):

        description["Agent"] = Agent
        if public:
            self._publicDB.addDocument(projectName=projectName, resource="", type="Agent", dataFormat="string", desc=description)
        else:
            self._userDB.addDocument(projectName=projectName, resource="", type="Agent", dataFormat="string", desc=description)

    def getAgent(self, Agent, projectName="AgentsCollection"):
        if len(self._userDB.getDocuments(projectName=projectName, Agent=Agent)) > 0:
            newAgent = self._userDB.getDocuments(projectName=projectName, Agent=Agent)[0].asDict()["desc"]
        elif len(self._publicDB.getDocuments(projectName=projectName, Agent=Agent)) > 0:
            newAgent = self._publicDB.getDocuments(projectName=projectName, Agent=Agent)[0].asDict()["desc"]
        else:
            raise KeyError(f"Agent {agent} was not found.")
        return newAgent