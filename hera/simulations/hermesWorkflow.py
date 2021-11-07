from ..utils.logging.loggedObject import loggedObject
import jsonpath_ng as jsonpath
from ..utils.jsonutils import loadJSON
import json
import os

class hermesWorkflow(loggedObject):
    """
        This class serves as an interface to edit the hermes workflow.

        We assume that the workflow is expanded.

    """
    _workflow = None # The workflow JSON.

    @property
    def _workflowJSON(self):
        return self._workflow["workflow"]

    @property
    def nodeList(self):
        return self._workflowJSON['nodeList']

    @property
    def nodes(self):
        return self._workflowJSON['nodes']


    def __init__(self,workflow):
        """
            Initialize with workflow
        Parameters
        ----------
        workflow: str, dict, file type.
            Can be the file name, a string of the JSON, or a file type.
        """
        super().__init__(loggerName=None)

        self.logger.info("Initializing")

        self._workflow = loadJSON(workflow)


    def __getitem__(self, item):
        """
            Returns a node.
        Parameters
        ----------
        item: str
            The node name

        Returns
        -------
            A node object of the requested node.

        """
        nodeJSON = self._workflowJSON['nodes'][item]
        return hermesNode(item,nodeJSON)


    def __delitem__(self, key):
        """
            Removes a node from the workflow.
            Raises ValueError if node not found.

        Parameters
        ----------
        key: The name of the node

        Returns
        -------
            None

        """

        # 1. Remove the node from the nodelist in key: "workflow.nodeList"
        try:
            self.nodeList.remove(key)

            # 2. Remove the node from the nodes. "workflow.nodes"
            del self._workflowJSON['nodes'][key]

        except ValueError:
            raise ValueError(f"{key} node is not found. Found nodes: {','.join(self.nodeList)}")


    def getNodeValue(self,jsonpath):
        """
            Returns a value from the JSON path.
            The search is relative to the 'nodes' node in the workflow.

        Parameters
        ----------
        jsonpath: str
            The path to obtain.

        Returns
        -------
            List
            jsonpath DatumInContext object with the query results.
        """
        jsonexpr = jsonpath.parse(jsonpath)
        return jsonexpr.find(self._workflowJSON['nodes'])


    def updateNodeValue(self,jsonpath,value):
        """

            Returns a value from the JSON path.
            The search is relative to the 'nodes' node in the workflow.

        Parameters
        ----------
        jsonpath: str
            The path to obtain.

        value: str
            The new value.

        Returns
        -------
            None
        """
        jsonexpr = jsonpath.parse(jsonpath)
        jsonexpr.update(self._workflowJSON['nodes'],value)


    def write(self,filename,overwrite=False):
        """
            writes the new workflow to the file.
        Parameters
        ----------
        filename : str
            The file name

        overwrite: bool
            If true, the writes over existing file. Otherwise raises an exception.

        Returns
        -------

        """
        self.logger.info("-- Start --")
        if not overwrite:
            if os.path.exists(filename):
                err = f"{filename} alread exists. Use overwrite=True to overwite it"
                self.logger.error(err)
                raise FileExistsError(err)

        with open(filename,'w') as outputFile:
            json.dump(self._workflow,outputFile)

        self.logger.info("-- End --")

class hermesNode:
    """
        An interface to the JSON of an hermes workflow.
    """
    _nodeJSON= None
    _nodeName = None

    def __init__(self,nodeName,nodeJSON):
        self._nodeJSON =nodeJSON
        self._nodeName = nodeName

    @property
    def name(self):
        return self._nodeName

    @property
    def parameters(self):
        return self._nodeJSON['Execution']['input_parameters']


    def __getitem__(self, item):
        return self._nodeJSON['Execution']['input_parameters'][item]

    def keys(self):
        return self._nodeJSON['Execution']['input_parameters'].keys()

    def values(self):
        return self._nodeJSON['Execution']['input_parameters'].values()

    def items(self):
        return self._nodeJSON['Execution']['input_parameters'].items()