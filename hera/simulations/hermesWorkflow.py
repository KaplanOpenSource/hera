from ..utils.logging.loggedObject import loggedObject
import jsonpath_ng as jsonpath
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

        if hasattr(workflow,'read'):
            self.logger.debug("loading from file-like")
            self._workflow = json.load(workflow)
        elif isinstance(workflow,str):

            if os.path.exists(workflow):
                self.logger.debug("loading from file")
                with open(workflow) as jsonFile:
                    self._workflow = json.load(jsonFile)
            else:
                self.logger.debug("loading from str")
                self._workflow = json.loads(workflow)
        elif isinstance(workflow,dict):
            self.logger.debug("using dict")
            self._workflow = workflow
        else:
            err = f"workflow type: {type(workflow)} is unknonw. Must be str, file-like or dict. "
            self.logger.error(err)
            raise ValueError(err)

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