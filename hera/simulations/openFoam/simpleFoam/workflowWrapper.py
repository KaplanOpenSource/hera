from ..OFWorkflow import Workflow_Eulerian

class Workflow_simpleFoam(Workflow_Eulerian):

    def __init__(self ,workflowJSON,workflowHeraDocument=None,name=None):
        super().__init__(workflowJSON=workflowJSON, workflowHeraDocument=workflowHeraDocument,name=name)
