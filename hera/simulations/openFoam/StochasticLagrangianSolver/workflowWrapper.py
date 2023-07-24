from ..OFWorkflow import Workflow_Lagrangian

class Workflow_StochasticLagrangianSolver(Workflow_Lagrangian):

    def __init__(self ,workflowJSON,workflowHeraDocument=None,name=None):
        super().__init__(workflowJSON=workflowJSON, workflowHeraDocument=workflowHeraDocument,name=name)

        # Make sure that the
        # dispersionFlowField exists
        if 'dispersionFlowField' not in self.parameters.parameters:
            err = "The StochasticLagrangianSolver must have a dispersionFlowField specification in the parameters node"
            self.logger.error(err)
            raise ValueError(err)

    @property
    def dispersionFlowField(self):
        return

    @property
    def dispersionFlowField(self):
        return self.parameters.parameters['dispersionFlowField']

    @dispersionFlowField.setter
    def dispersionFlowField(self, value):
        self.parameters.parameters['dispersionFlowField'] = str(value)
