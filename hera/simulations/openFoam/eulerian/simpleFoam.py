from .abstractEulerianSolver import absractEulerianSolver_toolkitExtension


class simpleFoam_toolkitExtension(absractEulerianSolver_toolkitExtension):

    def __init__(self,toolkit):
        super().__init__(toolkit=toolkit,solverName="simpleFOAM",incompressible=True)


    def buildSpecializedField(self,workflow,turbulenceType,):
        """
            Updating the field of the workflow.

            - The solver name
            - The field names (depends on the turbulence models).
            - turbulence node.

        Parameters
        ----------
        workflow : hermesWorkflow
            The workflow object
        turbulenceType : str
            Can be :
                - laminar : laminar turbulence model.
                - RAS model : The name RAS model.
                - LES model : The name of RAS model.

        Returns
        -------

        """
        pass

