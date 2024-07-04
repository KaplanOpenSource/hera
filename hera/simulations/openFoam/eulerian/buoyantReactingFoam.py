from .abstractEulerianSolver import absractEulerianSolver_toolkitExtension

class buoyantReactingFoam_toolkitExtension(absractEulerianSolver_toolkitExtension):

    def __init__(self,toolkit):
        super().__init__(toolkit=toolkit,solverName="buoyantReactingFoam",incompressible=False)

