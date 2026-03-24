from .abstractEulerianSolver import absractEulerianSolver_toolkitExtension

class buoyantReactingFoam_toolkitExtension(absractEulerianSolver_toolkitExtension):
    """Toolkit extension for the compressible buoyantReactingFoam solver."""

    def __init__(self,toolkit):
        """Initialize the buoyantReactingFoam extension as compressible."""
        super().__init__(toolkit=toolkit,solverName="buoyantReactingFoam",incompressible=False)

