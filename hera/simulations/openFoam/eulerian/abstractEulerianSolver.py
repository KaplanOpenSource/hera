
class absractEulerianSolver_toolkitExtension:

    toolkit = None

    analysis = None # link to the analysis of the StochasticLagrnagian
    presentation = None # link to the presentation of the stochasticLagrangian

    DOCTYPE_EULERIAN_CACHE = "eulerianCacheDocType"

    solverName = None
    incompressible = None

    def __init__(self,toolkit,solverName,incompressible):
        self.toolkit = toolkit
        self.solverName = solverName
        self.incompressible = incompressible
        #self.analysis = analysis(self)

    def specialize_snappyHexMeshFromMesh(self):
        """
            Change the snappy hex mesh according to the
            flow field.
        Returns
        -------

        """
        pass