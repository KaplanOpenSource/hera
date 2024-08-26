from .abstractLagrangianSolver import absractStochasticLagrangianSolver_toolkitExtension
from .. import FLOWTYPE_INCOMPRESSIBLE


class StochasticLagrangianSolver_toolkitExtension(absractStochasticLagrangianSolver_toolkitExtension):

    def __init__(self, toolkit):
        super().__init__(toolkit=toolkit)

    def createDispersionFlowField(self, flowName, flowData, OriginalFlowField, dispersionDuration,flowType=FLOWTYPE_INCOMPRESSIBLE,
                                  overwrite: bool = False, useDBSupport: bool = True):
        """
            Building the dispersion field,
        Parameters
        ----------
        flowName
        flowData : JSON
        A JSON with the format:
        {
                            "originalFlow": {
                                "time": {
                                    "temporalType": "steadyState",
                                    "timestep": 400
                                },
                                "linkMeshSymbolically": true
                            },
                            "dispersionFields": {
                                "ustar": 0.25,
                                "Hmix": 1000
                            }
                        }

            The dispersion fields are defines in the dispersionFields key.

        OriginalFlowField
        overwrite
        useDBSupport

        Returns
        -------

        """
        for requiredField,defaultValue in  [('ustar',0),('distanceFromWalls',0)]:
            flowData['dispersionFields'].setdefault(requiredField,defaultValue)

        super().createDispersionFlowField(flowName=flowName,
                                          flowData=flowData,
                                          OriginalFlowField=OriginalFlowField,
                                          dispersionDuration=dispersionDuration,
                                          flowType=flowType,
                                          overwrite=overwrite,
                                          useDBSupport=useDBSupport)
