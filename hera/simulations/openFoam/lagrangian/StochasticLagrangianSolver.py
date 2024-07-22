from .abstractLagrangianSolver import absractStochasticLagrangianSolver_toolkitExtension
from .. import FLOWTYPE_INCOMPRESSIBLE


class StochasticLagrangianSolver_toolkitExtension(absractStochasticLagrangianSolver_toolkitExtension):

    def __init__(self, toolkit):
        super().__init__(toolkit=toolkit)

    def createDispersionFlowField(self, flowName, flowData, OriginalFlowField, flowType=FLOWTYPE_INCOMPRESSIBLE,
                                  overwrite: bool = False, useDBSupport: bool = True):
        """
            Building the dispersion field,
        Parameters
        ----------
        flowName
        flowData
        OriginalFlowField
        overwrite
        useDBSupport

        Returns
        -------

        """
        dispersionFieldList = ['Ustar,cellHeights']

        super().createDispersionFlowField(flowName=flowName,
                                          flowData=flowData,
                                          OriginalFlowField=OriginalFlowField,
                                          dispersionFieldList=dispersionFieldList,
                                          flowType=flowType,
                                          overwrite=overwrite,
                                          useDBSupport=useDBSupport)
