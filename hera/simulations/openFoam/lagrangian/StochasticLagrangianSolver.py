from .abstractLagrangianSolver import absractStochasticLagrangianSolver_toolkitExtension
from .. import FLOWTYPE_INCOMPRESSIBLE


class StochasticLagrangianSolver_toolkitExtension(absractStochasticLagrangianSolver_toolkitExtension):

    def __init__(self, toolkit):
        super().__init__(toolkit=toolkit)

    def createDispersionFlowField(self, flowName, flowData, OriginalFlowField, dispersionDuration,flowType=FLOWTYPE_INCOMPRESSIBLE,
                                  overwrite: bool = False, useDBSupport: bool = True,dispersionFieldList=[]):
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
        for requiredField in  ['ustar','distanceFromWalls']:
            if requiredField not in dispersionFieldList:
                dispersionFieldList.append(requiredField)

        super().createDispersionFlowField(flowName=flowName,
                                          flowData=flowData,
                                          OriginalFlowField=OriginalFlowField,
                                          dispersionFieldList=dispersionFieldList,
                                          dispersionDuration=dispersionDuration,
                                          flowType=flowType,
                                          overwrite=overwrite,
                                          useDBSupport=useDBSupport)
