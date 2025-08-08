from hera.utils.logging import with_logger, get_classMethod_logger
from hera.toolkit import abstractToolkit

from hera.simulations.machineLearningDeepLearning.torch.torchModels import torchLightingModel
from hera.utils import dictToMongoQuery
from hera.utils.jsonutils import compareJSONS

class machineLearningDeepLearningToolkit(abstractToolkit):
    """
        The class handles machine/deep learning models.

        It helps saving hyper parameters and it provide simple
        tools (like batch/train splitting).

        Notes:
            * Torch models it requires pytorch installed.
            * SkiLearn requires scikitlearn installed.

    """

    def __init__(self, projectName: str, filesDirectory: str = None):
        """
            Initializes the machineLearning/deepLearning toolkit.

        Parameters
        ----------
        projectName: str
            The project where the models are stored.

        filesDirectory : str
            The directory to write all the Workflow and the outputs. default is current directory.
        """
        super().__init__(projectName=projectName,
                         filesDirectory=filesDirectory,
                         toolkitName= "machineLearningDeepLearningToolkit")

    def getEmptyTorchModel(self):
        return torchLightingModel(self)

    def listTorchModels(self,longFormat=True,**qry,):
        qryMongo = dictToMongoQuery(qry)
        docList = self.getSimulationsDocuments(type=torchLightingModel.MODEL,**qryMongo)

        return compareJSONS(**dict([(f"M{mdl.desc['modelID']}", mdl.desc['model']) for mdl in docList]),
                            longFormat=longFormat,changeDotToUnderscore=True)

    def getTorchModelByID(self,modelID,**qry):
        qryModngo = dictToMongoQuery(qry,prefix="model")
        docList = self.getSimulationsDocuments(type=torchLightingModel.MODEL, modelID=modelID,**qryModngo)
        if len(docList)>0:
            mdlDesc = self.getEmptyTorchModel()
            mdlDesc.modelJSON = docList[0].desc['model']
            mdlDesc.load()
        else:
            mdlDesc= None
        return mdlDesc