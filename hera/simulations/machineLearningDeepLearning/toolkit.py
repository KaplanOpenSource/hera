import inspect
import os
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

    def listTorchModels(self, modelObjectOrName=None, longFormat=True, **qry):
        qryMongo = dictToMongoQuery(qry)
        if modelObjectOrName is not None:
            if isinstance(modelObjectOrName, str):
                qryMongo["model__model__classpath"] = modelObjectOrName
            else:
                qryMongo["model__model__classpath"] = self.get_model_fullname(modelObjectOrName)
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


    ## ====================================================================================================
    ## ====================================================================================================
    ## ===================================== CLASS METHODS ================================================
    ## ====================================================================================================
    ## ====================================================================================================

    @classmethod
    def get_model_fullname(cls,modelCls):
        name,data = cls.get_class_info(modelCls)
        return data['classpath']

    @classmethod
    def get_class_info(cls,modelCls):
        module = modelCls.__module__
        name = modelCls.__name__
        file_path = inspect.getfile(modelCls)
        file_path = os.path.dirname(os.path.abspath(file_path))

        full_path = f"{module}.{name}"
        patList = file_path.split(os.path.sep)
        moduleNameIndex = patList.index(full_path.split(".")[0])
        patList[0] = '/'
        module_file_path = os.path.join(*patList[:moduleNameIndex])

        return name, dict(classpath=full_path, filepath=module_file_path)
