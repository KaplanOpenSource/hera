import inspect
import os
from hera.utils.logging import with_logger, get_classMethod_logger
from hera.toolkit import abstractToolkit

from hera.simulations.machineLearningDeepLearning.torch.modelContainer import torchLightingModelContainer
from hera.utils import dictToMongoQuery
from hera.utils.jsonutils import compareJSONS

try:
    import SALib
    import SALib.sample as SA_sample
    import SALib.analyze as SA_analyze
except ImportError:
    print("SALib not installed, cannot support sensitivity analysis")


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

    def getEmptyTorchModelContainer(self):
        return torchLightingModelContainer(self)

    def listTorchModels(self, modelObjectOrName=None, longFormat=True, **qry):
        qryMongo = dictToMongoQuery(qry)
        if modelObjectOrName is not None:
            if isinstance(modelObjectOrName, str):
                qryMongo["model__model__classpath"] = modelObjectOrName
            else:
                qryMongo["model__model__classpath"] = self.get_model_fullname(modelObjectOrName)
        docList = self.getSimulationsDocuments(type=torchLightingModelContainer.MODEL, **qryMongo)
        return compareJSONS(**dict([(f"M{mdl.desc['modelID']}", mdl.desc['model']) for mdl in docList]),
                            longFormat=longFormat,changeDotToUnderscore=True)

    def getTorchModelContainerByID(self, modelID, **qry):
        qryModngo = dictToMongoQuery(qry,prefix="model")
        docList = self.getSimulationsDocuments(type=torchLightingModelContainer.MODEL, modelID=modelID, **qryModngo)
        if len(docList)>0:
            mdlDesc = self.getEmptyTorchModelContainer()
            mdlDesc.modelJSON = docList[0].desc['model']
            mdlDesc.load()
        else:
            mdlDesc= None
        return mdlDesc



    def sensitivityAnalysis_morris(self,modelContainer,SALibProblem,categoryParameters=[],sampleParameters=dict()):
        """
            Performs the sensitivity analysis of Morris to identify the
        Parameters
        ----------
        modelContainer : The torch container wrapper model.
            This is the hera container that contains all the data required for the
            dataset, dataloaders and initializing the model.

        SALibProblem : JSON
            The JSON that defines the parameters and the values for the SALib.

        parametersType: dict
            parameter Name -> 'category' | 'int' | 'log' | <add more types like gaussian in the future>

            A dict to determine the parameter types.

            Since sample is uniform, we need a transformer to get either cate

            The names of the parameters that are catory parameters.
            The values of these parameters are rounded and then casted to int.

        Returns
        -------

        """
        morris_sample_parameters = {
            'N': 100,  # Number of trajectories (you can change as needed)
            'num_levels': 4,  # Number of levels in the grid
            'grid_jump': None,  # Defaults to num_levels // 2 if None
            'optimal_trajectories': None,  # No trajectory optimization by default
            'local_optimization': False  # No local optimization by default
        }

        baseJson = modelContainer.modelJSON
        morris_sample_parameters.update(sampleParameters)

        param_values = sobol.sample(SALibProblem,**morris_sample_parameters)

        Y = np.zeros(param_values.shape[0])

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
