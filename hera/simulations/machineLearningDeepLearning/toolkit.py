import numpy
import inspect
import os
from hera.utils.logging import with_logger, get_classMethod_logger
from hera.toolkit import abstractToolkit

from hera.simulations.machineLearningDeepLearning.torch.modelContainer import torchLightingModelContainer
from hera.utils import dictToMongoQuery
from hera.utils.jsonutils import compareJSONS
from hera.utils.SALibUtils import SALibUtils
from hera.utils.jsonutils import setJSONPath

try:
    import SALib
    from SALib.sample import morris as morris_sample
    from SALib.analyze import morris as morris_analyze
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
        return compareJSONS(**dict([(f"{mdl.desc['modelID']}", mdl.desc['model']) for mdl in docList]),
                            longFormat=longFormat,changeDotToUnderscore=True).rename(columns=dict(datasetName="modelID"))

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



    def sensitivityAnalysis_morris(self,modelContainer,problemContainer,maxEpoch=500,sampleParameters=dict(),analysisParameters=dict()):
        """
            Performs the sensitivity analysis of Morris to identify the
        Parameters
        ----------
        modelContainer : The torch container wrapper model.
            This is the hera container that contains all the data required for the
            dataset, dataloaders and initializing the model.

        problemContainer : JSON
            The JSON that defines the parameters, values and their type for the SALib.
            Created with the SALibUtils.buildSAProblem

        Returns
        -------

        """
        morris_sample_parameters = {
            'N': 100,  # Number of trajectories (you can change as needed)
            'num_levels': 4,  # Number of levels in the grid
            'optimal_trajectories': None,  # No trajectory optimization by default
            'local_optimization': False  # No local optimization by default
        }
        morris_sample_parameters.update(sampleParameters)


        morris_analyze_parameters = {
            'num_levels': morris_sample_parameters['num_levels'],  # must match the sample design
            'grid_jump':  morris_sample_parameters['num_levels'] // 2,
            'conf_level': 0.95,  # default confidence level
            'print_to_console': False,  # suppress printing by default
            'num_resamples': 1000,  # bootstrap iterations for CIs
            'seed': None  # random seed (can be set for reproducibility)
        }
        morris_analyze_parameters.update(analysisParameters)

        baseJson = modelContainer.modelJSON

        raw_param_values = morris_sample.sample(problemContainer['problem'],**morris_sample_parameters)
        samples = SALibUtils.transformSample(batchList=raw_param_values,problemContainer=problemContainer)

        Y = numpy.zeros(len(samples))

        for i,sample in enumerate(samples):
            # Transfer to a dict of param name -> real value.
            paramDict = dict([(name,value) for name,value in zip(problemContainer['problem']['names'],sample)])
            sampleJSON = setJSONPath(base=baseJson,valuesDict = paramDict,inPlace=False)
            emptyContainer = self.getEmptyTorchModelContainer()


            emptyContainer.modelJSON = sampleJSON
            emptyContainer.fit(maxEpoch)
            stats = emptyContainer.getStatistics()

            result = stats.loc[stats.groupby("tag")["step"].idxmax(), ["tag", "value"]].set_index("tag")
            Y[i] = result.loc["val_loss_epoch"].item()

        Si = morris_analyze.analyze(problemContainer['problem'], samples, Y, **morris_analyze_parameters)
        return Si

    #def sensitivityAnalysisExecute_morris

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
