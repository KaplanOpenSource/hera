import json
import pydoc

import torch
import numpy
import inspect
import os
import pandas
from collections.abc import Iterable

from hera import toolkitHome
from hera.toolkit import abstractToolkit
from hera.simulations.machineLearningDeepLearning.torch.modelContainer import torchLightingModelContainer
from hera.utils.logging import with_logger, get_classMethod_logger
from hera.utils import dictToMongoQuery
from hera.utils.jsonutils import compareJSONS
from hera.utils.SALibUtils import SALibUtils
from hera.utils.jsonutils import setJSONPath
from hera.utils.logging import get_classMethod_logger
from hera.utils.zipUtils import zip_items,list_json_files_in_zip

try:
    import SALib
    from SALib.sample import morris as morris_sample
    from SALib.analyze import morris as morris_analyze
except ImportError:
    print("SALib not installed, cannot support sensitivity analysis")

try:
    from joblib import Parallel, delayed
except ImportError:
    print("joblib not installed, cannot support parallel sensitivity analysis")


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


    def getTorchModelFromJSON(self,modelJSON):
        newModel = torchLightingModelContainer(self)
        newModel.modelJSON = modelJSON
        return newModel


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

    def packTorchModelByID(self,modelIDorListID,packFileName):
        """
            Pack one or more models, checkpoint and
        Parameters
        ----------
        modelIDorListID : int/list of it
                ID or list of model id
        packFileName: string
                The name of the output file.

        Returns
        -------

        """
        if isinstance(modelIDorListID,Iterable):
            modelContainer = [self.getTorchModelContainerByID(x) for x in modelIDorListID]
            if len(modelContainer):
                raise ValueError(f"Non of the Models ID  {modelIDorListID} not found project {self.projectName}. If this is not the project you ment, make sure caseConfiguration.json exists or that you initialized the toolkit with the desired project name")
        else:
            modelContainer = self.getTorchModelContainerByID(modelIDorListID)

            if modelContainer is None:
                raise ValueError(f"Model {modelIDorListID} not found in project {self.projectName}. If this is not the project you ment, make sure caseConfiguration.json exists or that you initialized the toolkit with the desired project name")

        self.packTorchModel(modelContainer,packFileName)

    def packTorchModel(self,modelContainerOrListContainers,packFileName):
        """

        Parameters
        ----------
        modelContainer

        Returns
        -------

        """
        logger = get_classMethod_logger(self,"packModel")
        logger.info("Packing models")
        if not isinstance(modelContainerOrListContainers,Iterable):
            modelContainerOrListContainers = [modelContainerOrListContainers] # make it iterable.

        itemsToZip = []
        modelsDict = dict()
        for modelContainer in modelContainerOrListContainers:
            logger.info(f"Pack model {modelContainer.modelID}")
            modelDoc = modelContainer.getModelDocument()
            logger.debug("Remove all the locations of the classes. They will be found again when loading")
            descWithoutFilePaths = remove_key_recursively_new(modelDoc.desc,"filepath")
            logger.debug(f"Remove all the prefix of the dataset. We assume that it is relative to the files diction of the project (={self.filesDirectory}")
            descWithoutFilePaths = remove_prefix_from_values(descWithoutFilePaths, "pathToData", self.filesDirectory)
            modelDescStr = json.dumps(descWithoutFilePaths,indent=4)
            modelName    = f"{modelContainer.modelName}_{modelContainer.modelID}.json"
            modelsDict[modelName] = modelDescStr
            logger.debug("Pack items")
            itemsToZip.append(modelDoc.getData())

        itemsToZip.append(modelsDict)
        zip_items(packFileName,itemsToZip)

    def loadPackedModel(self,archiveFile,overwrite=False):
        """
            Loads the model to the database, and extracts the runtime data to the directory.

            If the model is in the database, we will overwrite the data if the overwrite is True.
            Otherwise, it will skip.

        Parameters
        ----------
        archiveFile
        overwrite

        Returns
        -------

        """
        logger = get_classMethod_logger(self,"packModel")
        logger.info("Unpacking models")

        models  = list_json_files_in_zip(archiveFile)
        for model in models:
            modelName = model['name']
            logger.info(f"Loading the model {modelName}")
            modelJSON = model['content']

            modelJSON = self.update_classes_filepath(modelJSON)
            modelJSON = self.append_filesDirectory_to_pathToData(modelJSON)

            modelContainer = self.getTorchModelFromJSON(modelJSON['model'])
            modelContainerDoc = modelContainer.getModelDocument()

            # Check if the directory exists.

            # it it does, overwrite only with overwrite=True.




    def append_filesDirectory_to_pathToData(self, modelJSON):
        """

        Parameters
        ----------
        modelJSON

        Returns
        -------

        """
        if isinstance(modelJSON, dict):
            new_dict = {}
            for k, v in modelJSON.items():
                if k == "pathToData" and isinstance(v, str):
                    new_dict[k] = os.path.join(self.filesDirectory,v)
                elif isinstance(v, dict):
                    new_dict[k] = self.append_filesDirectory_to_pathToData(v)
                elif isinstance(v, list):
                    new_dict[k] = [self.append_filesDirectory_to_pathToData(v) if isinstance(item, dict) else item for item in v]
                else:
                    new_dict[k] = v
            return new_dict
        else:
            return modelJSON



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
    def get_class_info(cls, modelClsOrName):
        if isinstance(modelClsOrName, str):
            modelCls = pydoc.locate(modelClsOrName)
        else:
            modelCls = modelClsOrName
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

    @classmethod
    def update_classes_filepath(cls,modelJSON):
        """
        Recursively adds the correct files path to all the classpath of the machine.

        Parameters
        ----------
        modelJSON : dictionary to process

        Returns
        -------

        """

        if isinstance(modelJSON, dict):
            new_dict = {}
            for k, v in modelJSON.items():
                if k == "classpath" and isinstance(v, str):
                    _, fileData = cls.get_class_info(v)
                    new_dict["classpath"] = v
                    new_dict['filepath'] = fileData['filepath']
                elif isinstance(v, dict):
                    new_dict[k] = cls.update_classes_filepath(v)
                elif isinstance(v, list):
                    new_dict[k] = [cls.update_classes_filepath(v) if isinstance(item, dict) else item for item in v]
                else:
                    new_dict[k] = v
            return new_dict
        else:
            return modelJSON



## ====================================================================================================
##
## =====================================     ANALYSIS  ================================================
##
## ====================================================================================================


class analysis:

    datalayer = None
    def __init__(self,datalayer):
        self.datalayer = datalayer


    def sensitivityAnalysis_morris(self,modelContainer,problemContainer,maxEpoch,sampleParameters=dict(),analysisParameters=dict(),parallel=True):
        """
            Performs the sensitivity analysis of Morris to identify the important parameters.

        Parameters
        ----------
        modelContainer : The torch container wrapper model.
            This is the hera container that contains all the data required for the
            dataset, dataloaders and initializing the model.

        problemContainer : JSON
            The JSON that defines the parameters, values and their type for the SALib.
            Created with the SALibUtils.buildSAProblem
        maxEpoch : int
            The number of epochs to train
        sampleParameters : dict
            parameters for the morris sample

        analysisParameters : dict
            parameters for the morris analyze.
        parallel : bool [default=True]
            If tue, use the parallel option if there are more than 2 gpgpue.

        Returns
        -------
            The morris SI as pandas.DataFrame.
        """
        logger = get_classMethod_logger(self,"sensitivityAnalysis_morris")
        morris_sample_parameters = {
            'N': 4,  # Number of trajectories (you can change as needed)
            'num_levels': 4,  # Number of levels in the grid
            'optimal_trajectories': None,  # No trajectory optimization by default
            'local_optimization': False  # No local optimization by default
        }
        morris_sample_parameters.update(sampleParameters)
        morris_analyze_parameters = {
            'num_levels': morris_sample_parameters['num_levels'],  # must match the sample design
            'conf_level': 0.95,  # default confidence level
            'print_to_console': False,  # suppress printing by default
            'num_resamples': 1000,  # bootstrap iterations for CIs
            'seed': None  # random seed (can be set for reproducibility)
        }
        morris_analyze_parameters.update(analysisParameters)

        raw_samples = morris_sample.sample(problemContainer['problem'],**morris_sample_parameters)
        samples = SALibUtils.transformSample(batchList=raw_samples,problemContainer=problemContainer)

        num_gpus = torch.cuda.device_count()
        gpu_ids = [x for x in range(num_gpus)]  # Use both GPUs

        if num_gpus > 1 and parallel:
            logger.info("Preparing the models JSON")

            modelJSONList = []
            for i, sample in enumerate(samples):
                gpu_id = gpu_ids[i % num_gpus]
                modelContainer.setTrainer(devices=[gpu_id])
                baseJson = modelContainer.modelJSON

                paramDict = dict([(name, value) for name, value in zip(problemContainer['problem']['names'], sample)])
                sampleJSON = setJSONPath(base=baseJson, valuesDict=paramDict, inPlace=False)
                modelJSONList.append(sampleJSON)

            Y = Parallel(n_jobs=num_gpus)(
                delayed(_evaluateSample)(i, sample,problemContainer['problem'],modelJSONList[i],maxEpoch) for i, sample in enumerate(samples)
            ) # ,modelContainer
        else:
            Y = numpy.zeros(len(samples))
            for i,sample in enumerate(samples):
                logger.debug(f"Running iteration {i} with sample {sample}")
                # Transfer to a dict of param name -> real value.
                baseJson = modelContainer.modelJSON
                paramDict = dict([(name,value) for name,value in zip(problemContainer['problem']['names'],sample)])
                sampleJSON = setJSONPath(base=baseJson,valuesDict = paramDict,inPlace=False)
                emptyContainer = self.getEmptyTorchModelContainer()
                logger.debug("Fitting the model")

                emptyContainer.modelJSON = sampleJSON
                emptyContainer.fit(maxEpoch,continueTraining=False)
                stats = emptyContainer.getStatistics()

                result = stats.loc[stats.groupby("tag")["step"].idxmax(), ["tag", "value"]].set_index("tag")
                Y[i] = result.loc["val_loss_epoch"].item()

        Si = morris_analyze.analyze(problemContainer['problem'], numpy.array(raw_samples), numpy.array(Y), **morris_analyze_parameters)
        return pandas.DataFrame(Si)

## ====================================================================================================
##
## =====================================     GLOBAL  ================================================
##
## ====================================================================================================

def _evaluateSample(i,sample,problem,modelJSON,maxEpoch): #,modelContainer):
    """
        Used to created the JSON and run the fit for the parallel running over couple of GPU.
    Parameters
    ----------
    gpu_id
    sample
    problem
    modelContainer

    Returns
    -------

    """
    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.MACHINELEARNING_DEEPLEARNING)
    modelContainer = tk.getEmptyTorchModelContainer()
    modelContainer.modelJSON = modelJSON
    modelContainer.fit(maxEpoch, continueTraining=False)
    stats = modelContainer.getStatistics()

    result = stats.loc[stats.groupby("tag")["step"].idxmax(), ["tag", "value"]].set_index("tag")
    return result.loc["val_loss_epoch"].item()



def remove_key_recursively_new(d, key_to_remove):
    """
    Return a new dict with key `key_to_remove` removed from all nested dicts.

    :param d: original dictionary
    :param key_to_remove: key to remove
    :return: new dictionary with the key removed
    """
    if isinstance(d, dict):
        new_dict = {}
        for k, v in d.items():
            if k != key_to_remove:
                if isinstance(v, dict):
                    new_dict[k] = remove_key_recursively_new(v, key_to_remove)
                elif isinstance(v, list):
                    new_dict[k] = [
                        remove_key_recursively_new(item, key_to_remove) if isinstance(item, dict) else item
                        for item in v
                    ]
                else:
                    new_dict[k] = v
        return new_dict
    else:
        return d  # non-dict objects are returned as-is


def remove_prefix_from_values(d, target_key, prefix):
    """
    Recursively remove a prefix from the value of all occurrences of `target_key` in a dict.

    Does not take into accoutn JSON files that have list in their root.

    :param d: dictionary to process
    :param target_key: key whose values will have the prefix removed
    :param prefix: prefix to remove from the value
    :return: new dictionary with updated values
    """
    if isinstance(d, dict):
        new_dict = {}
        for k, v in d.items():
            if k == target_key and isinstance(v, str) and v.startswith(prefix):
                new_dict[k] = v[len(prefix)+1:]
            elif isinstance(v, dict):
                new_dict[k] = remove_prefix_from_values(v, target_key, prefix)
            elif isinstance(v, list):
                new_dict[k] = [
                    remove_prefix_from_values(item, target_key, prefix) if isinstance(item, dict) else item
                    for item in v
                ]
            else:
                new_dict[k] = v
        return new_dict
    else:
        return d


