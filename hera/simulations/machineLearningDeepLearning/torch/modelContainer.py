import pandas
import pydoc
import copy
from hera.simulations.machineLearningDeepLearning.LightingModule import LightningModuleHera

try:
    import torch
    from pytorch_lightning.utilities.exceptions import MisconfigurationException
except ImportError:
    print("pyTorch or pytorch lighting are not Found. Cannot uyse this toolkit")

import inspect
import os
import sys
import pydoc
from hera.datalayer import Project
from hera.utils.logging import get_classMethod_logger
from torch.utils.data import DataLoader
from pytorch_lightning.callbacks import ModelCheckpoint
from tensorboard.backend.event_processing import event_accumulator
from glob import glob

from hera.toolkit import abstractToolkit
import pytorch_lightning as pl
from pytorch_lightning.loggers import TensorBoardLogger
from hera.utils import dictToMongoQuery

class torchLightingModelContainer(Project):
    MODEL= "torchModel"

    machineLearningDeepLearning =None

    _modelDict = None
    modelID = None
    modelResource = None

    initFromState = None

    @property
    def modelDict(self):
        return self._modelDict

    @modelDict.setter
    def modelDict(self, value):
        self._modelDict = value
        self.modelID = None

    @property
    def modelJSON(self):
        return self._modelDict

    @modelJSON.setter
    def modelJSON(self, value):
        self._modelDict = value
        self.modelID = None

    @property
    def modelName(self):
        return self.modelJSON['model']['classpath'].replace(".", "_")

    def __init__(self,mldlModels, connectionName=None):

        self.machineLearningDeepLearning = mldlModels
        super().__init__(projectName=mldlModels.projectName,filesDirectory=mldlModels.filesDirectory, connectionName=connectionName)
        self.state_other_models = []
        self.initModel()

    def initModel(self):
        self.modelJSON = dict()
        self.modelJSON['dataset'] = dict()
        self.modelJSON['trainDataset'] = dict()
        self.modelJSON['validateDataset'] = dict()
        self.modelJSON['trainer'] = dict()
        self.modelJSON['model'] = dict()
        self.modelJSON['checkpoint'] = dict()

    def loadOtherComponent(self,componentName,otherModelID,nameOnOtherModel):
        self.state_other_models.append(dict(componentName=componentName,nameOnOtherModel=nameOnOtherModel,otherModelID=otherModelID))


    def setDataset(self,datasetName,datasetClass,**kwargs):
        """
            Adds the dataset name and class to the

        Parameters
        ----------
        datasetName : str`
        datasetClass : obj`
        kwargs : the constuctor parameters of class.

        Returns
        -------

        """
        name,info = self.machineLearningDeepLearning.get_class_info(datasetClass)
        params = self.get_init_params(datasetClass)
        if 'kwargs' in params:
            del params['kwargs']
        if datasetName in self.modelJSON['dataset']:
            params.update(self.modelJSON['dataset'][datasetName]['parameters'])

        params.update(kwargs)
        info['parameters'] = params
        self.modelJSON['dataset'][datasetName] = info

    def setTrainDataLoader(self, datasetName,**kwargs):
        name,info = self.machineLearningDeepLearning.get_class_info(DataLoader)
        params = self.get_init_params(DataLoader)
        info['dataset'] = datasetName
        if  'parameters' in self.modelJSON['trainDataset']:
            params.update(self.modelJSON['trainDataset']['parameters'])

        params.update(kwargs)
        info['parameters'] = params
        self.modelJSON['trainDataset'] = info

    def setTValidateDataLoader(self,datasetName,**kwargs):
        name, info = self.machineLearningDeepLearning.get_class_info(DataLoader)
        params = self.get_init_params(DataLoader)
        info['dataset'] = datasetName
        if  'parameters' in self.modelJSON['validateDataset']:
            params.update(self.modelJSON['validateDataset']['parameters'])

        params.update(kwargs)
        info['parameters'] = params
        self.modelJSON['validateDataset'] = info

    def setModel(self,modelClass,**kwargs):
        name,info = self.machineLearningDeepLearning.get_class_info(modelClass)
        params = self.get_init_params(modelClass)
        if 'kwargs' in params:
            del params['kwargs']

        if 'parameters'  in self.modelJSON['model']:
            params.update(self.modelJSON['model']['parameters'])


        params.update(**kwargs)
        info['parameters'] = params
        self.modelJSON['model'] = info

    def setTrainer(self,val_check_interval=1,**kwargs):
        name, info = self.machineLearningDeepLearning.get_class_info(pl.Trainer)
        params = self.get_init_params(pl.Trainer)
        if 'parameters' in self.modelJSON['trainer']:
            params.update(self.modelJSON['trainer']['parameters'])

        params.update(kwargs)
        info['parameters'] = params
        self.modelJSON['trainer'] = info

    def setCheckPoint(self,**kwargs):
        name, info = self.machineLearningDeepLearning.get_class_info(ModelCheckpoint)
        params = self.get_init_params(ModelCheckpoint)
        params.update(kwargs)
        info['parameters'] = params
        self.modelJSON['checkpoint'] = info

    def load(self):
        logger = get_classMethod_logger(self,"load")
        logger.info(f"Loadin model. Current ID is {self.modelID}")
        if self.modelID is None:
            logger.debug(f"Loading or creating a new model in DB")
            doc = self.getModelDocument()
            self.modelJSON = copy.deepcopy(doc.desc['model'])
            self.modelID = doc.desc['modelID']
            self.modelResource = doc.getData()
            logger.debug(f"Set up the new model {self.modelID}")

    @property
    def max_epoch(self):
        ckpt_path = self.checkpoint_path
        ckpt = torch.load(ckpt_path, map_location="cpu")
        return ckpt["epoch"]+1

    @property
    def checkpoint_path(self):
        return os.path.join(self.modelResource,f"{self.modelName}.ckpt")

    def fit(self,max_epochs,continueTraining=True):
        """
            Initializes all the object and returns the trainer.
        Returns
        -------

        """
        logger = get_classMethod_logger(self, "fit")
        logger.info("Loading the model ")
        self.load()
        logger.info(f"Loaded  model ID {self.modelID}")


        # 1. Initialize the dataloaders .
        logger.info("Getting validation dataset")
        validateDatasetLoader = self.getValidateDatasetLoader()

        logger.info("Initializing model")
        model = self.initClass(self.modelJSON['model'])

        logger.info("Registering model container if necessary")
        if isinstance(model,LightningModuleHera):
            model.setModelContainer(self)

        if hasattr(model,"train_dataloader"):
            trainDatasetLoader = None
        else:
            trainDatasetLoader = self.getTrainDatasetLoader()

        trainer = self.getTrainer(max_epochs=max_epochs)

        ckpt_path_param = dict()
        ckpt_path = None

        if len(self.state_other_models) == 0:
            if continueTraining:
                ckpt_path = self.checkpoint_path
                if not os.path.exists(ckpt_path):
                    ckpt_path = None
            elif os.path.exists(self.checkpoint_path):
                    os.remove(self.checkpoint_path)

            ckpt_path_param['ckpt_path'] = ckpt_path

        else:
            # load the weights from another model.
            for otherState in self.state_other_models:
                if otherState['componentName'] == self.MODEL:
                    component = model
                else:
                    component = getattr(model,otherState['componentName'])
                prefix   = otherState['nameOnOtherModel']
                otherChkpnt = self.machineLearningDeepLearning.getTorchModelContainerByID(otherState['otherModelID']).checkpoint_path
                self.load_submodule_from_ckpt(component,otherChkpnt,prefix)

        logger.info("Starting to train")
        trainer.fit(model,val_dataloaders=validateDatasetLoader,train_dataloaders=trainDatasetLoader,**ckpt_path_param)

    def only_validate(self,continueTraining=True):
        """
            Initializes all the object and returns the trainer.
        Returns
        -------

        """
        logger = get_classMethod_logger(self, "fit")
        logger.info("Loading the model ")
        self.load()
        logger.info(f"Loaded  model ID {self.modelID}")

        # 1. Initialize the dataloaders .

        logger.info("Getting validation dataset")
        validateDatasetLoader = self.getValidateDatasetLoader()

        logger.info("Initializing model")
        model = self.initClass(self.modelJSON['model'])

        if isinstance(model,LightningModuleHera):
            model.setModelContainer(self)

        # if hasattr(model,"train_dataloader"):
        #     trainDatasetLoader = None
        # else:
        #     trainDatasetLoader = self.getTrainDatasetLoader()

        trainer = self.getTrainer(max_epochs=1)

        ckpt_path_param = dict()
        ckpt_path = None

        if len(self.state_other_models) == 0:
            if continueTraining:
                ckpt_path = self.checkpoint_path
                if not os.path.exists(ckpt_path):
                    ckpt_path = None
            elif os.path.exists(self.checkpoint_path):
                    os.remove(self.checkpoint_path)

            ckpt_path_param['ckpt_path'] = ckpt_path

        else:
            # load the weights from another model.
            for otherState in self.state_other_models:
                if otherState['componentName'] == self.MODEL:
                    component = model
                else:
                    component = getattr(model,otherState['componentName'])
                prefix   = otherState['nameOnOtherModel']
                otherChkpnt = self.machineLearningDeepLearning.getTorchModelContainerByID(otherState['otherModelID']).checkpoint_path
                self.load_submodule_from_ckpt(component,otherChkpnt,prefix)

        logger.info("Starting to train")
        trainer.validate(model, dataloaders=validateDatasetLoader)

    def getStatistics(self):
        if self.modelJSON is None:
            doc = self.getModelDocument()
        event_files_path = os.path.join(self.modelResource,"version_0", "events.out.tfevents.*")
        event_files = sorted(glob(event_files_path))
        if not event_files:
            raise FileNotFoundError(f"No TensorBoard event files found in {event_files_path} ")

        df_list = []
        for event_file in event_files:
            ea = event_accumulator.EventAccumulator(event_file)
            ea.Reload()

            tags = ea.Tags()["scalars"]

            for tag in tags:
                events = ea.Scalars(tag)
                temp_df = pandas.DataFrame({
                    "tag": tag,
                    "step": [e.step for e in events],
                    "value": [e.value for e in events],
                    "wall_time": [e.wall_time for e in events],
                    "max_epoch": self.max_epoch
                })
                df_list.append(temp_df)

        df = pandas.concat(df_list, ignore_index=True).sort_values("step")
        return df

    def getModelDocument(self):
        logger = get_classMethod_logger(self, "getModelDocument")
        qry = dictToMongoQuery(self.modelJSON,prefix="model")
        docList = self.getSimulationsDocuments(type=self.MODEL, **qry)
        if len(docList) == 0:
            modelID = self.getCounterAndAdd(self.MODEL)
            logger.debug(f"Model not found. Creating new with ID: {modelID} ")
            resource = os.path.join(self.filesDirectory, "modelData", f"{self.modelName}_{modelID}")
            doc = self.addSimulationsDocument(type=self.MODEL,
                                              resource=resource,
                                              dataFormat=self.datatypes.STRING,
                                              desc=dict(model=self.modelJSON,modelID=modelID))
            logger.debug(f"Model added with ID: {doc.desc['modelID']} ")
        else:
            doc = docList[0]
            logger.debug(f"Model found with ID: {doc.desc['modelID']} ")
        return doc


    def getTrainer(self,**kwargs):

        if self.modelResource is None:
            self.load()
        else:
            savedir = self.modelResource

        filename = self.modelName
        trainerJSON = self.modelJSON['trainer']
        logger = TensorBoardLogger(savedir, name=None, version=0)
        checkpoint_callback = self.initClass(self.modelJSON['checkpoint'],dirpath=savedir,filename=filename)

        trainer    = self.getClass(trainerJSON)
        params     = trainerJSON['parameters']
        params['enable_model_summary'] = True
        params.update(**kwargs)
        params['logger'] = logger
        params['callbacks'] = [checkpoint_callback]
        return trainer(**params)

    def getModel(self):
        model =  self.initClass(self.modelJSON['model'])
        if os.path.exists(self.checkpoint_path):
            checkpoint = torch.load(self.checkpoint_path)
            model.load_state_dict(checkpoint['state_dict'])
            model.eval()
        return model
        # now loading the state.

    def getTrainDatasetLoader(self):
        """
            Return the trainer class with the relevant dataset.
        Returns
        -------

        """
        return self._getDatasetLoader(self.modelJSON['trainDataset'])


    def getValidateDatasetLoader(self):
        """
            Returns the loader class with the relevant dataaset.
        Returns
        -------

        """
        return self._getDatasetLoader(self.modelJSON['validateDataset'])


    def getDatasetNames(self):
        return [x for x in self.modelJSON['dataset'].keys()]

    def getDataSet(self,datasetName,**kwargs):
        """
            Returns the dataset
        Parameters
        ----------
        datasetName

        Returns
        -------

        """
        datasetDict = dict(self.modelJSON['dataset'][datasetName])
        datasetDict['parameters'].update(kwargs)

        return self.initClass(datasetDict)


    def _getDatasetLoader(self, JSONdesc):
        """
            The
        Parameters
        ----------
        JSONdesc

        Returns
        -------

        """
        datasetName = JSONdesc['dataset']
        dataset = self.getDataSet(datasetName)
        datasetLoader = self.initClass(JSONdesc,dataset=dataset)
        return datasetLoader


    def getClass(self,JSONdesc):
        """
            Loads the class from the JSON description.
        Parameters
        ----------
        JSONdesc

        Returns
        -------

        """
        filePath = JSONdesc['filepath']
        classPath = JSONdesc['classpath']
        if filePath not in sys.path:
            sys.path.append(filePath)

        clss = pydoc.locate(classPath)
        return clss

    def initClass(self,JSONdesc,**kwargs):
        """
            load the class and initialize it initialization of the class.
        Parameters
        ----------
        JSONdesc

        Returns
        -------

        """
        clss = self.getClass(JSONdesc)
        params = JSONdesc['parameters'].copy()
        params.update(**kwargs)
        return clss(**params)


    def get_init_params(self,cls):
        init = cls.__init__
        sig = inspect.signature(init)

        params = {}
        for name, param in sig.parameters.items():
            if name == 'self':
                continue
            if param.default is inspect.Parameter.empty:
                params[name] = None
            else:
                params[name] = param.default
        return params


    def load_submodule_from_ckpt(self, submodule, ckpt_path, prefix=None, strict=False):
        """
        Load weights from a Lightning checkpoint into a submodule.

        Args:
            submodule (nn.Module): the part of your model (e.g. model.encoder).
            ckpt_path (str): path to the Lightning checkpoint (.ckpt).
            prefix (str or None): prefix of the submodule inside the checkpoint
                                  (e.g. "encoder." or "decoder.").
                                  If None, tries to match directly.
            strict (bool): whether to enforce that all keys match.
            map_location (str): device mapping for torch.load.

        Returns:
            missing_keys, unexpected_keys: lists from load_state_dict
        """

        map_location = "cpu"
        ckpt = torch.load(ckpt_path, map_location=map_location)
        state_dict = ckpt["state_dict"]

        if prefix is not None:
            if prefix[-1] != ".":
                prefix += '.'
            # keep only weights belonging to this submodule
            sub_state = {k.replace(prefix, "", 1): v
                         for k, v in state_dict.items()
                         if k.startswith(prefix)}
        else:
            sub_state = state_dict

        missing, unexpected = submodule.load_state_dict(sub_state, strict=strict)
        return missing, unexpected





