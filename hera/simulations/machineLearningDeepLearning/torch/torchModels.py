import pandas
import pydoc

try:
    import torch
except ImportError:
    print("pyTorch not Found. Cannot uyse this toolkit")

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

class torchLightingModel(Project):
    MODEL= "torchModel"

    machineLearningDeepLearning =None

    modelJSON = None


    @property
    def modelName(self):
        return self.modelJSON['model']['classpath'].replace(".", "_")

    def __init__(self,mldlModels):

        self.machineLearningDeepLearning = mldlModels
        super().__init__(projectName=mldlModels.projectName,filesDirectory=mldlModels.filesDirectory)
        self.initModel()

    def initModel(self):
        self.modelJSON = dict()
        self.modelJSON['dataset'] = dict()
        self.modelJSON['trainDataset'] = dict()
        self.modelJSON['validateDataset'] = dict()
        self.modelJSON['trainer'] = dict()
        self.modelJSON['model'] = dict()
        self.modelJSON['checkpoint'] = dict()

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
        name,info = self.get_class_info(datasetClass)
        params = self.get_init_params(datasetClass)
        if 'kwargs' in params:
            del params['kwargs']

        params.update(kwargs)
        info['parameters'] = params
        self.modelJSON['dataset'][datasetName] = info

    def setTrainDataLoader(self, datasetName,**kwargs):
        name,info = self.get_class_info(DataLoader)
        params = self.get_init_params(DataLoader)
        info['dataset'] = datasetName
        params.update(kwargs)
        info['parameters'] = params
        self.modelJSON['trainDataset'] = info

    def setTValidateDataLoader(self,datasetName,**kwargs):
        name, info = self.get_class_info(DataLoader)
        params = self.get_init_params(DataLoader)
        info['dataset'] = datasetName
        params.update(kwargs)
        info['parameters'] = params
        self.modelJSON['validateDataset'] = info

    def setModel(self,modelClass,**kwargs):
        name,info = self.get_class_info(modelClass)
        params = self.get_init_params(modelClass)
        if 'kwargs' in params:
            del params['kwargs']
        params.update(**kwargs)
        info['parameters'] = params
        self.modelJSON['model'] = info

    def setTrainer(self,val_check_interval=1,**kwargs):
        name, info = self.get_class_info(pl.Trainer)
        params = self.get_init_params(pl.Trainer)
        params.update(kwargs)
        info['parameters'] = params
        self.modelJSON['trainer'] = trainer=info


    def setCheckPoint(self,**kwargs):
        name, info = self.get_class_info(ModelCheckpoint)
        params = self.get_init_params(ModelCheckpoint)
        params.update(kwargs)
        info['parameters'] = params
        self.modelJSON['checkpoint'] = info


    def fit(self,max_epochs,continueTraining=True):
        """
            Initializes all the object and returns the trainer.
        Returns
        -------

        """
        doc = self.getModelDocument()

        # 1. Initialize the dataloaders .
        trainDatasetLoader = self.getDatasetLoader(self.modelJSON['trainDataset'])
        validateDatasetLoader = self.getDatasetLoader(self.modelJSON['validateDataset'])


        model = self.initClass(self.modelJSON['model'])
        trainer = self.getTrainer(max_epochs=max_epochs,doc=doc)

        ckpt_path = os.path.join(doc.getData(),f"{self.modelName}.ckpt")

        # Remove it before training starts
        if continueTraining:
            if not os.path.exists(ckpt_path):
                ckpt_path = None
        elif os.path.exists(ckpt_path):
            os.remove(ckpt_path)
            ckpt_path = None

        trainer.fit(model,trainDatasetLoader,validateDatasetLoader,ckpt_path=ckpt_path)


    def getStatistics(self):
        doc = self.getModelDocument()

        event_files = sorted(glob(os.path.join(doc.getData(),"version_0", "events.out.tfevents.*")))
        if not event_files:
            raise FileNotFoundError("No TensorBoard event files found.")

        event_file = event_files[0]
        ea = event_accumulator.EventAccumulator(event_file)
        ea.Reload()

        tags = ea.Tags()["scalars"]
        df_list = []

        for tag in tags:
            events = ea.Scalars(tag)
            temp_df = pandas.DataFrame({
                "tag": tag,
                "step": [e.step for e in events],
                "value": [e.value for e in events],
                "wall_time": [e.wall_time for e in events]
            })
            df_list.append(temp_df)

        df = pandas.concat(df_list, ignore_index=True)
        return df

    def getModelDocument(self):
        qry = dictToMongoQuery(self.modelJSON)
        docList = self.getSimulationsDocuments(type=self.MODEL, **qry)

        if len(docList) == 0:

            modelID = self.getCounterAndAdd(self.MODEL)
            resource = os.path.join(self.filesDirectory, "modelData", f"{self.modelName}_{modelID}")
            desc = self.modelJSON
            desc['modelID'] = modelID
            doc = self.addSimulationsDocument(type=self.MODEL,
                                              resource=resource,
                                              dataFormat=self.datatypes.STRING,
                                              desc=self.modelJSON)
        else:
            doc = docList[0]
        return doc


    def getTrainer(self,doc=None,**kwargs):

        if doc is None:
            doc = self.getModelDocument()

        savedir = doc.getData()
        filename = self.modelName

        trainerJSON = self.modelJSON['trainer']
        logger = TensorBoardLogger(savedir, name=None, version=0)
        checkpoint_callback = self.initClass(self.modelJSON['checkpoint'],dirpath=savedir,filename=filename,monitor = "val_loss")

        trainer    = self.getClass(trainerJSON)
        params     = trainerJSON['parameters']
        params['enable_model_summary'] = True
        params['limit_val_batches']= 1.0
        params.update(**kwargs)
        params['logger'] = logger
        params['callbacks'] = [checkpoint_callback]

        return trainer(**params)



    def getDatasetLoader(self, JSONdesc):
        """
            The
        Parameters
        ----------
        JSONdesc

        Returns
        -------

        """
        datasetName = JSONdesc['dataset']
        dataset = self.initClass(self.modelJSON['dataset'][datasetName])
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

    def getModel(self,modelName,**hyperParameters):
        """
            Get the model and load it.
        Parameters
        ----------
        modelName
        hyperParameters

        Returns
        -------

        """
        doc = self.getModelDocument(modelName,**hyperParameters)
        clsPath = doc.desc['modelPath']
        modelName = doc.desc['modelName']
        hyperParameters  = doc.desc['hyperParameters']

        os.path.append(clsPath)

        mdlCls = pydoc.locate(modelName)

        return modelName(**hyperParameters)



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

    def get_class_info(self,cls):
        module = cls.__module__
        name = cls.__name__
        file_path = inspect.getfile(cls)
        file_path = os.path.dirname(os.path.abspath(file_path))

        full_path = f"{module}.{name}"
        return name, dict(classpath=full_path, filepath=file_path)






# class torchModels(Project):
#
#     DEVICE_GPGPU="cuda"
#     DEVICE_CPU = "cpu"
#
#     TYPE_MODEL = "torchModel_results"
#
#     def __init__(self,projectName,filesDirectory):
#         """
#             Initializes the model class.
#
#         Parametersmonth
#         ----------
#         projectName : str
#             The name of the project that this class belongs to.
#
#         filesDirectory : str
#             The directory to write the files in.
#         """
#         super().__init__(projectName=projectName,filesDirectory=filesDirectory)
#         os.makedirs(os.path.join(self.filesDirectory,"torchModels","modelWeights"),exist_ok=True)
#
#
#     def getTrainingStrategy(self,name,dataLoader,optimizer,**kwargs):
#         """
#             Creates an instance of the requested training stragety.
#
#         Parameters
#         ----------
#         dataLoader : obj
#             The dataloader to use.
#         optimizer : obj
#             The optimizer to use
#         kwargs:
#             Specific parameters for the strategy.
#
#             Will allow to set the evaulation function in the reinforcement learning strategy.
#
#         Returns
#         -------
#             A trainingStrategy object.
#         """
#
#     def trainModel(self,modelInstance,trainingStrategy,totalEpoch,startFromCheckpoint=True,deviceName=None):
#         """
#             Save a model to the DB.
#
#             The procedure saves:
#                 1. The value of the hyper parameters.
#                     We assume that all were specified in the constructor of the class.
#                 2. The code of the model.
#                 3. Checkpoints are saved automatically, with the current epoch.
#
#             The procedure will attempt to train in parallel on all the gpus available.
#
#         Parameters
#         ----------
#         modelInstance : obj
#             A class of pyTorch nn.Model.
#
#         trainingStrategy: obj
#             A class that implements the training strategy. Whether just plain vanilla, autoregressive (feeds the
#             output step by step to the model), or Reinforcement learning (that would require the evaluation function).
#
#             The trainingStrategy also includes the batch to work on.
#
#         epoch : int
#             The total epochs to train
#         startFromCheckpoint : bool [default = True]
#             If true, try to load the existing
#
#         deviceName: str
#             cuda or cpu.
#             if None, choose cuda if available.
#             Use the DEVICE_GPGPU or DEVICE_CPU constants.
#
#         Returns
#         -------
#
#         """
#         logger = get_classMethod_logger(self,"trainModel")
#         deviceName = "cuda" if torch.cuda.is_available() else "cpu" if deviceName is None else deviceName
#         logger.info(f"training the model using device {deviceName}")
#
#         device = torch.device(deviceName)
#         parallelTrain = False
#         if deviceName == self.DEVICE_GPGPU:
#             gpu_count = torch.cuda.device_count()
#             if gpu_count >1:
#                 logger.info(f"Found {gpu_count} gpus, training in parallel.")
#                 parallelTrain = True
#
#
#         signatureMap = dictToMongoQuery(self.get_constructor_args(model),prefix=hyperParameters)
#         modelName,_ = self.get_class_info(model)
#
#
#         docList = self.getSimulationsDocuments(type=TYPE_MODEL,modelName=modelName,**signatureMap)
#         if len(docList) == 0:
#             logger.debug(f"Not found. Training with the requested epoch and adding to the database")
#             cntr = self.getCounterAndAdd(modelName)
#             newFileName = os.path.join(self.filesDirectory,f"{modelName.replace(".","_")}_{cntr}.pth")
#             doc = self.addSimulationsDocument(type=self.TYPE_MODEL,
#                                               resource=newFileName,
#                                               dataFormat=self.datatypes.STRING,
#                                               desc=dict(
#                                                   modelName=modelName,
#                                                   hyperParameters=signatureMap,
#                                                   training=dict(state=trainingStrategy.getState(),
#                                                                 epoch = epoch,
#                                                                 loss  = 1e10
#                                                                 )
#                                             )
#             )
#             current_epoch = 0
#         else:
#             logger.debug(f"Found. Loading the checkpoint and continue from there")
#             doc = docList[0]
#             current_epoch = doc.desc['training']['epoch']
#             trainingStrategy.loadState(doc.desc['training'])
#             model.load_state_dict(doc.getData())
#
#         for epoch in range(current_epoch,totalEpoch):
#
#             avg_loss = trainingStrategy.train(model)
#             print(f"Epoch {epoch + 1}/{num_epochs}, Loss: {avg_loss:.6f}")
#
#             # Save checkpoint after every epoch
#             trainingStrategy.storeState(doc,epoch)
#
#         return avg_loss
#
#     def listModels(self,modelName=None,**hyperParameters):
#         """
#             Returns the list of models and hyper parameters that are in the project.
#         Parameters
#         ----------
#         modelName
#         hyperParameters
#
#         Returns
#         -------
#
#         """
#         signatureMap = dictToMongoQuery(hyperParameters, prefix=hyperParameters)
#
#         docList = self.getSimulationsDocuments(type=TYPE_MODEL, modelName=modelName, **signatureMap)
#         modelList = []
#         for doc in docList:
#             item  = convertJSONtoPandas(doc.desc['hyperParameters'])
#             item.assign(doc.desc['modelName'])
#             modelList,append(item)
#
#         return pandas.concat(modelList,ignore_index=True)
#
#     def getModel(self,modelName,**hyperParameters):
#         """
#             Gets a model from the DB and loads it.
#         Parameters
#         ----------
#         modelName
#         hyperParameters
#
#         Returns
#         -------
#
#         """
#         signatureMap = dictToMongoQuery(hyperParameters, prefix=hyperParameters)
#
#         docList = self.getSimulationsDocuments(type=TYPE_MODEL, modelName=modelName, **signatureMap)
#         if len(docList)==0:
#             return None
#         return docList[0]
#
#
#     ##########################################################
#     ##
#     ##              Utility functions.
#     ##
#     ##########################################################
#     def save_checkpoint(self,model, optimizer, epoch, loss, path):
#         logger = get_classMethod_logger(self,"save_checkpoint")
#         torch.save({
#             'epoch': epoch,
#             'model_state_dict': model.state_dict(),
#             'optimizer_state_dict': optimizer.state_dict(),
#             'loss': loss
#         }, path)
#         logger.info(f"Saved the checkpoint to {path}")
#
#     def load_checkpoint(self,path, model, optimizer):
#         logger = get_classMethod_logger(self,"load_checkpoint")
#         if os.path.exists(path):
#             checkpoint = torch.load(path)
#             model.load_state_dict(checkpoint['model_state_dict'])
#             optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
#             epoch = checkpoint['epoch']
#             loss = checkpoint['loss']
#             logger.info(f"Loaded checkpoint from epoch {epoch + 1} with loss {loss:.6f}")
#             return epoch + 1  # resume from next epoch
#         else:
#             print("No checkpoint found. Starting from scratch.")
#             return 0
