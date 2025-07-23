import pandas

try:
    import torch
except ImportError:
    print("pyTorch not Found. Cannot uyse this toolkit")

import inspect
import os
from hera.datalayer import Project
from hera.utils.logging import get_classMethod_logger
from hera.utils import dictToMongoQuery
from hera.utils import convertJSONtoPandas


class torchModels(Project):

    DEVICE_GPGPU="cuda"
    DEVICE_CPU = "cpu"

    TYPE_MODEL = "torchModel_results"

    def __init__(self,projectName,filesDirectory):
        """
            Initializes the model class.

        Parametersmonth
        ----------
        projectName : str
            The name of the project that this class belongs to.

        filesDirectory : str
            The directory to write the files in.
        """
        super().__init__(projectName=projectName,filesDirectory=filesDirectory)
        os.makedirs(os.path.join(self.filesDirectory,"torchModels","modelWeights"),exist_ok=True)


    def getTrainingStrategy(self,name,dataLoader,optimizer,**kwargs):
        """
            Creates an instance of the requested training stragety.

        Parameters
        ----------
        dataLoader : obj
            The dataloader to use.
        optimizer : obj
            The optimizer to use
        kwargs:
            Specific parameters for the strategy.

            Will allow to set the evaulation function in the reinforcement learning strategy.

        Returns
        -------
            A trainingStrategy object.
        """

    def trainModel(self,modelInstance,trainingStrategy,totalEpoch,startFromCheckpoint=True,deviceName=None):
        """
            Save a model to the DB.

            The procedure saves:
                1. The value of the hyper parameters.
                    We assume that all were specified in the constructor of the class.
                2. The code of the model.
                3. Checkpoints are saved automatically, with the current epoch.

            The procedure will attempt to train in parallel on all the gpus available.

        Parameters
        ----------
        modelInstance : obj
            A class of pyTorch nn.Model.

        trainingStrategy: obj
            A class that implements the training strategy. Whether just plain vanilla, autoregressive (feeds the
            output step by step to the model), or Reinforcement learning (that would require the evaluation function).

            The trainingStrategy also includes the batch to work on.

        epoch : int
            The total epochs to train
        startFromCheckpoint : bool [default = True]
            If true, try to load the existing

        deviceName: str
            cuda or cpu.
            if None, choose cuda if available.
            Use the DEVICE_GPGPU or DEVICE_CPU constants.

        Returns
        -------

        """
        logger = get_classMethod_logger(self,"trainModel")
        deviceName = "cuda" if torch.cuda.is_available() else "cpu" if deviceName is None else deviceName
        logger.info(f"training the model using device {deviceName}")

        device = torch.device(deviceName)
        parallelTrain = False
        if deviceName == self.DEVICE_GPGPU:
            gpu_count = torch.cuda.device_count()
            if gpu_count >1:
                logger.info(f"Found {gpu_count} gpus, training in parallel.")
                parallelTrain = True


        signatureMap = dictToMongoQuery(self.get_constructor_args(model),prefix=hyperParameters)
        modelName,_ = self.get_class_info(model)


        docList = self.getSimulationsDocuments(type=TYPE_MODEL,modelName=modelName,**signatureMap)
        if len(docList) == 0:
            logger.debug(f"Not found. Training with the requested epoch and adding to the database")
            cntr = self.getCounterAndAdd(modelName)
            newFileName = os.path.join(self.filesDirectory,f"{modelName.replace(".","_")}_{cntr}.pth")
            doc = self.addSimulationsDocument(type=self.TYPE_MODEL,
                                              resource=newFileName,
                                              dataFormat=self.datatypes.STRING,
                                              desc=dict(
                                                  modelName=modelName,
                                                  hyperParameters=signatureMap,
                                                  training=dict(state=trainingStrategy.getState(),
                                                                epoch = epoch,
                                                                loss  = 1e10
                                                                )
                                            )
            )
            current_epoch = 0
        else:
            logger.debug(f"Found. Loading the checkpoint and continue from there")
            doc = docList[0]
            current_epoch = doc.desc['training']['epoch']
            trainingStrategy.loadState(doc.desc['training'])
            model.load_state_dict(doc.getData())

        for epoch in range(current_epoch,totalEpoch):

            avg_loss = trainingStrategy.train(model)
            print(f"Epoch {epoch + 1}/{num_epochs}, Loss: {avg_loss:.6f}")

            # Save checkpoint after every epoch
            trainingStrategy.storeState(doc,epoch)

        return avg_loss

    def listModels(self,modelName=None,**hyperParameters):
        """
            Returns the list of models and hyper parameters that are in the project.
        Parameters
        ----------
        modelName
        hyperParameters

        Returns
        -------

        """
        signatureMap = dictToMongoQuery(hyperParameters, prefix=hyperParameters)

        docList = self.getSimulationsDocuments(type=TYPE_MODEL, modelName=modelName, **signatureMap)
        modelList = []
        for doc in docList:
            item  = convertJSONtoPandas(doc.desc['hyperParameters'])
            item.assign(doc.desc['modelName'])
            modelList,append(item)

        return pandas.concat(modelList,ignore_index=True)

    def getModel(self,modelName,**hyperParameters):
        """
            Gets a model from the DB and loads it.
        Parameters
        ----------
        modelName
        hyperParameters

        Returns
        -------

        """
        signatureMap = dictToMongoQuery(hyperParameters, prefix=hyperParameters)

        docList = self.getSimulationsDocuments(type=TYPE_MODEL, modelName=modelName, **signatureMap)
        if len(docList)==0:
            return None
        return docList[0]


    ##########################################################
    ##
    ##              Utility functions.
    ##
    ##########################################################
    def save_checkpoint(self,model, optimizer, epoch, loss, path):
        logger = get_classMethod_logger(self,"save_checkpoint")
        torch.save({
            'epoch': epoch,
            'model_state_dict': model.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'loss': loss
        }, path)
        logger.info(f"Saved the checkpoint to {path}")

    def load_checkpoint(self,path, model, optimizer):
        logger = get_classMethod_logger(self,"load_checkpoint")
        if os.path.exists(path):
            checkpoint = torch.load(path)
            model.load_state_dict(checkpoint['model_state_dict'])
            optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
            epoch = checkpoint['epoch']
            loss = checkpoint['loss']
            logger.info(f"Loaded checkpoint from epoch {epoch + 1} with loss {loss:.6f}")
            return epoch + 1  # resume from next epoch
        else:
            print("No checkpoint found. Starting from scratch.")
            return 0

    def get_constructor_args(self,instance):
        cls = instance.__class__
        try:
            sig = inspect.signature(cls.__init__)
        except (TypeError, ValueError):
            return {}

        bound_args = {}
        for name, param in sig.parameters.items():
            if name == "self":
                continue
            if hasattr(instance, name):
                bound_args[name] = getattr(instance, name)
        return bound_args


    def get_class_info(self,instance):
        cls = instance.__class__
        module = cls.__module__
        name = cls.__name__
        full_path = f"{module}.{name}"
        return full_path, name
