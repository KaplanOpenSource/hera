import json
import os
import pandas
import inspect
from promise.utils import deprecated

from hera.datalayer.datahandler import datatypes
from hera.utils.logging import get_classMethod_logger
from hera.utils import loadJSON
from hera import toolkit
from hera.datalayer.collection import AbstractCollection,\
    Cache_Collection,\
    Measurements_Collection,\
    Simulations_Collection

from hera.utils import ConfigurationToJSON

def getProjectList(connectionName=None):
    """
        Return the list with the names of the existing projects .

    :param connectionName: str
        The name of the database.

    Returns
    -------
        list.
    """
    return list(set(AbstractCollection(connectionName=connectionName).getProjectList()))

def createProjectDirectory(outputPath,projectName=None):
    """
        Creates a basic caseConfiguration file
        with the requested project name.

    Parameters
    ----------
    outputPath : str
        The path to create the configuration file in.
        Create if does not exist.

    projectName : str
        The nme of the project.

    Returns
    -------
        None.
    """
    os.makedirs(os.path.abspath(outputPath),exist_ok=True)
    basicOut = dict(projectName=projectName)
    with open(os.path.join(os.path.abspath(outputPath),"caseConfiguration.json"),'w') as outFile:
        json.dump(basicOut,outFile,indent=4)


class Project:
    """
        Provides a simple interface to the data of a specific project.

        The class has all the following functions for the measurements.old, simulations.old and cache.

        **Measurements**


        -  getMeasurementsDocumentsAsDict"
        -  getMeasurementsDocuments
        -  addMeasurementsDocument
        -  deleteMeasurementsDocuments

        **Simulations**

        -  getSimulationsDocumentsAsDict"
        -  getSimulationsDocuments
        -  addSimulationsDocument
        -  deleteSimulationsDocuments

        **Cache**

        -  getCacheDocumentsAsDict"
        -  getCacheDocuments
        -  addCacheDocument
        -  deleteCacheDocuments

    """
    _projectName = None
    _all          = None
    _measurements = None
    _cache     = None
    _simulations  = None

    datatypes = datatypes

    DEFAULTPROJECT = "defaultProject"

    _allowWritingToDefaultProject = None # A flag to allow the update of the default project. Used to add the datasources to it.

    _FilesDirectory = None

    @property
    def FilesDirectory(self):
        """
            The directory to save files (when creating files).
        :return:
        """
        return self._FilesDirectory

    @property
    def filesDirectory(self):
        """
            The directory to save files (when creating files).
        :return:
        """
        return self._FilesDirectory



    @property
    def measurements(self) -> Measurements_Collection:
        """
            Access the measurement type documents.

    Returns
    -------
            hera.datalayer.collection.Measurements_Collection
        """
        return self._measurements

    @property
    def cache(self) -> Cache_Collection:
        """
            Access the Cache type documents.

    Returns
    -------
            hera.datalayer.collection.Cache_Collection

        """
        return self._cache

    @property
    def all(self) -> AbstractCollection:
        """
            Access to all document types.

    Returns
    -------
            hera.datalayer.collection.AbstractCollection

        """
        return self._all

    @property
    def projectName(self):
        return self._projectName


    def setConfig(self,keep_old_values=True, **kwargs):
        """
            Create a config document or updates an existing config document.
        """
        if self._projectName == self.DEFAULTPROJECT:
            raise ValueError("Default project cannot use configuration")

        doc = self._getConfigDocument()
        if keep_old_values:
            doc.desc.update(kwargs)
        else:
            doc['desc'] = kwargs
        doc.save()

    def __init__(self, projectName=None, connectionName=None, configurationPath=None,filesDirectory=None):
        """
            Initialize the project class.

        Parameters
        ----------
        projectName: str default None
                The name of the project.
                If projectName is None, try to load it from the [configurationPath]/caseConfiguration.json
                if configurationPath is None, load it from the current directory.
                the structure of this file is

                ```
                    {
                        "projectName" : [project Name]
                    }
                ```

        connectionName: str
                The name of the DB connection. If not specified, use the
                connection with the linux user name.

        :param configurationPath: str
                The path to the caseConfiguration.json file.
                If not supplied, use the current directory.
        """
        logger = get_classMethod_logger(self,"init")
        self._allowWritingToDefaultProject = False

        if projectName is None:
            configurationPath = os.getcwd() if configurationPath is None else configurationPath
            confFile = os.path.join(configurationPath, "caseConfiguration.json")
            logger.debug(f"projectName is None, try to load file {confFile}")
            if os.path.exists(confFile):
                logger.debug(f"Load as JSON")
                configuration = loadJSON(confFile)
                if 'projectName' not in configuration:
                    err = f"Got projectName=None and the key 'projectName' does not exist in the JSON. "
                    err += """conifguration should be :
{
    'projectName' : [project name]
}                                        
"""
                    logger.error(err)
                    raise ValueError(err)
                else:
                    projectName   = configuration['projectName']
            else:
                logger.debug(f"configuration file {confFile} is not found. Using the default project: DEFAULTPROJECT={self.DEFAULTPROJECT}")
                projectName = self.DEFAULTPROJECT

        logger.info(f"Initializing with logger {projectName}")
        self._projectName = projectName

        self._measurements  = Measurements_Collection(connectionName=connectionName)
        self._cache      = Cache_Collection(connectionName=connectionName)
        self._simulations   = Simulations_Collection(connectionName=connectionName)
        self._all           =   AbstractCollection(connectionName=connectionName)

        savedFilesDirectory = self.getConfig().get("filesDirectory", None)

        if savedFilesDirectory is None:
            if filesDirectory is None:
                filesDirectory = os.path.abspath(os.getcwd())
            else:
                filesDirectory= os.path.abspath(filesDirectory)

            logger.info(f"Files directory is not saved for the project, using {filesDirectory}")
            self.setConfig(filesDirectory=filesDirectory)
        else:
            filesDirectory = savedFilesDirectory

        os.makedirs(os.path.abspath(filesDirectory),exist_ok=True)
        self._FilesDirectory = filesDirectory

    @property
    def simulations(self) -> Simulations_Collection:
        """
            Access the simulation type documents.

    Returns
    -------
            hera.datalayer.collection.Simulation_Collection

        """
        return self._simulations

    def _getConfigDocument(self):
        """
        Returns the document of the config.
        If there is no config document, return empty dictionary.

        Returns
        -------

         dict
                The configuration of the toolkit.
        """
        config_type = f"{self.projectName}__config__"
        documents = self.getCacheDocuments(type=config_type)
        if len(documents) == 0:
            documents = self.addCacheDocument(type=config_type,
                                  resource="",
                                  dataFormat=datatypes.STRING,
                                  desc={})
            ret = documents
        else:
            ret =documents[0]

        return ret

    def setCounter(self,counterName,defaultValue=0):
        """
            Defines a counter in the config of the project.
            The counter is specific to this project.
        Parameters
        ----------
        counterName :  str
            The counter name.

        Returns
        -------

        """
        cnfg = self.getConfig().copy()
        coutnerDict = cnfg.setdefault("counters",{}).copy()
        coutnerDict[counterName] =defaultValue
        cnfg["counters"] = coutnerDict
        self.setConfig(**cnfg)
        return cnfg

    def defineCounter(self,counterName,defaultValue=0):
        """
            Defines a counter in the config of the project, if it does not exist
            The counter is specific to this project.
        Parameters
        ----------
        counterName :  str
            The counter name.

        Returns
        -------

        """
        cnfg = self.getConfig().copy()
        coutnerDict = cnfg.get("counters", {}).copy()
        coutnerDict.setdefault(counterName,defaultValue)
        cnfg["counters"] = coutnerDict
        self.setConfig(**cnfg)
        return cnfg


    def getCounter(self,counterName):
        """
            Return the value of the counter and add [addition].
        Parameters
        ----------
        counterName :  str
            The name of the counter.

        addition : int
            The amount to add to the counter. The default is 1

        Returns
        -------

        """
        cnfg = self.getConfig().copy()
        coutnerDict = cnfg.setdefault("counters", {})
        ret = coutnerDict[counterName]
        return ret

    def getCounterAndAdd(self, counterName, addition=1):
            """
                Return the value of the counter and add [addition].
                If the counter is not defined it is initialized to 0.
            Parameters
            ----------
            counterName :  str
                The name of the counter.

            addition : int
                The amount to add to the counter. The default is 1

            Returns
            -------

            """
            cnfg =self.getConfig()
            counterDict = cnfg.get("counters",{}).copy()
            ret = counterDict.setdefault(counterName,0)
            counterDict[counterName] += addition
            cnfg["counters"] =counterDict
            self.setConfig(**cnfg)
            return ret

    @deprecated(reason="Use getCounterAndAdd instead")
    def addCounter(self, counterName, addition=1):
        return self.getCounterAndAdd(counterName,addition)

    def getConfig(self):
        """
        Returns the config document's description.
        If there is no config document, return empty dictionary.

        Returns
        -------
        dict
                The configuration of the toolkit.
        """
        if self._projectName == self.DEFAULTPROJECT:
            raise ValueError("Default project cannot use configuration")
        doc = self._getConfigDocument()
        return dict(doc["desc"])


    def initConfig(self,**kwargs):
        """
            Sets the value of the config, if the keys does not exist. If they exist, leave the old value.
        Parameters
        ----------
        kwargs

        Returns
        -------

        """
        if self._projectName == self.DEFAULTPROJECT:
            raise ValueError("Default project cannot use configuration")

        doc = self._getConfigDocument()
        for key,value in doc['desc'].items():
            doc['desc'].setdefault(key,value)
        doc.save()



    def getMetadata(self):
        """
        Return the description of all ot the documents in the current project.
        It assumes that description is a key-value format (does not support more complex structures).

    Returns
    -------
             pandas.DataFrame
        """
        descList = [doc.desc for doc in AbstractCollection().getDocuments(projectName=self._projectName)]
        return pandas.DataFrame(descList)

    def getDocumentByID(self,id):
        return self._all.getDocumentByID(id)

    def getMeasurementsDocumentsAsDict(self, with_id=False, **kwargs):
        """
            Querying the DB for measurements.old documents and return the results as a list of dict

        Parameters
        ----------


        with_id : bool, optional, default False
            rather or not should the 'id' key be in the documents.

        kwargs: parameters
            Filters for the query

        Returns
        -------
            List of dicts
        """
        return self.measurements.getDocumentsAsDict(projectName=self._projectName, with_id=with_id, **kwargs)

    def getMeasurementsDocuments(self,  resource=None, dataFormat=None, type=None, **desc):
        """
            Query measurements.old documents.

        Parameters
        ----------
        resource: str
            query by resource, optional.

        dataFormat: str
            query by data format, optional.

        type: str
            query by type

        desc: dict
            query by the measurement document

        Returns
        -------
            List of documents.
        """
        return self.measurements.getDocuments(projectName=self._projectName, resource=resource, dataFormat=dataFormat, type=type, **desc)

    def addDocumentFromDict(self,documentDict):
        """
            Load the document to the project.
            The structure of the dict is:
            {
                _cls : "Metadata.<Cache|Measurements|Simulations>
                desc : {...},
                type : "...",
                resource : ""
                dataFormat : ""
            }

        Parameters
        ----------
        documentDict :  dict
            The dictionary with the data of the document

        Returns
        -------

        """
        addingDict = dict(documentDict)
        if 'projectName' in addingDict:
            del addingDict['projectName']
        docType = addingDict['_cls'].split(".")[1]
        del addingDict['_cls']

        addingFunc = getattr(self,f"add{docType}Document")
        if addingFunc is None:
            raise ValueError(f"The document type {docType} is not found")

        addingFunc(**addingDict)

    def addMeasurementsDocument(self, resource="", dataFormat="string", type="", desc={}):
        """
            Adds a new measurment document.

        Parameters
        ----------
        resource: str
            query by resource, optional.

        dataFormat: str
            query by data format, optional.

        type: str
            query by type

        desc: dict
            query by the measurement document

        Returns
        -------
            The new document
        """
        logger = get_classMethod_logger(self, "init")
        if self.projectName == self.DEFAULTPROJECT and not self._allowWritingToDefaultProject:
            err = f"project {self.projectName} is read-only. "
            logger.error(err)
            raise RuntimeError(err)

        return self.measurements.addDocument(projectName=self._projectName, resource=resource, dataFormat=dataFormat, type=type, desc=desc)

    def deleteMeasurementsDocuments(self, **kwargs):
        """
            Delete the measurements.old documents that fit the query.

        Parameters
        ----------
        kwargs: query dicts.

        Returns
        -------
            The list of documents that was deleted.
        """
        return self.measurements.deleteDocuments(projectName=self._projectName, **kwargs)

    def getSimulationsDocumentsAsDict(self, with_id=False, **kwargs):
        """
            Querying the DB for simulation documents and return the results as a list of dict

        Parameters
        ----------
        with_id : bool, optional, default False
            rather or not should the 'id' key be in the documents.

        kwargs: parameters
            Filters for the query

        Returns
        -------
            List of dicts
        """

        return self.simulations.getDocumentsAsDict(projectName=self._projectName, with_id=with_id, **kwargs)

    def getSimulationsDocuments(self, resource=None, dataFormat=None, type=None, **desc):
        """
            Query simulation documents.

        Parameters
        ----------
        resource: str
            query by resource, optional.

        dataFormat: str
            query by data format, optional.

        type: str
            query by type

        desc: dict
            query by the measurement document

        Returns
        -------
            List of documents.
        """
        return self.simulations.getDocuments(projectName=self._projectName, resource=resource, dataFormat=dataFormat, type=type,
                                **desc)

    def addSimulationsDocument(self, resource="", dataFormat="string", type="", desc={}):
        """
            Adds a new simulations.old document.

        Parameters
        ----------
        resource: str
            query by resource, optional.

        dataFormat: str
            query by data format, optional.

        type: str
            query by type

        desc: dict
            query by the measurement document

        Returns
        -------
            The new document
        """
        logger = get_classMethod_logger(self, "init")
        if self.projectName == self.DEFAULTPROJECT and not self._allowWritingToDefaultProject:
            err = f"project {self.projectName} is read-only. "
            logger.error(err)
            raise RuntimeError(err)

        return self.simulations.addDocument(projectName=self._projectName, resource=resource, dataFormat=dataFormat, type=type,
                               desc=desc)

    def deleteSimulationsDocuments(self, **kwargs):
        """
            Delete the simulations.old documents that fit the query.

        Parameters
        ----------
        kwargs: query dicts.

        Returns
        -------
            The list of documents that was deleted.
        """

        return self.simulations.deleteDocuments(projectName=self._projectName, **kwargs)

    def getCacheDocumentsAsDict(self,  with_id=False, **kwargs):
        """
            Querying the DB for cache documents and return the results as a list of dict

        Parameters
        ----------


        with_id : bool, optional, default False
            rather or not should the 'id' key be in the documents.

        kwargs: parameters
            Filters for the query

        Returns
        -------
            List of dicts
        """

        return self.cache.getDocumentsAsDict(projectName=self._projectName, with_id=with_id, **kwargs)

    def getCacheDocuments(self, resource=None, dataFormat=None, type=None, **desc):
        """
            Querying the DB for cache documents and return the results as a list of dict

        Parameters
        ----------
        with_id : bool, optional, default False
            rather or not should the 'id' key be in the documents.

        kwargs: parameters
            Filters for the query

        Returns
        -------
            List of dicts
        """

        return self.cache.getDocuments(projectName=self._projectName, resource=resource, dataFormat=dataFormat, type=type,
                                       **desc)

    def addCacheDocument(self, resource="", dataFormat="string", type="", desc={}):
        """
            Adds a new cache document.

        Parameters
        ----------
        resource: str
            query by resource, optional.

        dataFormat: str
            query by data format, optional.

        type: str
            query by type

        desc: dict
            query by the measurement document

        Returns
        -------
            The new document
        """
        logger = get_classMethod_logger(self, "init")
        if self.projectName == self.DEFAULTPROJECT and not self._allowWritingToDefaultProject:
            err = f"project {self.projectName} is read-only. "
            logger.error(err)
            raise RuntimeError(err)

        return self.cache.addDocument(projectName=self._projectName, resource=resource, dataFormat=dataFormat, type=type,desc=desc)

    def deleteCacheDocuments(self, **kwargs):
        """
            Delete the cache documents that fit the query.

        Parameters
        ----------
        kwargs: query dicts.

        Returns
        -------
            The list of documents that was deleted.
        """

        return self.cache.deleteDocuments(projectName=self._projectName, **kwargs)

    def saveData(self,name,data,desc,kind,type=None,dataFormat=None,**kwargs):
        """
            Adds a cache document with the data.
            Estimates the dataFormat from the data type.

            The type is
        Parameters
        ----------
        data : a data that can be dataframe, xarray, numpy ....

        desc : dict
            A dict with the meatadata to save
        type : str
            If None, then set the type to be the name of the function that called this method.

        Returns
        -------
            The new document
        """
        guessedDataFormat = self.datatypes.getDataFormatName(data) if dataFormat is None else dataFormat

        handler = self.datatypes.getHandler(guessedDataFormat)
        file_extension = self.datatypes.getDataFormatExtension(data)

        cacheDirectory = os.path.join(self.filesDirectory, "cache")
        os.makedirs(cacheDirectory, exist_ok=True)
        fileID = self.getCounterAndAdd(name)
        fileName = os.path.join(cacheDirectory, f"{name}_{fileID}.{file_extension}")

        saveParamsUpdatedDict = handler.saveData(data, fileName,**kwargs)

        funcName = getattr(self,f"add{kind}Document")
        fullType = type if type is not None else name

        qry = ConfigurationToJSON(desc, standardize=True, splitUnits=True, keepOriginalUnits=True)
        storeParamsDict = qry.get("storeParameters",{})
        storeParamsDict.update(saveParamsUpdatedDict)
        qry["storeParameters"] = storeParamsDict

        doc = funcName(type=fullType, dataFormat=guessedDataFormat, resource=fileName, desc=qry)
        return doc

    def saveMeasurementData(self,name,data,desc,type=None,dataFormat=None,**kwargs):
        self.saveData(name=name,data=data,desc=desc,kind="Measurement",type=type,dataFormat=dataFormat,**kwargs)

    def saveCacheData(self,name,data,desc,type=None,dataFormat=None,**kwargs):
        self.saveData(name=name,data=data,desc=desc,kind="Cache",type=type,dataFormat=dataFormat,**kwargs)

    def saveSimulationData(self,name,data,desc,type=None,dataFormat=None,**kwargs):
        self.saveData(name=name,data=data,desc=desc,kind="Simulation",type=type,dataFormat=dataFormat,**kwargs)

    def _get_full_func_name(self,func):
        """Returns the full qualified path: module.[class.]function_name"""
        if not callable(func):
            raise TypeError("Provided object is not callable.")

        # Handle bound methods by unwrapping them
        if inspect.ismethod(func):
            # Get the original function and its class
            cls = func.__self__.__class__
            method_name = func.__name__
            class_qualname = cls.__qualname__
            module = func.__module__
            if module == "__main__":
                ret = f"{class_qualname}.{method_name}"
            else:
                ret = f"{module}.{class_qualname}.{method_name}"

        elif inspect.isfunction(func):
            ret = func.__qualname__
        else:

            # Handle unbound class or static methods and plain functions
            qualname = func.__qualname__
            module = func.__module__
            ret = f"{module}.{qualname}"
        return ret



    @staticmethod
    def getProjectList(cls,user=None):
        """
            Return the list with the names of the existing projects .

        :param user: str
            The name of the database.

        Returns
        -------
            list.
        """
        return list(set(AbstractCollection(connectionName=user).getProjectList()))

