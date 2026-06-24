import itertools
import json
import os
import pickle
import zipfile
import pandas
import inspect

from tqdm import tqdm
from deprecated import deprecated

from hera.datalayer.datahandler import datatypes
from hera.utils.logging import get_classMethod_logger
from hera.datalayer.collection import AbstractCollection,\
    Cache_Collection,\
    Measurements_Collection,\
    Simulations_Collection

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
    from hera.utils import get_logger
    logger = get_logger(None,"hera.datalayer.project.createProjectDirectory")
    if projectName == "defaultProject":
        logger.warning("Cannot create default project")
        return
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
        """
        The name of the current project.

        Returns
        -------
        str
        """
        return self._projectName


    def setConfig(self,keep_old_values=True, **kwargs):
        """
            Create a config document or updates an existing config document.
        """
        if self._projectName == self.DEFAULTPROJECT:
            err = """
                    Default project was initiated and it is cannot use configuration.
                    
                    This error os obtained because you have initiated the project or the toolkit 
                    using projectName=None, and you don't case a caseConfiguration.json in the directory that 
                    specifies the project name. 
                    
                    To solve this project either 
                    - Create a new project in the directory: 
                      run  
                        >> hera-project project create <the directory name>
                      in the parent directory 
                    - Create manually the file caseConfiguration.json in the directory with the key 'projectName' that holds 
                      the name of the project as a string.     
            """
            raise ValueError(err)

        doc = self._getConfigDocument()
        if keep_old_values:
            update_kwargs = {
                f"set__desc__{k}": v
                for k, v in kwargs.items()
            }
            doc.update(**update_kwargs)
        else:
            doc.update(set__desc=kwargs)
        

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
                from hera.utils.jsonutils import loadJSON
                configuration = loadJSON(confFile)
                if 'projectName' not in configuration:
                    err = f"Got projectName=None and the key 'projectName' does not exist in the JSON. "
                    err += """configuration should be :
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

        if self.projectName != self.DEFAULTPROJECT:
            logger.info(f"Attempting to get default directory from the disk")
            savedFilesDirectory = self.getConfig().get("filesDirectory", None)
        else:
            logger.info(f"Default project, setting to current directory(not saving on disk)")
            savedFilesDirectory = None

        if filesDirectory is not None:
            # Caller explicitly passed a directory — always honor it and persist it.
            filesDirectory = os.path.abspath(filesDirectory)
            if filesDirectory != savedFilesDirectory and self.projectName != self.DEFAULTPROJECT:
                self.setConfig(filesDirectory=filesDirectory)
        elif savedFilesDirectory is not None:
            filesDirectory = savedFilesDirectory
        elif self.projectName != self.DEFAULTPROJECT:
            # No saved path and none passed — default to ~/.hera/<projectName>
            filesDirectory = os.path.join(os.path.expanduser("~"), ".hera", projectName)
            logger.info(f"Files directory is not saved for the project, using {filesDirectory}")
            self.setConfig(filesDirectory=filesDirectory)

        if self.projectName != self.DEFAULTPROJECT:
            os.makedirs(os.path.abspath(filesDirectory),exist_ok=True)
        self._FilesDirectory = filesDirectory

    @staticmethod
    def _batched_cursor(cur, export_chunk_size):# -> Generator[list, Any, None]:
        """
            turns an iterator to an iterator that returns batches of the original.
            if you cur is a documents list iterator and export_chunk_size=5 then each iteration of  
            _batched_cursor will return a list with 5 documents

        Parameters
        ----------
        cur: Iterator

        export_chunk_size: int

        """
        while True:
            batch = list(itertools.islice(cur, export_chunk_size))
            if not batch:
                break
            else:
                yield batch
    
    def export(self, path, export_chunk_size=1024, show_progressbar=True):
        """
            exports the project in chunks contained to one zip file

        Parameters
        ----------
        path: str
                the path to the zip file to be created, overrides file if already exists

        export_chunk_size: int
                the number of documents in each chunk

        show_progressbar: bool
        """
        logger = get_classMethod_logger(self,"export")
        if self.projectName == self.DEFAULTPROJECT:
            logger.warning("Can't export default project")
            return
        docs_cursor = self._all._metadataCol._get_collection().find({"projectName":self.projectName})
        with zipfile.ZipFile(path, mode='w', compression=zipfile.ZIP_DEFLATED) as zf:
            i = 0
            if show_progressbar:
                docs_iterator= tqdm(self._batched_cursor(docs_cursor, export_chunk_size),
                                    desc="Exporting documents", unit="batchedDocs", unit_scale=True) 
            else:
                docs_iterator = docs_cursor
            for docs_batch in docs_iterator:
                filename = f"chunk_{i}"
                with zf.open(filename, 'w') as zf_archive:
                    pickle.dump(docs_batch, zf_archive, protocol=pickle.HIGHEST_PROTOCOL)  # nosec B301 — internal format, not from untrusted input; migration to JSON planned
                i+=1

    @staticmethod
    def _iter_pickled_docs(zf, return_batched):
        """
            creates an iterator over the documents in an exported project,
            can either returned the batched data or single documents 

        Parameters
        ----------
        zf: zipfile.ZipFile
                handler for the exported project

        return_batched: str
                should each iteration return an entire batch(with the size it was exported as)
        """
        for name in zf.namelist():
            with zf.open(name) as f:
                depickled_docs_batch=pickle.load(f)  # nosec B301 — loads from project's own export; migration to JSON planned in Part 2
                if return_batched:
                    yield depickled_docs_batch
                else:
                    for depickled_doc in depickled_docs_batch:
                        yield depickled_doc 

    def updateProjectNameOnDoc(self, doc_son):
        """
            updates the projectName field of a document to be assigned to this project 

        Parameters
        ----------
        doc_son: str
                the document to be updated
        """
        doc_son.update({"projectName": self.projectName})
        return doc_son
    
    @staticmethod
    def load(proj, exported_project_path, is_hard_import, show_progressbar=True):
        """
            loads an exported project's documents to a project, either hard import with the original ids or not 

        Parameters
        ----------
        proj: str or Project
                either name of the project to import or an existing Project instance,
                creates a new project if it doesn't exist

        is_hard_import: bool
                should the original ids be used when importing or generate new ones per document(using original ids might fail due to duplicate ids)

        show_progressbar: bool 
        """
        if isinstance(proj, str):
            proj = Project(proj)
        
        with zipfile.ZipFile(exported_project_path, 'r') as zf:
            if is_hard_import:
                pickled_docs_iterator = Project._iter_pickled_docs(zf, return_batched=True)
                if show_progressbar:
                    pickled_docs_iterator = tqdm(pickled_docs_iterator, desc="Loading documents", unit="docsBatch", unit_scale=True)
                for pickled_docs_batch in pickled_docs_iterator:
                    proj._all._metadataCol._get_collection().insert_many([proj.updateProjectNameOnDoc(doc_son) for doc_son in pickled_docs_batch])
            else:
                pickled_docs_iterator = Project._iter_pickled_docs(zf, return_batched=False)
                if show_progressbar:
                    pickled_docs_iterator = tqdm(pickled_docs_iterator, desc="Loading documents", unit="docs", unit_scale=True)
                for pickled_doc in pickled_docs_iterator:
                    proj._all._metadataCol.objects.insert(proj._all._metadataCol(**(proj.updateProjectNameOnDoc(pickled_doc))), load_bulk=False)

    @property
    def simulations(self) -> Simulations_Collection:
        """
            Access the simulation type documents.

    Returns
    -------
            hera.datalayer.collection.Simulation_Collection

        """
        return self._simulations

    def _getConfigDocument(self, evaluate=False):
        """
        Returns the document of the config.
        If there is no config document, return empty dictionary.

        Parameters
        -------
        evaluate :  bool
            if False returns a QuerySet otherwise evaluates to a document

        Returns
        -------

         dict
                The configuration of the toolkit.
        """
        config_type = f"{self.projectName}__config__"
        documents = self.getCacheDocuments(type=config_type)
        # keeping lazy so we just return documents which is a QuerySet
        if documents.first() is None:
            self.addCacheDocument(type=config_type,
                                  resource="",
                                  dataFormat=datatypes.STRING,
                                  desc={})

        return documents[0] if evaluate else documents

    def setCounter(self,counterName:str,defaultValue=0):
        """
            Defines a counter in the config of the project.
            The counter is specific to this project.
        Parameters
        ----------
        counterName :  str
            The counter name.

        Returns
        -------
        None
        """
        counterName = self._normalizeCounterName(counterName)
        cnfg_doc = self._getConfigDocument()
        self._enforce_counter_field(cnfg_doc)
        # defineCounter only creates if not exists; then we force-set the value.
        self.defineCounter(counterName, 0)
        cnfg_doc.update_one(**{f"set__desc__counters__{counterName}": defaultValue})

    def _normalizeCounterName(self, counterName):
        """Normalize a counter name by replacing dots with underscores."""
        counterName=counterName.replace('.','_')
        # avoiding conflicts with mongodb
        if '__' in counterName:
            raise RuntimeError("a counter's name cannot contain either of __ / ._ / _. / ..")
        return counterName

    def defineCounter(self,counterName,defaultValue=0) -> None:
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
        counterName = self._normalizeCounterName(counterName)
        cnfg_doc = self._getConfigDocument()
        self._enforce_counter_field(cnfg_doc)
        if cnfg_doc.filter(**{f"desc__counters__{counterName}__exists": True}).first() is None: 
            cnfg_doc.update(**{f"set__desc__counters__{counterName}": defaultValue})
            return True
        return False


    def getCounter(self,counterName):
        """
            Return the value of the counter if doesn't exist returning None
        Parameters
        ----------
        counterName :  str
            The name of the counter.

        Returns
        -------

        """
        counterName = self._normalizeCounterName(counterName)
        cnfg_doc = self._getConfigDocument()
        self._enforce_counter_field(cnfg_doc)
        #if it doesn't exist there is nothing to get, returning None
        if cnfg_doc.filter(**{f"desc__counters__{counterName}__exists": True}).first() is None: 
            return None
        doc = cnfg_doc[0]
        return doc['desc']['counters'][counterName]

    def _enforce_counter_field(self, cnfg_doc):
        """Ensure the ``counters`` field exists in the config document."""
        if cnfg_doc.filter(desc__counters__exists=True).first() is None:
            cnfg_doc.update(**{f"set__desc__counters": {}})


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
        int
            The updated counter value, or 0 if the counter was newly created.
        """
        counterName = self._normalizeCounterName(counterName)
        isNew = self.defineCounter(counterName,0)
        if isNew:
            return 0

        cnfg_doc = self._getConfigDocument()
        doc = cnfg_doc.modify(
            upsert=True,
            new=True,
            **{
                f"inc__desc__counters__{counterName}":addition
            })
        return doc['desc']['counters'][counterName]

    @deprecated(reason="Use getCounterAndAdd instead")
    def addCounter(self, counterName, addition=1):
        """Deprecated. Use ``getCounterAndAdd`` instead."""
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
        doc = self._getConfigDocument(evaluate=True)
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
        for key in kwargs.keys():
            doc.filter(**{f"desc__{key}__exists":False}).update(**{f"set__desc__{key}":kwargs[key]})



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

    def getDocumentByID(self, id):
        """
        Returns a document by its ID from any collection.

        Parameters
        ----------
        id : str
            The document ID.

        Returns
        -------
        document
            The document with the given ID.
        """
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

    def getMeasurementsDocuments(self, resource=None, dataFormat=None, type=None, **desc):
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

    def getAllDocuments(self, resource=None, dataFormat=None, type=None, **desc):
        """
            Runs getXDocuments for measurements, cache and simulations aggregating all to one list

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
        docs = []
        docs.extend(self.getSimulationsDocuments(resource=resource, dataFormat=dataFormat, type=type, **desc))
        docs.extend(self.getMeasurementsDocuments(resource=resource, dataFormat=dataFormat, type=type, **desc))
        docs.extend(self.getCacheDocuments(resource=resource, dataFormat=dataFormat, type=type, **desc))
        return docs
    
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

    def addMeasurementsDocument(self, resource="", dataFormat="string", type="", desc=None):
        """
            Adds a new measurement document.

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
        if desc is None:
            desc = {}
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

    def addSimulationsDocument(self, resource="", dataFormat="string", type="", desc=None):
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
        if desc is None:
            desc = {}
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

    def addCacheDocument(self, resource="", dataFormat="string", type="", desc=None):
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
        if desc is None:
            desc = {}
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
            A dict with the metadata to save
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

        from hera.utils.jsonutils import ConfigurationToJSON
        qry = ConfigurationToJSON(desc, standardize=True, splitUnits=True, keepOriginalUnits=True)
        storeParamsDict = qry.get("storeParameters",{})
        storeParamsDict.update(saveParamsUpdatedDict)
        qry["storeParameters"] = storeParamsDict

        doc = funcName(type=fullType, dataFormat=guessedDataFormat, resource=fileName, desc=qry)
        return doc

    def saveMeasurementData(self, name, data, desc, type=None, dataFormat=None, **kwargs):
        """
        Save data as a measurement document. See ``saveData`` for parameter details.
        """
        return self.saveData(name=name, data=data, desc=desc, kind="Measurements", type=type, dataFormat=dataFormat, **kwargs)

    def saveCacheData(self, name, data, desc, type=None, dataFormat=None, **kwargs):
        """
        Save data as a cache document. See ``saveData`` for parameter details.
        """
        return self.saveData(name=name, data=data, desc=desc, kind="Cache", type=type, dataFormat=dataFormat, **kwargs)

    def saveSimulationData(self, name, data, desc, type=None, dataFormat=None, **kwargs):
        """
        Save data as a simulation document. See ``saveData`` for parameter details.
        """
        return self.saveData(name=name, data=data, desc=desc, kind="Simulations", type=type, dataFormat=dataFormat, **kwargs)

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

