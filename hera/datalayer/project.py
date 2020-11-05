import pandas
import getpass
import logging
import numpy
from . import getDBNamesFromJSON

from .collection import AbstractCollection,\
    Cache_Collection,\
    Measurements_Collection,\
    Simulations_Collection


def getProjectList(user=None):
    """
        Return the list with the names of the existing projects .

    :param user: str
        The name of the database.

    :return:
        list.
    """
    return list(set(AbstractCollection(user=user).getProjectList()))



class Project(object):
    """
        Provides a simple interface to the data of a specific project.

        The class has all the following functions for the measurements, simulations and Cache:

        -  getMeasurementsDocumentsAsDict"
        -  getMeasurementsDocuments
        -  addMeasurementsDocument
        -  deleteMeasurementsDocuments


        -  getSimulationsDocumentsAsDict"
        -  getSimulationsDocuments
        -  addSimulationsDocument
        -  deleteSimulationsDocuments

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

    _logger     = None

    @property
    def projectName(self):
        return self._projectName

    def getProjectName(self):
        return self._projectName

    @property
    def logger(self):
        return self._logger

    @property
    def measurements(self):
        """
            Access the measurement type documents.

        :return:
            hera.datalayer.collection.Measurements_Collection
        """
        return self._measurements

    @property
    def cache(self):
        """
            Access the Cache type documents.

        :return:
            hera.datalayer.collection.Cache_Collection

        """
        return self._cache

    @property
    def all(self):
        """
            Access the Cache type documents.

        :return:
            hera.datalayer.collection.Cache_Collection

        """
        return self._all

    @property
    def simulations(self):
        """
            Access the simulation type documents.

        :return:
            hera.datalayer.collection.Simulation_Collection

        """
        return self._simulations

    def __init__(self, projectName, databaseName=None):
        """
            Initialize the project.

        :param projectName: str
                The name of the project.

        :param databaseName: str
                the name of the database to use. If None, use the default database (the name of the current databaseName).

        """
        self._projectName = projectName

        self._measurements  = Measurements_Collection(user=databaseName)
        self._cache      = Cache_Collection(user=databaseName)
        self._simulations   = Simulations_Collection(user=databaseName)

        self._setLogger()

    def getConfig(self):
        """
        Returns the config document's description.
        If there is no config document, return None.
        """
        documents = self.getCacheDocumentsAsDict(type="__config__")
        return dict()  if len(documents) == 0 else documents[0].desc


    def setConfig(self, config):
        """
        Create a config documnet or updates an existing config document.
        """
        documents = self.getCacheDocuments(type="__config__")
        if len(documents) == 0:
            self.addCacheDocument(type="__config__",desc=config)
        else:
            documents[0].update(desc=config)

    def _setLogger(self):
        """
            Set the logger with the name of the class.
        :return:
        """
        name =  ".".join(str(self.__class__)[8:-2].split(".")[1:])
        self._logger = logging.getLogger(name)

    def getMetadata(self):
        """
        Returns a pandas dataframe which contains all the description of all ot the documents in the current project.

        :return: pandas
        """
        descList = [doc.desc for doc in AbstractCollection().getDocuments(projectName=self._projectName)]
        return pandas.DataFrame(descList)

    def getMeasurementsDocumentsAsDict(self, with_id=False, **kwargs):
        return self.measurements.getDocumentsAsDict(projectName=self._projectName, with_id=with_id, **kwargs)

    def getMeasurementsDocuments(self,  resource=None, dataFormat=None, type=None, **desc):
        return self.measurements.getDocuments(projectName=self._projectName, resource=resource, dataFormat=dataFormat, type=type, **desc)

    def addMeasurementsDocument(self, resource="", dataFormat="string", type="", desc={}):
        return self.measurements.addDocument(projectName=self._projectName, resource=resource, dataFormat=dataFormat, type=type, desc=desc)

    def deleteMeasurementsDocuments(self, **kwargs):
        return self.measurements.deleteDocuments(projectName=self._projectName, **kwargs)

    def getSimulationsDocumentsAsDict(self, with_id=False, **kwargs):
        return self.simulations.getDocumentsAsDict(projectName=self._projectName, with_id=with_id, **kwargs)

    def getSimulationsDocuments(self, resource=None, dataFormat=None, type=None, **desc):
        return self.simulations.getDocuments(projectName=self._projectName, resource=resource, dataFormat=dataFormat, type=type,
                                **desc)

    def addSimulationsDocument(self, resource="", dataFormat="string", type="", desc={}):
        return self.simulations.addDocument(projectName=self._projectName, resource=resource, dataFormat=dataFormat, type=type,
                               desc=desc)

    def deleteSimulationsDocuments(self, **kwargs):
        return self.simulations.deleteDocuments(projectName=self._projectName, **kwargs)

    def getCacheDocumentsAsDict(self,  with_id=False, **kwargs):
        return self.cache.getDocumentsAsDict(projectName=self._projectName, with_id=with_id, **kwargs)

    def getCacheDocuments(self, resource=None, dataFormat=None, type=None, **desc):
        return self.cache.getDocuments(projectName=self._projectName, resource=resource, dataFormat=dataFormat, type=type,
                                       **desc)

    def addCacheDocument(self, resource="", dataFormat="string", type="", desc={}):
        return self.cache.addDocument(projectName=self._projectName, resource=resource, dataFormat=dataFormat, type=type,
                                      desc=desc)

    def deleteCacheDocuments(self, **kwargs):
        return self.cache.deleteDocuments(projectName=self._projectName, **kwargs)


class ProjectMultiDB:
    """
        Provides a simple interface to the data of a specific project.

        The class has all the following functions for the measurements, simulations and Cache:

        -  getMeasurementsDocumentsAsDict"
        -  getMeasurementsDocuments
        -  addMeasurementsDocument
        -  deleteMeasurementsDocuments


        -  getSimulationsDocumentsAsDict"
        -  getSimulationsDocuments
        -  addSimulationsDocument
        -  deleteSimulationsDocuments

        -  getCacheDocumentsAsDict"
        -  getCacheDocuments
        -  addCacheDocument
        -  deleteCacheDocuments

    """
    _all = None
    _measurements = None
    _cache     = None
    _simulations  = None
    _useAll = None

    _logger     = None

    _projectNameDict = None

    @property
    def logger(self):
        return self._logger

    def _setLogger(self):
        """
            Set the logger with the name of the class.
        :return:
        """
        name =  ".".join(str(self.__class__)[8:-2].split(".")[1:])
        self._logger = logging.getLogger(name)



    @property
    def measurements(self):
        """
            Access the measurement type documents.

        :return:
            hera.datalayer.collection.Measurements_Collection
        """
        return self._measurements

    @property
    def cache(self):
        """
            Access the Cache type documents.

        :return:
            hera.datalayer.collection.Cache_Collection

        """
        return self._cache

    @property
    def simulations(self):
        """
            Access the simulation type documents.

        :return:
            hera.datalayer.collection.Simulation_Collection

        """
        return self._simulations

    @property
    def projectNameDict(self):
        """
            return the DB->project name map
        :return:
        """
        return self._projectNameDict

    @projectNameDict.setter
    def projectNameDict(self, value):
        """
            Creates the DB->project name map.

        :param value: str or dict.

                if str, define a map with the default DB (user name)
                else just use the map passed.
        :return:
        """
        if isinstance(value,str):
            self._projectNameDict ={getpass.getuser() : value}
        elif isinstance(value, dict):
            self._projectNameDict = value
        else:
            raise ValueError("project name dict must be str or dict. ")

    @property
    def databaseNames(self):
        return [x for x in self.projectNameDict.keys()]


    @property
    def useAll(self):
        return self._useAll

    @useAll.setter
    def useAll(self,newUseAll):
        self._useAll=newUseAll


    def getProjectName(self, databaseName=None):
        """
            Return the project name of the relevant database.

        :param databaseName: str
                    The name of the database. Return the name of the default database (the user name)
                    if None.
        :return:
        """
        if isinstance(self._projectNameDict,str):
            projectName = self._projectNameDict
        else:
            databaseName  = getpass.getuser() if databaseName is None else databaseName
            projectName = self._projectNameDict[databaseName]

        return  projectName


    def __init__(self, projectNameDict, useAll=False):
        """
            Initialize the project.

        Parameters
        ----------

        projectNameDict: str, dict .
                The name of the project.

                if str, use the project name only in the default DB.

                if dict, defines the DB to look in and the corresponding project name.

        databaseNameList: str,list
                the name of the database to use.
                If None, use the default database (the name of the current user).


        """

        self._setLogger()
        self.projectNameDict = projectNameDict
        self._useAll = useAll
        self._measurements  = dict([(user,Measurements_Collection(user=user)) for user in self.databaseNames])
        self._cache         = dict([(user,Cache_Collection(user=user)) for user in self.databaseNames])
        self._simulations   = dict([(user,Simulations_Collection(user=user)) for user in self.databaseNames])
        self._all           = dict([(user,AbstractCollection(user=user)) for user in self.databaseNames])



    def getConfig(self):
        """
        Returns the config document's description.
        If there is no config document, return None.
        """
        documents = self.getCacheDocuments(type="__config__")
        return dict() if len(documents) == 0 else documents[0].desc


    def setConfig(self, config):
        """
        Create a config documnet or updates an existing config document.
        """
        documents = self.getCacheDocuments(type="__config__")
        if len(documents) == 0:
            self.addCacheDocument(type="__config__",desc=config)
        else:
            doc = documents[0]
            doc.desc.update(desc=config)
            documents[0].save()

    def getMetadata(self):
        """
        Returns a pandas dataframe which contains all the description of all ot the documents in the current project.

        :return: pandas
        """
        descList = []
        for userName,allDB in self._all.item():
            projectName = self.getProjectName(userName)
            descList = [doc.desc for doc in self._all.getDocuments(projectName=projectName)]

        return pandas.DataFrame(descList)

    def _getSomeTypeDocumentsAsDict(self, searchtype, with_id, users=None, **kwargs):
        returnData = []
        searchtype = searchtype if users is None else dict([(user, searchtype[user]) for user in users])
        for userName, searched in searchtype.items():
            projectName = self.getProjectName(userName)
            data = searched.getDocumentsAsDict(projectName=projectName, with_id=with_id, **kwargs)
            if len(data["documents"]) != 0:
                if self._useAll:
                    returnData.append(data)
                else:
                    returnData = data
                    break

        return returnData

    def _getSomeTypeDocuments(self, searchtype, resource, dataFormat, type, **desc):
        returnData = []
        for userName, searched in searchtype.items():
            projectName = self.getProjectName(userName)
            data = searched.getDocuments(projectName=projectName, resource=resource, dataFormat=dataFormat,type=type, **desc)
            if len(data) != 0:
                if self._useAll:
                    returnData += data
                else:
                    returnData = data
                    break
        return returnData

    def _addSomeTypeDocuments(self, searchtype, resource, dataFormat, type, users=None, **desc):
        if users is None:
            userName = getpass.getuser()
            projectName = self.getProjectName(userName)
            ret = searchtype[userName].addDocument(projectName=projectName, resource=resource, dataFormat=dataFormat, type=type, **desc)
        else:
            for user in numpy.atleast_1d(users):
                projectName = self.getProjectName(user)
                ret = searchtype[user].addDocument(projectName=projectName, resource=resource, dataFormat=dataFormat, type=type, **desc)
        return ret

    def _deleteSomeTypeDocuments(self, searchtype, users=None, **kwargs):
        if users is None:
            userName = getpass.getuser()
            projectName = self.getProjectName(userName)
            searchtype[userName ].deleteDocuments(projectName=projectName, **kwargs)
        else:
            for user in numpy.atleast_1d(users):
                projectName = self.getProjectName(user)
                searchtype[user].deleteDocuments(projectName=projectName, **kwargs)

    def getMeasurementsDocumentsAsDict(self, with_id=False, users=None, **kwargs):
        return self._getSomeTypeDocumentsAsDict(searchtype=self._measurements, with_id=with_id, users=users, **kwargs)

    def getMeasurementsDocuments(self, resource=None, dataFormat=None, type=None, **desc):
        return self._getSomeTypeDocuments(searchtype=self._measurements, resource=resource, dataFormat=dataFormat, type=type, **desc)

    def addMeasurementsDocument(self, resource="", dataFormat="string", type="", desc={}, users=None):
        return self._addSomeTypeDocuments(searchtype=self._measurements, resource=resource, dataFormat=dataFormat, type=type, desc=desc, users=users)

    def deleteMeasurementsDocuments(self, users=None, **kwargs):
        return self._deleteSomeTypeDocuments(searchType=self._measurements, users=users, **kwargs)

    def getSimulationsDocumentsAsDict(self, with_id=False, users=None, **kwargs):
        return self._getSomeTypeDocumentsAsDict(searchtype=self._simulations, with_id=with_id,users=users, **kwargs)

    def getSimulationsDocuments(self, resource=None, dataFormat=None, type=None, **desc):
        return self._getSomeTypeDocuments(searchtype=self._simulations, resource=resource, dataFormat=dataFormat, type=type, **desc)

    def addSimulationsDocument(self, resource="", dataFormat="string", type="", desc={}, users=None):
        return self._addSomeTypeDocuments(searchtype=self._simulations, resource=resource, dataFormat=dataFormat, type=type, desc=desc, users=users)

    def deleteSimulationsDocuments(self, users=None, **kwargs):
        return self._deleteSomeTypeDocuments(searchtype=self._simulations, users=users, **kwargs)

    def getCacheDocumentsAsDict(self,  with_id=False, users=None, **kwargs):
        return self._getSomeTypeDocumentsAsDict(searchtype=self._cache, with_id=with_id, users=users, **kwargs)

    def getCacheDocuments(self, resource=None, dataFormat=None, type=None, **desc):
        return self._getSomeTypeDocuments(searchtype=self._cache, resource=resource, dataFormat=dataFormat, type=type, **desc)

    def addCacheDocument(self, resource="", dataFormat="string", type="", desc={}, users=None):
        return self._addSomeTypeDocuments(searchtype=self._cache, resource=resource, dataFormat=dataFormat, type=type, desc=desc, users=users)

    def deleteCacheDocuments(self, users=None, **kwargs):
        return self._deleteSomeTypeDocuments(searchtype=self._cache, users=users, **kwargs)


class ProjectMultiDBPublic(ProjectMultiDB):
    """
        A multi-project db, but adds the Public (or public) to the search db list.

        The class accepts the default public project name.

    """
    def __init__(self, projectName, publicProjectName, useAll=False):
        """
            Initializing a multi-database project that has in data in a Public database.



        Parameters:
        -----------

         projectName: str, dict
            if Str, look only in the default database (the name of the user)
            using that project name.

            if dict, the map of project name -> DB. will define the DB to look in.

         publicProjectName: str
                The project name in the public DB.

         useAll: bool
                If true, return a union of all the results from all the DB.

        """
        if isinstance(projectName,str):
            projectNamesDict = {getpass.getuser() : projectName}
        elif isinstance(projectName,dict):
            projectNamesDict = projectName
        else:
            raise ValueError(f"projectName must be str or dict, not {type(projectName)}")


        if 'public' in getDBNamesFromJSON():
            projectNamesDict['public'] = publicProjectName
        if ('Public' in getDBNamesFromJSON()):
            projectNamesDict['Public'] = publicProjectName

        super().__init__(projectNamesDict, useAll=useAll)