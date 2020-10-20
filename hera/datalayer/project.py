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

        self._measurements  = Measurements_Collection(user=user)
        self._cache      = Cache_Collection(user=user)
        self._simulations   = Simulations_Collection(user=user)
        self._all           =   AbstractCollection(user=user)

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
    _projectName = None


    _all = None
    _measurements = None
    _cache     = None
    _simulations  = None
    _databaseNameList = None
    _useAll = None

    _logger     = None

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
    def databaseName(self):
        return self._databaseNameList

    @databaseName.setter
    def databaseName(self, newDatabaseList):
        self._databaseNameList = newDatabaseList
        self._measurements  = [Measurements_Collection(user=user) for user in newDatabaseList]
        self._cache         = [Cache_Collection(user=user) for user in newDatabaseList]
        self._simulations   = [Simulations_Collection(user=user) for user in newDatabaseList]
        self._all           = [AbstractCollection(user=user) for user in newDatabaseList]

    @property
    def useAll(self):
        return self._useAll

    @useAll.setter
    def useAll(self,newUseAll):
        self._useAll=newUseAll


    def getProjectName(self, databaseName=None):

        if databaseName is None:
            projectName = self._projectName
        if isinstance(self._projectName,str):
            projectName = self._projectName
        else:
            projectName = self._projectName[databaseName]

        return  projectName


    def __init__(self, projectName, databaseNameList=None, useAll=False):
        """
            Initialize the project.

        Parameters
        ----------

        projectName: str, dict .
                The name of the project.
                if dict, the project name depends on the database.

        databaseNameList: str,list
                the name of the database to use.
                If None, use the default database (the name of the current user).


        """
        self._projectName = projectName
        self._databaseNameList = numpy.atleast_1d(databaseNameList)
        self._useAll = useAll
        self._measurements  = dict([(user,Measurements_Collection(user=user)) for user in self._databaseNameList])
        self._cache         = dict([(user,Cache_Collection(user=user)) for user in self._databaseNameList])
        self._simulations   = dict([(user,Simulations_Collection(user=user)) for user in self._databaseNameList])
        self._all           = dict([(user,AbstractCollection(user=user)) for user in self._databaseNameList])

        self._setLogger()

    def getConfig(self):
        """
        Returns the config document's description.
        If there is no config document, return None.
        """
        documents = self.getCacheDocumentsAsDict(type="__config__")
        if len(documents) == 0:
            raise KeyError("There is no config document.")
        else:
            if type(documents)==list:
                desc = documents[0]["documents"][0]["desc"]
            else:
                desc = documents["documents"][0]["desc"]
        return desc

    def setConfig(self, config, user=None):
        """
        Create a config documnet or updates an existing config document.
        """
        documents = self.getCacheDocuments(type="__config__",user=user)
        if len(documents) == 0:
            if self._databaseNameList[0] == "public" or self._databaseNameList[0] == "Public":
                if len(self._databaseNameList) == 1:
                    raise KeyError("Can't set config document in public, choose aditional user/s.")
                else:
                    user = self._databaseNameList[1] if user is None else user
            else:
                user = self._databaseNameList[0] if user is None else user
            self.addCacheDocument(type="__config__",desc=config,users=[user])
        else:
            documents[0].update(desc=config)

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
                    returnData.append(data)
                else:
                    returnData = data
                    break
        return returnData

    def _addSomeTypeDocuments(self, searchtype, resource, dataFormat, type, users=None, **desc):
        if users is None:
            if self._databaseNameList[0] == "public" or self._databaseNameList[0] == "Public" and len(self._databaseNameList) > 1:
                userName = self._databaseNameList[1]
            else:
                userName = self._databaseNameList[0]
            projectName = self.getProjectName(userName)
            searchtype[userName].addDocument(projectName=projectName, resource=resource, dataFormat=dataFormat, type=type, **desc)
        else:
            for user in numpy.atleast_1d(users):
                projectName = self.getProjectName(user)
                searchtype[user].addDocument(projectName=projectName, resource=resource, dataFormat=dataFormat, type=type, **desc)

    def _deleteSomeTypeDocuments(self, searchtype, users=None, **kwargs):
        if users is None:
            userName = self._databaseNameList[0]
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
    def __init__(self, projectName, publicProjectName, databaseNameList=None, useAll=False):
        """
            Initializes the search list of the DB.

            The class is initiated with the default project name for the public DB
            and the list of DB's and project names to look for.

            The public is initiated as the first DB to look in.

        Parameters:
        -----------

         projectName: str, dict
            The project name (if str).
            if dict, the map of project name for a DB.

         publicProjectName: str
                The project name in the public DB.
         databaseNameList: str, list of str
                The name of the DB to look in (except for public).
                Can be a str or a list.

         useAll: bool
                If true, return a union of all the results from all the DB.

        """
        projectNamesDict = dict()
        dbListNames = []
        if ('public' in getDBNamesFromJSON()):
            dbListNames = ['public']
            projectNamesDict['public'] = publicProjectName
        if ('Public' in getDBNamesFromJSON()):
            dbListNames = ['Public']
            projectNamesDict['Public'] = publicProjectName

        elif isinstance(projectName,dict):
                projectNamesDict.update(projectName)

        if databaseNameList is None:
            users = [getpass.getuser()]
            databaseNameList_full = dbListNames + users
            if isinstance(projectName, str):
                for user in numpy.atleast_1d(users):
                    projectNamesDict[user] = projectName
        else:
            if isinstance(projectName, str):
                for user in numpy.atleast_1d(databaseNameList):
                    projectNamesDict[user] = projectName
            databaseNameList_full = dbListNames + list(numpy.atleast_1d(databaseNameList))
        super().__init__(projectNamesDict,databaseNameList=databaseNameList_full, useAll=useAll)