from .collection import AbstractCollection
from .collection import Cache_Collection,Measurements_Collection,Simulations_Collection
import pandas
from functools import partial
from itertools import product



def getProjectList(user=None):
    """
        Return the list with the names of the existing projects .

    :param user: str
        The name of the database.

    :return:
        list.
    """
    return list(set(AbstractCollection(user=user).getProjectList()))



class ProjectMultiDB(object):
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

    _measurements = None
    _cache     = None
    _simulations  = None
    _users = None
    _useAll = None

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
    def users(self):
        return self._users

    @users.setter
    def users(self, newUsers):
        self._users = newUsers
        self._measurements  = [Measurements_Collection(user=user) for user in newUsers]
        self._cache      = [Cache_Collection(user=user) for user in newUsers]
        self._simulations   = [Simulations_Collection(user=user) for user in newUsers]

    @property
    def useAll(self):
        return self._useAll

    @useAll.setter
    def useAll(self,newUseAll):
        self._useAll=newUseAll

    @property
    def projectName(self):
        return self._projectName

    @projectName.setter
    def projectName(self, newProjectName):
        self._projectName = newProjectName

    def __init__(self, projectName,users=None, useAll=False):
        """
            Initialize the project.

        :param projectName: str
                The name of the project.

        :param user: str
                the name of the database to use. If None, use the default database (the name of the current user).

        """
        self._projectName = projectName
        self._users = users
        self._useAll = useAll
        self._measurements  = [Measurements_Collection(user=user) for user in users]
        self._cache      = [Cache_Collection(user=user) for user in users]
        self._simulations   = [Simulations_Collection(user=user) for user in users]

    def getMetadata(self):
        """
        Returns a pandas dataframe which contains all the description of all ot the documents in the current project.

        :return: pandas
        """
        descList = [doc.desc for doc in AbstractCollection().getDocuments(projectName=self._projectName)]
        return pandas.DataFrame(descList)

    def getSomeTypeDocumentsAsDict(self, searchtype, with_id, **kwargs):
        returnData = []
        for searched in searchtype:
            data = searched.getDocumentsAsDict(projectName=self._projectName, with_id=with_id, **kwargs)
            if len(data["documents"]) != 0:
                if self._useAll:
                    returnData.append(data)
                else:
                    returnData = data
                    break
        return returnData

    def getSomeTypeDocuments(self, searchtype, resource, dataFormat, type, **desc):
        returnData = []
        for searched in searchtype:
            data = searched.getDocuments(projectName=self._projectName, resource=resource, dataFormat=dataFormat,type=type, **desc)
            if len(data) != 0:
                if self._useAll:
                    returnData.append(data)
                else:
                    returnData = data
                    break
        return returnData

    def addSomeTypeDocuments(self, searchType, resource, dataFormat, type, users=None, **desc):
        if users is None:
            searchType[0].addDocument(projectName=self._projectName, resource=resource, dataFormat=dataFormat, type=type, desc=desc)
        else:
            for user in users:
                searchType[self._users.index(user)].addDocument(projectName=self._projectName, resource=resource, dataFormat=dataFormat, type=type, desc=desc)

    def deleteSomeTypeDocuments(self, searchType, users=None, **kwargs):
        if users is None:
            searchType[0].deleteDocuments(projectName=self._projectName, **kwargs)
        else:
            for user in users:
                searchType[self._users.index(user)].deleteDocuments(projectName=self._projectName, **kwargs)

    def getMeasurementsDocumentsAsDict(self, with_id=False, **kwargs):
        return self.getSomeTypeDocumentsAsDict(searchtype=self._measurements, with_id=with_id, **kwargs)

    def getMeasurementsDocuments(self, resource=None, dataFormat=None, type=None, **desc):
        return self.getSomeTypeDocuments(searchtype=self._measurements, resource=resource, dataFormat=dataFormat, type=type, **desc)

    def addMeasurementsDocument(self, resource="", dataFormat="string", type="", desc={}, users=None):
        return self.addSomeTypeDocuments(searchtype=self._measurements, resource=resource, dataFormat=dataFormat, type=type, desc=desc, users=users)

    def deleteMeasurementsDocuments(self, users=None, **kwargs):
        return self.deleteSomeTypeDocuments(searchType=self._measurements, users=users, **kwargs)

    def getSimulationsDocumentsAsDict(self, with_id=False, **kwargs):
        return self.getSomeTypeDocumentsAsDict(searchtype=self._simulations, with_id=with_id, **kwargs)

    def getSimulationsDocuments(self, resource=None, dataFormat=None, type=None, **desc):
        return self.getSomeTypeDocuments(searchtype=self._simulations, resource=resource, dataFormat=dataFormat, type=type, **desc)

    def addSimulationsDocument(self, resource="", dataFormat="string", type="", desc={}, users=None):
        return self.addSomeTypeDocuments(searchtype=self._simulations, resource=resource, dataFormat=dataFormat, type=type, desc=desc, users=users)

    def deleteSimulationsDocuments(self, users=None, **kwargs):
        return self.deleteSomeTypeDocuments(searchType=self._simulations,  users=users, **kwargs)

    def getCacheDocumentsAsDict(self,  with_id=False, **kwargs):
        return self.getSomeTypeDocumentsAsDict(searchtype=self._cache, with_id=with_id, **kwargs)

    def getCacheDocuments(self, resource=None, dataFormat=None, type=None, **desc):
        return self.getSomeTypeDocuments(searchtype=self._cache, resource=resource, dataFormat=dataFormat, type=type,**desc)

    def addCacheDocument(self, resource="", dataFormat="string", type="", desc={}, users=None):
        return self.addSomeTypeDocuments(searchtype=self._cache, resource=resource, dataFormat=dataFormat, type=type, desc=desc, users=users)

    def deleteCacheDocuments(self, users=None, **kwargs):
        return self.deleteSomeTypeDocuments(searchType=self._cache,  users=users, **kwargs)
