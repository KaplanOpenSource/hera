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

    def __init__(self, projectName,user=None):
        """
            Initialize the project.

        :param projectName: str
                The name of the project.

        :param user: str
                the name of the database to use. If None, use the default database (the name of the current user).

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
