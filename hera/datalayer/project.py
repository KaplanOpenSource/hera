from .collection import AbstractCollection
from .collection import Analysis_Collection,Measurements_Collection,Simulations_Collection
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




def getDocumentsAsDict(self, obj, with_id=False, **kwargs):
    return obj.getDocumentsAsDict(projectName=self._projectName, with_id=with_id, **kwargs)


def getDocuments(self, obj, resource=None, dataFormat=None, type=None, **desc):
    return obj.getDocuments(projectName=self._projectName, resource=resource, dataFormat=dataFormat, type=type, **desc)

def addDocument(self, obj, resource="", dataFormat="string", type="", desc={}):
    return obj.addDocument(projectName=self._projectName, resource=resource, dataFormat=dataFormat, type=type, desc=desc)

def deleteDocuments(self, obj,  **kwargs):
    return obj.deleteDocuments(projectName=self._projectName, **kwargs)


class Project(object):
    """
        Provides a simple interface to the data of a specific project.

        The class has all the following functions for the measurements, simulations and analysis:

        -  getMeasurementsDocumentsAsDict"
        -  getMeasurementsDocuments
        -  addMeasurementsDocument
        -  deleteMeasurementsDocuments


        -  getSimulationsDocumentsAsDict"
        -  getSimulationsDocuments
        -  addSimulationsDocument
        -  deleteSimulationsDocuments

        -  getAnalysisDocumentsAsDict"
        -  getAnalysisDocuments
        -  addAnalysisDocument
        -  deleteAnalysisDocuments

    """
    _projectName = None

    _measurements = None
    _analysis     = None
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
    def analysis(self):
        """
            Access the analysis type documents.

        :return:
            hera.datalayer.collection.Analysis_Collection

        """
        return self._analysis

    @property
    def simulations(self):
        """
            Access the analysis type documents.

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
        self._analysis      = Analysis_Collection(user=user)
        self._simulations   = Simulations_Collection(user=user)

        # registering the all the functions for measurements, simulations and analysis.
        #
        funcListName = ["get{func}DocumentsAsDict","get{func}Documents","add{func}Document","delete{func}Documents"]
        funcList = [getDocumentsAsDict,getDocuments,addDocument,deleteDocuments]
        func_class = ["measurements","simulations","analysis"]

        for fclass,(func,funcName) in product(func_class,zip(funcList,funcListName)):
            func_class_obj = getattr(self,fclass)
            partial_func = partial(func,self,func_class_obj)

            setattr(self,funcName.format(func=fclass.title()),partial_func)




    def getMetadata(self):
        """
        Returns a pandas dataframe which contains all the description of all ot the documents in the current project.

        :return: pandas
        """
        descList = [doc.desc for doc in AbstractCollection().getDocuments(projectName=self._projectName)]
        return pandas.DataFrame(descList)
