# from .document.metadataDocument import Metadata,Measurements,Simulations,Analysis,Projects
from . import getDBObject
from mongoengine import ValidationError, MultipleObjectsReturned, DoesNotExist
import warnings


class AbstractCollection(object):
    _metadataCol = None
    _type = None

    @property
    def type(self):
        return self._type

    def __init__(self, ctype=None, user=None):
        self._type = ctype
        self._metadataCol = getDBObject('Metadata', user) if self.type is None else getDBObject(ctype, user)

    def getDocumentsAsDict(self, projectName, with_id=False, **kwargs):
        """
        Returns a dict with a 'documents' key and list of documents in a dict formats as value.
        The list of the documents are the result of your query.

        :param projectName: The projectName.
        :param with_id: rather or not should the id key be in the documents.
        :param kwargs: query arguments.
        :return: dict with 'documents' key and list of dicts in its value.
        """
        dictList = [doc.asDict(with_id=with_id) for doc in self.getDocuments(projectName=projectName, **kwargs)]
        return dict(documents=dictList)

    def getDocuments(self, projectName, resource=None, dataFormat=None, type=None, **desc):
        """
        Get the documents that satisfy the given query.
        If projectName is None search over all projects.

        :param projectName: The project name.
        :param resource: The data resource.
        :param dataFormat: The data format.
        :param type: The type which the data belongs to.
        :param desc: Other metadata arguments.
        :return:
        """
        query = {}
        if resource is not None:
            query['resource'] = resource
        if dataFormat is not None:
            query['dataFormat'] = dataFormat
        if type is not None:
            query['type'] = type
        if projectName is not None:
            query['projectName'] = projectName
        for key, value in desc.items():
                query['desc__%s' % key] = value
        return self._metadataCol.objects(**query)

    def _getAllValueByKey(self, key, **query):
        return [doc[key] for doc in self.getDocuments(projectName=None, **query)]

    def getProjectList(self):
        return self._getAllValueByKey(key='projectName')

    def getDocumentByID(self, id):
        """
        Returns a document by its ID.

        :param id: The document ID.
        :return: document
        """
        return self._metadataCol.objects.get(id=id)

    def addDocument(self,projectName,resource="",dataFormat="string",type="",desc={}):
        """
            Adds a document to the database.


        :param projectName: str
                The project to add the document
        :param resource:
                The data of the document.
        :param dataFormat: str
                The type of the dataformat.
                See datahandler for the available types.

        :param desc: dict
                Holds any additional fields that describe the
        :param type:

        :param kwargs:
        :return:
            mongoDB.documents.
        """
        try:
            obj = self._metadataCol(projectName=projectName,resource=resource,dataFormat=dataFormat,type=type,desc=desc).save()
        except ValidationError as e:
            raise ValidationError("Not all of the required fields are delivered "
                                  "or one of the fields type is not proper. %s " % str(e))
        return obj

    def addDocumentFromJSON(self, json_data):
        self._metadataCol.from_json(json_data).save()

    def deleteDocuments(self, projectName, **kwargs):
        """
        Deletes documents that satisfy the given query.

        :param projectName: The project name.
        :param kwargs: Other query arguments.
        :return:
        """
        for doc in self.getDocuments(projectName=projectName, **kwargs):
            doc.delete()

    def deleteDocumentByID(self, id):
        """
        Deletes a documents by its ID.

        :param id: The document ID.
        :return:
        """
        self.getDocumentByID(id=id).delete()

    def getData(self, projectName, usePandas=None, **kwargs):
        """
        Returns the data by the given parameters.

        This is obeselete. Use getDocuments and the get the data from the desired document.

        :param projectName: The name of the project.
        :param usePandas: Return the data as pandas if True, dask if False. Default is None, and returns the data depends on the data format.
        :param kwargs: Other properties of the data.
        :return: pandas/dask dataframe.
        """
        warnings.warn('getData is going to be deprecated in the next version.'
                      'use getDocuments to get documents and use getData of the document object' ,DeprecationWarning)
        docList = self.getDocuments(projectName=projectName, **kwargs)
        if usePandas is None:
            return [doc.getData() for doc in docList]
        else:
            return [doc.getData(usePandas=usePandas) for doc in docList]


class Measurements_Collection(AbstractCollection):

    def __init__(self, user=None):
        super().__init__(ctype='Measurements', user=user)
    #
    # def meta(self):
    #     return self._metadataCol


class Simulations_Collection(AbstractCollection):

    def __init__(self, user=None):
        super().__init__(ctype='Simulations', user=user)


class Analysis_Collection(AbstractCollection):

    def __init__(self, user=None):
        super().__init__(ctype='Analysis', user=user)
