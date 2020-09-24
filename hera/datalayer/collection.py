from . import getDBObject
from mongoengine import ValidationError, MultipleObjectsReturned, DoesNotExist
import warnings
import sys
version = sys.version_info[0]

class AbstractCollection(object):
    _metadataCol = None
    _type = None

    @property
    def type(self):
        return self._type

    def __init__(self, ctype=None, user=None):
        self._type = ctype
        self._metadataCol = getDBObject('Metadata', user) if self.type is None else getDBObject(ctype, user)

    def getDocumentsAsDict(self, projectName, with_id=False, **query):
        """
        Returns a dict with a 'documents' key and list of documents in a dict formats as value.
        The list of the documents are the result of your query.

        Parameters
        ----------
        projectName : str
            The projectName.

        with_id : bool, optional, default False
            rather or not should the 'id' key be in the documents.

        query :
            query arguments.

        Returns
        -------
        dict
            A dict with 'documents' key and the value is a list of dicts that represent the documents that fulfills the query.
        """
        dictList = [doc.asDict(with_id=with_id) for doc in self.getDocuments(projectName=projectName, **query)]
        return dict(documents=dictList)

    def getDocuments(self, projectName, resource=None, dataFormat=None, type=None, **desc):
        """
        Get the documents that satisfy the given query.
        If projectName is None search over all projects.

        Parameters
        ----------
        projectName : str
            The project name.

        resource :
            The data resource.

        dataFormat : str
            The data format.
        type : str
            The type which the data belongs to.
        desc :
            Other metadata arguments.

        Returns
        -------
        list
            List of documents that fulfill the query.
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
        return list(set([doc[key] for doc in self.getDocuments(projectName=None, **query)]))

    def getProjectList(self):
        return self._getAllValueByKey(key='projectName')

    def getDocumentByID(self, id):
        """
        Returns a document by its ID.

        Parameters
        ----------
        id : str
            The document ID.

        Returns
        -------
        document
            The document with the relevant ID.
        """
        return self._metadataCol.objects.get(id=id)

    def addDocument(self,projectName,resource="",dataFormat="string",type="",desc={}):
        """
            Adds a document to the database.

        Parameters
        ----------
        projectName : str
            The project to add the document

        resource :
            The data of the document.

        dataFormat : str
            The type of the dataformat.
            See datahandler for the available types.

        desc : dict
            Holds any additional fields that describe the

        type : str
            The type of the data

        Returns
        -------
        mongoengine document
        """
        try:
            obj = self._metadataCol(projectName=projectName,resource=resource,dataFormat=dataFormat,type=type,desc=desc).save()
        except ValidationError as e:
            raise ValidationError("Not all of the required fields are delivered "
                                  "or one of the fields type is not proper. %s " % str(e))
        return obj

    def addDocumentFromJSON(self, json_data):
        self._metadataCol.from_json(json_data).save()

    def deleteDocuments(self, projectName, **query):
        """
        Deletes documents that satisfy the given query.

        Parameters
        ----------
        projectName : str
            The project name.

        query :
            Other query arguments.

        Returns
        -------

        """
        for doc in self.getDocuments(projectName=projectName, **query):
            doc.delete()

    def deleteDocumentByID(self, id):
        """
        Deletes a documents by its ID.

        Parameters
        ----------
        id : str
            The document ID.

        Returns
        -------

        """
        doc = self.getDocumentByID(id=id)
        doc.delete()
        return doc

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
            return [doc.getDocFromDB() for doc in docList]
        else:
            return [doc.getDocFromDB(usePandas=usePandas) for doc in docList]


class Measurements_Collection(AbstractCollection):

    def __init__(self, user=None):
        if version == 2:
            super(Measurements_Collection, self).__init__(ctype='Measurements', user=user)
        elif version == 3:
            super().__init__(ctype='Measurements', user=user)
    #
    # def meta(self):
    #     return self._metadataCol


class Simulations_Collection(AbstractCollection):

    def __init__(self, user=None):
        if version == 2:
            super(Simulations_Collection, self).__init__(ctype='Simulations', user=user)
        elif version == 3:
            super().__init__(ctype='Simulations', user=user)

class Cache_Collection(AbstractCollection):

    def __init__(self, user=None):
        if version == 2:
            super(Cache_Collection, self).__init__(ctype='Cache', user=user)
        elif version == 3:
            super().__init__(ctype='Cache', user=user)
