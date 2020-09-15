from mongoengine import *
import json
from ..datahandler import getHandler

class MetadataFrame(object):
    """
        A basic structure for a document.
    """
    projectName = StringField(required=True)
    desc = DictField(required=True)
    type = StringField(required=True)
    resource = DynamicField(required=True)
    dataFormat = StringField(required=True)

    def asDict(self, with_id=False):
        docDict = json.loads(self.to_json())
        if not with_id:
            docDict.pop('_id')
        # docDict.pop('_cls')
        return docDict

    def getData(self, **kwargs):
        """
        Returns the data of the document.

        the kwargs passed to the datahandler.
        See the datahandler class for your specific datatype.

        Parameters
        ----------
        kwargs : dict
        
        Returns
        -------
            object according to the datahandler. 
        """
        return getHandler(self.dataFormat).getDocFromDB(self.resource, **kwargs)


class nonDBMetadataFrame(object):
    """
        A wrapper class to use when the data is not loaded into the
        DB.

        This class will be used when getting data from local files.
    """
    _projectName = None
    _desc = None
    _type = None
    _resource = None
    _dataFormat = None

    _data = None

    def __init__(self,data, projectName=None, type=None, resource=None, dataFormat=None,**desc):
        self._projectName = projectName
        self._type = type
        self._resource = resource
        self._dataFormat = dataFormat
        self._desc = desc

        self._data = data

    @property
    def projectName(self):
        return self._projectName

    @property
    def type(self):
        return self._type

    @property
    def resource(self):
        return self._resource

    @property
    def dataFormat(self):
        return self._dataFormat

    @property
    def desc(self):
        return self._desc


    def getData(self, **kwargs):
        return self._data


