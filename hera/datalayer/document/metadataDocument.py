from mongoengine import *
import json
from hera.datalayer.datahandler import getHandler

class MetadataFrame(object):
    """
        A basic structure for a document.

        Each document is related to a project and described by the following fields:

        - type : str : The type of the document.
                       This is an helper attribute that is used to query the data.

        - resource: str: The resource that the document represents.
                         This can be either path to a file on the disk or the data itself.

        - dataFormat : str: The format of the data. Taken from ::class:`..datatypes.datatypes`

        - desc: dict: A dictionary of arbitrary format that holds the metadata of the record.

        - id : str : The id of the record in the DB.

    """
    projectName = StringField(required=True)
    desc = DictField(required=False)
    type = StringField(required=True)
    resource = DynamicField(required=True)
    dataFormat = StringField(required=True)

    def asDict(self, with_id=False):
        """
        Convert the document to a plain dictionary.

        Parameters
        ----------
        with_id : bool, optional
            If True, include the ``_id`` key. Default is False.

        Returns
        -------
        dict
            Dictionary representation of the document.
        """
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
        storeParametersDict = self.desc.get("storeParameters",{})
        storeParametersDict.update(kwargs)
        return getHandler(self.dataFormat).getData(resource=self.resource,desc=self.desc, **storeParametersDict)

    def __str__(self):
        return json.dumps(self.asDict(with_id=False),indent=4)


class nonDBMetadataFrame(object):
    """
        A wrapper class to use when the data is not loaded into the
        DB.

        This class will be used when getting data from local files.
    """
    _data = None

    def __init__(self, data, projectName=None, type=None, resource=None, dataFormat=None, **desc):
        """
        Initialize a non-database metadata frame.

        Parameters
        ----------
        data : object
            The data to wrap.
        projectName : str, optional
            The project name.
        type : str, optional
            The document type.
        resource : str, optional
            The resource path or identifier.
        dataFormat : str, optional
            The data format name.
        desc : dict
            Additional metadata fields.
        """
        self.projectName = projectName
        self.type = type
        self.resource = resource
        self.dataFormat = dataFormat
        self.desc = desc

        self._data = data

    def getData(self, **kwargs):
        """
        Return the wrapped data object.

        Returns
        -------
        object
            The data passed at initialization.
        """
        return self._data

    def __getitem__(self, item):
        """
        Access document attributes by key.

        Parameters
        ----------
        item : str
            The attribute name.

        Returns
        -------
        object
        """
        return self.__dict__[item]

