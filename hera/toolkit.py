from .datalayer import project,datatypes

class toolkit(project):
    """
        A base class for toolkits.

        *  Like project, it is initialized with a project name.
           If the toolkit works on data, it should be present in that project.

        *  Iherits from project and therefore exposes all the datalayer functions.

        *  Holds the toolkit name, and references to the analysis and presentation layers.

        *  Adds a mechanism (setConfig,getConfig) for saving configuration in the DB. the settings are specific for a project.

    """
    _toolkitname = None

    _analysis = None # holds the analysis layer.
    _presentation = None # holds the presentation layer


    @property
    def presentation(self):
        return self._presentation

    @property
    def analysis(self):
        return self._analysis

    @property
    def name(self):
        return self._toolkitname


    @classmethod
    def listSources(cls):
        """
            Lists the sources that are associated with this class.

        :return:
        """
        pass

    @classmethod
    def getSourceSignature(cls):
        """
            Return the source signature name (will be the type of the document).


        :return:
        """
        pass


    def __init__(self,toolkitName,projectName):
        """
            Initializes a new toolkit.

        :param toolkitName: str
                    The name of the toolkit

        :param projectName: str
                    The project that the toolkit works in.
        """
        super().__init__(projectName=projectName)
        self._toolkitname = toolkitName

    def _getConfigDocument(self):
        """
        Returns the document of the config.
        If there is no config document, return empty dictionary.

        :return: dict
                The configuration of the toolkit.
        """
        """

        """
        documents = self.getCacheDocumentsAsDict(type=f"{self.projectName}__{self.name}__config__")
        if len(documents) == 0:
            self.addCacheDocument(type=f"{self.projectName}__{self.name}__config__",
                                  resource="",
                                  dataFormat=datatypes.STRING,
                                  desc={})

        return documents[0]

    def getConfig(self):
        """
        Returns the config document's description.
        If there is no config document, return empty dictionary.

        :return: dict
                The configuration of the toolkit.
        """
        """
        
        """
        doc = self._getConfigDocument()
        return doc['desc']


    def setConfig(self, config):
        """
        Create a config documnet or updates an existing config document.
        """
        doc = self._getConfigDocument()
        doc['desc'].update(config)
        doc.save()

    def getSource(self,name=None):
        """
            Return the document of the signature.

        :param name: str
                The name of the source
        :return:
                the document of the source. (None if not found)
        """
        source_signature = self.getSourceSignature()

        docList = self.getMeasurementsDocuments(type=source_signature,name=name)
        if len(docList) ==0:
            return None
        else:
            return docList[0]

