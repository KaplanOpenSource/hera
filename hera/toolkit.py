from .datalayer import project,datatypes
import pandas
import numpy

TOOLKIT_DATASOURCE_TYPE = "ToolkitDataSource"
TOOLKIT_TOOLKITNAME_FIELD       = "toolkit"
TOOLKIT_DATASOURCE_NAME = "datasourceName"
TOOLKIT_DATASOURCE_VERSION = "version"

TOOLKIT_SAVEMODE_NOSAVE = None
TOOLKIT_SAVEMODE_ONLYFILE = "File"
TOOLKIT_SAVEMODE_ONLYFILE_REPLACE = "File_overwrite"
TOOLKIT_SAVEMODE_FILEANDDB = "DB"
TOOLKIT_SAVEMODE_FILEANDDB_REPLACE = "DB_overwrite"


class abstractToolkit(project):
    """
        A base class for Toolkits.

        *  Like project, it is initialized with a project name.
           If the toolkit works on data, it should be present in that project.

        *  Inherits from project and therefore exposes all the datalayer functions.

        *  Holds the toolkit name, and references to the analysis and presentation layers.

        *  Adds a mechanism (setConfig,getConfig) for saving configuration in the DB. the settings are specific for a project.

        *  Adds a mechanism to list, get and add data sources.
            - A data source will always be saved as a measurement document.
            - Each source has the following properties in the description (except for the other properties):
                    * name : str
                    * toolkit : str
                    * projectName :str
                    * version : tuple (major version, minor varsion, bug fix).
                    * the type is TOOLKIT_DATASOURCE_TYPE.
                    * metadata: dict with additional metadata of the datasource.

            - The toolkit can have a default source for the project.
                    A default data source is defined with its name and version
                    If the version is not supplied, takes the latest version.
            -

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


    def __init__(self,toolkitName,projectName):
        """
            Initializes a new toolkit.

        Params

        toolkitName: str
                    The name of the toolkit

        projectName: str
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

        Returns
        -------
        dict
                The configuration of the toolkit.
        """
        doc = self._getConfigDocument()
        return doc.desc

    def setConfig(self, **kwargs):
        """
            Create a config document or updates an existing config document.
        """
        doc = self._getConfigDocument()
        doc.desc.update(kwargs)
        doc.save()

    def setDefaultDataSource(self, datasourceName, version=None):
        """
            Set the default datasource.
            If the version is not supplied, take the latest.

            Raise ValueError if source datasourceName is not present in the DB.

        Parameters
        ----------
        datasourceName: str
            The datasourceName of the data source to set as default.
        version: tuple
            If not present, take the latest version in the project DB.

        Returns
        -------
            None
        """
        sourceDoc = self.getDatasourceDocument(datasourceName=datasourceName, version=version)

        if sourceDoc is None:
            raise ValueError(f"Data Source {datasourceName} is not found in the project {self.projectName}")

        if isinstance(sourceDoc,list):
            # There is more than one version, take the latest.
            versionsList =[ doc.desc[TOOLKIT_DATASOURCE_VERSION] for doc in sourceDoc]
            latestVersion = numpy.max(versionsList,axis=0)
        else:
            latestVersion = version

        self.setConfig(toolkit=dict(name=datasourceName, version=latestVersion))

    def getSourceList(self,asPandas=True):
        """
            Return the list of sources and their versions that are related to this toolkit

        Parameters
        ----------
            asPandas: bool
                If true, convert to pandas.

        Returns
        -------
            list of dicts or pandas
        """
        docList = self.getMeasurementsDocuments(type=TOOLKIT_DATASOURCE_TYPE,
                                                toolkit=self.name)

        ret = []
        for doc in docList:
            ret.append(dict(datasourceName=doc.desc[TOOLKIT_DATASOURCE_NAME],
                            version=f"({','.join(doc.desc[TOOLKIT_DATASOURCE_VERSION])}",
                            toolkit = doc.desc['toolkit'],
                            dataFormat=doc['dataFormat'],
                            resource=doc['resource']
                            )
                       )

        return pandas.DataFrame(ret) if asPandas else ret

    def getDatasourceDocument(self, datasourceName, version=None,**kwargs):
        """
            Return the document of the datasource.
            Returns a single document.

        Parameters
        ------
        datasourceName: str
            The datasourceName of the source
            if None, return the default source (if set).

        version: tuple
            The version of the source.
            if not found, return the latest source


        Returns
        -------
                The document of the source. (None if not found)
        """
        docList = self.getMeasurementsDocuments(type=TOOLKIT_DATASOURCE_TYPE,
                                                name=datasourceName,
                                                toolkit=self.name,
                                                version=version,**kwargs)

        if len(docList) ==0:
            ret =  None
        elif len(docList)>1:
            versionsList =[ doc['desc']['version'] for doc in docList]
            latestVersion = numpy.max(versionsList,axis=0)
            docList = [doc for doc in docList if doc['desc']['version']==latestVersion]
            ret =docList[0]
        return ret

    def _getDataDescriptor(self):
        """
            Returns a map with the description of the toolkit.
            Used to save the analysis products of the toolkit such that they will be considered as a
            datasources.

            the dict holds the keys:
                type : the const TOOLKIT_DATASOURCE_TYPE
                toolkit: the name of the toolkit


        :return: dict

        """
        return dict(type=TOOLKIT_DATASOURCE_TYPE, toolkit=self.name)

    def getDatasourceData(self, datasourceName=None, version=None):
        """
            Returns the data from the datasource.

        Parameters
        ------
        datasourceName: str
            The datasourceName of the source
            if None, return the default source (if set).

        version: tuple
            The version of the source.
            if not found, return the latest source


        Returns
        -------
                The data of the source. (None if not found)
        """
        doc = self.getDatasourceDocument(datasourceName=datasourceName, version=version)
        return None if doc is None else doc.getData()


    def addDataSource(self,dataSourceName,resource,dataFormat,version=None,**kwargs):
        """
            Adds a resource to the toolkit.
            The type is always TOOLKIT_DATASOURCE_TYPE.
            The toolkit name is added to the description.

        :param dataSourceName:
        :param version:
        :param resource:
        :param dataFormat:
        :param kwargs:
        :return:
        """

        kwargs[TOOLKIT_TOOLKITNAME_FIELD] = self.name
        kwargs[TOOLKIT_DATASOURCE_NAME] = dataSourceName
        kwargs[TOOLKIT_DATASOURCE_VERSION] = version

        doc  = self.addMeasurementsDocument(type=TOOLKIT_DATASOURCE_TYPE,
                                     resource=resource,
                                     dataFormat=dataFormat,
                                     desc=kwargs)
        return doc

