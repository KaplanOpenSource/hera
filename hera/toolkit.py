from .datalayer import Project,datatypes
import os
import pandas
import numpy
import pydoc

TOOLKIT_DATASOURCE_TYPE = "ToolkitDataSource"
TOOLKIT_TOOLKITNAME_FIELD       = "toolkit"
TOOLKIT_DATASOURCE_NAME = "datasourceName"
TOOLKIT_DATASOURCE_VERSION = "version"

TOOLKIT_SAVEMODE_NOSAVE = None
TOOLKIT_SAVEMODE_ONLYFILE = "File"
TOOLKIT_SAVEMODE_ONLYFILE_REPLACE = "File_overwrite"
TOOLKIT_SAVEMODE_FILEANDDB = "DB"
TOOLKIT_SAVEMODE_FILEANDDB_REPLACE = "DB_overwrite"

class ToolkitHome:

    ############### Duplicates of the above for comfortability
    TOOLKIT_SAVEMODE_NOSAVE = None
    TOOLKIT_SAVEMODE_ONLYFILE = "File"
    TOOLKIT_SAVEMODE_ONLYFILE_REPLACE = "File_overwrite"
    TOOLKIT_SAVEMODE_FILEANDDB = "DB"
    TOOLKIT_SAVEMODE_FILEANDDB_REPLACE = "DB_overwrite"
    ############################################################


    GIS_BUILDINGS  = "GIS_Buildings"
    GIS_TILES       = "GIS_Tiles"
    #GIS_RASTER     = "GIS_Raster"
    GIS_TOPOGRAPHY = "GIS_Topography"
    GIS_DEMOGRAPHY = "GIS_Demography"
    GIS_SHAPES     = "GIS_Shapes"
    RISKASSESSMENT = "RiskAssessment"
    LSM            = "LSM.old"
    OF_LSM         =  "OF_LSM"

    SIMULATIONS_WORKFLOWS = "hermesWorkflows"
    SIMULATIONS_OPENFOAM = "OpenFOAM"

    METEOROLOGY_HIGHFREQ = "MeteoHighFreq"
    METEOROLOGY_LOWFREQ = "MeteoLowFreq"

    EXPERIMENT = "experiment"

    _toolkits = None

    def __init__(self):
        self._toolkits = dict(

            Logging        = dict(cls = "hera.utils.logging.toolkit.loggingToolkit"),

            GIS_Vector=dict(cls="hera.measurements.GIS.vector.vector.VectorToolkit",desc=None),

            GIS_Buildings  = dict(cls = "hera.measurements.GIS.vector.buildings.toolkit.BuildingsToolkit",desc=None),

            #GIS_Raster     = dict(cls = "hera.measurements.GIS.raster.toolkit.RasterToolkit",desc=None),

            GIS_Tiles      =  dict(cls = "hera.measurements.GIS.raster.tiles.TilesToolkit",desc=None),

            GIS_Topography = dict(cls = "hera.measurements.GIS.vector.topography.TopographyToolkit",desc=None),

            GIS_Demography = dict(cls = "hera.measurements.GIS.vector.demography.DemographyToolkit",desc=None),

            GIS_Shapes     = dict(cls = "hera.measurements.GIS.shapes.ShapesToolKit",desc=None),

            RiskAssessment = dict(cls = "hera.riskassessment.riskToolkit.RiskToolkit",desc=None),
            LSM            = dict(cls = "hera.simulations.LSM.old.toolkit.LSMToolkit",desc=None),

            OF_LSM         = dict(cls="hera.simulations.openFoam.LSM.old.toolkit.OFLSMToolkit"),

            MeteoHighFreq  = dict(cls="hera.measurements.meteorology.highfreqdata.datalayer.HighFreqToolKit"),

            MeteoLowFreq = dict(cls="hera.measurements.meteorology.lowfreqdata.datalayer.lowFreqToolKit"),

            experiment =dict(cls="hera.measurements.experiment.experiment.experimentHome"),

            hermesWorkflows = dict(cls="hera.simulations.hermesWorkflowToolkit.workflowToolkit"),
            OpenFOAM = dict(cls="hera.simulations.openFoam.toolkit.OFToolkit")
        )


    def getToolkit(self, toolkitName, projectName, filesDirectory=None, **kwargs):
        """
            Returns a toolkit for the requested project.

        Parameters
        ----------
        projectName: str
            The name of the project


        filesDirectory: str
            The directory to save file (if necessary).
            If None, use the current directory.

        kwargs: dict
            The parameters for the toolkit.
            See the specific toolkit documentation for further details.

        Returns
        -------
            The tookit
        """
        if toolkitName not in self._toolkits.keys():
            raise ValueError(f"Toolkit name must be one of [{','.join(self._toolkits.keys())}]. Got {toolkitName} instead")
        clsName = self._toolkits[toolkitName]['cls']

        tookit = pydoc.locate(clsName)(projectName, filesDirectory=filesDirectory, **kwargs)
        return tookit



class abstractToolkit(Project):
    """
        A base class for Toolkits.

        *  Like project, it is initialized with a project name.
           If the toolkit works on data, it should be present in that project.

        *  Inherits from project and therefore exposes all the datalayer functions.

        *  Holds the toolkit name, and references to the datalayer and presentation layers.

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
    _projectName = None

    _analysis = None # holds the datalayer layer.
    _presentation = None # holds the presentation layer

    _FilesDirectory = None

    @property
    def FilesDirectory(self):
        """
            The directory to save files (when creating files).
        :return:
        """
        return self._FilesDirectory

    @property
    def filesDirectory(self):
        """
            The directory to save files (when creating files).
        :return:
        """
        return self._FilesDirectory


    @property
    def presentation(self):
        """
            Access to the presentation layer
        :return:
        """
        return self._presentation

    @property
    def analysis(self):
        """
            Access to the datalayer layer
        :return:
        """
        return self._analysis

    @property
    def toolkitName(self):
        """
            The name of the toolkit name
        :return:
        """
        return self._toolkitname

    @property
    def projectName(self):
        """
            The name of the project
        :return:
        """
        return self._projectName

    def __init__(self, toolkitName, projectName, filesDirectory=None):
        """
            Initializes a new toolkit.

        Parameters
        ----------

        toolkitName: str
            The name of the toolkit

        projectName: str
            The project that the toolkit works in.

        filesDirectory: str
            The directory to save datasource

        """
        super().__init__(projectName=projectName)
        self._toolkitname = toolkitName
        self._projectName = projectName

        if filesDirectory is None:
            self.logger.execution("Directory is not given, tries to load from default or using the current directory")
            self._FilesDirectory = self.getConfig().get("filesDirectory",os.getcwd())
            self.logger.execution(f"Using {self._FilesDirectory}")
        else:
            self.logger.execution(f"Using {os.path.abspath(filesDirectory)}. Creating if does not exist")
            os.system("mkdir -p %s" % os.path.abspath(filesDirectory))
            self._FilesDirectory = filesDirectory



    def getDataSourceMap(self,**filters):
        """
            Return the list of all data sources and their versions that are related to this toolkit

        Parameters
        ----------
            asPandas: bool
                If true, convert to pandas.

            filters: parameters
                Additional parameters to query the templates

        Returns
        -------
            list of dicts or pandas
        """
        docList = self.getMeasurementsDocuments(type=TOOLKIT_DATASOURCE_TYPE,
                                                toolkit=self.toolkitName,
                                                **filters)

        ret = []
        for doc in docList:
            ret.append(dict(datasourceName=doc.desc[TOOLKIT_DATASOURCE_NAME],
                            version=  str(doc.desc[TOOLKIT_DATASOURCE_VERSION]),
                            toolkit = doc.desc['toolkit'],
                            dataFormat=doc['dataFormat'],
                            resource=doc['resource']
                            )
                       )
        return ret

    def getDataSourceTable(self, **filters):

        Table = []
        for sourceMap in self.getDataSourceMap(**filters):
            table = pandas.json_normalize(sourceMap)
            Table.append(table)

        if len(Table) == 0:
            return pandas.DataFrame()
        else:
            return pandas.concat((Table),ignore_index=True)

    def getDatasourceDocumentsList(self, **kwargs):
        """
            Return all the datasources associated with this toolkit.

        Returns
        -------
            List of docs.
        """
        queryDict = {"type": TOOLKIT_DATASOURCE_TYPE,
                     TOOLKIT_TOOLKITNAME_FIELD : self.toolkitName}
        queryDict.update(**kwargs)
        return self.getMeasurementsDocuments(**queryDict)

    def getDatasourceDocument(self, datasourceName, version=None, **filters):
        """
            Return the document of the datasource.
            If version is not specified, return the latest version.

            Returns a single document.

        Parameters
        ----------
        datasourceName: str
            The datasourceName of the source
            if None, return the default source (if set).

        version: tuple
            The version of the source.
            if not found, return the latest source


        filters:
            Additional parameters to the query.

        Returns
        -------
                The document of the source. (None if not found)
        """
        if datasourceName is not None:
            filters[TOOLKIT_DATASOURCE_NAME] = datasourceName
        if version is not None:
            filters[TOOLKIT_DATASOURCE_VERSION] = version

        filters[TOOLKIT_TOOLKITNAME_FIELD] = self.toolkitName  # {'toolkit' : self.toolkitName}

        docList = self.getMeasurementsDocuments(type=TOOLKIT_DATASOURCE_TYPE, **filters)

        if len(docList) ==0:
            ret =  None

        elif len(docList) ==1:
            ret = docList[0]

        elif len(docList)>1:
            versionsList =[ doc['desc']['version'] for doc in docList]
            latestVersion = max(versionsList)
            docList = [doc for doc in docList if doc['desc']['version']==latestVersion]
            ret =docList[0]
        return ret

    def getDatasourceData(self, datasourceName=None, version=None,**filters):
        """
            Returns the data from the datasource.

        Parameters
        ----------

        datasourceName: str
            The datasourceName of the source
            if None, return the default source (if set).

        version: tuple
            The version of the source.
            if not found, return the latest source


        filters: dict
                additional filters to the query.

        Returns
        -------
                The data of the source. (None if not found)
        """
        filters[TOOLKIT_TOOLKITNAME_FIELD] = self.toolkitName  # {'toolkit' : self.toolkitName}
        doc = self.getDatasourceDocument(datasourceName=datasourceName, version=version,**filters)

        return None if doc is None else doc.getData()

    def addDataSource(self,dataSourceName,resource,dataFormat,version=(0,0,1),overwrite=False,**kwargs):
        """
            Adds a resource to the toolkit.
            The type is always TOOLKIT_DATASOURCE_TYPE.
            The toolkit name is added to the description.

        Parameters
        ----------
        dataSourceName: str
                The name of the data source

        version: tuple (of int)
                A 3-tuple of the version

        resource: str
                The resource

        dataFormat: str
                A string of a datatypes.

        kwargs: dict
                The parameters

        Returns
        -------
            The document of the datasource.
        """

        kwargs[TOOLKIT_TOOLKITNAME_FIELD] = self.toolkitName
        kwargs[TOOLKIT_DATASOURCE_NAME] = dataSourceName
        kwargs[TOOLKIT_DATASOURCE_VERSION] = version

        if (self.getDatasourceDocument(dataSourceName,version=version) is None) or overwrite:
            if self.getDatasourceDocument(dataSourceName,version=version) is not None:  # not None = Exist
                # print("Delete existing, and add new data source.")
                delargs = {TOOLKIT_DATASOURCE_NAME : dataSourceName,
                           TOOLKIT_DATASOURCE_VERSION :  version}

                self.deleteDataSourceDocuments(**delargs)
            #else:
                # print("Does not exist: add data source.")

            doc = self.addMeasurementsDocument(type=TOOLKIT_DATASOURCE_TYPE,
                                               resource=resource,
                                               dataFormat=dataFormat,
                                               desc=kwargs)
        else:
            raise ValueError(f"Record {dataSourceName} (version {version}) already exists in project {self.projectName}. use overwrite=True to overwrite on the existing document")
            print("exist: Raise exception (ValueError) that the record with the name that was given in the input already exists")

        return doc


    def deleteDataSourceDocuments(self,**filters):

        filters[TOOLKIT_TOOLKITNAME_FIELD] = self.toolkitName # {'toolkit' : self.toolkitName}

        return self.deleteMeasurementsDocuments(type=TOOLKIT_DATASOURCE_TYPE,**filters)

