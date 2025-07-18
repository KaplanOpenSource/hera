import json


from .datalayer import Project
import os
import pandas
import numpy
import pydoc
from .utils.logging import get_classMethod_logger

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
    GIS_LANDCOVER   = "GIS_LandCover"
    #GIS_RASTER     = "GIS_Raster"
    GIS_VECTOR_TOPOGRAPHY = "GIS_Vector_Topography"
    GIS_RASTER_TOPOGRAPHY = "GIS_Raster_Topography"
    GIS_DEMOGRAPHY = "GIS_Demography"
    GIS_SHAPES     = "GIS_Shapes"
    RISKASSESSMENT = "RiskAssessment"
    LSM            = "LSM"

    DATA           = "heraData"

    SIMULATIONS_WORKFLOWS = "hermesWorkflows"
    SIMULATIONS_OPENFOAM = "OpenFOAM"

    METEOROLOGY_HIGHFREQ = "MeteoHighFreq"
    METEOROLOGY_LOWFREQ = "MeteoLowFreq"

    EXPERIMENT = "experiment"

    WINDPROFILE = "WindProfile"
    GAUSSIANDISPERSION = "GaussianDispersion"

    _toolkits = None

    def __init__(self):
        self._toolkits = dict(
            GIS_Buildings  = dict(cls = "hera.measurements.GIS.vector.buildings.toolkit.BuildingsToolkit",desc=None),
            GIS_Tiles      =  dict(cls = "hera.measurements.GIS.raster.tiles.TilesToolkit",desc=None),
            GIS_Vector_Topography = dict(cls = "hera.measurements.GIS.vector.topography.TopographyToolkit",desc=None),
            GIS_Raster_Topography = dict(cls = "hera.measurements.GIS.raster.topography.TopographyToolkit",desc=None),
            GIS_Demography = dict(cls = "hera.measurements.GIS.vector.demography.DemographyToolkit",desc=None),
            GIS_LandCover     = dict(cls = "hera.measurements.GIS.raster.landcover.LandCoverToolkit",desc=None),

            RiskAssessment = dict(cls = "hera.riskassessment.riskToolkit.RiskToolkit",desc=None),
            LSM            = dict(cls = "hera.simulations.LSM.toolkit.LSMToolkit",desc=None),

            OF_LSM         = dict(cls="hera.simulations.openFoam.LSM.toolkit.OFLSMToolkit"),

            MeteoHighFreq  = dict(cls="hera.measurements.meteorology.highfreqdata.toolkit.HighFreqToolKit"),

            MeteoLowFreq = dict(cls="hera.measurements.meteorology.lowfreqdata.toolkit.lowFreqToolKit"),

            experiment =dict(cls="hera.measurements.experiment.experiment.experimentHome"),

            hermesWorkflows = dict(cls="hera.simulations.hermesWorkflowToolkit.hermesWorkflowToolkit"),
            OpenFOAM = dict(cls="hera.simulations.openFoam.toolkit.OFToolkit"),

            WindProfile = dict(cls="hera.simulations.windProfile.toolkit.WindProfileToolkit"),
            GaussianDispersion = dict(cls="hera.simulations.gaussian.toolkit.gaussianToolkit")

        )


    def getToolkit(self, toolkitName, projectName=None, filesDirectory=None, **kwargs):
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
        super().__init__(projectName=projectName,filesDirectory=filesDirectory)
        logger = get_classMethod_logger(self,"init")
        self._toolkitname = toolkitName

    @property
    def classLoggerName(self):
        return str(get_classMethod_logger(self,"{the_function_name}")).split(" ")[1]

    def getDataSourceList(self,**filters):
        """
            Returns a list of the data source names
        Parameters
        ----------
        filters

        Returns
        -------

        """
        docList = self.getMeasurementsDocuments(type=TOOLKIT_DATASOURCE_TYPE,
                                                toolkit=self.toolkitName,
                                                **filters)

        ret = []
        for doc in docList:
            ret.append(doc['desc']['datasourceName'])

        return ret

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
            dta = dict(dataFormat=doc['dataFormat'],
                 resource=doc['resource'])
            dta.update(doc.desc)
            ret.append(dta)
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

    def getDataSourceDocumentsList(self, **kwargs):
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

    def getDataSourceDocument(self, datasourceName, version=None, **filters):
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
        else:
            try:
                defaultVersion = self.getConfig()[f"{datasourceName}_defaultVersion"]
                filters[TOOLKIT_DATASOURCE_VERSION] = defaultVersion
            except:
                pass


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

    def getDataSourceDocuments(self, datasourceName, version=None, **filters):
        """
            Returns a list with the datasource. This is for the complteness of the interface.
            That is, making it similar to the Measurement, Cache and Simulation document retrieval.

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
                A list that containes the document of the source. (empty list  if not found)
        """
        doc = self.getDataSourceDocument(datasourceName=datasourceName, version=version, **filters)
        return [] if doc is None else [doc]

    def getDataSourceData(self, datasourceName=None, version=None, **filters):
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
        doc = self.getDataSourceDocument(datasourceName=datasourceName, version=version, **filters)
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
        if (self.getDataSourceDocument(dataSourceName, version=version) is None) or overwrite:
            if self.getDataSourceDocument(dataSourceName, version=version) is not None:  # not None = Exist
                # print("Delete existing, and add new data source.")
                delargs = {TOOLKIT_DATASOURCE_NAME : dataSourceName,
                           TOOLKIT_DATASOURCE_VERSION :  version}

                self.deleteDataSource(**delargs)
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


    def deleteDataSource(self, datasourceName, version=None, **filters):

        doc = self.getDataSourceDocument(datasourceName=datasourceName, version=version, **filters)
        doc.delete()

        return doc


    def setDataSourceDefaultVersion(self,datasourceName:str,version:tuple):
        if len(self.getMeasurementsDocuments(type="ToolkitDataSource", **{"datasourceName": datasourceName ,
                                                                            "version": version}))==0:
            raise ValueError(f"No DataSource with name={datasourceName} and version={version}.")

        self.setConfig(**{f"{datasourceName}_defaultVersion": version})
        print(f"{version} for dataSource {datasourceName} is now set to default.")

