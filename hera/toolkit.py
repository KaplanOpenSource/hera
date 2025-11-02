from hera.datalayer import Project
from hera.datalayer.datahandler import datatypes  # for datatypes.CLASS
from hera.datalayer.datahandler import DataHandler_Class  # הוסף אם לא קיים

import inspect
import os
import pandas
import pydoc
from hera.utils.logging import get_classMethod_logger

TOOLKIT_DATASOURCE_TYPE = "ToolkitDataSource"
TOOLKIT_TOOLKITNAME_FIELD = "toolkit"
TOOLKIT_DATASOURCE_NAME = "datasourceName"
TOOLKIT_DATASOURCE_VERSION = "version"

TOOLKIT_SAVEMODE_NOSAVE = None
TOOLKIT_SAVEMODE_ONLYFILE = "File"
TOOLKIT_SAVEMODE_ONLYFILE_REPLACE = "File_overwrite"
TOOLKIT_SAVEMODE_FILEANDDB = "DB"
TOOLKIT_SAVEMODE_FILEANDDB_REPLACE = "DB_overwrite"

import pydoc
import pandas as pd
from hera.utils.data.toolkit_repository import ToolkitRepository  # new import for DB integration


class ToolkitHome:
    """
    Central registry for available toolkits (static + dynamic).
    Provides:
      - getToolkit(toolkitName, ...): locate & instantiate a toolkit class
      - getToolkitTable(projectName): table of all toolkits (static + DB)
      - registerToolkit(toolkitclass, ...): register a class into project datasources (dataFormat=Class)
    """

    # -------- Save modes (kept for compatibility) --------
    TOOLKIT_SAVEMODE_NOSAVE = None
    TOOLKIT_SAVEMODE_ONLYFILE = "File"
    TOOLKIT_SAVEMODE_ONLYFILE_REPLACE = "File_overwrite"
    TOOLKIT_SAVEMODE_FILEANDDB = "DB"
    TOOLKIT_SAVEMODE_FILEANDDB_REPLACE = "DB_overwrite"

    # -------- Static toolkit identifiers --------
    GIS_BUILDINGS = "GIS_Buildings"
    GIS_TILES = "GIS_Tiles"
    GIS_LANDCOVER = "GIS_LandCover"
    GIS_VECTOR_TOPOGRAPHY = "GIS_Vector_Topography"
    GIS_RASTER_TOPOGRAPHY = "GIS_Raster_Topography"
    GIS_DEMOGRAPHY = "GIS_Demography"
    GIS_SHAPES = "GIS_Shapes"
    RISKASSESSMENT = "RiskAssessment"
    LSM = "LSM"

    DATA = "heraData"

    SIMULATIONS_WORKFLOWS = "hermesWorkflows"
    SIMULATIONS_OPENFOAM = "OpenFOAM"

    METEOROLOGY_HIGHFREQ = "MeteoHighFreq"
    METEOROLOGY_LOWFREQ = "MeteoLowFreq"

    EXPERIMENT = "experiment"

    WINDPROFILE = "WindProfile"
    GAUSSIANDISPERSION = "GaussianDispersion"
    MACHINELEARNING_DEEPLEARNING = "machine_deep_learning"

    _toolkits = None

    def __init__(self):
        # Static built-in toolkits (internal source)
        self._toolkits = dict(
            GIS_Buildings=dict(
                cls="hera.measurements.GIS.vector.buildings.toolkit.BuildingsToolkit",
                desc=None,
                type="measurements"
            ),
            GIS_Tiles=dict(
                cls="hera.measurements.GIS.raster.tiles.TilesToolkit",
                desc=None,
                type="measurements"
            ),
            GIS_Vector_Topography=dict(
                cls="hera.measurements.GIS.vector.topography.TopographyToolkit",
                desc=None,
                type="measurements"
            ),
            GIS_Raster_Topography=dict(
                cls="hera.measurements.GIS.raster.topography.TopographyToolkit",
                desc=None,
                type="measurements"
            ),
            GIS_Demography=dict(
                cls="hera.measurements.GIS.vector.demography.DemographyToolkit",
                desc=None,
                type="measurements"
            ),
            GIS_LandCover=dict(
                cls="hera.measurements.GIS.raster.landcover.LandCoverToolkit",
                desc=None,
                type="measurements"
            ),
            RiskAssessment=dict(
                cls="hera.riskassessment.riskToolkit.RiskToolkit",
                desc=None,
                type="riskassessment"
            ),
            LSM=dict(
                cls="hera.simulations.LSM.toolkit.LSMToolkit",
                desc=None,
                type="simulations"
            ),
            OF_LSM=dict(
                cls="hera.simulations.openFoam.LSM.toolkit.OFLSMToolkit",
                desc=None,
                type="simulations"
            ),
            MeteoHighFreq=dict(
                cls="hera.measurements.meteorology.highfreqdata.toolkit.HighFreqToolKit",
                desc=None,
                type="measurements"
            ),
            MeteoLowFreq=dict(
                cls="hera.measurements.meteorology.lowfreqdata.toolkit.lowFreqToolKit",
                desc=None,
                type="measurements"
            ),
            hermesWorkflows=dict(
                cls="hera.simulations.hermesWorkflowToolkit.hermesWorkflowToolkit",
                desc=None,
                type="simulations"
            ),
            OpenFOAM=dict(
                cls="hera.simulations.openFoam.toolkit.OFToolkit",
                desc=None,
                type="simulations"
            ),
            WindProfile=dict(
                cls="hera.simulations.windProfile.toolkit.WindProfileToolkit",
                desc=None,
                type="simulations"
            ),
            GaussianDispersion=dict(
                cls="hera.simulations.gaussian.toolkit.gaussianToolkit",
                desc=None,
                type="simulations"
            ),
            machine_deep_learning=dict(
                cls="hera.simulations.machineLearningDeepLearning.toolkit.machineLearningDeepLearningToolkit",
                desc=None,
                type="simulations"
            ),

        )

    # --- Place this near the top of the file imports if needed ---
    from hera.datalayer import Project

    # --- Inside class ToolkitHome, replace ONLY the "not found" branch in getToolkit(...) ---

    def getToolkit(self, toolkitName, projectName=None, filesDirectory=None, **kwargs):
        """
        Locate a toolkit class by name (static registry or DB), then instantiate it.
        If not found anywhere, return a lightweight fallback that wraps Project so that
        repository JSON loading can still proceed without a concrete Toolkit class.
        """
        # 1) Static registry (unchanged)
        if toolkitName in self._toolkits:
            clsName = self._toolkits[toolkitName]['cls']
            toolkitClass = pydoc.locate(clsName)
            if toolkitClass is None:
                raise ImportError(f"Cannot locate class: {clsName}")
            return toolkitClass(projectName, filesDirectory=filesDirectory, **kwargs)

        # 2) Dynamic registry via DB (unchanged)
        repo = ToolkitRepository(projectName or "DefaultProject")
        doc = repo.getToolkitDocument(toolkitName)
        if doc:
            desc = getattr(doc, "desc", None) or (doc.get("desc", {}) if isinstance(doc, dict) else {})
            resource = getattr(doc, "resource", None) or (doc.get("resource", "") if isinstance(doc, dict) else "")
            classpath = desc.get("classpath") or desc.get("cls")
            if classpath:
                norm_desc = dict(desc)
                norm_desc["classpath"] = classpath
                norm_desc.pop("cls", None)
                # Use the dynamic Class loader path when classpath exists
                return DataHandler_Class.getData(resource=resource, desc=norm_desc)
            # If there is a dynamic doc but no classpath, we'll fall through to the shim

        # 3) Fallback SHIM: no static and no usable dynamic class found -> wrap Project
        class _FallbackToolkit(Project):
            """
            Minimal shim so repository JSON loader can operate on the project
            even without a concrete Toolkit class.
            """

            def __init__(self, toolkitName, projectName, filesDirectory=None):
                super().__init__(projectName=projectName)
                self._toolkitname = toolkitName
                self._projectName = projectName
                self._filesDirectory = filesDirectory

            # Keep parity with expected attributes in loaders
            @property
            def toolkitName(self):
                return self._toolkitname

            @property
            def projectName(self):
                return self._projectName

            # Map the methods used by loaders to Project APIs

            # DataSource layer
            def getDataSourceDocuments(self, **qry):
                return super().getDataSourceDocuments(**qry)

            def deleteDataSource(self, *, datasourceName):
                # remove by name if exists
                try:
                    docs = super().getDataSourceDocuments(datasourceName=datasourceName)
                    for d in docs:
                        d.delete()
                except Exception:
                    pass

            def addDataSource(self, **item):
                # Repository JSON usually passes resource, dataFormat, dataSourceName, etc.
                return super().addDataSource(**item)

            # Measurements / Cache / Simulations generic helpers used by _DocumentHandler
            def getMeasurementsDocuments(self, **qry):
                return super().getMeasurementsDocuments(**qry)

            def addMeasurementsDocument(self, **kwargs):
                return super().addMeasurementsDocument(**kwargs)

            def getCacheDocuments(self, **qry):
                return super().getCacheDocuments(**qry)

            def addCacheDocument(self, **kwargs):
                return super().addCacheDocument(**kwargs)

            def getSimulationsDocuments(self, **qry):
                return super().getSimulationsDocuments(**qry)

            def addSimulationsDocument(self, **kwargs):
                return super().addSimulationsDocument(**kwargs)

            # Config used by _handle_Config
            def setConfig(self, **cfg):
                current = super().getConfig() or {}
                current.update(cfg)
                super().setConfig(**current)

        # Return the shim instance
        return _FallbackToolkit(toolkitName=toolkitName, projectName=projectName, filesDirectory=filesDirectory)

    # hera/toolkit.py  (inside class ToolkitHome)
    # -----------------------------------------------------------------------------
    # Auto-register a missing toolkit using classpath hints (from repository JSON
    # or from the Toolkit document in DB) and then return an instance via getToolkit.
    # -----------------------------------------------------------------------------
    def auto_register_and_get(self,
                              toolkitName: str,
                              projectName: str,
                              repositoryJSON: dict = None,
                              repositoryName: str = None,
                              params: dict = None,
                              version: tuple = (0, 0, 1)):
        """
        Attempts to auto-register a missing toolkit and return an instance.
        1) Try to find a classpath hint in the repositoryJSON (if provided).
           We look for keys like: repositoryJSON[toolkitName]["Registry"]["classpath"]
           or ...["Registry"]["cls"].
        2) If not found, try the DB-backed Toolkit document (ToolkitRepository).
        3) Import the class, choose a repository to register into:
             - repositoryName argument if provided,
             - else the project's default repository (must exist).
        4) Register via registerToolkit(...), then getToolkit(...) and return it.
        """
        params = params or {}
        classpath_hint = None

        # 1) Classpath hint in the repository JSON
        if repositoryJSON:
            try:
                tk_section = repositoryJSON.get(toolkitName, {})
                reg = tk_section.get("Registry", {})
                classpath_hint = reg.get("classpath") or reg.get("cls")
            except Exception:
                pass

        # 2) If still not found, try DB Toolkit document
        if not classpath_hint:
            from hera.utils.data.toolkit_repository import ToolkitRepository
            repo = ToolkitRepository(projectName)
            doc = repo.getToolkitDocument(toolkitName)
            if doc and getattr(doc, "desc", None):
                classpath_hint = doc.desc.get("classpath") or doc.desc.get("cls")

        if not classpath_hint:
            raise ValueError(
                f"auto_register_and_get: no classpath hint found for toolkit '{toolkitName}'. "
                f"Provide a 'Registry.classpath'/'cls' in repository JSON or seed a Toolkit document."
            )

        # Import the class
        mod_name, _, cls_name = classpath_hint.rpartition(".")
        if not mod_name or not cls_name:
            raise ValueError(f"Invalid classpath hint: '{classpath_hint}'")
        from importlib import import_module
        try:
            mod = import_module(mod_name)
            toolkit_cls = getattr(mod, cls_name)
        except Exception as exc:
            raise ImportError(f"Failed to import '{classpath_hint}' for toolkit '{toolkitName}'") from exc

        # Decide target repository for registration
        repo_to_use = repositoryName or self.getDefaultRepository(projectName=projectName)
        if not repo_to_use:
            raise ValueError(
                f"auto_register_and_get: no target repository for project '{projectName}'. "
                f"Set a default repository or pass repositoryName explicitly."
            )

        # Register (idempotent if overwrite=True)
        self.registerToolkit(
            toolkitclass=toolkit_cls,
            datasource_name=toolkitName,
            params=params,
            version=version,
            overwrite=True,
            projectName=projectName,
            repositoryName=repo_to_use,
        )

        # Return an instance
        return self.getToolkit(toolkitName=toolkitName, projectName=projectName)

    def getToolkitTable(self, projectName):
        """
        Return a DataFrame that merges static toolkits with dynamic toolkits from DB.
        """
        repo = ToolkitRepository(projectName)
        dynamic = repo.getToolkitTable()

        static = []
        for name, info in self._toolkits.items():
            static.append({
                "toolkit": name,
                "cls": info["cls"],
                "source": "internal",
                "type": info.get("type", "measurements"),
                "description": info.get("desc", "")
            })

        df_static = pd.DataFrame(static)
        all_toolkits = pd.concat([df_static, dynamic], ignore_index=True).drop_duplicates("toolkit")
        return all_toolkits

    # בתוך class ToolkitHome (בקובץ hera/toolkit.py)

    def registerToolkit(
            self,
            toolkitclass,
            *,
            projectName,
            repositoryName,  # <<< חדש: דרישת רפוזיטורי
            datasource_name=None,
            params=None,
            version=(0, 0, 1),
            overwrite=False,
    ):
        """
        Register a toolkit class as a datasource document in the given project & repository.

        It stores:
          - resource: the directory that contains the module file (DataHandler_Class adds to sys.path)
          - dataFormat: datatypes.CLASS
          - desc: {
                'toolkit'       : <datasource_name>,
                'datasourceName': <datasource_name>,
                'repository'    : <repositoryName>,   # <<< נשמר במסמך
                'version'       : (major, minor, patch),
                'classpath'     : '<module.Class>',
                'parameters'    : { ... }
            }
        """
        if projectName is None:
            raise ValueError("registerToolkit: 'projectName' is required")
        if not repositoryName:
            raise ValueError("registerToolkit: 'repositoryName' is required")

        import inspect, os
        module_path = inspect.getfile(toolkitclass)
        resource_dir = os.path.dirname(os.path.abspath(module_path))
        classpath = f"{toolkitclass.__module__}.{toolkitclass.__qualname__}"

        ds_name = datasource_name or toolkitclass.__name__
        params = params or {}

        desc = {
            "toolkit": ds_name,
            "datasourceName": ds_name,
            "repository": repositoryName,  # <<< שדה רפוזיטורי
            "version": tuple(version),
            "classpath": classpath,
            "parameters": params,
        }

        proj = Project(projectName=projectName)

        # בדיקת קיום לפי (type, repository, datasourceName, version)
        existing = proj.getMeasurementsDocuments(
            type="ToolkitDataSource",
            repository=repositoryName,  # <<< סינון לפי רפוזיטורי
            datasourceName=ds_name,
            version=tuple(version),
        )
        if existing:
            if not overwrite:
                raise ValueError(
                    f"Toolkit datasource '{ds_name}' (version {version}) already exists in "
                    f"repository '{repositoryName}' of project '{projectName}'. "
                    f"Use overwrite=True to replace."
                )
            for doc in existing:
                doc.delete()

        doc = proj.addMeasurementsDocument(
            type="ToolkitDataSource",
            resource=resource_dir,
            dataFormat=datatypes.CLASS,
            desc=desc,
        )
        return doc

    def setDefaultRepository(self, *, projectName: str, repositoryName: str, overwrite: bool = True):
        """
        Persist default repository name for a project so future calls can omit --repository.
        We store it as a tiny Project document under type='RepositoryConfig'.
        """
        if not projectName:
            raise ValueError("setDefaultRepository: 'projectName' is required")
        if not repositoryName:
            raise ValueError("setDefaultRepository: 'repositoryName' is required")

        proj = Project(projectName=projectName)
        # delete previous config if exists (by type)
        if overwrite:
            old = proj.getMeasurementsDocuments(type="RepositoryConfig")
            for d in old:
                d.delete()

        desc = {"defaultRepository": repositoryName}

        # Try to pick a dataFormat constant if available. Fallback: omit the arg.
        df_arg = {}
        try:
            from hera.datalayer import datatypes as _dt
            dfmt = getattr(_dt, "JSON", None) or getattr(_dt, "json", None) or getattr(_dt, "TEXT", None)
            if dfmt is not None:
                df_arg["dataFormat"] = dfmt
        except Exception:
            pass

        return proj.addMeasurementsDocument(
            type="RepositoryConfig",
            resource=".",  # trivial
            desc=desc,
            **df_arg,
        )

    def getDefaultRepository(self, *, projectName: str) -> str:
        """
        Read the saved default repository name (if exists). Returns '' if missing.
        """
        if not projectName:
            raise ValueError("getDefaultRepository: 'projectName' is required")
        proj = Project(projectName=projectName)
        docs = proj.getMeasurementsDocuments(type="RepositoryConfig")
        if not docs:
            return ""
        # Take the newest (or first)
        return docs[0].desc.get("defaultRepository", "") or ""

    def getDatasourceDocument(
            self,
            *,
            projectName: str,
            datasourceName: str,
            repositoryName: str = None,
            version=None,  # tuple like (0,0,1) or None
    ):
        """
        Fetch a ToolkitDataSource by (repository, datasourceName [, version]).
        If repositoryName is None or '', fallback to the project's defaultRepository.
        """
        if not projectName:
            raise ValueError("getDatasourceDocument: 'projectName' is required")
        if not datasourceName:
            raise ValueError("getDatasourceDocument: 'datasourceName' is required")

        repo = (repositoryName or "").strip()
        if not repo:
            repo = self.getDefaultRepository(projectName=projectName)
            if not repo:
                raise ValueError(
                    "Repository name is not provided and no defaultRepository is set for the project. "
                    "Call setDefaultRepository(...) first, or pass repositoryName explicitly."
                )

        proj = Project(projectName=projectName)

        q = {
            "type": "ToolkitDataSource",
            "repository": repo,
            "datasourceName": datasourceName,
        }
        if version is not None:
            q["version"] = tuple(version)

        docs = proj.getMeasurementsDocuments(**q)
        return docs[0] if docs else None


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

    _analysis = None  # holds the datalayer layer.
    _presentation = None  # holds the presentation layer

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
        super().__init__(projectName=projectName, filesDirectory=filesDirectory)
        logger = get_classMethod_logger(self, "init")
        self._toolkitname = toolkitName

    @property
    def classLoggerName(self):
        return str(get_classMethod_logger(self, "{the_function_name}")).split(" ")[1]

    def getDataSourceList(self, **filters):
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

    def getDataSourceMap(self, **filters):
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
            return pandas.concat((Table), ignore_index=True)

    def getDataSourceDocumentsList(self, **kwargs):
        """
            Return all the datasources associated with this toolkit.

        Returns
        -------
            List of docs.
        """
        queryDict = {"type": TOOLKIT_DATASOURCE_TYPE,
                     TOOLKIT_TOOLKITNAME_FIELD: self.toolkitName}
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

        if len(docList) == 0:
            ret = None

        elif len(docList) == 1:
            ret = docList[0]

        elif len(docList) > 1:
            versionsList = [doc['desc']['version'] for doc in docList]
            latestVersion = max(versionsList)
            docList = [doc for doc in docList if doc['desc']['version'] == latestVersion]
            ret = docList[0]
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

    def addDataSource(self, dataSourceName, resource, dataFormat, version=(0, 0, 1), overwrite=False, **kwargs):
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
                delargs = {TOOLKIT_DATASOURCE_NAME: dataSourceName,
                           TOOLKIT_DATASOURCE_VERSION: version}

                self.deleteDataSource(**delargs)
            # else:
            # print("Does not exist: add data source.")

            doc = self.addMeasurementsDocument(type=TOOLKIT_DATASOURCE_TYPE,
                                               resource=resource,
                                               dataFormat=dataFormat,
                                               desc=kwargs)
        else:
            raise ValueError(
                f"Record {dataSourceName} (version {version}) already exists in project {self.projectName}. use overwrite=True to overwrite on the existing document")
            print(
                "exist: Raise exception (ValueError) that the record with the name that was given in the input already exists")

        return doc

    def deleteDataSource(self, datasourceName, version=None, **filters):

        doc = self.getDataSourceDocument(datasourceName=datasourceName, version=version, **filters)
        doc.delete()

        return doc

    def setDataSourceDefaultVersion(self, datasourceName: str, version: tuple):
        if len(self.getMeasurementsDocuments(type="ToolkitDataSource", **{"datasourceName": datasourceName,
                                                                          "version": version})) == 0:
            raise ValueError(f"No DataSource with name={datasourceName} and version={version}.")

        self.setConfig(**{f"{datasourceName}_defaultVersion": version})
        print(f"{version} for dataSource {datasourceName} is now set to default.")
