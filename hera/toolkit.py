from hera.datalayer import Project
from hera.datalayer.datahandler import datatypes  # for datatypes.CLASS

import os
import inspect
import pydoc
import pandas as pd
from typing import Optional, List, Dict

from hera.utils.logging import get_classMethod_logger

# Optional: keep these constants if they are used elsewhere
TOOLKIT_DATASOURCE_TYPE = "ToolkitDataSource"
TOOLKIT_TOOLKITNAME_FIELD = "toolkit"
TOOLKIT_DATASOURCE_NAME = "datasourceName"
TOOLKIT_DATASOURCE_VERSION = "version"

TOOLKIT_SAVEMODE_NOSAVE = None
TOOLKIT_SAVEMODE_ONLYFILE = "File"
TOOLKIT_SAVEMODE_ONLYFILE_REPLACE = "File_overwrite"
TOOLKIT_SAVEMODE_FILEANDDB = "DB"
TOOLKIT_SAVEMODE_FILEANDDB_REPLACE = "DB_overwrite"



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
            experiment=dict(
                cls="hera.measurements.experiment.experiment.experimentHome",
                desc=None,
                type="measurements"
            ),
        )

    def _get_data_toolkit(self, projectName: str = None):
        """
        Helper that returns a dataToolkit instance.

        We import dataToolkit lazily here to avoid circular imports
        between hera.toolkit and hera.utils.data.toolkit.
        dataToolkit itself works on the DEFAULT project internally.
        """
        from hera.utils.data.toolkit import dataToolkit
        return dataToolkit()



    def getToolkit(self, toolkitName, projectName=None, filesDirectory=None, **kwargs):
        """
        Locate and instantiate a toolkit by name.

        Resolution order:
          1) Static registry (self._toolkits).
          2) Dynamic ToolkitDataSource document in the project (type='ToolkitDataSource'),
             then use doc.getData(), which is the standard Hera mechanism to load a class.
        """
        # 1) Static built-in toolkits
        if toolkitName in (self._toolkits or {}):
            info = self._toolkits[toolkitName]
            cls_path = info.get("cls")
            toolkit_cls = pydoc.locate(cls_path)
            if toolkit_cls is None:
                raise ImportError(f"Cannot locate class: {cls_path}")
            return toolkit_cls(projectName, filesDirectory=filesDirectory, **kwargs)

        # 2) Dynamic: look for a ToolkitDataSource in the project
        if not projectName:
            raise ValueError(
                f"Toolkit '{toolkitName}' not found in static registry and no projectName was provided "
                f"to search for a dynamic ToolkitDataSource."
            )

        proj = Project(projectName=projectName)

        # We expect toolkits to be stored as ToolkitDataSource measurements
        docs = proj.getMeasurementsDocuments(
            type="ToolkitDataSource",
            datasourceName=toolkitName,
        )

        if not docs:
            raise ValueError(
                f"Toolkit '{toolkitName}' not found in static registry or as ToolkitDataSource "
                f"in project '{projectName}'."
            )

        # Take the first matching document (or later you can add version/default logic)
        doc = docs[0]

        # Standard Hera way: the document knows how to load itself via datalayer handlers
        tk = doc.getData()

        # Optionally, pass filesDirectory or other kwargs to the toolkit if it supports it
        if hasattr(tk, "setFilesDirectory") and filesDirectory is not None:
            tk.setFilesDirectory(filesDirectory)

        return tk

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

    def getToolkitTable(self, projectName: Optional[str]):
        """
        Build a DataFrame from getToolkitDocuments(...).
        This avoids duplicated logic and guarantees both static + DB rows are represented.
        """
        docs = self.getToolkitDocuments(name=None, projectName=projectName) or []
        rows = []
        for d in docs:
            desc = d.get("desc", {})
            rows.append({
                "toolkit": d.get("toolkit", ""),
                "cls": desc.get("classpath", ""),
                "source": desc.get("source", ""),
                "type": desc.get("type", ""),
                "repositoryName": desc.get("repositoryName", ""),
                "version": desc.get("version", ""),
            })

        if not rows:
            return pd.DataFrame(columns=["toolkit", "cls", "source", "type", "repositoryName", "version"])

        # Drop duplicates by (toolkit, source) to avoid double rows for the same name/source
        df = pd.DataFrame(rows).drop_duplicates(subset=["toolkit", "source"], keep="first")
        return df


    def registerToolkit(
            self,
            toolkitclass,
            *,
            projectName,
            repositoryName,  # required repository
            datasource_name=None,
            params=None,
            version=(0, 0, 1),
            overwrite=False,
    ):
        """
        Register a toolkit class as a datasource document in the given project & repository.

        It stores:
          - resource: the directory that contains the module file
          - dataFormat: datatypes.CLASS
          - desc: {
                'toolkit'       : <datasource_name>,
                'datasourceName': <datasource_name>,
                'repository'    : <repositoryName>,
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
            "repository": repositoryName,
            "version": tuple(version),
            "classpath": classpath,
            "parameters": params,
        }

        proj = Project(projectName=projectName)

        # Check existence by (type, repository, datasourceName, version)
        existing = proj.getMeasurementsDocuments(
            type="ToolkitDataSource",
            repository=repositoryName,
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
        The configuration is stored as a measurement document with type='RepositoryConfig'.
        """
        if not projectName:
            raise ValueError("setDefaultRepository: 'projectName' is required")
        if not repositoryName:
            raise ValueError("setDefaultRepository: 'repositoryName' is required")

        dt = self._get_data_toolkit(projectName=projectName)

        # delete previous config if exists (by type)
        if overwrite:
            old = dt.getMeasurementsDocuments(type="RepositoryConfig")
            for d in old:
                d.delete()

        desc = {"defaultRepository": repositoryName}

        df_arg = {}
        try:
            from hera.datalayer import datatypes as _dt
            dfmt = getattr(_dt, "JSON", None) or getattr(_dt, "json", None) or getattr(_dt, "TEXT", None)
            if dfmt is not None:
                df_arg["dataFormat"] = dfmt
        except Exception:
            pass

        return dt.addMeasurementsDocument(
            type="RepositoryConfig",
            resource=".",
            desc=desc,
            **df_arg,
        )

    def getDefaultRepository(self, *, projectName: str) -> str:
        """
        Read the saved default repository name (if exists). Returns '' if missing.
        """
        if not projectName:
            raise ValueError("getDefaultRepository: 'projectName' is required")

        dt = self._get_data_toolkit(projectName=projectName)
        docs = dt.getMeasurementsDocuments(type="RepositoryConfig")
        if not docs:
            return ""
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
        If repositoryName is None or empty, fallback to the project's defaultRepository.
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

        dt = self._get_data_toolkit(projectName=projectName)

        query = {
            "type": "ToolkitDataSource",
            "repository": repo,
            "datasourceName": datasourceName,
        }
        if version is not None:
            query["version"] = tuple(version)

        docs = dt.getMeasurementsDocuments(**query)
        return docs[0] if docs else None

    # --- inside class ToolkitHome (toolkit.py) ---
    from typing import Optional, List, Dict

    def getToolkitDocuments(self, name: Optional[str] = None, projectName: Optional[str] = None) -> List[Dict]:
        """
        Single source of truth for listing toolkits.
        Returns a uniform list of "document-like" dicts coming from:
          1) Static registry (self._toolkits)
          2) Dynamic DB documents (type='ToolkitDataSource') of the given project.
        """
        docs: List[Dict] = []

        # 1) Static: normalized entries from the in-memory registry
        for tk_name, info in (self._toolkits or {}).items():
            if name and tk_name != name:
                continue
            docs.append({
                "toolkit": tk_name,
                "desc": {
                    "classpath": info.get("cls", ""),
                    "type": info.get("type", "measurements"),
                    "source": "internal",
                    "repositoryName": "",
                    "version": "",
                }
            })

        # 2) Dynamic (DB): query the project's measurements for type='ToolkitDataSource'
        if projectName:
            try:
                dt = self._get_data_toolkit(projectName=projectName)
                dyn_docs = dt.getMeasurementsDocuments(type="ToolkitDataSource") or []
                for d in dyn_docs:
                    try:
                        desc = getattr(d, "desc", {}) or {}
                        tk_name = desc.get("datasourceName") or getattr(d, "datasourceName", None)
                        if not tk_name:
                            continue
                        if name and tk_name != name:
                            continue

                        docs.append({
                            "toolkit": tk_name,
                            "desc": {
                                "classpath": desc.get("classpath", ""),
                                "type": desc.get("type", "") or "measurements",
                                "source": "db",
                                "repositoryName": desc.get("repository", "") or getattr(d, "repository", ""),
                                "version": tuple(desc.get("version", ())) or getattr(d, "version", ""),
                            }
                        })
                    except Exception:
                        # Be forgiving for partially-formed rows
                        pass
            except Exception:
                # No project/DB available: static list is still valuable
                pass

        return docs

    # --- Add inside class ToolkitHome (toolkit.py) ---

    def import_toolkits_from_json(self, *, projectName: str, json_path: str,
                                  default_repository: str = None, overwrite: bool = True) -> list:
        """
        Read a JSON file and register all Toolkits it declares into the given project.
        """
        import json
        from pydoc import locate

        if not projectName:
            raise ValueError("import_toolkits_from_json: projectName is required")

        with open(json_path, "r") as f:
            payload = json.load(f) or {}

        repo_from_json = (payload.get("repository") or "").strip()
        repo_to_use = repo_from_json or (default_repository or self.getDefaultRepository(projectName=projectName))
        if not repo_to_use:
            raise ValueError(
                f"No repository provided in JSON and no default repository set for project '{projectName}'. "
                f"Set one via ToolkitHome.setDefaultRepository(...) or pass default_repository."
            )

        toolkits = payload.get("toolkits") or []
        if not isinstance(toolkits, list):
            raise ValueError("import_toolkits_from_json: 'toolkits' must be a list")

        registered = []
        for item in toolkits:
            name = item.get("name")
            classpath = item.get("classpath")
            version = item.get("version", [0, 0, 1])
            params = item.get("parameters", {})

            if not name or not classpath:
                raise ValueError(f"Toolkit entry missing name/classpath: {item}")

            tk_class = locate(classpath)
            if tk_class is None:
                raise ImportError(f"Cannot locate class by classpath: {classpath}")

            self.registerToolkit(
                toolkitclass=tk_class,
                projectName=projectName,
                repositoryName=repo_to_use,
                datasource_name=name,
                params=params,
                version=tuple(int(x) for x in version),
                overwrite=overwrite,
            )
            registered.append(name)

        return registered


    def import_experiments_from_json(self, *, projectName: str, json_path: str) -> list:
        """
        Load experiments from the same JSON file into the given project.

        Each entry under 'experiments' is of the form:
          {
            "name": "Haifa2014",
            "data": [
              {
                "type": "Experiment_rawData",
                "resource": "data/experiment/data/Sonic.parquet",
                "dataFormat": "parquet",
                "desc": {...},
                "isRelativePath": true
              },
              ...
            ]
          }

        We write them as measurement documents into the *projectName* project.
        """
        import json
        import os
        from hera.datalayer import Project

        if not projectName:
            raise ValueError("import_experiments_from_json: projectName is required")

        # Load JSON payload
        with open(json_path, "r") as f:
            payload = json.load(f) or {}

        experiments = payload.get("experiments") or []
        if not isinstance(experiments, list):
            raise ValueError("'experiments' must be a list in the JSON file")

        # Work on the target project (NOT defaultProject)
        proj = Project(projectName=projectName)

        loaded = []
        base_dir = os.path.dirname(os.path.abspath(json_path))

        for exp in experiments:
            exp_name = exp.get("name")
            data_items = exp.get("data", [])

            for di in data_items:
                typ = di.get("type")
                resource = di.get("resource")
                data_fmt = di.get("dataFormat")
                desc = di.get("desc", {})
                is_rel = bool(di.get("isRelativePath", False))

                if not typ or not resource or not data_fmt:
                    # Skip invalid rows quietly
                    continue

                # Resolve relative path against JSON location
                res_path = resource
                if is_rel:
                    res_path = os.path.abspath(os.path.join(base_dir, resource))

                proj.addMeasurementsDocument(
                    type=typ,
                    resource=res_path,
                    dataFormat=data_fmt,
                    desc=desc,
                )

            if exp_name and exp_name not in loaded:
                loaded.append(exp_name)

        return loaded


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
        """
        Build a pandas DataFrame from all data sources of this toolkit.
        """
        tables = []

        for sourceMap in self.getDataSourceMap(**filters):
            table = pd.json_normalize(sourceMap)
            tables.append(table)

        if not tables:
            return pd.DataFrame()
        else:
            return pd.concat(tables, ignore_index=True)

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