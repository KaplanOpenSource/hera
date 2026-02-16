from hera.datalayer import Project
from hera.datalayer.datahandler import datatypes  # for datatypes.CLASS

import os
import sys
import logging
import inspect
import pydoc
import pandas as pd
from typing import Optional, List, Dict, Any

from hera.utils.logging import get_classMethod_logger
from sympy.physics.units import second

# ---------------------------------------------------------------------------
# Constants for Toolkit data sources
# ---------------------------------------------------------------------------

TOOLKIT_DATASOURCE_TYPE = "ToolkitDataSource"
TOOLKIT_TOOLKITNAME_FIELD = "toolkit"
TOOLKIT_DATASOURCE_NAME = "datasourceName"
TOOLKIT_DATASOURCE_VERSION = "version"

TOOLKIT_SAVEMODE_NOSAVE = None
TOOLKIT_SAVEMODE_ONLYFILE = "File"
TOOLKIT_SAVEMODE_ONLYFILE_REPLACE = "File_overwrite"
TOOLKIT_SAVEMODE_FILEANDDB = "DB"
TOOLKIT_SAVEMODE_FILEANDDB_REPLACE = "DB_overwrite"


# ======================================================================
#  abstractToolkit
# ======================================================================

class abstractToolkit(Project):
    """
    Base class for Toolkits.

    * Inherits from Project – ולכן יש גישה לכל פונקציות ה־datalayer.
    * מחזיק toolkitName ו־projectName.
    * מוסיף מנגנון data sources שנשמרים כ-measurement documents מסוג
      TOOLKIT_DATASOURCE_TYPE.
    """

    _toolkitname = None
    _projectName = None

    _analysis = None  # holds the datalayer layer.
    _presentation = None  # holds the presentation layer

    @property
    def presentation(self):
        """Access to the presentation layer."""
        return self._presentation

    @property
    def analysis(self):
        """Access to the datalayer layer."""
        return self._analysis

    @property
    def toolkitName(self):
        """The name of the toolkit."""
        return self._toolkitname

    @property
    def projectName(self):
        """The name of the project."""
        return self._projectName

    def __init__(self, toolkitName: str, projectName: Optional[str] = None, filesDirectory: Optional[str] = None):
        """
        Initializes a new toolkit.

        Parameters
        ----------
        toolkitName : str
            The name of the toolkit.

        projectName : str or None
            The project that the toolkit works in.
            If None – Project's automatic project-name mechanism is used.

        filesDirectory : str
            Directory to save datasource files.
        """
        super().__init__(projectName=projectName, filesDirectory=filesDirectory)
        logger = get_classMethod_logger(self, "init")
        self._toolkitname = toolkitName
        self._projectName = projectName

    @property
    def classLoggerName(self):
        return str(get_classMethod_logger(self, "{the_function_name}")).split(" ")[1]

    # ------------------------------------------------------------------
    # Data sources API
    # ------------------------------------------------------------------

    def getDataSourceList(self, **filters) -> List[str]:
        """Returns a list of data source names for this toolkit."""
        docList = self.getMeasurementsDocuments(
            type=TOOLKIT_DATASOURCE_TYPE,
            toolkit=self.toolkitName,
            **filters,
        )

        ret = []
        for doc in docList:
            ret.append(doc["desc"]["datasourceName"])
        return ret

    def getDataSourceMap(self, **filters) -> List[Dict[str, Any]]:
        """
        Return list of all data sources and their versions related to this toolkit.
        """
        docList = self.getMeasurementsDocuments(
            type=TOOLKIT_DATASOURCE_TYPE,
            toolkit=self.toolkitName,
            **filters,
        )

        ret = []
        for doc in docList:
            dta = dict(dataFormat=doc["dataFormat"], resource=doc["resource"])
            dta.update(doc.desc)
            ret.append(dta)
        return ret

    def getDataSourceTable(self, **filters) -> pd.DataFrame:
        """Build a pandas DataFrame from all data sources of this toolkit."""
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
        Return all the data source documents associated with this toolkit.
        """
        queryDict = {
            "type": TOOLKIT_DATASOURCE_TYPE,
            TOOLKIT_TOOLKITNAME_FIELD: self.toolkitName,
        }
        queryDict.update(**kwargs)
        return self.getMeasurementsDocuments(**queryDict)

    def getDataSourceDocument(self, datasourceName: Optional[str], version=None, **filters):
        """
        Return the document of the datasource.
        If version is not specified, return the latest version or default version (if configured).
        """
        if datasourceName is not None:
            filters[TOOLKIT_DATASOURCE_NAME] = datasourceName
        if version is not None:
            filters[TOOLKIT_DATASOURCE_VERSION] = version
        else:
            try:
                defaultVersion = self.getConfig()[f"{datasourceName}_defaultVersion"]
                filters[TOOLKIT_DATASOURCE_VERSION] = defaultVersion
            except Exception:
                pass

        filters[TOOLKIT_TOOLKITNAME_FIELD] = self.toolkitName

        docList = self.getMeasurementsDocuments(
            type=TOOLKIT_DATASOURCE_TYPE,
            **filters,
        )

        if len(docList) == 0:
            ret = None
        elif len(docList) == 1:
            ret = docList[0]
        else:
            versionsList = [doc["desc"]["version"] for doc in docList]
            latestVersion = max(versionsList)
            docList = [doc for doc in docList if doc["desc"]["version"] == latestVersion]
            ret = docList[0]
        return ret

    def getToolkitDocument(self, toolkit_name: str):
        """
        Find a dynamic toolkit document by name (either desc.datasourceName or desc.toolkit).
        Returns the mongoengine document or None.
        """
        # First: direct filter on datasourceName (works on most implementations)
        try:
            q = self.getMeasurementsDocuments(
                type="ToolkitDataSource", datasourceName=toolkit_name
            )
            if q and len(q) > 0:
                return q[0]
        except Exception:
            # fall through to broader search below
            pass

        # Second: scan all ToolkitDataSource docs and match by desc fields
        try:
            q = self.getMeasurementsDocuments(type="ToolkitDataSource")
            for d in q:
                desc = d.desc or {}
                if desc.get("datasourceName") == toolkit_name or desc.get("toolkit") == toolkit_name:
                    return d
        except Exception:
            pass

        # Optional: also look in DataSource collection if your project uses it
        try:
            q = self.getDataSourceDocuments(datasourceName=toolkit_name)
            if q and len(q) > 0:
                return q[0]
        except Exception:
            pass

        return None


    def getDataSourceDocuments(self, datasourceName, version=None, **filters):
        """
        Returns a list with the datasource (for API symmetry with measurements/cache).
        """
        doc = self.getDataSourceDocument(datasourceName=datasourceName, version=version, **filters)
        return [] if doc is None else [doc]

    def getDataSourceData(self, datasourceName=None, version=None, **filters):
        """
            Return the data payload of a toolkit data source.

            Args:
                datasourceName (Optional[str]): Name of the datasource.
                version (Optional[tuple]): Specific version to retrieve.
                **filters: Additional query filters.

            Returns:
                Any: The loaded datasource data, or None if not found.
        """
        filters[TOOLKIT_TOOLKITNAME_FIELD] = self.toolkitName
        doc = self.getDataSourceDocument(datasourceName=datasourceName, version=version, **filters)
        return None if doc is None else doc.getData()

    def addDataSource(
        self,
        dataSourceName: str,
        resource: str,
        dataFormat: str,
        version=(0, 0, 1),
        overwrite: bool = False,
        **kwargs,
    ):
        """
        Adds a resource to the toolkit.
        The type is always TOOLKIT_DATASOURCE_TYPE.
        The toolkit name is added automatically to the description.
        """
        kwargs[TOOLKIT_TOOLKITNAME_FIELD] = self.toolkitName
        kwargs[TOOLKIT_DATASOURCE_NAME] = dataSourceName
        kwargs[TOOLKIT_DATASOURCE_VERSION] = version

        if (self.getDataSourceDocument(dataSourceName, version=version) is None) or overwrite:
            if self.getDataSourceDocument(dataSourceName, version=version) is not None:
                # Exists -> delete and re-add
                delargs = {
                    TOOLKIT_DATASOURCE_NAME: dataSourceName,
                    TOOLKIT_DATASOURCE_VERSION: version,
                }
                self.deleteDataSource(**delargs)

            doc = self.addMeasurementsDocument(
                type=TOOLKIT_DATASOURCE_TYPE,
                resource=resource,
                dataFormat=dataFormat,
                desc=kwargs,
            )
        else:
            raise ValueError(
                f"Record {dataSourceName} (version {version}) already exists in project {self.projectName}. "
                f"use overwrite=True to overwrite the existing document"
            )

        return doc

    def deleteDataSource(self, datasourceName, version=None, **filters):
        """Delete a data source document."""
        doc = self.getDataSourceDocument(datasourceName=datasourceName, version=version, **filters)
        if doc is not None:
            doc.delete()
        return doc

    def setDataSourceDefaultVersion(self, datasourceName: str, version: tuple):
        """Set the default version for a given data source."""
        if len(
            self.getMeasurementsDocuments(
                type=TOOLKIT_DATASOURCE_TYPE,
                **{"datasourceName": datasourceName, "version": version},
            )
        ) == 0:
            raise ValueError(f"No DataSource with name={datasourceName} and version={version}.")

        self.setConfig(**{f"{datasourceName}_defaultVersion": version})
        print(f"{version} for dataSource {datasourceName} is now set to default.")


# ======================================================================
#  ToolkitHome
# ======================================================================

class ToolkitHome(abstractToolkit):
    """
    Central registry for available toolkits (static + dynamic).

    Responsibilities:
      - getToolkit(toolkitName, ...): locate & instantiate a toolkit class.
      - getToolkitTable(projectName): table of all toolkits (static + DB).
      - registerToolkit(...): register a toolkit class as a ToolkitDataSource
        using the abstractToolkit data source interface.
    """

    # Save modes (kept for compatibility)
    TOOLKIT_SAVEMODE_NOSAVE = TOOLKIT_SAVEMODE_NOSAVE
    TOOLKIT_SAVEMODE_ONLYFILE = TOOLKIT_SAVEMODE_ONLYFILE
    TOOLKIT_SAVEMODE_ONLYFILE_REPLACE = TOOLKIT_SAVEMODE_ONLYFILE_REPLACE
    TOOLKIT_SAVEMODE_FILEANDDB = TOOLKIT_SAVEMODE_FILEANDDB
    TOOLKIT_SAVEMODE_FILEANDDB_REPLACE = TOOLKIT_SAVEMODE_FILEANDDB_REPLACE

    # Static toolkit identifiers
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

    _toolkits: Dict[str, Dict[str, Any]] = None

    def __init__(self, projectName: Optional[str] = None, filesDirectory: Optional[str] = None):
        """
        ToolkitHome itself הוא Toolkit (abstractToolkit):
        - projectName יטען אוטומטית מ־caseConfiguration אם None.
        - כל ה־ToolkitDataSource הדינמיים נרשמים תחת toolkitName="ToolkitHome".
        """
        super().__init__(toolkitName="ToolkitHome", projectName=projectName, filesDirectory=filesDirectory)

        # Static built-in toolkits (internal source)
        self._toolkits = dict(
            GIS_Buildings=dict(
                cls="hera.measurements.GIS.vector.buildings.toolkit.BuildingsToolkit",
                desc=None,
                type="measurements",
            ),
            GIS_Tiles=dict(
                cls="hera.measurements.GIS.raster.tiles.TilesToolkit",
                desc=None,
                type="measurements",
            ),
            GIS_Vector_Topography=dict(
                cls="hera.measurements.GIS.vector.topography.TopographyToolkit",
                desc=None,
                type="measurements",
            ),
            GIS_Raster_Topography=dict(
                cls="hera.measurements.GIS.raster.topography.TopographyToolkit",
                desc=None,
                type="measurements",
            ),
            GIS_Demography=dict(
                cls="hera.measurements.GIS.vector.demography.DemographyToolkit",
                desc=None,
                type="measurements",
            ),
            GIS_LandCover=dict(
                cls="hera.measurements.GIS.raster.landcover.LandCoverToolkit",
                desc=None,
                type="measurements",
            ),
            RiskAssessment=dict(
                cls="hera.riskassessment.riskToolkit.RiskToolkit",
                desc=None,
                type="riskassessment",
            ),
            LSM=dict(
                cls="hera.simulations.LSM.toolkit.LSMToolkit",
                desc=None,
                type="simulations",
            ),
            OF_LSM=dict(
                cls="hera.simulations.openFoam.LSM.toolkit.OFLSMToolkit",
                desc=None,
                type="simulations",
            ),
            MeteoHighFreq=dict(
                cls="hera.measurements.meteorology.highfreqdata.toolkit.HighFreqToolKit",
                desc=None,
                type="measurements",
            ),
            MeteoLowFreq=dict(
                cls="hera.measurements.meteorology.lowfreqdata.toolkit.lowFreqToolKit",
                desc=None,
                type="measurements",
            ),
            hermesWorkflows=dict(
                cls="hera.simulations.hermesWorkflowToolkit.hermesWorkflowToolkit",
                desc=None,
                type="simulations",
            ),
            OpenFOAM=dict(
                cls="hera.simulations.openFoam.toolkit.OFToolkit",
                desc=None,
                type="simulations",
            ),
            WindProfile=dict(
                cls="hera.simulations.windProfile.toolkit.WindProfileToolkit",
                desc=None,
                type="simulations",
            ),
            GaussianDispersion=dict(
                cls="hera.simulations.gaussian.toolkit.gaussianToolkit",
                desc=None,
                type="simulations",
            ),
            machine_deep_learning=dict(
                cls="hera.simulations.machineLearningDeepLearning.toolkit.machineLearningDeepLearningToolkit",
                desc=None,
                type="simulations",
            ),
            experiment=dict(
                cls="hera.measurements.experiment.experiment.experimentHome",
                desc=None,
                type="measurements",
            ),
        )

        # Optional: keep a handle to the experiment toolkit (if available)
        self.experimentTK = None
        try:
            self.experimentTK = self.getToolkit(self.EXPERIMENT)
        except Exception:
            self.experimentTK = None

    # ------------------------------------------------------------------
    # Internal helper for repository config (uses generic dataToolkit)
    # ------------------------------------------------------------------

    def _get_data_toolkit(self, projectName: str = None):
        """
        Helper that returns a dataToolkit instance.

        We import dataToolkit lazily here to avoid circular imports
        between hera.toolkit and hera.utils.data.toolkit.
        dataToolkit itself works on the DEFAULT project internally.
        """
        from hera.utils.data.toolkit import dataToolkit
        return dataToolkit()

    # ------------------------------------------------------------------
    # Main API: getToolkit
    # ------------------------------------------------------------------

    def getToolkit(self, toolkitName: str, filesDirectory: Optional[str] = None, **kwargs):
        """
        Locate and instantiate a toolkit by name.

        This is the main public entry point for accessing toolkits in Hera.
        The method resolves the requested toolkit using the following order:
          1. Static built-in toolkit registry.
          2. Dynamically registered ToolkitDataSource documents (DB-backed).
          3. Experiment toolkits exposed via the experiment toolkit.

        Args:
            toolkitName (str): Logical name of the toolkit to load.
            filesDirectory (Optional[str]): Optional directory for toolkit file outputs.
            **kwargs: Additional keyword arguments forwarded to the toolkit constructor.
                      Commonly includes `projectName`.

        Returns:
            abstractToolkit: An initialized toolkit instance.

        Raises:
            ImportError: If the toolkit class cannot be imported.
            ValueError: If the toolkit cannot be found in any registry.
        """
        from hera.utils.data.toolkit import dataToolkit

        projectName = kwargs.pop("projectName", None) or self.projectName

        # 1) Static built-in toolkits
        if toolkitName in (self._toolkits or {}):
            info = self._toolkits[toolkitName]
            cls_path = info.get("cls")
            toolkit_cls = pydoc.locate(cls_path)
            if toolkit_cls is None:
                raise ImportError(f"Cannot locate class: {cls_path}")

            # Instantiate the toolkit in the resolved project context
            return toolkit_cls(
                projectName=projectName,
                filesDirectory=filesDirectory,
                **kwargs,
            )

        # 2) Dynamic toolkits registered as ToolkitDataSource of ToolkitHome
        dtk = dataToolkit()
        doc = dtk.getDataSourceDocument(datasourceName=toolkitName)
        if doc is not None:
            # self.logger.info(f"Found dynamic toolkit {toolkitName}. Loading")

            toolkitPath_raw = doc.getData()

            # ------------------------------------------------------------
            # Robust path resolution (same spirit as getExperiment)
            # ------------------------------------------------------------
            candidates = []

            # 1) Raw value
            candidates.append(toolkitPath_raw)

            # 2) Absolute relative to CWD
            try:
                candidates.append(os.path.abspath(toolkitPath_raw))
            except Exception:
                pass

            # 3) Relative to repo root if provided
            repo_root = os.environ.get("HERA_REPO_ROOT")
            if repo_root:
                candidates.append(os.path.join(repo_root, toolkitPath_raw))

            toolkitPath = None
            for cand in candidates:
                try:
                    cand_abs = os.path.abspath(cand)

                    # Valid if directory exists
                    if os.path.isdir(cand_abs):
                        toolkitPath = cand_abs
                        break

                    # Or if a single-file toolkit exists: <toolkitName>.py
                    py_file = os.path.join(cand_abs, f"{toolkitName}.py")
                    if os.path.isfile(py_file):
                        toolkitPath = cand_abs
                        break
                except Exception:
                    continue

            # Best-effort fallback (keeps error messages informative)
            if toolkitPath is None:
                toolkitPath = os.path.abspath(toolkitPath_raw)

            # ------------------------------------------------------------
            # Add toolkit path to sys.path (highest priority)
            # ------------------------------------------------------------
            if toolkitPath in sys.path:
                try:
                    sys.path.remove(toolkitPath)
                except ValueError:
                    pass
            sys.path.insert(0, toolkitPath)

            # self.logger.debug(f"Toolkit path (raw): {toolkitPath_raw}")
            # self.logger.debug(f"Toolkit path (resolved): {toolkitPath}")
            # self.logger.debug(f"Adding path {toolkitPath} to sys.path (priority)")

            # ------------------------------------------------------------
            # Import convention: <toolkitName>.<toolkitName>
            # ------------------------------------------------------------
            import_path = f"{toolkitName}.{toolkitName}"
            # self.logger.debug(f"Loading toolkit class: {import_path}")

            toolkit_cls = pydoc.locate(import_path)
            if toolkit_cls is None:
                raise ImportError(
                    f"Cannot find toolkit class {import_path} in {toolkitPath}"
                )

            return toolkit_cls(
                projectName=projectName,
                filesDirectory=filesDirectory,
                **kwargs,
            )
        # 3) Experiment toolkits fallback (experimentHome)
        if self.experimentTK is not None:
            try:
                return self.experimentTK.getExperiment(
                    experimentName=toolkitName,
                    filesDirectory=filesDirectory,
                )
            except Exception:
                # experimentHome does not recognize this experiment name
                pass

        # Nothing found in any registry
        raise ValueError(
            f"Toolkit '{toolkitName}' not found in static registry, ToolkitDataSource, "
            f"or experiment toolkit in project '{projectName}'."
        )

    # ------------------------------------------------------------------
    # Auto-register + get (kept for compatibility – but uses datasource API)
    # ------------------------------------------------------------------

    def auto_register_and_get(
        self,
        toolkitName: str,
        repositoryJSON: dict = None,
        repositoryName: Optional[str] = None,
        params: Optional[dict] = None,
        version: tuple = (0, 0, 1),
    ):
        """
           Automatically register a toolkit (if missing) and return an instance.

           This helper attempts to locate a toolkit class using repository metadata
           or JSON configuration, registers it as a ToolkitDataSource, and then
           returns an initialized toolkit instance.

           Args:
               toolkitName (str): Name of the toolkit to load.
               repositoryJSON (Optional[dict]): Optional repository metadata.
               repositoryName (Optional[str]): Target repository name.
               params (Optional[dict]): Toolkit initialization parameters.
               version (tuple): Toolkit version.

           Returns:
               abstractToolkit: An initialized toolkit instance.

           Raises:
               ValueError: If registration information is incomplete.
               ImportError: If the toolkit class cannot be imported.
           """
        from importlib import import_module

        params = params or {}
        classpath_hint = None
        projectName = self.projectName

        # 1) Classpath hint in the repository JSON
        if repositoryJSON:
            try:
                tk_section = repositoryJSON.get(toolkitName, {})
                reg = tk_section.get("Registry", {})
                classpath_hint = reg.get("classpath") or reg.get("cls")
            except Exception:
                pass

        if not classpath_hint:
            raise ValueError(
                f"auto_register_and_get: no classpath hint found for toolkit '{toolkitName}'. "
                f"Provide a 'Registry.classpath'/'cls' in repository JSON or seed a Toolkit document."
            )

        # Import the class
        mod_name, _, cls_name = classpath_hint.rpartition(".")
        if not mod_name or not cls_name:
            raise ValueError(f"Invalid classpath hint: '{classpath_hint}'")
        try:
            mod = import_module(mod_name)
            toolkit_cls = getattr(mod, cls_name)
        except Exception as exc:
            raise ImportError(f"Failed to import '{classpath_hint}' for toolkit '{toolkitName}'") from exc

        # Decide target repository for registration
        repo_to_use = (repositoryName or "").strip()
        if not repo_to_use:
            if projectName is None:
                raise ValueError(
                    "auto_register_and_get: projectName is None and no repositoryName provided. "
                    "Cannot resolve default repository."
                )
            repo_to_use = self.getDefaultRepository(projectName=projectName)
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
            repositoryName=repo_to_use,
        )

        # Return an instance
        return self.getToolkit(toolkitName=toolkitName)

    # ------------------------------------------------------------------
    # Listing toolkits (static + dynamic)
    # ------------------------------------------------------------------

    from typing import Optional, List, Dict

    def getToolkitDocuments(self, name: Optional[str] = None) -> List[Dict]:
        """
        Single source of truth for listing toolkits.

        This method returns a uniform list of "document-like" dictionaries coming from:
          1) The static in-memory registry (self._toolkits).
          2) Dynamic DB-backed toolkit records (type='ToolkitDataSource').
          3) Experiments that are exposed via the 'experiment' toolkit
             (so that experiments appear as toolkits in the same view).

        Each returned dict has the general shape:
            {
                "toolkit": <name>,
                "desc": {
                    "classpath": <string>,
                    "type": <string>,
                    "source": <string>,
                    "repositoryName": <string>,
                    "version": <any>
                }
            }
        """
        from hera.utils.data.toolkit import dataToolkit

        docs: List[Dict] = []

        # ------------------------------------------------------------------
        # 1) Static toolkits from the in-memory registry
        # ------------------------------------------------------------------
        for tk_name, info in (self._toolkits or {}).items():
            if name and tk_name != name:
                continue

            docs.append(
                {
                    "toolkit": tk_name,
                    "desc": {
                        # Fully-qualified classpath of the toolkit implementation
                        "classpath": info.get("cls", ""),
                        # Logical type of the toolkit (e.g. 'measurements', 'simulations', ...)
                        "type": info.get("type", "measurements"),
                        # Static entries are considered 'internal'
                        "source": "internal",
                        # Static toolkits do not come from a specific repository
                        "repositoryName": "",
                        # No explicit version for static entries
                        "version": "",
                    },
                }
            )
        # ------------------------------------------------------------------
        # 2) Dynamic toolkits stored in the DB as ToolkitDataSource documents
        # ------------------------------------------------------------------
        dtk = dataToolkit()
        try:
            # The dataToolkit is used as a generic interface to measurements
            # and to the underlying MongoDB-backed documents.
            dyn_docs = dtk.getMeasurementsDocuments(type="ToolkitDataSource") or []
            for d in dyn_docs:
                try:
                    desc = getattr(d, "desc", {}) or {}
                    tk_name = desc.get("datasourceName") or getattr(d, "datasourceName", None)
                    if not tk_name:
                        continue
                    if name and tk_name != name:
                        continue
                    docs.append(
                        {
                            "toolkit": tk_name,
                            "desc": {
                                # Dynamic entries may carry a classpath for dynamic import
                                "classpath": getattr(d, "resource", ""),
                                # Toolkit logical type; default to 'measurements' if missing
                                "type": desc.get("type", "") or "measurements",
                                # DB-backed entries are marked as coming from the DB
                                "source": desc.get("source", "") or "db",
                                # Repository is taken from the desc or from the document itself
                                "repositoryName": desc.get("repository", "") or getattr(d, "repository", ""),
                                # Version may be saved as a list or any other structure
                                "version": tuple(desc.get("version", ())) or getattr(d, "version", ""),
                            },
                        }
                    )
                except Exception:
                    # Be forgiving in case some records are partially formed
                    continue
        except Exception:
            # If the project or DB is not available, we still return the static list.
            pass

        # ------------------------------------------------------------------
        # 3) Experiments exposed as toolkits (via the 'experiment' toolkit)
        # ------------------------------------------------------------------
        docs.extend(self.getExperimentToolkitDocuments(name=name))

        return docs

    def getExperimentToolkitDocuments(self, name: Optional[str] = None) -> List[Dict]:
        """
        Return experiment definitions as "toolkit-like" documents.

        The 'experiment' toolkit (experimentHome) exposes all experiments
        of the project via getExperimentsMap(). This helper converts them
        into the same normalized "document-like" shape used by
        getToolkitDocuments(), so that experiments appear in the unified
        toolkits table and can be discovered via the same CLI.

        Notes:
        - Experiments do not have a direct classpath here; they are
          instantiated via experimentHome.getExperiment(...), so the
          'classpath' field is left empty.
        - The 'type' is marked as 'experiment'.
        - The 'source' is marked as 'experiment' to distinguish them from
          static or DB-backed toolkits.
        """
        # If the experiment toolkit is not available, return an empty list
        if self.experimentTK is None:
            return []

        try:
            # experimentHome.getExperimentsMap() returns a dictionary where:
            #   keys   = experiment names
            #   values = experiment metadata / configuration
            exp_map = self.experimentTK.getExperimentsMap()
        except Exception:
            # If anything goes wrong while querying experiments, do not
            # break the unified toolkit listing.
            return []

        docs: List[Dict] = []
        for exp_name in exp_map.keys():
            if name and exp_name != name:
                continue

            docs.append(
                {
                    "toolkit": exp_name,
                    "desc": {
                        # Experiments are not imported via a direct classpath
                        "classpath": "",
                        # Logical type to mark this as an experiment
                        "type": "experiment",
                        # Source tag to differentiate experiments from static/DB toolkits
                        "source": "experiment",
                        # Experiments are not associated with a repository name here
                        "repositoryName": "",
                        # No explicit version is tracked at this layer
                        "version": "",
                    },
                }
            )

        return docs


    def getToolkitTable(self, projectName: Optional[str]):
        """
        Build a DataFrame from getToolkitDocuments(...).
        """
        docs = self.getToolkitDocuments(name=None, projectName=projectName) or []
        rows = []
        for d in docs:
            desc = d.get("desc", {})
            rows.append(
                {
                    "toolkit": d.get("toolkit", ""),
                    "cls": desc.get("classpath", ""),
                    "source": desc.get("source", ""),
                    "type": desc.get("type", ""),
                    "repositoryName": desc.get("repositoryName", ""),
                    "version": desc.get("version", ""),
                }
            )

        if not rows:
            return pd.DataFrame(
                columns=["toolkit", "cls", "source", "type", "repositoryName", "version"]
            )

        df = pd.DataFrame(rows).drop_duplicates(subset=["toolkit", "source"], keep="first")
        return df

    # ------------------------------------------------------------------
    # Registration helpers (default repository config)
    # ------------------------------------------------------------------

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
            dfmt = getattr(_dt, "JSON", None) or getattr(_dt, "json", None) or getattr(
                _dt, "TEXT", None
            )
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

    # ------------------------------------------------------------------
    # Registration of toolkits using datasource interface
    # ------------------------------------------------------------------

    def registerToolkit(
            self,
            toolkit_name: str,
            toolkit_path: str,
            params: Optional[dict] = None,
            version = (0, 0, 1),
            overwrite: bool = False,
            **kwargs
    ):
        from hera.utils.data.toolkit import dataToolkit

        if not toolkit_name:
            raise ValueError("toolkit_name must be provided")

        if not toolkit_path:
            raise ValueError("toolkit_path must be provided")

        # If a file is given, store its parent directory
        toolkit_path = os.path.abspath(toolkit_path)
        if os.path.isfile(toolkit_path):
            toolkit_path = os.path.dirname(toolkit_path)

        params = params or {}
        dtk = dataToolkit()
        dtk._allowWritingToDefaultProject = True
        return dtk.addDataSource(
            dataSourceName=toolkit_name,
            resource = toolkit_path,  # directory only
            dataFormat = "string",
            version=version,
            overwrite = overwrite,
            params = params
        )
    # ------------------------------------------------------------------
    # JSON import helpers (unchanged, still valid עם הממשק החדש)
    # ------------------------------------------------------------------

    def import_toolkits_from_json(
        self,
        *,
        projectName: str,
        json_path: str,
        default_repository: str = None,
        overwrite: bool = True,
    ) -> list:
        """
        Read a JSON file and register all Toolkits it declares into the given project.
        (Uses registerToolkit -> datasource interface.)
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
        Load experiments from JSON into the given project as measurements documents.
        (לוגיקה קיימת – לא קשורה ישירות ל-datasource של toolkits.)
        """
        import json
        from hera.datalayer import Project

        if not projectName:
            raise ValueError("import_experiments_from_json: projectName is required")

        with open(json_path, "r") as f:
            payload = json.load(f) or {}

        experiments = payload.get("experiments") or []
        if not isinstance(experiments, list):
            raise ValueError("'experiments' must be a list in the JSON file")

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
                    continue

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
