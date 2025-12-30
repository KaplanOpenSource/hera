from hera.datalayer import Project
from hera.datalayer.datahandler import datatypes  # for datatypes.CLASS

import os
import inspect
import pydoc
import pandas as pd
from typing import Optional, List, Dict, Any

from hera.utils.logging import get_classMethod_logger

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
    Base class for all Hera toolkits.

    Notes
    -----
    - Inherits from :class:`hera.datalayer.Project`, so toolkits can access the full datalayer API.
    - Holds the logical toolkit name (toolkitName) and the project context (projectName).
    - Provides a unified "data sources" mechanism backed by measurement documents of type
      ``TOOLKIT_DATASOURCE_TYPE`` (i.e., ``"ToolkitDataSource"``).
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
        Initialize a new toolkit instance.

        This base class binds a toolkit to a Hera project context by inheriting
        from `Project`, and stores the toolkit logical name for later queries
        (e.g., when listing/adding ToolkitDataSource documents).

        Args:
            toolkitName (str): Logical name of the toolkit.
            projectName (Optional[str]): Project that the toolkit operates in.
                If None, the Project default/auto project mechanism is used.
            filesDirectory (Optional[str]): Optional directory for storing files
                created by this toolkit (datasource resources, outputs, etc.).
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
        """
        List datasource names registered for this toolkit.

        This queries measurement documents of type `TOOLKIT_DATASOURCE_TYPE`
        that are associated with this toolkit name.

        Args:
            **filters: Additional filters forwarded to `getMeasurementsDocuments(...)`
                (e.g., datasourceName, version, repository, etc.).

        Returns:
            List[str]: Datasource names (desc["datasourceName"]) for matching documents.
        """
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
        Return metadata maps for datasource documents of this toolkit.

        Each returned dict is a merged view of:
          - document fields (dataFormat, resource)
          - document description (`doc.desc`) containing datasource metadata.

        Args:
            **filters: Additional query filters forwarded to `getMeasurementsDocuments(...)`.

        Returns:
            List[Dict[str, Any]]: One dict per datasource document.
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
        """
        Build a pandas table of datasource metadata for this toolkit.

        This converts `getDataSourceMap()` into a normalized DataFrame.

        Args:
            **filters: Additional filters forwarded to `getDataSourceMap()`.

        Returns:
            pd.DataFrame: A table of datasource metadata (empty if no datasources found).
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
        Return all datasource documents associated with this toolkit.

        This is a thin convenience wrapper over `getMeasurementsDocuments(...)`
        that injects:
          - type=TOOLKIT_DATASOURCE_TYPE
          - toolkit=<this toolkit name>

        Args:
            **kwargs: Additional query fields for `getMeasurementsDocuments(...)`.

        Returns:
            List[Any]: List of MeasurementDocument objects (or document-like objects)
            returned by the datalayer.
        """
        queryDict = {
            "type": TOOLKIT_DATASOURCE_TYPE,
            TOOLKIT_TOOLKITNAME_FIELD: self.toolkitName,
        }
        queryDict.update(**kwargs)
        return self.getMeasurementsDocuments(**queryDict)

    def getDataSourceDocument(self, datasourceName: Optional[str], version=None, **filters):
        """
        Resolve a single datasource document for this toolkit.

        Resolution rules
        ----------------
        - If ``version`` is provided: return the matching document (if exists).
        - If ``version`` is not provided:
            1) If a default version is configured in project config under
               ``<datasourceName>_defaultVersion``, prefer it.
            2) Otherwise, if multiple versions exist, return the latest
               (by ``max(version)`` comparison).
        - If no matching documents exist: return ``None``.

        Parameters
        ----------
        datasourceName : Optional[str]
            Logical datasource name. If None, only ``filters`` are used.
        version : Optional[Any]
            Specific version to fetch (commonly a tuple like ``(0, 0, 1)``).
        **filters
            Additional filters forwarded to ``getMeasurementsDocuments(...)``.

        Returns
        -------
        Optional[Any]
            The resolved datasource document, or None if not found.
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

    def getDataSourceDocuments(self, datasourceName, version=None, **filters):
        """
        Return datasource documents as a list (API symmetry helper).

        This mirrors other APIs that return lists (measurements/cache) even when
        only one document is expected.

        Args:
            datasourceName (str): Datasource name.
            version (Optional[Any]): Specific version to retrieve.
            **filters: Additional query filters.

        Returns:
            List[Any]: [] if not found, otherwise [document].
        """
        doc = self.getDataSourceDocument(datasourceName=datasourceName, version=version, **filters)
        return [] if doc is None else [doc]

    def getDataSourceData(self, datasourceName=None, version=None, **filters):
        """
        Load the data payload of a toolkit datasource.

        This resolves the datasource document using `getDataSourceDocument(...)`
        and returns `doc.getData()`.

        Args:
            datasourceName (Optional[str]): Name of the datasource.
            version (Optional[Any]): Specific version to retrieve.
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
        Register a new datasource document under this toolkit.

        This creates a measurement document with:
          - type = TOOLKIT_DATASOURCE_TYPE
          - desc.toolkit = <this toolkit name>
          - desc.datasourceName = dataSourceName
          - desc.version = version
          - resource/dataFormat as provided

        If the datasource already exists and `overwrite=False`, raises ValueError.
        If `overwrite=True`, the old document is deleted and a new one is created.

        Args:
            dataSourceName (str): Logical datasource name.
            resource (str): Resource path (directory or file) associated with the datasource.
            dataFormat (str): Data format identifier (e.g., datatypes.CLASS).
            version (tuple): Datasource version identifier.
            overwrite (bool): Whether to overwrite an existing datasource with same name/version.
            **kwargs: Additional fields added into the document desc.

        Returns:
            Any: The created MeasurementDocument.

        Raises:
            ValueError: If the datasource already exists and overwrite=False.
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
        """
        Delete a datasource document from this toolkit.

        Args:
            datasourceName (str): Datasource name to delete.
            version (Optional[Any]): Version to delete. If None, the resolved document
                is determined by `getDataSourceDocument(...)` rules.
            **filters: Additional filters for locating the datasource document.

        Returns:
            Optional[Any]: The deleted document (or the resolved document), or None if not found.
        """
        doc = self.getDataSourceDocument(datasourceName=datasourceName, version=version, **filters)
        if doc is not None:
            doc.delete()
        return doc

    def setDataSourceDefaultVersion(self, datasourceName: str, version: tuple):
        """
        Set the default version for a datasource in the project configuration.

        This writes a config entry:
            <datasourceName>_defaultVersion = version

        The default is later used by `getDataSourceDocument(...)` when `version`
        is not explicitly provided.

        Args:
            datasourceName (str): Datasource name.
            version (tuple): Version tuple to set as default.

        Raises:
            ValueError: If no datasource document exists with the given name/version.
        """
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
        Initialize ToolkitHome (the central toolkit registry).

        ToolkitHome is itself a toolkit (inherits from abstractToolkit) and acts as:
          - A static registry of built-in toolkits (in-memory mapping).
          - A dynamic registry backed by measurement documents (ToolkitDataSource).
          - A gateway to experiments via the experiment toolkit (if available).

        Args:
            projectName (Optional[str]): Project context for registry operations.
                If None, Project default/auto mechanism is used.
            filesDirectory (Optional[str]): Optional directory for toolkit file outputs.
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
        Create a generic dataToolkit helper for repository/measurement queries.

        Notes:
            Imported lazily to avoid circular imports between `hera.toolkit`
            and `hera.utils.data.toolkit`.

        Args:
            projectName (Optional[str]): Currently unused; kept for compatibility.

        Returns:
            dataToolkit: An instance of `hera.utils.data.toolkit.dataToolkit`.
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
        doc = self.getDataSourceDocument(datasourceName=toolkitName)
        if doc is not None:
            tk = doc.getData()

            # Best-effort: enforce filesDirectory if the toolkit supports it
            if hasattr(tk, "setFilesDirectory") and filesDirectory is not None:
                tk.setFilesDirectory(filesDirectory)

            # Note: dynamic toolkits are returned as-is (they may already be bound to a project)
            return tk

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
        Auto-register a missing toolkit and return an initialized instance.

        This helper is used when a toolkit is not found in the static registry.
        It attempts to locate a classpath hint (typically from repository JSON),
        registers the toolkit as a dynamic ToolkitDataSource, and then calls
        `getToolkit(...)` to return an instance.

        Args:
            toolkitName (str): Name of the toolkit to load/register.
            repositoryJSON (Optional[dict]): Repository metadata that may contain
                a classpath hint under `toolkitName -> Registry -> classpath/cls`.
            repositoryName (Optional[str]): Target repository name. If not provided,
                the project's default repository is used.
            params (Optional[dict]): Toolkit initialization parameters to store in the registration.
            version (tuple): Version tuple to store in the registration.

        Returns:
            abstractToolkit: An initialized toolkit instance.

        Raises:
            ValueError: If classpath hint / repository resolution fails.
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

    def getToolkitDocuments(self, name: Optional[str] = None, projectName: Optional[str] = None) -> List[Dict]:
        """
        List all known toolkits (static, dynamic, and experiments) as normalized records.

        The returned records are "document-like" dictionaries that unify:
          1) Static in-memory registry (built-in toolkits).
          2) Dynamic DB-backed ToolkitDataSource documents.
          3) Experiments exposed via the experiment toolkit.

        Args:
            name (Optional[str]): If provided, only return records for this toolkit name.
            projectName (Optional[str]): If provided, include DB-backed toolkits
                from this project.

        Returns:
            List[Dict]: Normalized toolkit records of the form:
                {"toolkit": <name>, "desc": {...}}.
        """
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
        if projectName:
            try:
                # The dataToolkit is used as a generic interface to measurements
                # and to the underlying MongoDB-backed documents.
                dt = self._get_data_toolkit(projectName=projectName)
                dyn_docs = dt.getMeasurementsDocuments(type=TOOLKIT_DATASOURCE_TYPE) or []

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
                                    "classpath": desc.get("classpath", ""),
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
        List experiments as "toolkit-like" records for unified discovery.

        Experiments are provided by the experiment toolkit (experimentHome) via
        `getExperimentsMap()`. This method normalizes them into the same record
        shape returned by `getToolkitDocuments()`.

        Args:
            name (Optional[str]): If provided, only return this experiment.

        Returns:
            List[Dict]: Normalized experiment records (toolkit-like).
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
        Build a pandas table of toolkits for display/CLI usage.

        Args:
            projectName (Optional[str]): If provided, include DB-backed toolkits
                from this project.

        Returns:
            pd.DataFrame: Columns: toolkit, cls, source, type, repositoryName, version.
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
        Persist the default repository for a project.

        This stores a measurement document with:
          - type = "RepositoryConfig"
          - desc.defaultRepository = repositoryName

        Future operations can omit an explicit repository name and rely on this default.

        Args:
            projectName (str): Project name to store the default repository for.
            repositoryName (str): Default repository name.
            overwrite (bool): If True, remove previous RepositoryConfig documents first.

        Returns:
            Any: The created RepositoryConfig document.

        Raises:
            ValueError: If projectName or repositoryName is missing.
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
        Read the default repository for a project.

        Args:
            projectName (str): Project name.

        Returns:
            str: Default repository name, or "" if not configured.

        Raises:
            ValueError: If projectName is missing.
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
            toolkitclass,
            *,
            repositoryName: str,
            datasource_name: Optional[str] = None,
            params: Optional[dict] = None,
            version=(0, 0, 1),
            overwrite: bool = False,
    ):
        """
        Register a toolkit class as a dynamic ToolkitDataSource.

        This stores a measurement document under ToolkitHome with:
          - type = TOOLKIT_DATASOURCE_TYPE
          - resource = directory containing the toolkit class
          - dataFormat = datatypes.CLASS
          - desc fields:
                datasourceName, version, toolkit="ToolkitHome",
                repository, classpath, parameters, source="db"

        Args:
            toolkitclass: The toolkit class to register.
            repositoryName (str): Repository that owns the toolkit definition.
            datasource_name (Optional[str]): Logical name used to load the toolkit.
                Defaults to the class name.
            params (Optional[dict]): Initialization parameters saved in the registration.
            version (tuple): Toolkit version tuple.
            overwrite (bool): Whether to overwrite an existing registration.

        Returns:
            Any: The created/updated MeasurementDocument.

        Raises:
            ValueError: If repositoryName is missing.
        """
        if not repositoryName:
            raise ValueError("registerToolkit: 'repositoryName' is required")

        module_path = inspect.getfile(toolkitclass)
        resource_dir = os.path.dirname(os.path.abspath(module_path))
        classpath = f"{toolkitclass.__module__}.{toolkitclass.__qualname__}"

        ds_name = datasource_name or toolkitclass.__name__
        params = params or {}

        extra_desc = {
            "repository": repositoryName,
            "classpath": classpath,
            "parameters": params,
            "type": "ToolkitDataSource",
            "source": "db",
        }

        return self.addDataSource(
            dataSourceName=ds_name,
            resource=resource_dir,
            dataFormat=datatypes.CLASS,
            version=tuple(version),
            overwrite=overwrite,
            **extra_desc,
        )

    def import_toolkits_from_json(
        self,
        *,
        projectName: str,
        json_path: str,
        default_repository: str = None,
        overwrite: bool = True,
    ) -> list:
        """
        Import and register toolkits from a JSON declaration.

        The JSON is expected to contain:
          - optional "repository"
          - "toolkits": list of items with {name, classpath, version, parameters}

        For each entry:
          1) locate(classpath)
          2) registerToolkit(...)

        Args:
            projectName (str): Project name (used to resolve default repository if needed).
            json_path (str): Path to JSON file containing toolkit declarations.
            default_repository (Optional[str]): Fallback repository name if JSON does not contain one.
            overwrite (bool): Whether to overwrite existing registrations.

        Returns:
            list: Names of registered toolkits.

        Raises:
            ValueError: If JSON structure is invalid or repository cannot be resolved.
            ImportError: If a classpath cannot be located/imported.
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
        Import experiments from a JSON file into a project.

        This loads measurement documents described in the JSON under the given project.
        Relative resource paths (if marked) are resolved relative to the JSON file location.

        Args:
            projectName (str): Target project to load experiments into.
            json_path (str): Path to JSON file containing experiment declarations.

        Returns:
            list: Names of loaded experiments (unique).

        Raises:
            ValueError: If projectName is missing, or JSON has invalid structure.
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