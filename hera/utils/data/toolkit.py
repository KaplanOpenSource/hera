import json
import argparse
import pathlib
import os

from hera.utils import loadJSON, dictToMongoQuery
from hera.utils.logging import get_classMethod_logger
from hera.toolkit import abstractToolkit, ToolkitHome


class dataToolkit(abstractToolkit):
    """
    Toolkit for managing data repositories (replacing the old hera-data).

    It is initialized only with the DEFAULT project.

    The structure of a datasource file is:

        {
            "<toolkit name>": {
                "<datasource name>": {
                    "resource": "<location of datasource>",
                    "dataFormat": "<type of data source>",
                    "desc": {
                        ... metadata ...
                    }
                },
                ...
            },
            ...
        }
    """

    def __init__(self, connectionName=None):
        """
        Initialize the dataToolkit on the default project.

        Parameters
        ----------
        connectionName : str, optional
            The DB connection name. If None, uses the current OS username.
        """
        super().__init__(toolkitName="heradata", projectName=self.DEFAULTPROJECT, filesDirectory=None, connectionName=connectionName)

    def addRepository(self, repositoryName, repositoryPath, overwrite=False):
        """
        Register a repository JSON file as a data source.

        Parameters
        ----------
        repositoryName : str
            The name to register the repository under.
        repositoryPath : str
            Path to the repository JSON file. ``.json`` extension is appended if missing.
        overwrite : bool
            If True, overwrite an existing repository with the same name.
        """
        self._allowWritingToDefaultProject = True  # allows the addition of datasource to the Default project.

        repositoryPath = f"{repositoryPath}.json" if "json" not in repositoryPath else repositoryPath
        self.addDataSource(dataSourceName=repositoryName, resource=os.path.abspath(repositoryPath),
                           dataFormat=self.datatypes.JSON_DICT, overwrite=overwrite)
        self._allowWritingToDefaultProject = False

    def getRepositoryTable(self):
        """
        Return a DataFrame listing all registered repositories.

        Returns
        -------
        pandas.DataFrame
        """
        return self.getDataSourceTable()

    def getRepository(self, repositoryName):
        """
        Load and return a repository's JSON content by name.

        Parameters
        ----------
        repositoryName : str
            The name of the registered repository.

        Returns
        -------
        dict
            The parsed repository JSON.
        """
        logger = get_classMethod_logger(self, "getRepository")
        logger.info(f"Trying to find repository {repositoryName} in project {self.DEFAULTPROJECT}")
        repo = self.getDataSourceData(datasourceName=repositoryName)

        return loadJSON(repo)

    def loadAllDatasourcesInAllRepositoriesToProject(self, projectName, overwrite=False):
        """
        Load all data sources from all registered repositories into a project.

        Parameters
        ----------
        projectName : str
            The target project name.
        overwrite : bool
            If True, overwrite existing data sources.
        """
        logger = get_classMethod_logger(self, "loadAllDatasourcesInAllRepositoriesToProject")
        for repository in self.getDataSourceList():
            try:
                logger.info(f"Loading the repository {repository} to project {projectName}")
                self.loadAllDatasourcesInRepositoryToProject(projectName, repositoryName=repository,
                                                             overwrite=overwrite)
            except ValueError as e:
                logger.info(
                    f"Did not loaded repository: {repository}, since an error occured when tried to load it.\n The error message: {e}")

    def loadAllDatasourcesInRepositoryToProject(self, projectName, repositoryName, overwrite=False):
        """
        Load all data sources from a specific repository into a project.

        Parameters
        ----------
        projectName : str
            The target project name.
        repositoryName : str
            The name of the registered repository to load from.
        overwrite : bool
            If True, overwrite existing data sources.
        """
        logger = get_classMethod_logger(self, "loadAllDatasourcesInRepositoryToProject")
        logger.info(f"Loading repository {repositoryName}")
        repdoc = self.getDataSourceDocument(repositoryName)
        conf = repdoc.getData()
        logger.info(f"Data: {conf}")
        basedir = os.path.dirname(repdoc.resource)
        logger.info(f"basedir: {basedir}")
        logger.info(f"Loading the items in {repositoryName} repository to the {projectName}")
        self.loadAllDatasourcesInRepositoryJSONToProject(projectName=projectName,
                                                         repositoryJSON=conf,
                                                         basedir=basedir,
                                                         overwrite=overwrite)

    # hera/utils/data/toolkit.py  (inside class dataToolkit)
    # -----------------------------------------------------------------------------
    # Load all datasources from a repository JSON into a project.
    # If a toolkit is missing, try to auto-register it using classpath hints.
    # -----------------------------------------------------------------------------
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


    def loadAllDatasourcesInRepositoryJSONToProject(self,
                                                    projectName: str,
                                                    repositoryJSON: dict,
                                                    basedir: str = "",
                                                    overwrite: bool = False,
                                                    auto_register_missing: bool = True):
        """
        Iterate through the repository JSON and for each toolkit:
        - Try to get an instance via ToolkitHome.getToolkit.
        - If missing and auto_register_missing=True, attempt auto-register ONLY if there is
          a clear classpath hint in the JSON (Registry.classpath or Registry.cls).
        - After we have a valid instance, dispatch to the appropriate handler per section.
        """
        logger = get_classMethod_logger(self, "loadAllDatasourcesInRepositoryJSONToProject")
        if isinstance(repositoryJSON, str):
            if  repositoryJSON.startswith('/'): # if there is no data
                logger.info("skipping dynamic toolkit")
                return
            try:
                repositoryJSON = json.loads(repositoryJSON)
            except json.JSONDecodeError:
                logger.error("repositoryJSON is a string but not a valid JSON format.")
                return
        if not isinstance(repositoryJSON, dict):
            logger.warning(f"Expected dict for repositoryJSON, got {type(repositoryJSON)}. Skipping.")
            return
        if not repositoryJSON:
            logger.info("repositoryJSON is empty. Nothing to load.")
            return
        handlerDict = dict(
            Config=self._handle_Config,
            Datasource=self._handle_DataSource,
            Measurements=lambda toolkit, itemName, docTypeDict, overwrite, basedir: self._DocumentHandler(
                toolkit, itemName, docTypeDict, overwrite, "Measurements", basedir
            ),
            Simulations=lambda toolkit, itemName, docTypeDict, overwrite, basedir: self._DocumentHandler(
                toolkit, itemName, docTypeDict, overwrite, "Simulations", basedir
            ),
            Cache=lambda toolkit, itemName, itemDesc, overwrite, basedir: self._DocumentHandler(
                toolkit, itemName, itemDesc, overwrite, "Cache", basedir
            ),
            Function=self._handle_Function,
        )

        tk_home = ToolkitHome(projectName=projectName)

        for toolkitName, toolkitDict in (repositoryJSON or {}).items():
            # 1) Try static/dynamic resolution via ToolkitHome.getToolkit
            try:
                toolkit = tk_home.getToolkit(toolkitName=toolkitName)

            except Exception as e:
                logger.info(f"Toolkit '{toolkitName}' not found via getToolkit: {e}")
                toolkit = None



            # 3) If we still do not have a toolkit instance, skip this key quietly
            if toolkit is None:
                logger.info(
                    f"Skipping key '{toolkitName}' in repository JSON – "
                    f"no matching toolkit and no auto-registration performed."
                )
                continue

            # 4) Dispatch sections (Config, Datasource, Measurements, Simulations, Cache, Function)
            for key, docTypeDict in toolkitDict.items():
                logger.info(f"Loading document type {key} to toolkit {toolkitName}")
                handler = handlerDict.get(key.title(), None)

                if handler is None:
                    err = (
                        f"Unkonw Handler {key.title()}. "
                        f"The handler must be {', '.join(handlerDict.keys())}. "
                    )
                    logger.error(err)
                    raise ValueError(err)

                try:
                    handler(
                        toolkit=toolkit,
                        itemName=key,
                        docTypeDict=docTypeDict,
                        overwrite=overwrite,
                        basedir=basedir,
                    )
                except Exception as e:
                    err = (
                        f"The error {e} occured while adding *{key}* to toolkit {toolkitName}... skipping!!!"
                    )
                    logger.error(err)


    def _handle_Config(self, toolkit, itemName, docTypeDict, overwrite, basedir):
        """
        Handle a Config section from a repository JSON by calling ``toolkit.setConfig``.

        Parameters
        ----------
        toolkit : abstractToolkit
            The toolkit instance to configure.
        itemName : str
            The section name (unused, always 'Config').
        docTypeDict : dict
            Key-value pairs to set as configuration.
        overwrite : bool
            Whether to overwrite existing values.
        basedir : str
            Base directory for resolving relative paths (unused for Config).
        """
        toolkit.setConfig(**docTypeDict)

    def _handle_DataSource(self, toolkit, itemName, docTypeDict, overwrite, basedir):
        """
        Handle a DataSource section from a repository JSON by adding data sources to the toolkit.

        Parameters
        ----------
        toolkit : abstractToolkit
            The toolkit instance to add data sources to.
        itemName : str
            The section name.
        docTypeDict : dict
            Dictionary mapping data source names to their descriptions.
        overwrite : bool
            If True, overwrite existing data sources.
        basedir : str
            Base directory for resolving relative resource paths.
        """
        logger = get_classMethod_logger(self, "_handle_DataSource")

        for itemName, itemDesc in docTypeDict.items():
            theItem = itemDesc["item"]

            isRelativePath = itemDesc.get("isRelativePath")
            assert (isRelativePath=='True' or isRelativePath=='False') or isinstance(isRelativePath,bool), "isRelativePath must be defined as 'True' or 'False'. "


            if 'resource' in theItem and "resourceFilePath" in theItem:
                logger.warning(f"both resource and resourceFilePath are defined for datasource {itemName}, using just resource")
                theItem.pop("resourceFilePath")

            if 'resource' not in theItem and "resourceFilePath" in theItem:
                if isRelativePath=='True' or isRelativePath is True:
                    logger.debug(
                        f"The input is not absolute (it is relative). Adding the path {basedir} to the resource {theItem['resourceFilePath']}")
                    theItem["resourceFilePath"] = os.path.join(basedir, theItem["resourceFilePath"])
                
                logger.info("detected dataSource resource specified using file's contents")
                try:
                    with open(theItem.pop("resourceFilePath")) as dataSourceResourceFile:
                        theItem['resource'] = json.load(dataSourceResourceFile)
                        logger.info("extracted resource from file successfully")
                except Exception as e:
                    logger.error(f"failed reading resource from file, {e}")
            else:
                # logger.debug(f"Checking if {itemName} resource is a path {isRelativePath}, is it absolute? {isAbsolute}")
                if isRelativePath=='True' or isRelativePath is True:
                    logger.debug(
                        f"The input is not absolute (it is relative). Adding the path {basedir} to the resource {theItem['resource']}")
                    theItem["resource"] = os.path.join(basedir, theItem["resource"])




            logger.debug(f"Checking if the data item {itemName} is already in project {toolkit.projectName}")
            datasource = toolkit.getDataSourceDocuments(datasourceName=itemName)
            if len(datasource) == 0 or overwrite:

                if len(datasource) == 1:
                    logger.debug("Remove the old datasource")
                    toolkit.deleteDataSource(datasourceName=itemName)

                logger.debug("Adding a new datasource")
                theItem['dataSourceName'] = itemName
                theItem['overwrite'] = overwrite
                toolkit.addDataSource(**theItem)
                logger.info(f"Added source {itemName} to tool {toolkit.toolkitName} in project {toolkit.projectName}")
            else:
                logger.error(f"Source {itemName} already exists in {toolkit.projectName}. Use --overwrite to force update")

    def _DocumentHandler(self, toolkit, itemName, docTypeDict, overwrite, documentType, basedir):
        """
        Handle a Measurements, Simulations, or Cache section from a repository JSON.

        Parameters
        ----------
        toolkit : abstractToolkit
            The toolkit instance to add documents to.
        itemName : str
            The section name.
        docTypeDict : dict
            Dictionary mapping document names to their descriptions.
        overwrite : bool
            If True, overwrite existing documents.
        documentType : str
            One of 'Measurements', 'Simulations', or 'Cache'.
        basedir : str
            Base directory for resolving relative resource paths.
        """
        logger = get_classMethod_logger(self, "_handle_Document")
        logger.info(f"Loading {itemName} to toolkit {toolkit.toolkitName} (ProjectName {toolkit.projectName}")
        for itemName, itemDesc in docTypeDict.items():
            theItem = itemDesc["item"]
            theItem["resource"] = self._makeItemPathAbsolute(theItem,basedir)

            logger.debug(f"Checking if the data item {itemName} is already in the project")
            retrieveFuncName = f"get{documentType}Documents"
            retrieveFunc = getattr(toolkit, retrieveFuncName)
            if retrieveFunc is None:
                raise ValueError(
                    f"function {retrieveFuncName} not found. Key {documentType} must be : DataSource, Measurement, Cache, or Simulation")
            qrydict = dict(theItem)
            del qrydict['resource']
            del qrydict['dataFormat']
            itemQry = dictToMongoQuery(qrydict)
            datasource = retrieveFunc(**itemQry)
            logger.debug(f"Found {len(datasource)} documents")

            if len(datasource) == 0:
                funcName = f"add{documentType}Document"

                logger.debug(f"Adding the document of type {documentType} using the function {funcName}")
                func = getattr(toolkit, funcName)

                func(**theItem)
                logger.info(f"Added source {itemName} to tool {toolkit.toolkitName} in project {toolkit.projectName}")

            elif overwrite:
                logger.debug("Updating an existing document")
                dataitem = datasource[0]
                dataitem['resource'] = theItem["resource"]
                dataitem['dataFormat'] = theItem['dataFormat']
                curDesc = theItem.get("desc", {})
                curDesc.update(dataitem['desc'])
                dataitem['desc'] = curDesc
                dataitem.save()
                logger.info(f"Updated source {itemName} in tool {toolkit.toolkitName} in project {toolkit.projectName}")
            else:
                logger.error(
                    f"Source {itemName} already exists in {toolkit.projectName}. Use --overwrite to force update")

    def _handle_Function(self, toolkit, itemName, docTypeDict, overwrite, basedir):
        """
        Handle a Function section by calling named methods on the toolkit.

        Each key in ``docTypeDict`` is a method name on ``self``. The value can be:
        - A dict: passed as keyword arguments to a single call.
        - A list of dicts: each dict triggers a separate call.

        The called method must accept an ``overwrite`` keyword argument.

        Parameters
        ----------
        toolkit : abstractToolkit
            The toolkit instance (unused directly; methods are called on ``self``).
        itemName : str
            The section name.
        docTypeDict : dict
            Maps method names to their argument(s).
        overwrite : bool
            Passed to each method call.
        basedir : str
            Base directory (unused for Function).
        """
        logger = get_classMethod_logger(self, "_handle_GeneralFunction")
        for itemName, itemDesc in docTypeDict.items():
            retrieveFunc = getattr(self,itemName)

            if isinstance(itemDesc,dict):
                retrieveFunc(**itemDesc,overwrite=overwrite)
            elif isinstance(itemDesc,list):
                for imt in itemDesc:
                    if isinstance(imt,dict):
                        retrieveFunc(**imt, overwrite=overwrite)
                    else:
                        err = f"{itemName} has a non dict item in the list : {imt}... ignoring."
                        logger.error(err)
            else:
                err = f"{itemName} value must be dict of a list of dicts. "
                logger.error(err)
                raise ValueError(err)


    def _makeItemPathAbsolute(self, theItem, basedir):
        """
        Convert a resource path to absolute if the ``isRelativePath`` flag is set.

        Parameters
        ----------
        theItem : dict
            The item data containing ``resource`` and optionally ``isRelativePath``.
        basedir : str
            Base directory to resolve relative paths against.

        Returns
        -------
        str
            The absolute resource path.
        """
        logger = get_classMethod_logger(self, "_makeItemPathAbsolute")
        isRelativePath = bool(theItem.get("isRelativePath", True))
        # logger.debug(f"Checking if {itemName} resource is a path {isRelativePath}, is it absolute? {isAbsolute}")

        if isRelativePath:
            logger.debug(
                f"The input is not absolute (it is relative). Adding the path {basedir} to the resource {theItem['resource']}")

        return os.path.join(basedir, theItem["resource"]) if isRelativePath else theItem["resource"]

    # -------------------------------------------------------------------------
    # Direct-load helpers (no MongoDB round-trip required)
    # -------------------------------------------------------------------------

    @staticmethod
    def resolveDataSourcePaths(repositoryJSON, basedir=""):
        """
        Walk a repository JSON dict and resolve every ``resource`` field to an
        absolute path, respecting the ``isRelativePath`` flag on each entry.

        Parameters
        ----------
        repositoryJSON : dict
            The parsed repository JSON (toolkit-name -> section dict).
        basedir : str
            The base directory against which relative paths are resolved.
            Typically the directory that contains the repository JSON file.

        Returns
        -------
        dict
            A *deep copy* of ``repositoryJSON`` with all ``resource`` fields
            converted to absolute paths.
        """
        import copy
        resolved = copy.deepcopy(repositoryJSON)

        for _toolkitName, toolkitDict in resolved.items():
            if not isinstance(toolkitDict, dict):
                continue
            for sectionKey, sectionDict in toolkitDict.items():
                if not isinstance(sectionDict, dict):
                    continue
                for itemName, itemDesc in sectionDict.items():
                    if not isinstance(itemDesc, dict):
                        continue
                    # Handle entries that have an "item" wrapper
                    item = itemDesc.get("item", itemDesc)
                    if "resource" not in item:
                        continue
                    is_rel = itemDesc.get("isRelativePath", item.get("isRelativePath"))
                    if is_rel == "True" or is_rel is True:
                        item["resource"] = os.path.abspath(
                            os.path.join(basedir, item["resource"])
                        )
        return resolved

    @staticmethod
    def loadRepositoryFromPath(json_path):
        """
        Read a repository JSON file directly from disk, resolve all relative
        ``resource`` paths to absolute paths based on the JSON file's directory,
        and return the resulting dict.

        This allows tests (and lightweight scripts) to work with repository
        data without going through ``addRepository`` + MongoDB storage.

        Parameters
        ----------
        json_path : str
            Path to the repository JSON file.

        Returns
        -------
        dict
            The repository dict with all resource paths resolved to absolute.

        Raises
        ------
        FileNotFoundError
            If *json_path* does not exist.
        """
        json_path = os.path.abspath(json_path)
        if not os.path.isfile(json_path):
            raise FileNotFoundError(f"Repository JSON not found: {json_path}")

        with open(json_path, "r", encoding="utf-8") as fh:
            repo_json = json.load(fh)

        basedir = os.path.dirname(json_path)
        return dataToolkit.resolveDataSourcePaths(repo_json, basedir=basedir)
