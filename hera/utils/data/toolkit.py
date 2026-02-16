import json
import argparse
from hera import toolkitHome,toolkit
from hera.utils import loadJSON, dictToMongoQuery
from hera.utils.logging import get_classMethod_logger
import pathlib
import os


class dataToolkit(toolkit.abstractToolkit):
    """
        A toolkit to handle the data (replacing the function of hera-data).

        It is initialized only with the DEFAULT project.

        The structure of a datasource file is:
        {
            <toolkit name>: {
              <datasource name>: {
                "resource": <location of datasource>,
                "dataFormat": <type of data source>,
                "desc": {
                  < any other parameters/ metadata descriptions of the datasource>
                },
                .
                .

            },
            .
            .

       }
    """

    def __init__(self,connectionName=None):
        super().__init__(toolkitName="heradata", projectName=self.DEFAULTPROJECT, filesDirectory=None,connectionName=connectionName)

    def addRepository(self, repositoryName, repositoryPath, overwrite=False):
        """
            A path to the repository
        Parameters
        ----------
        repositoryName
        repositoryPath

        Returns
        -------

        """
        self._allowWritingToDefaultProject = True  # allows the addition of datasource to the Default project.

        repositoryPath = f"{repositoryPath}.json" if "json" not in repositoryPath else repositoryPath
        self.addDataSource(dataSourceName=repositoryName, resource=os.path.abspath(repositoryPath),
                           dataFormat=self.datatypes.JSON_DICT, overwrite=overwrite)
        self._allowWritingToDefaultProject = False

    def getRepositoryTable(self):
        return self.getDataSourceTable()

    def getRepository(self, repositoryName):
        logger = get_classMethod_logger(self, "getRepository")
        logger.info(f"Trying to find repository {repositoryName} in project {self.DEFAULTPROJECT}")
        repo = self.getDataSourceData(datasourceName=repositoryName)

        return loadJSON(repo)

    def loadAllDatasourcesInAllRepositoriesToProject(self, projectName, overwrite=False):
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
            Loads all the datasets from the requested repository
        Parameters
        ----------
        projectName
        repositoryName
        overwrite

        Returns
        -------

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
    def loadAllDatasourcesInRepositoryJSONToProject(self,
                                                    projectName: str,
                                                    repositoryJSON: dict,
                                                    basedir: str = "",
                                                    overwrite: bool = False,
                                                    auto_register_missing: bool = True):
        """
        Iterate through the repository JSON and for each toolkit:
        - Ensure we can instantiate it (ToolkitHome.getToolkit).
        - If missing and auto_register_missing=True, attempt auto-register using:
            * Registry.classpath or Registry.cls in the JSON, or
            * Toolkit document from DB.
        - After we have the instance, proceed with per-toolkit loading logic.
        """
        logger = get_classMethod_logger(self, "loadAllDatasourcesInRepositoryJSONToProject")

        handlerDict = dict(Config = self._handle_Config,
                           Datasource = self._handle_DataSource,
                           Measurements = lambda toolkit, itemName, docTypeDict, overwrite,basedir: self._DocumentHandler(toolkit, itemName, docTypeDict, overwrite,"Measurements",basedir),
                           Simulations = lambda toolkit, itemName, docTypeDict, overwrite,basedir: self._DocumentHandler(toolkit, itemName, docTypeDict, overwrite,"Simulations",basedir),
                           Cache = lambda toolkit, itemName, itemDesc, overwrite,basedir: self._DocumentHandler(toolkit, itemName, itemDesc, overwrite,"Cache",basedir),
                           Function = self._handle_Function)

        # repositoryJSON is expected to be a dict mapping: toolkitName -> section
        for toolkitName, toolkitDict in (repositoryJSON or {}).items():
            toolkit = toolkitHome.getToolkit(toolkitName=toolkitName, projectName=projectName)

            for key, docTypeDict in toolkitDict.items():
                logger.info(f"Loading document type {key} to toolkit {toolkitName}")
                handler = handlerDict.get(key.title(),None)

                if handler is None:
                    err = f"Unkonw Handler {key.title()}. The handler must be {','.join(handlerDict.keys())}. "
                    logger.error(err)
                    raise ValueError(err)
                try:
                    handler(toolkit=toolkit, itemName=key, docTypeDict=docTypeDict, overwrite=overwrite,basedir=basedir)
                except Exception as e:
                    err = f"The error {e} occured while adding *{key}* to toolkit {toolkitName}... skipping!!!"
                    logger.error(err)


    def _handle_Config(self,toolkit,itemName,docTypeDict,overwrite,basedir):
        """
            The procedure to handle the config node.
        Parameters
        ----------
        node

        Returns
        -------

        """
        toolkit.setConfig(**docTypeDict)

    def _handle_DataSource(self,toolkit,itemName,docTypeDict,overwrite,basedir):
        """
            Adds a datasource to the toolkit.
        Parameters
        ----------
        itemName : string
            The name of the item to add
        itemDesc : dict
            The JSON description

        overwrite : bool
            If true, overwrite.

        Returns
        -------

        """
        logger = get_classMethod_logger(self, "_handle_DataSource")

        for itemName, itemDesc in docTypeDict.items():
            theItem = itemDesc["item"]

            isRelativePath = itemDesc.get("isRelativePath")
            assert (isRelativePath=='True' or isRelativePath=='False') or isinstance(isRelativePath,bool), "isRelativePath must be defined as 'True' or 'False'. "
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

    def _DocumentHandler(self, toolkit, itemName, docTypeDict, overwrite, documentType,basedir):
        """
            Handles measurement/cache or simulation document
        Parameters
        ----------
        node
        overwrite

        Returns
        -------
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

    def _handle_Function(self,toolkit,itemName,docTypeDict,overwrite,basedir):
        """
            A general function.

            The key of docTypeDict is
             - dict: params
             - list - list of dicts that will be the params to multiple calls.

            Must have the signature :

            - overwrite: bool

            - [other paraeters] - passed from the item description.

        Parameters
        ----------
        node
        overwrite

        Returns
        -------
        """
        logger = get_classMethod_logger(self, "_handle_GeneralFunction")
        for itemName, itemDesc in docTypeDict.items():
            retrieveFunc = getattr(self,itemName)

            if isinstance(itemDesc,dict):
                retrieveFunc(**itemDesc,overwrite=overwrite)
            elif isinstance(itemDesc,list):
                for imt in itemDesc:
                    if isintance(imt,dict):
                        retrieveFunc(**imt, overwrite=overwrite)
                    else:
                        err = f"{itemName} has a non dict item in the list : {imt}... ignoring."
                        logger.error(err)
            else:
                err = f"{itemName} value must be dict of a list of dicts. "
                logger.error(err)
                raise ValueError(err)


    def _makeItemPathAbsolute(self, theItem,basedir):
        """
            Change the path of the item to be absolution be chechin the is relative flag.

            An internal function.

        Parameters
        ----------
        theItem : dict
            The item data to be added.
            Has the fields:
                 - isRelative field : bool
                - resource : string.
                        Either absolute or relative path.

        Returns
        -------

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
