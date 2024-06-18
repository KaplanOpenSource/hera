import json
import argparse
from ... import toolkitHome
from ... import toolkit
from ...utils import loadJSON, dictToMongoQuery
from ..logging import get_classMethod_logger
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

    def __init__(self):
        super().__init__(toolkitName="heradata", projectName=self.DEFAULTPROJECT, filesDirectory=None)

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
                logger.info(f"Loading the repository {repository}")
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

    def loadAllDatasourcesInRepositoryJSONToProject(self, projectName, repositoryJSON,basedir, overwrite=False):
        logger = get_classMethod_logger(self, "loadAllDatasourcesInRepositoryJSONToProject")
        for toolkitName, toolkitDict in repositoryJSON.items():
            logger.info(f"Loading into toolkit  {toolkitName}")
            try:
                toolkit = toolkitHome.getToolkit(toolkitName=toolkitName, projectName=projectName)

                for key, docTypeDict in toolkitDict.items():
                    logger.info(f"Loading document type {key} to toolkit {toolkitName}")
                    print(toolkitName)
                    if key.title()=='Config':
                        toolkit.setConfig(**docTypeDict)
                    else:
                        for itemName, itemDesc in docTypeDict.items():
                            logger.info(f"Loading {itemName} to toolkit {toolkitName} (ProjectName {toolkit.projectName}")

                            theItem = itemDesc["item"]

                            isRelativePath = bool(itemDesc.get("isRelativePath", True))
                            # logger.debug(f"Checking if {itemName} resource is a path {isRelativePath}, is it absolute? {isAbsolute}")

                            if isRelativePath:
                                logger.debug(
                                    f"The input is not absolute (it is relative). Adding the path {basedir} to the resource {theItem['resource']}")
                                theItem["resource"] = os.path.join(basedir, theItem["resource"])

                            logger.debug(f"Checking if the data item {itemName} is already in the project")
                            if key == 'DataSource':
                                datasource = toolkit.getDataSourceDocuments(datasourceName=itemName)
                            else:
                                retrieveFuncName = f"get{key}Documents"
                                retrieveFunc = getattr(toolkit, retrieveFuncName)
                                if retrieveFunc is None:
                                    raise ValueError(
                                        f"function {retrieveFuncName} not found. Key {key} must be : DataSource, Measurement, Cache, or Simulation")
                                qrydict = dict(theItem)
                                del qrydict['resource']
                                del qrydict['dataFormat']

                                itemQry = dictToMongoQuery(qrydict)
                                datasource = retrieveFunc(**itemQry)

                            logger.debug(f"Found {len(datasource)} documents")

                            if len(datasource) == 0:
                                logger.debug("Adding a new document")
                                if key == 'DataSource':
                                    funcName = f"add{key}"
                                    theItem['dataSourceName'] = itemName
                                    theItem['overwrite'] = overwrite
                                else:
                                    funcName = f"add{key}Document"

                                logger.debug(f"Adding the document of type {key} using the function {funcName}")
                                func = getattr(toolkit, funcName)

                                func(**theItem)
                                logger.info(f"Added source {itemName} to tool {toolkitName} in project {projectName}")

                            elif overwrite:
                                logger.debug("Updating an existing document")
                                dataitem = datasource[0]
                                dataitem['resource'] = theItem["resource"]
                                dataitem['dataFormat'] = theItem['dataFormat']
                                curDesc = theItem.get("desc", {})
                                curDesc.update(dataitem['desc'])
                                dataitem['desc'] = curDesc
                                dataitem.save()
                                logger.info(f"Updated source {itemName} in tool {toolkitName} in project {projectName}")
                            else:
                                logger.error(
                                    f"Source {itemName} already exists in {projectName}. Use --overwrite to force update")
            except Exception as e:
                err = f"The error {e} occured while adding toolkit {toolkitName}... skipping!!!"
                logger.error(err)