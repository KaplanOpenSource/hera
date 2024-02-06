import json
import argparse
from hera import toolkitHome
from utils.jsonutils import loadJSON
from ..logging import get_classMethod_logger
import pathlib
import os




class dataToolkit(toolkitHome.abstractToolkit):
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


    def addRepository(self,repositoryName,repositoryPath):
        """
            A path to the repository
        Parameters
        ----------
        repositoryName
        repositoryPath

        Returns
        -------

        """
        self._allowWritingToDefaultProject = True # allows the addition of datasource to the Default project.
        self.addDataSource(dataSourceName=repositoryName, resource=os.path.abspath(repositoryPath))
        self._allowWritingToDefaultProject = False

    def getRepositoryTable(self):
        return self.getDataSourceTable()

    def getRepository(self,repositoryName):
        repo = self.getDatasourceData(datasourceName=repositoryName)
        return loadJSON(repo)

    def loadDataSource(self,projectName,dataItem,toolkitName,repositoryName):
        """
            Loads a specific datasource to the toolkit of the requested project.
        Parameters
        ----------
        projectName
        dataItem
        toolkitName
        repositoryName

        Returns
        -------

        """
        logger = get_classMethod_logger(self,"loadDataSource")
        conf = self.getRepository(repositoryName)
        if toolkitName not in conf:
            err =f"There is not data for toolkit {toolkitName}. Existing toolkits are are {''.join(conf.keys())} "
            logger.error(err)
            raise ValueError(err)

        datasourcesDict = conf[dataItem]

        if dataSourceName not in datasourcesDict:
            err = f"The datasource {dataSourceName} does not exist. Existing datasource are: {','.join(datasourcesDict.keys())}"
            logger.error(err)
            raise ValueError(err)

        dataSourceDesc = datasourcesDict[dataSourceName]
        # dataSourceDesc["resource"] = os.path.join(pathlib.Path(__file__).parent.absolute(),dataSourceDesc["resource"])

        toolkit = toolkitHome.getToolkit(toolkitName=toolkitName, projectName=projectName)

        datasource = toolkit.getDatasourceDocument(datasourceName=dataSourceName)

        if datasource is None:
            toolkit.addDataSource(dataSourceName=dataSourceName, **dataSourceDesc)
            logger.info(f"Added source {dataSourceName} to toolkit {toolkitName} in project {projectName}")
        else:
            logger.info(f"Source {dataSourceName} already exists in {projectName}")

    def loadAllRepositoriesInList(self,projectName,overwrite=False):
        for repository in self.getDataSourceMap().keys():
            logger.info(f"Loading the repository {repository}")
            self.loadAllRepository(projectName,overwrite=overwrite)

    def loadAllRepository(self, projectName, repositoryName,overwrite=False):
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
        conf = self.getRepository(repositoryName)

        for toolkit, datasourcesDict in conf.items():

            for dataSourceName in datasourcesDict.keys():
                print(f"Loading {dataSourceName} to toolkit {toolkit}")

                dataSourceDesc = datasourcesDict[dataSourceName]
                dataSourceDesc["resource"] = os.path.join(pathlib.Path(__file__).parent.absolute(),
                                                          dataSourceDesc["resource"])

                toolkit = toolkitHome.getToolkit(toolkitName=toolkit, projectName=projectName)
                datasource = toolkit.getDatasourceDocument(datasourceName=dataSourceName)

                if datasource is None:
                    toolkit.addDataSource(dataSourceName=dataSourceName, **dataSourceDesc)
                    print(f"Added source {dataSourceName} to tool {toolkit} in project {projectName}")
                elif args.overwrite:
                    datasource = toolkit.getDatasourceDocument(datasourceName=dataSourceName)
                    datasource['resource'] = os.path.abspath(dataSourceDesc["resource"])
                    datasource['desc']['desc'] = dataSourceDesc.get("desc", {})
                    datasource.save()
                    print(f"Updated source {dataSourceName} in tool {toolkit} in project {projectName}")
                else:
                    print(f"Source {dataSourceName} already exists in {projectName}. Use --overwrite to force update")
