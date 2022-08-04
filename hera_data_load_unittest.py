import json
import pathlib
from hera import toolkitHome
import os
import sys

sys.path.append(str(pathlib.Path(__file__).parent.parent))


def loader(projectName, dataType, dataSourceName):

    with open(os.path.join(pathlib.Path(__file__).parent,"datasources.json"), "r") as datasourceConfig:
        conf = json.load(datasourceConfig)

    if dataType not in conf['dataType']:
        raise ValueError(f"Data type {dataType} does not exist. Datasources are {''.join(conf['dataType'].keys())} ")

    datasourcesDict = conf['dataType'][dataType]

    if dataSourceName not in datasourcesDict:
        raise ValueError(f"The datasource {dataSourceName} does not exist. Existing datasource are: {','.join(datasourcesDict.keys())}")

    dataSourceDesc = datasourcesDict[dataSourceName]

       
    dataSourceDesc["resource"] = os.path.join(pathlib.Path(__file__).parent.absolute(),dataSourceDesc["resource"])

    toolkit = toolkitHome.getToolkit(toolkitName=dataType, projectName=projectName)

    datasource = toolkit.getDatasourceDocument(datasourceName=dataSourceName)

    if datasource is None:
        toolkit.addDataSource(dataSourceName=dataSourceName,**dataSourceDesc)
        print(f"Added source {dataSourceName} to tool {dataType} in project {projectName}")
        return(f"Added source {dataSourceName} to tool {dataType} in project {projectName}")
    else:
        print(f"Source {dataSourceName} already exists in {projectName}")
        return(f"Source {dataSourceName} already exists in {projectName}")