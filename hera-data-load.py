#! /usr/bin/env python
import json
import argparse
import pathlib
from hera import toolkitHome
import os
import sys

# sys.path.append(str(pathlib.Path(__file__).parent.parent))


with open(os.path.join(pathlib.Path(__file__).parent,"datasources.json"), "r") as datasourceConfig:
    conf = json.load(datasourceConfig)

def load(projectName, dataType, dataSourceName):


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
    else:
        print(f"Source {dataSourceName} already exists in {projectName}")

if __name__ =="__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subparser_name")
    list_parser = subparsers.add_parser('list')
    load_parser = subparsers.add_parser('load')

    load_parser.add_argument('projectName', type=str,help='The project to load the data')
    load_parser.add_argument('dataType', type=str,help='The data type to load')
    load_parser.add_argument('source', type=str,help='The name of the datasource to load')

    list_parser.add_argument('dataType',type=str,default=None,nargs='*')
    args = parser.parse_args()

    if (args.subparser_name == 'list'):
        if len(args.dataType)==0:

            print("The available data types are:")
            print("\n".join(conf['dataType'].keys()))
        else:
            dataType = args.dataType[0]
            if dataType not in conf['dataType'].keys():
                print("The available data types are:")
                print("\n".join(conf['dataType'].keys()))
            else:
                print("The available data sources are:")
                print("\n".join(conf['dataType'][dataType].keys()))


    elif (args.subparser_name == 'load'):
        load(args.projectName, args.dataType, args.source)
