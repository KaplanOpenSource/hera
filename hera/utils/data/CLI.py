import os
from ... import toolkitHome
import getpass
import json
import logging
from ...datalayer import getProjectList,Project,createProjectDirectory,removeConnection,addOrUpdateDatabase,getMongoJSON
from ...datalayer import All as datalayer_All
from .. import loadJSON
from .toolkit import dataToolkit
import pandas

def project_list(arguments):
    """
        List all the projects of the user.
    Parameters
    ----------
    arguments

    Returns
    -------

    """
    connectionName = getpass.getuser() if arguments.connectionName is None else arguments.connectionName

    projectList = getProjectList(connectionName=connectionName)
    ttl = f"Projects in the connection {connectionName}"
    print("\n")
    print(ttl)
    print("-"*len(ttl))
    projList = []
    for projName in projectList:
        projDesct = {"Project Name" : projName}
        if arguments.fulldetails:
            proj = Project(projectName=projName, connectionName=connectionName)
            cacheCount = len(proj.getCacheDocuments())
            measureCount = len(proj.getMeasurementsDocuments())
            simCount = len(proj.getSimulationsDocuments())

            projDesct['Cache documents'] = cacheCount
            projDesct['Measurements documents'] = measureCount
            projDesct['Simulation documents'] = simCount

        projList.append(projDesct)

    df = pandas.DataFrame(projList).sort_values("Project Name")

    with pandas.option_context('display.max_rows', None,
                           'display.max_columns', None,
                           'display.width', 1000,
                           'display.precision', 3,
                           'display.colheader_justify', 'center'):
        print(df)


    print("-"*len(ttl))

def project_create(arguments):
    """
        Creating a directory and a project.

        The project is a caseConfiguration file with the configuration name.

    Parameters
    ----------
    arguments :
        -- directory: the directory to use
        -- database: the name of the DB to use

    Returns
    -------

    """
    if arguments.directory is None:
        directory = os.path.join(os.getcwd(),arguments.projectName)
    else:
        directory = arguments.directory

    createProjectDirectory(outputPath=directory,projectName=arguments.projectName)
    print(f"Created project {arguments.projectName} in directory {directory}")

    if arguments.loadRepositories:
        dtk = dataToolkit()
        dtk.loadAllDatasourcesInAllRepositoriesToProject(projectName=arguments.projectName,overwrite=arguments.overwrite)

def project_dump(arguments):

    fullQuery=dict(projectName = arguments.projectName)

    for queryElement in arguments.query:
        fullQuery[queryElement.split('=')[0]] = eval(queryElement.split('=')[1])


    docList = []
    for doc in datalayer_All.getDocuments(**fullQuery):
        docDict = doc.asDict()
        if ('docid' not in docDict['desc']):
            docDict['desc']['docid'] = str(doc.id)

        docList.append(docDict)

    outStr = json.dumps(docList,indent=4)

    outputFileName = arguments.fileName
    if outputFileName is not None:
        with open(outputFileName,"w") as outputFile:
            outputFile.write(outStr)

    if arguments.outputFormat=='json':
        print(outStr)
    elif arguments.outputFormat=='table':
        df = pandas.DataFrame(docList)
        with pandas.option_context('display.max_rows', None,
                                   'display.max_columns', None,
                                   'display.width', 1000,
                                   'display.precision', 3,
                                   'display.colheader_justify', 'center'):
            print(df)
    else:
        print("The outputFormat must be 'json' or 'table'")


def project_load(arguments):

    docsDict = loadJSON(arguments.fileName)
    proj     = Project(projectName=arguments.projectName)

    for indx,doc in enumerate(docsDict):
        print(f"Loading document {indx}/{len(docsDict)}")
        proj.addDocumentFromDict(doc)


def repository_list(argumets):
    dtk = dataToolkit()

    repDataframe = dtk.getRepositoryTable()
    if len(repDataframe) ==0:
        print("No repositories loaded")
    else:
        with pandas.option_context('display.max_rows', None,
                                   'display.max_columns', None,
                                   'display.width', 1000,
                                   'display.precision', 3,
                                   'display.colheader_justify', 'center'):
            print(repDataframe)

def repository_add(argumets):
    logger = logging.getLogger("hera.bin.repository_add")
    dtk = dataToolkit()

    repositoryName = os.path.basename(argumets.repositoryName).split(".")[0]

    logger.info(f"Adding the repository {argumets.repositoryName} as name {repositoryName}")
    dtk.addRepository(repositoryName=repositoryName,
                      repositoryPath=argumets.repositoryName,
                      overwrite=argumets.overwrite)

def repository_remove(arguments):
    logger = logging.getLogger("hera.bin.repository_remove")
    dtk = dataToolkit()

    datasourceName = arguments.repositoryName
    logger.info(f"Removing the datasource {datasourceName}")
    dtk.deleteDataSource(datasourceName=datasourceName)


def repository_show(arguments):
    logger = logging.getLogger("hera.bin.repository_remove")
    dtk = dataToolkit()

    datasourceName = arguments.repositoryName
    logger.info(f"Listing the datasource {datasourceName}")
    repositoryData = dtk.getDataSourceData(datasourceName=datasourceName)
    for toolkitName, toolDesc in repositoryData.items():
        ttl = f"\t\t\033[1mToolkit:\033[0m {toolkitName}"
        print("#"*(2*len(ttl.expandtabs())))
        print(ttl)
        print("#"*(2*len(ttl.expandtabs())))

        print(f"DataSource")
        print("===========")
        for repName,repItems in toolDesc.get("dataSource",{}).items():
            ttl = f"\033[1mRepository name:\033[0m {repName}"
            print(f"\t{ttl}")
            print("\t"+"-"*len(ttl))
            with pandas.option_context('display.max_rows', None,
                                       'display.max_columns', None,
                                       'display.width', 1000,
                                       'display.precision', 3,
                                       'display.colheader_justify', 'center'):
                print(pandas.DataFrame.from_dict(repItems,orient='index',columns=['Value']))

        for additionalType in ['Measuerments','Cache','Simulations']:
            if additionalType in toolDesc:
                print(additionalType)
                print("="*len(additionalType))

                for repName,repItems in toolDesc.get(additionalType,{}).items():
                    ttl = f"\tRepository name: {repName}"
                    print(ttl)
                    print("-" * (2 * len(ttl.expandtabs())))
                    print(f"Is path? {repItems['isPath']}")
                    print(f"Is path? {repItems['isPath']}")
                    if repItems['isPath']:
                        print(f"\tAbsolute path{repItems.get('absolutePath',False)}")

                    with pandas.option_context('display.max_rows', None,
                                               'display.max_columns', None,
                                               'display.width', 1000,
                                               'display.precision', 3,
                                               'display.colheader_justify', 'center'):
                        print(pandas.DataFrame.from_dict(repItems['item'],orient='index',columns=['Value']))



def db_list(arguments):
    """
        List the databases in the
    Parameters
    ----------
    arguments

    Returns
    -------

    """
    dbconfig = getMongoJSON()
    conList = []
    for connectionName,connectionData in dbconfig.items():
        condict = {"Connection Name" : connectionName}
        if arguments.fulldetails:
            condict.update(connectionData)

        conList.append(condict)

    df = pandas.DataFrame(conList).rename(columns=dict(dbIP="IP",dbName="databaseName"))

    with pandas.option_context('display.max_rows', None,
                           'display.max_columns', None,
                           'display.width', 1000,
                           'display.precision', 3,
                           'display.colheader_justify', 'center'):
        print(df)

def db_create(arguments):
    addOrUpdateDatabase(connectionName=arguments.connectionName,
                        username=arguments.username,
                        password=arguments.password,
                        databaseIP=arguments.IP,
                        databaseName=arguments.databaseName)

def db_remove(arguments):
    removeConnection(arguments.connectionName)
