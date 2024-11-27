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
from tabulate import tabulate

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
        print("The user does not have repositories.")
    else:
        with pandas.option_context('display.max_rows', None,
                                   'display.max_columns', None,
                                   'display.max_colwidth', None,
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
    dataTypeList = ['DataSource','Measurements','Cache','Simulations']

    for toolkitName, toolDesc in repositoryData.items():
        ttl = f"\t\t\033[1mToolkit:\033[0m {toolkitName}"
        print("#"*(2*len(ttl.expandtabs())))
        print(ttl)
        print("#"*(2*len(ttl.expandtabs())))

        for datatype in dataTypeList:
            print("="*len(datatype))
            print(f"{datatype}")
            print("="*len(datatype))

            for repName,repItems in toolDesc.get(datatype,{}).items():
                ttl = f"\033[1mName:\033[0m {repName}"
                print(f"\t{ttl}")
                print("-" * (2 * len(ttl.expandtabs())))
                print(f"Is it relative path? {repItems['isRelativePath']}")
                with pandas.option_context('display.max_rows', None,
                                           'display.max_columns', None,
                                           'display.width', 1000,
                                           'display.max_colwidth', None,
                                           'display.precision', 3,
                                           'display.colheader_justify', 'center'):
                    print(pandas.DataFrame.from_dict(repItems['item'],orient='index',columns=['Value']))
                    print("\n")

def repository_load(arguments):
    logger = logging.getLogger("hera.bin.repository_load")
    dtk = dataToolkit()

    repositoryFile = arguments.repositoryName
    if 'projectName' in arguments:
        projectName = arguments.projectName
    else:
        projectName = None

    logger.info(f"Loading the repository {repositoryFile} to the project {projectName if projectName is not None else 'default project'}")
    repositoryJSON= loadJSON(repositoryFile)
    dtk.loadAllDatasourcesInRepositoryJSONToProject(projectName=projectName,
                                                    repositoryJSON=repositoryJSON,
                                                    basedir=os.path.dirname(os.path.abspath(arguments.repositoryName)),
                                                    overwrite=arguments.overwrite)

def display_datasource_versions(arguments):
    proj = Project(projectName=arguments.projectName)
    datasources = []
    for document in proj.getMeasurementsDocumentsAsDict()['documents']:
        try:
            d = {}
            d['toolkit'] = document['desc']['toolkit']
            d['datasourceName'] = document['desc']['datasourceName']
            d['version'] = document['desc']['version']
            if d not in datasources:
                if arguments.datasource:
                    if arguments.datasource==d['datasourceName']:
                        datasources.append(d)
                else:
                    datasources.append(d)
        except:
            pass

    if len(datasources)!=0:
        headers = datasources[0].keys()
        rows = [d.values() for d in datasources]
        print(tabulate(rows, headers=headers, tablefmt="grid"))
    else:
        if not arguments.datasource:
            print(f"No data to display. Are you sure project {arguments.projectName} exists?")
        else:
            print(f"No data to display. Are you sure datasource {arguments.datasource} and project {arguments.projectName} exists?")

def update_datasource_default_version(arguments):
    logger = logging.getLogger("hera.bin.update_datasource_version")
    arguments.version = tuple(int(item.strip()) for item in arguments.version.split(','))
    proj = Project(projectName=arguments.projectName)
    proj.setDataSourceDefaultVersion(datasourceName=arguments.datasource,version=arguments.version)

def update(arguments):
    logger = logging.getLogger("hera.bin.update")
    if not arguments.projectName:
        confFile = os.path.join(os.getcwd(), "caseConfiguration.json")
        if not os.path.exists(confFile):
            raise ValueError(f"If projectName is not provided, caseConfiguration json file must be in folder.")
        else:
            configuration = loadJSON(confFile)
            if 'projectName' not in configuration:
                err = f"Got projectName=None and the key 'projectName' does not exist in the JSON. "
                err += """conifguration should be :
                        {
                            'projectName' : [project name]
                        }                                        
                       """
                raise ValueError(err)
            else:
                arguments.projectName = configuration['projectName']

    if 'projectName' in arguments:
        projectName = arguments.projectName
    else:
        projectName = None

    dtk = dataToolkit()
    dtk.loadAllDatasourcesInAllRepositoriesToProject(projectName=projectName, overwrite=arguments.overwrite)

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
