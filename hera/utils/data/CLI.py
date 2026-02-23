import os
from ... import toolkitHome
import getpass
import json
import logging
from ...datalayer import getProjectList, Project, createProjectDirectory, removeConnection, addOrUpdateDatabase, getMongoJSON
from ...datalayer import All as datalayer_All
from .. import loadJSON
from .toolkit import dataToolkit
import pandas
from ...toolkit import ToolkitHome
from pydoc import locate  # for resolving classpath -> class
from tabulate import tabulate


def project_list(arguments):
    """
        List all the projects of the user.
    """
    connectionName = getpass.getuser() if arguments.connectionName is None else arguments.connectionName

    projectList = getProjectList(connectionName=connectionName)
    ttl = f"Projects in the connection {connectionName}"
    print("\n")
    print(ttl)
    print("-" * len(ttl))
    projList = []
    for projName in projectList:
        projDesct = {"Project Name": projName}
        if arguments.fulldetails:
            proj = Project(projectName=projName, connectionName=connectionName)
            cacheCount = len(proj.getCacheDocuments())
            measureCount = len(proj.getMeasurementsDocuments())
            simCount = len(proj.getSimulationsDocuments())

            projDesct['Cache documents'] = cacheCount
            projDesct['Measurements documents'] = measureCount
            projDesct['Simulation documents'] = simCount

        projList.append(projDesct)

    df = pandas.DataFrame(projList)
    if df.size > 0:
        df = df.sort_values("Project Name")

    with pandas.option_context('display.max_rows', None,
                               'display.max_columns', None,
                               'display.width', 1000,
                               'display.precision', 3,
                               'display.colheader_justify', 'center'):
        print(df)

    print("-" * len(ttl))


def project_create(arguments):
    """
        Creating a directory and a project.
    """
    if arguments.directory is None:
        directory = os.path.join(os.getcwd(), arguments.projectName)
    else:
        directory = arguments.directory

    createProjectDirectory(outputPath=directory, projectName=arguments.projectName)
    print(f"Created project {arguments.projectName} in directory {directory}")

    if arguments.loadRepositories:
        dtk = dataToolkit()
        dtk.loadAllDatasourcesInAllRepositoriesToProject(projectName=arguments.projectName,
                                                         overwrite=arguments.overwrite)


def project_dump(arguments):

    fullQuery = dict(projectName=arguments.projectName)

    for queryElement in arguments.query:
        fullQuery[queryElement.split('=')[0]] = eval(queryElement.split('=')[1])

    docList = []
    for doc in datalayer_All.getDocuments(**fullQuery):
        docDict = doc.asDict()
        if 'docid' not in docDict['desc']:
            docDict['desc']['docid'] = str(doc.id)

        docList.append(docDict)

    outStr = json.dumps(docList, indent=4)

    outputFileName = arguments.fileName
    if outputFileName is not None:
        with open(outputFileName, "w") as outputFile:
            outputFile.write(outStr)

    if arguments.outputFormat == 'json':
        print(outStr)
    elif arguments.outputFormat == 'table':
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

    docsDict = loadJSON(arguments.file)
    proj = Project(projectName=arguments.projectName)

    for indx, doc in enumerate(docsDict):
        print(f"Loading document {indx}/{len(docsDict)}")
        proj.addDocumentFromDict(docsDict.get(doc))


def repository_list(argumets):
    dtk = dataToolkit()

    repDataframe = dtk.getRepositoryTable()
    if len(repDataframe) == 0:
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
    dataTypeList = ['DataSource', 'Measurements', 'Cache', 'Simulations']

    for toolkitName, toolDesc in repositoryData.items():
        ttl = f"\t\t\033[1mToolkit:\033[0m {toolkitName}"
        print("#" * (2 * len(ttl.expandtabs())))
        print(ttl)
        print("#" * (2 * len(ttl.expandtabs())))

        for datatype in dataTypeList:
            print("=" * len(datatype))
            print(f"{datatype}")
            print("=" * len(datatype))

            for repName, repItems in toolDesc.get(datatype, {}).items():
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
                    print(pandas.DataFrame.from_dict(repItems['item'], orient='index', columns=['Value']))
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
    repositoryJSON = loadJSON(repositoryFile)
    dtk.loadAllDatasourcesInRepositoryJSONToProject(projectName=projectName,
                                                    repositoryJSON=repositoryJSON,
                                                    basedir=os.path.dirname(os.path.abspath(arguments.repositoryName)),
                                                    overwrite=arguments.overwrite)

def add_toolkit(arguments):
    """
    CLI entry point for:
        hera-project addToolkit <toolkit_name> <toolkit_path>
            [--params JSON]
            [--version d,d,d]
            [--overwrite]
    """
    logger = logging.getLogger("hera.bin.hera_projects")
    logger.debug(f"Adding toolkit with args: {arguments}")
    # -------- Parse version --------
    try:
        version_tuple = tuple(int(x.strip()) for x in arguments.version.split(","))
        if len(version_tuple) != 3:
            raise ValueError
    except Exception:
        raise ValueError(
            f"Invalid version format '{arguments.version}'. "
            "Expected format: d,d,d (example: 0,0,1)"
        )
    # -------- Parse params --------
    params_dict = None
    if arguments.params:
        try:
            params_dict = json.loads(arguments.params)
            if not isinstance(params_dict, dict):
                raise ValueError("Params must be a JSON object.")
        except Exception as e:
            raise ValueError(
                f"Invalid JSON passed to --params: {arguments.params}"
            ) from e
    # -------- Call registerToolkit --------
    try:
        tkh = ToolkitHome()
        tkh.registerToolkit(
            toolkit_name=arguments.toolkit_name,
            toolkit_path=arguments.toolkit_path,
            params=params_dict,
            version=version_tuple,
            overwrite=arguments.overwrite,
        )
        logger.info(
            f"Toolkit '{arguments.toolkit_name}' registered successfully "
            f"(version={version_tuple}, overwrite={arguments.overwrite})"
        )
    except Exception as e:
        logger.error(f"Failed to register toolkit '{arguments.toolkit_name}': {e}")
        raise

def display_datasource_versions(arguments):
    proj = Project(projectName=arguments.projectName)
    datasources = []

    if not arguments.default:
        for document in proj.getMeasurementsDocumentsAsDict()['documents']:
            try:
                d = {}
                d['toolkit'] = document['desc']['toolkit']
                d['datasourceName'] = document['desc']['datasourceName']
                d['version'] = document['desc']['version']

                if arguments.datasource:
                    if arguments.datasource == d['datasourceName']:
                        datasources.append(d)
                else:
                    datasources.append(d)
            except:
                pass
    else:
        config = proj.getConfig()
        for document in proj.getMeasurementsDocumentsAsDict()['documents']:
            try:
                d = {}
                d['toolkit'] = document['desc']['toolkit']
                d['datasourceName'] = document['desc']['datasourceName']

                if arguments.datasource:
                    if arguments.datasource == d['datasourceName']:
                        default_version = config.get(f"{arguments.datasource}_defaultVersion")
                    else:
                        default_version = None
                else:
                    default_version = config.get(f"{d['datasourceName']}_defaultVersion")

                if default_version:
                    d['DEFAULT_VERSION'] = default_version
                    datasources.append(d)

            except:
                pass

    if len(datasources) != 0:
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
    proj.setDataSourceDefaultVersion(datasourceName=arguments.datasource, version=arguments.version)


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
    """
    dbconfig = getMongoJSON()
    conList = []
    for connectionName, connectionData in dbconfig.items():
        condict = {"Connection Name": connectionName}
        if arguments.fulldetails:
            condict.update(connectionData)

        conList.append(condict)

    df = pandas.DataFrame(conList).rename(columns=dict(dbIP="IP", dbName="databaseName"))

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


# --- Toolkit related CLI ---

def toolkit_list(arguments):
    """
    Print a combined list of toolkits (static + dynamic from DB) for a project.
    Uses ToolkitHome.getToolkitDocuments(...) as the single source of truth.
    """
    logger = logging.getLogger("hera.utils.CLI.toolkit_list")
    project = arguments.project

    try:
        th = ToolkitHome(projectName=project)

        docs = th.getToolkitDocuments(name=None, projectName=project) or []
        rows = []
        for d in docs:
            desc = d.get("desc", {}) if isinstance(d, dict) else getattr(d, "desc", {}) or {}
            rows.append({
                "toolkit": d.get("toolkit") if isinstance(d, dict) else getattr(d, "toolkit", ""),
                "cls": desc.get("classpath", ""),
                "source": desc.get("source", "internal"),
                "type": desc.get("type", ""),
                "repository": desc.get("repositoryName", ""),
                "version": desc.get("version", ""),
            })

        if rows:
            print(tabulate(rows, headers="keys", tablefmt="github"))
        else:
            print("No toolkits found (static + DB are empty).")

    except Exception as e:
        logger.exception(e)
        print(f"[ERROR] {e}")


def toolkit_register(arguments):
    """
    Register a toolkit using a classpath and metadata.
    Delegates to ToolkitHome.registerToolkit (expects a class object).
    """
    logger = logging.getLogger("hera.utils.CLI.toolkit_register")
    project = arguments.project
    cls_path = arguments.cls
    name = arguments.name
    repository = arguments.repository
    version_txt = arguments.version

    # parse version "0,0,1" -> (0,0,1)
    try:
        version = tuple(int(x.strip()) for x in version_txt.split(","))
    except Exception:
        version = (0, 0, 1)

    try:
        th = ToolkitHome(projectName=project)

        # Resolve classpath -> class
        toolkit_cls = locate(cls_path)
        if toolkit_cls is None:
            raise ImportError(f"Cannot locate class by classpath: {cls_path}")

        # Call registerToolkit
        th.registerToolkit(
            toolkit_name = toolkit_name,
            toolkit_path = toolkit_path,
            params = None,
            version = version,
            overwrite = False,
        )
        print(f"Toolkit '{name}' registered (cls={cls_path}, repo={repository}, version={version}).")
    except Exception as e:
        logger.exception(e)
        print(f"[ERROR] {e}")

def load_project_name(config_file: str = "caseConfiguration.json") -> str:
    """
    Load projectName from caseConfiguration.json.
    """
    if not os.path.exists(config_file):
        raise FileNotFoundError(
            f"Could not find '{config_file}' in {os.getcwd()}"
        )

    with open(config_file, "r") as f:
        data = json.load(f)

    project_name = data.get("projectName")
    if not project_name:
        raise ValueError("caseConfiguration.json does not contain 'projectName'")

    return project_name


def toolkit_load(arguments):
    """
    Instantiate a toolkit by name.
    Delegates to ToolkitHome.getToolkit (static + dynamic + experiments).
    """
    logger = logging.getLogger("hera.utils.CLI.toolkit_load")
    project = arguments.project
    name = arguments.name
    try:
        # ToolkitHome itself is a toolkit; projectName is loaded automatically
        th = ToolkitHome(projectName=project)
        tk = th.getToolkit(toolkitName=name)
        print(f"Loaded toolkit: {getattr(tk, 'name', name)}")
    except Exception as e:
        logger.exception(e)
        print(f"[ERROR] {e}")


def toolkit_default_repo_show(arguments):
    """
    Show the project's default repository, via ToolkitHome.getDefaultRepository(projectName=...).
    """
    logger = logging.getLogger("hera.utils.CLI.toolkit_default_repo_show")
    project = getattr(arguments, "project", None) or "DefaultProject"
    try:
        th = ToolkitHome(projectName=project)
        repo = th.getDefaultRepository(projectName=project)
        print(repo if repo else "<no default repository>")
    except Exception as e:
        logger.exception(e)
        print(f"[ERROR] {e}")


def toolkit_default_repo_set(arguments):
    """
    Set the project's default repository, via ToolkitHome.setDefaultRepository(projectName=..., repositoryName=...).
    """
    logger = logging.getLogger("hera.utils.CLI.toolkit_default_repo_set")
    project = getattr(arguments, "project", None) or "DefaultProject"
    repo_name = getattr(arguments, "repository", None)
    if not repo_name:
        print("[ERROR] --repository is required")
        return
    try:
        th = ToolkitHome(projectName=project)
        th.setDefaultRepository(projectName=project, repositoryName=repo_name)
        print(f"Default repository set to '{repo_name}' for project '{project}'.")
    except Exception as e:
        logger.exception(e)
        print(f"[ERROR] {e}")


def toolkit_import_json(arguments):
    """
    Import a JSON repository that declares toolkits (and optionally experiments),
    and register them into the project.
    """
    logger = logging.getLogger("hera.utils.CLI.toolkit_import_json")
    project = getattr(arguments, "project", None)
    json_path = getattr(arguments, "file", None)
    no_experiments = getattr(arguments, "no_experiments", False)

    try:
        th = ToolkitHome(projectName=project)

        registered = th.import_toolkits_from_json(projectName=project, json_path=json_path)
        print(f"Registered toolkits: {registered}" if registered else "No toolkits in JSON.")

        if not no_experiments:
            exps = th.import_experiments_from_json(projectName=project, json_path=json_path)
            if exps:
                print(f"Loaded experiments data: {exps}")

    except Exception as e:
        logger.exception(e)
        print(f"[ERROR] {e}")


# שים לב: הקטע הבא נראה כאילו אמור להיות בתוך class, אבל אני משאיר כמו שהיה
def project_measurements_list(args):
    """
    Implementation for:
      hera-project project measurements list
      (used also by 'simulations list' and 'cache list' via shortcut)
    """
    from hera.datalayer import Project

    project = getattr(args, "project", None)
    mtype = getattr(args, "type", None)
    shortcut = getattr(args, "shortcut", None)
    contains = getattr(args, "contains", None)

    # Map shortcuts to measurement types
    shortcut_map = {
        "ds": "ToolkitDataSource",      # dynamic toolkits
        "exp": "Experiment_rawData",    # experiments
        "sim": "Simulation_rawData",    # (if any)
        "cache": "Cache_rawData",       # (if any)
        "all": None,                    # no type filter
    }

    if shortcut:
        if shortcut not in shortcut_map:
            print(f"Unknown shortcut '{shortcut}'. Valid: {', '.join(shortcut_map.keys())}")
            return
        mtype = shortcut_map[shortcut]

    # Open project (or default)
    proj = Project(projectName=project) if project else Project()

    query = {}
    if mtype:
        query["type"] = mtype

    docs = proj.getMeasurementsDocuments(**query)

    # Optional substring filter
    if contains:
        filtered = []
        for d in docs:
            name = getattr(d, "datasourceName", "") or d.desc.get("datasourceName", "")
            resource = getattr(d, "resource", "") or ""
            if contains in str(name) or contains in str(resource):
                filtered.append(d)
        docs = filtered

    if not docs:
        print("No measurements found for given filters.")
        return

    # Build rows for pretty table
    rows = []
    for d in docs:
        rows.append({
            "type": getattr(d, "type", ""),
            "datasourceName": getattr(d, "datasourceName", "") or d.desc.get("datasourceName", ""),
            "resource": getattr(d, "resource", ""),
            "dataFormat": getattr(d, "dataFormat", ""),
            "version": d.desc.get("version", ""),
            "repository": d.desc.get("repository", ""),
        })

    import pandas as pd
    df = pd.DataFrame(rows)
    print(df.to_markdown(index=False))
