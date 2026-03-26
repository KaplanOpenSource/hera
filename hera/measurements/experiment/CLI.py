import os
import logging
from ...utils.jsonutils import loadJSON
from ... import toolkitHome
from argos.experimentSetup.dataObjects import ExperimentZipFile
import pandas
import json
import shutil
from hera import datalayer
from ...utils.data.toolkit import dataToolkit
import glob
import requests
import zipfile
from hera import toolkitHome
from hera.utils.data import CLI as projectCLI

def _resolve_project_name(arguments):
    """
    Resolve project name:
    - Prefer --projectName if provided
    - Else fall back to caseConfiguration.json (legacy behavior)
    """
    projectName = getattr(arguments, "projectName", None)
    if projectName:
        return projectName

    configurationFile = getattr(arguments, "configurationFile", "caseConfiguration.json")
    configuration = loadJSON(configurationFile)
    return configuration["projectName"]


def experiments_list(arguments):
    """
    Print the list of experiment names in a project.

    Parameters
    ----------
    arguments : argparse.Namespace
        Must contain projectName or configurationFile.
    """
    logger = logging.getLogger("hera.bin.experiment_experiments_list")
    logger.debug(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")

    projectName = _resolve_project_name(arguments)

    # ✅ FIX: avoid toolkitHome.getToolkit(...) to prevent duplicate projectName injection
    from hera.measurements.experiment.experiment import experimentHome
    home = experimentHome(projectName=projectName)

    # keep old behavior: "print(tk.keys())" => print the experiment names
    print(list(home.getExperimentsMap().keys()))


def experiments_table(arguments):
    """
    Print a formatted table of experiments in a project.

    Parameters
    ----------
    arguments : argparse.Namespace
        Must contain projectName or configurationFile.
    """
    logger = logging.getLogger("hera.bin.experiment_experiments_table")
    logger.debug(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")

    projectName = _resolve_project_name(arguments)

    # ✅ FIX: instantiate directly
    from hera.measurements.experiment.experiment import experimentHome
    home = experimentHome(projectName=projectName)

    # If experimentHome has getExperimentsTable(), use it. Otherwise print a simple table.
    if hasattr(home, "getExperimentsTable"):
        print(home.getExperimentsTable())
        return

    names = list(home.getExperimentsMap().keys())
    print(f"\nExperiments in project '{projectName}':")
    if not names:
        print("(none)")
        return
    for i, n in enumerate(names, start=1):
        print(f"{i:>3}. {n}")


def get_experiment_data(arguments):
    """
    Retrieve and print experiment measurement data as a DataFrame.

    Parameters
    ----------
    arguments : argparse.Namespace
        Must contain projectName, experiment, deviceType, and optionally
        deviceName and perDevice.
    """
    logger = logging.getLogger("hera.bin.experiment_get_experiment_data")
    logger.debug(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")

    projectName = _resolve_project_name(arguments)

    if projectName not in datalayer.getProjectList():
        raise ValueError(f"Project '{projectName}' does not exists")

    # ✅ FIX: instantiate directly (no toolkitHome.getToolkit)
    from hera.measurements.experiment.experiment import experimentHome
    home = experimentHome(projectName=projectName)

    exp = home.getExperiment(arguments.experiment)
    parquet = exp.getExperimentData().getData(
        arguments.deviceType,
        deviceName=getattr(arguments, "deviceName", None),
        perDevice=getattr(arguments, "perDevice", None),
    )
    print(pandas.DataFrame(parquet))


def create_experiment(arguments):
    """
    Scaffold a new experiment directory with code, data, and repository files.

    Parameters
    ----------
    arguments : argparse.Namespace
        Must contain experimentName, path, zip, and relative.
    """
    logger = logging.getLogger("hera.bin.experiment_create_experiment")
    logger.debug(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")
    if arguments.path:
        experiment_path = arguments.path
    else:
        experiment_path = os.path.join(os.getcwd(),arguments.experimentName)

    logger.debug(f" creating code directory if not exist")
    os.makedirs(os.path.join(experiment_path,'code'),exist_ok=True)

    logger.debug(f" checking if class script already exists..")
    try:
        class_script = open(f"{os.path.join(experiment_path, 'code', arguments.experimentName)}.py", "x")
        logger.debug(f" creating the experiment class since it does not exists")
        _create_empty_class(experiment_path,arguments.experimentName)
    except:
        logger.debug(f" experiment class already exists ")

    logger.debug(f" creating data directory if not exist")
    os.makedirs(os.path.join(experiment_path, 'data'), exist_ok=True)

    arguments.projectName = arguments.experimentName
    arguments.overwrite   = True
    arguments.directory   = experiment_path
    arguments.loadRepositories = True

    logger.debug(f"Creating a project for local")
    projectCLI.project_create(arguments)

    _create_repository(arguments.zip,experiment_path,arguments.experimentName,arguments.relative)
    _make_runtimeExperimentData(arguments.zip,experiment_path,arguments.experimentName)

    logger.debug("Loading the experiment repository to the project")
    arguments.repositoryName = os.path.join(experiment_path,f"{arguments.experimentName}_repository.json")
    projectCLI.repository_load(arguments)

    with open("createNodeRedDeviceMap.sh","w") as noderedFile:
        noderedFile.writelines("argos-experiment-manager nodered createDeviceMap")

    os.chmod("createNodeRedDeviceMap.sh", 0o755)

def _create_empty_class(experiment_path,experimentName):
    """
    Generate a boilerplate Python class file for a new experiment.

    Parameters
    ----------
    experiment_path : str
        Root directory of the experiment.
    experimentName : str
        Name used for the class and file.
    """
    logger = logging.getLogger("hera.bin._create_empty_class")
    logger.debug(f" creating an empty class for implementation..")

    baseClass = f"""
from hera.measurements.experiment.experiment import experimentSetupWithData
from hera.measurements.experiment.analysis import experimentAnalysis
from hera.measurements.experiment.presentation import experimentPresentation


class {experimentName}(experimentSetupWithData):
    def __init__(self,projectName,pathToExperiment,filesDirectory):
        super().__init__(projectName=projectName,pathToExperiment=pathToExperiment,filesDirectory=filesDirectory)
        self._analysis = {experimentName}Analysis(self)
        self._presentation = {experimentName}Presentation(self,self.analysis)

class {experimentName}Analysis(experimentSetupWithData):
    datalayer = None 
    def __init__(self,datalayer):
        self.datalayer = datalayer


class {experimentName}Presentation(experimentSetupWithData):
    datalayer = None 
    analysis  = None 
    
    def __init__(self,datalayer,analysis):
        self.datalayer = datalayer
        self.analysis  = analysis

"""

    with open(f"{os.path.join(experiment_path, 'code', experimentName)}.py", "w") as class_script:
        class_script.write(baseClass)

    logger.debug(f" finished creating an empty class for implementation..")

def _create_repository(argos_zip,experiment_path,experimentName,relative):
    """
    Create a repository JSON file and datasource configuration from metadata.

    Parameters
    ----------
    argos_zip : str or None
        Path to an Argos experiment zip file, or None.
    experiment_path : str
        Root directory of the experiment.
    experimentName : str
        Name of the experiment.
    relative : bool
        If True, store resource paths as relative.
    """
    logger = logging.getLogger("hera.bin._create_repository")
    logger.info(f"Creating the repository")
    logger.debug(f" Since zip file is provided, creating a repository..")
    metadata = ExperimentZipFile(argos_zip) if argos_zip else None

    repo = {}
    perDevice = True         ##Will be defined by the updated zip format!

    if relative:
        is_relative = 'True'
        resource_path = '.'
    else:
        is_relative = 'False'
        resource_path = experiment_path

    repo['experiment'] = {}
    repo['experiment']['DataSource'] = {}
    repo['experiment']['DataSource'][experimentName] = {"isRelativePath": is_relative,
                                                                       "item":{
                                                                           "dataSourceName": experimentName,
                                                                           "resource": resource_path,
                                                                           "dataFormat": "string",
                                                                           "overwrite": "True"
                                                                       }

                                                                  }

    repo['experiment']['Measurements'] = {}
    if metadata:
        entities_dict_list = metadata.getExperimentEntities()
        for i,entity in enumerate(entities_dict_list):
            entityTypeName = entity['entityTypeName']
            entityName = entity['entityName']
            logger.debug(f" Creating record for {entityName} of type {entityTypeName}")

            parquet_name = entityName if metadata.entityType[entityTypeName][entityName].properties['StoreDataPerDevice'] else entityTypeName
            if parquet_name not in repo['experiment']['Measurements'].keys():
                repo['experiment']['Measurements'][parquet_name] = {"isRelativePath": "True",
                                                                      "item": {
                                                                      "type": "Experiment_rawData",
                                                                      "resource": os.path.join('data',f"{parquet_name}.parquet"),
                                                                      "dataFormat": "parquet",
                                                                      "desc": {
                                                                              "deviceType": entity['entityTypeName'],
                                                                              "experimentName": experimentName,
                                                                              "deviceName": entity['entityName']
                                                                            }
                                                                          }
                                                                      }

    with open(os.path.join(experiment_path,f'{experimentName}_repository.json'), "w") as f:
        json.dump(repo, f, indent=4)
    
    default_datasources = {"experimentName": experimentName}
    configurationFilePath = os.path.join(experiment_path, 'runtimeExperimentData', "Datasources_Configurations.json")
    with open(configurationFilePath, "w") as f:
        json.dump(default_datasources, f, indent=4)

    logger.debug(f" finished creating the repository json file")

def _make_runtimeExperimentData(argos_zip,experiment_path,experimentName):
    """
    Create the runtimeExperimentData directory and its configuration file.

    Parameters
    ----------
    argos_zip : str or None
        Path to an Argos zip file to copy, or None.
    experiment_path : str
        Root directory of the experiment.
    experimentName : str
        Name of the experiment.
    """
    logger = logging.getLogger("hera.bin._make_runtimeExperimentData")
    logger.debug(f" creating runtimeExperimentData directory if it does not exists")
    os.makedirs(os.path.join(experiment_path, 'runtimeExperimentData'), exist_ok=True)
    logger.debug(f" creating Datasources_Configurations json")
    config = {"experimentName": experimentName}
    with open(os.path.join(experiment_path, 'runtimeExperimentData','Datasources_Configurations.json'), "w") as f:
        json.dump(config, f, indent=2)
    logger.debug(f" saved Datasources_Configurations json")
    if argos_zip:
        shutil.copy(argos_zip, os.path.join(experiment_path, 'runtimeExperimentData',f'{experimentName}.zip'))

def load_experiment_to_project(arguments):
    """
    Load an experiment's repository into a project's data layer.

    Parameters
    ----------
    arguments : argparse.Namespace
        Must contain projectName, experiment (path), and overwrite flag.
    """
    logger = logging.getLogger("hera.bin.experiment_load_experiment_to_project")
    logger.debug(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")
    if arguments.experiment:
        experiment_path = arguments.experiment
    else:
        experiment_path = os.getcwd()

    repository = glob.glob(os.path.join(f"{experiment_path}","*_repository.json"))
    if len(repository)==0:
        raise ValueError(f"Can't find repository file in path directory: {experiment_path}. \n Make sure the path is an experiment directory")
    if len(repository)>1:
        raise ValueError(f" More than 1 repositories found in directory.")

    repository = repository[0]
    repository_name = os.path.split(repository)[-1]

    if arguments.projectName not in datalayer.getProjectList():
        logger.info(f" No project with name {arguments.projectName}, will create a new one.")

    data_tk = dataToolkit()
    data_tk.addRepository(repositoryName=repository_name,
                      repositoryPath=repository,
                      overwrite=False)
    try:
        data_tk.loadAllDatasourcesInRepositoryToProject(arguments.projectName,repositoryName=repository_name, overwrite=arguments.overwrite)
    except ValueError as e:
        logger.error(f"Couldn't load Repository. Error message: {e}")

    data_tk.deleteDataSource(datasourceName=repository_name)
