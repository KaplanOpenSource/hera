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

def experiments_list(arguments):
    logger = logging.getLogger("hera.bin.experiment_experiments_list")
    logger.debug(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")
    if 'projectName' not in arguments:
        configurationFile = arguments.configurationFile if 'configurationFile'  in arguments else "caseConfiguration.json"

        configuration = loadJSON(configurationFile)
        projectName = configuration['projectName']
    else:
        projectName = arguments.projectName

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.EXPERIMENT,projectName=projectName)
    print(tk.keys())


def experiments_table(arguments):
    logger = logging.getLogger("hera.bin.experiment_experiments_table")
    logger.debug(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")
    if 'projectName' not in arguments:
        configurationFile = arguments.configurationFile if 'configurationFile'  in arguments else "caseConfiguration.json"

        configuration = loadJSON(configurationFile)
        projectName = configuration['projectName']
    else:
        projectName = arguments.projectName

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.EXPERIMENT,projectName=projectName)
    print(tk.getExperimentsTable())


def get_experiment_data(arguments):
    logger = logging.getLogger("hera.bin.experiment_get_experiment_data")
    logger.debug(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")
    if 'projectName' not in arguments:
        configurationFile = arguments.configurationFile if 'configurationFile'  in arguments else "caseConfiguration.json"

        configuration = loadJSON(configurationFile)
        projectName = configuration['projectName']
    else:
        projectName = arguments.projectName

    if projectName not in datalayer.getProjectList():
        raise ValueError(f"Project '{projectName}' does not exists")

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.EXPERIMENT, projectName=projectName)
    parquet = tk.getExperiment(arguments.experiment).getExperimentData().getData(arguments.deviceType, deviceName=arguments.deviceName ,perDevice=arguments.perDevice)
    print(pandas.DataFrame(parquet))

def create_experiment(arguments):
    logger = logging.getLogger("hera.bin.experiment_create_experiment")
    logger.debug(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")
    if arguments.path:
        experiment_path = arguments.path
    else:
        experiment_path = os.getcwd()

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

    if arguments.zip:
        _create_repository(arguments.zip,experiment_path,arguments.experimentName,arguments.relative)
        _make_runtimeExperimentData(arguments.zip,experiment_path,arguments.experimentName)

    logger.debug("Loading the experiment repository to the project")
    arguments.repositoryName = os.path.join(experiment_path,f"{arguments.experimentName}_repository.json")
    projectCLI.repository_load(arguments)


def _create_empty_class(experiment_path,experimentName):
    logger = logging.getLogger("hera.bin._create_empty_class")
    logger.debug(f" creating an empty class for implementation..")
    with open(f"{os.path.join(experiment_path, 'code', experimentName)}.py", "w") as class_script:
        class_script.write(f"from hera.measurements.experiment.experiment import experimentSetupWithData")
        class_script.write("\n\n\n")
        class_script.write(f"class {experimentName}(experimentSetupWithData):")
        class_script.write("\n")
        class_script.write("\t###Implement your code here if you wish.\n")
        class_script.write("\tpass")
    logger.debug(f" finished creating an empty class for implementation..")



def _create_repository(zip,experiment_path,experimentName,relative):
    logger = logging.getLogger("hera.bin._create_repository")
    logger.info(f"Creating the repository")
    logger.debug(f" Since zip file is provided, creating a repository..")
    metadata = ExperimentZipFile(zip)

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

    logger.debug(f" finished creating the repository json file")

def _make_runtimeExperimentData(zip,experiment_path,experimentName):
    logger = logging.getLogger("hera.bin._make_runtimeExperimentData")
    logger.debug(f" creating runtimeExperimentData directory if it does not exists")
    os.makedirs(os.path.join(experiment_path, 'runtimeExperimentData'), exist_ok=True)
    logger.debug(f" creating Datasources_Configurations json")
    config = {"experimentName": experimentName}
    with open(os.path.join(experiment_path, 'runtimeExperimentData','Datasources_Configurations.json'), "w") as f:
        json.dump(config, f, indent=2)
    logger.debug(f" saved Datasources_Configurations json")
    shutil.copy(zip, os.path.join(experiment_path, 'runtimeExperimentData',f'{experimentName}.zip'))

def load_experiment_to_project(arguments):
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
