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

def experiments_list(arguments):
    logger = logging.getLogger("hera.bin.experiment_experiments_list")
    logger.execution(f"----- Start -----")
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
    logger.execution(f"----- Start -----")
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
    logger.execution(f"----- Start -----")
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
    logger.execution(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")
    if arguments.path:
        experiment_path = arguments.path
    else:
        experiment_path = os.getcwd()
    def create_empty_class():
        logger.debug(f" creating an empty class for implementation..")
        class_script = open(f"{os.path.join(experiment_path, 'code', arguments.experimentName)}.py", "w")
        class_script.write(f"from hera.measurements.experiment.experiment import experimentSetupWithData")
        class_script.write("\n\n\n")
        class_script.write(f"class {arguments.experimentName}(experimentSetupWithData):")
        class_script.write("\n")
        class_script.write("\t###Implement your code here if you wish.\n")
        class_script.write("\tpass")
        class_script.close()
        logger.debug(f" finished creating an empty class for implementation..")

    def create_repository():
        logger.debug(f" Since zip file is provided, creating a repository..")
        metadata = ExperimentZipFile(arguments.zip)

        repo = {}
        perDevice = False         ##Will be defined by the updated zip format!

        repo['experiment'] = {}
        repo['experiment']['DataSource'] = {}
        repo['experiment']['DataSource'][arguments.experimentName] = {"isRelativePath": "True",
                                                                           "item":{
                                                                               "dataSourceName": arguments.experimentName,
                                                                               "resource": "",
                                                                               "experimentPath": experiment_path,
                                                                               "dataFormat": "parquet",
                                                                               "overwrite": "True"
                                                                           }

                                                                      }

        repo['experiment']['Measurements'] = {}

        entities_dict_list = metadata.getExperimentEntities()

        for entity in entities_dict_list:
            if 'Station' != entity['entityTypeName']:
                if not perDevice:
                    parquet_name = entity['entityTypeName']
                else:
                    parquet_name = entity['entityName']

                if parquet_name not in repo['experiment']['Measurements'].keys():
                        repo['experiment']['Measurements'][entity['entityTypeName']] = {"isRelativePath": "True",
                                                                      "item": {
                                                                          "type": "Experiment_rawData",
                                                                          "resource": os.path.join(experiment_path,'data',parquet_name),
                                                                          "dataFormat": "parquet",
                                                                          "desc": {
                                                                              "deviceType": entity['entityTypeName'],
                                                                              "experimentName": arguments.experimentName,
                                                                            }
                                                                          }
                                                                      }

        with open(os.path.join(experiment_path,f'{arguments.experimentName}_repository.json'), "w") as f:
            json.dump(repo, f, indent=4)

        logger.debug(f" finished creating the repository json file")

    def make_runtimeExperimentData():
        logger.debug(f" creating runtimeExperimentData directory if it does not exists")
        os.makedirs(os.path.join(experiment_path, 'runtimeExperimentData'), exist_ok=True)
        logger.debug(f" creating Datasources_Configurations json")
        config = {"experimentName": arguments.experimentName}
        with open(os.path.join(experiment_path, 'runtimeExperimentData','Datasources_Configurations.json'), "w") as f:
            json.dump(config, f, indent=2)
        logger.debug(f" saved Datasources_Configurations json")
        shutil.copy(arguments.zip, os.path.join(experiment_path, 'runtimeExperimentData',f'{arguments.experimentName}.zip'))


    logger.debug(f" creating code directory if not exist")
    os.makedirs(os.path.join(experiment_path,'code'),exist_ok=True)

    logger.debug(f" checking if class script already exists..")
    try:
        class_script = open(f"{os.path.join(experiment_path, 'code', arguments.experimentName)}.py", "x")
        logger.debug(f" creating the experiment class since it does not exists")
        create_empty_class()
    except:
        logger.debug(f" experiment class already exists ")

    logger.debug(f" creating data directory if not exist")
    os.makedirs(os.path.join(experiment_path, 'data'), exist_ok=True)

    if arguments.zip:
        create_repository()
        make_runtimeExperimentData()

def load_experiment_to_project(arguments):
    logger = logging.getLogger("hera.bin.experiment_load_experiment_to_project")
    logger.execution(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")
    if arguments.experiment:
        experiment_path = arguments.experiment
    else:
        experiment_path = os.getcwd()

    repository = glob.glob(f"{experiment_path}/*_repository.json")
    if len(repository)==0:
        raise ValueError(f"Can't find repository file in path directory: {experiment_path}. \n Make sure the path is an experiment directory")
    if len(repository)>1:
        raise ValueError(f" More than 1 repositories found in directory.")

    repository = repository[0]
    repository_name = repository.split("/")[-1]

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

# def registerInProject(projectName, experimentName, experimentPath):
#
#     dataSourceDesc = dict()
#     dataSourceDesc['handlerPath'] = os.path.abspath(experimentPath)
#
#     dataSourceDesc['className'] = experimentName + 'Toolkit'
#
#     toolkit = toolkitHome.getToolkit(toolkitName=toolkitHome.EXPERIMENT, projectName=projectName)
#
#     datasource = toolkit.getDataSourceDocument(datasourceName=experimentName)
#
#     if datasource is not None:
#         toolkit.addDataSource(dataSourceName=experimentName, resource="t", dataFormat='str', **dataSourceDesc)
#         print(f"Added source {experimentName} to tool  in project {projectName}")
#     else:
#         toolkit.deleteDataSourceDocuments(datasourceName=experimentName)
#         toolkit.addDataSource(dataSourceName=experimentName, resource="t", dataFormat='str', **dataSourceDesc)
#         print(f"Source {experimentName} already exists in {projectName}, delete current source")
