import os
import logging
from ...utils.jsonutils import loadJSON
from ... import toolkitHome
from argos.experimentSetup.dataObjects import ExperimentZipFile
import pandas
import json
import shutil

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

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.EXPERIMENT, projectName=projectName)
    print(tk.getExperiment(arguments.experiment).getExperimentData().getData(arguments.deviceType, deviceName=arguments.deviceName ,perDevice=arguments.perDevice))


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
        deviceTable = metadata.trialSet['Measurements']['Measurement'].entitiesTable("design")

        repo = {}
        perDevice = False         ##Will be defined by the updated zip format!

        repo['experiment'] = {}
        repo['experiment']['DataSource'] = {}
        repo['experiment']['DataSource'][arguments.experimentName] = {"isRelativePath": "True",
                                                                           "item":{
                                                                               "dataSourceName": arguments.experimentName,
                                                                               "resource": "",
                                                                               "experimentPath": experiment_path,
                                                                               "dataFormat": "datatypes.PARQUET",
                                                                               "overwrite": "True"
                                                                           }
                                                                      }

        repo['experiment']['Measurements'] = {}
        for deviceType in list(set(deviceTable.entityType)):
            if not perDevice:
                repo['experiment']['Measurements'][deviceType] = {"isRelativePath": "True",
                                                                  "item": {
                                                                      "type": "Experiment_rawData",
                                                                      "resource": os.path.join(experiment_path,'data',deviceType),
                                                                      "dataFormat": "datatypes.PARQUET",
                                                                      "desc": {
                                                                          "deviceType": deviceType,
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

def registerInProject(projectName, experimentName, experimentPath):

    dataSourceDesc = dict()
    dataSourceDesc['handlerPath'] = os.path.abspath(experimentPath)

    dataSourceDesc['className'] = experimentName + 'Toolkit'

    toolkit = toolkitHome.getToolkit(toolkitName=toolkitHome.EXPERIMENT, projectName=projectName)

    datasource = toolkit.getDataSourceDocument(datasourceName=experimentName)

    if datasource is not None:
        toolkit.addDataSource(dataSourceName=experimentName, resource="t", dataFormat='str', **dataSourceDesc)
        print(f"Added source {experimentName} to tool  in project {projectName}")
    else:
        toolkit.deleteDataSourceDocuments(datasourceName=experimentName)
        toolkit.addDataSource(dataSourceName=experimentName, resource="t", dataFormat='str', **dataSourceDesc)
        print(f"Source {experimentName} already exists in {projectName}, delete current source")
