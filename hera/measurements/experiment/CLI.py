import os
import logging
from ...utils.jsonutils import loadJSON
from ... import toolkitHome


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

    os.makedirs(str(os.path.join(experiment_path,'code')),exist_ok=True)


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
