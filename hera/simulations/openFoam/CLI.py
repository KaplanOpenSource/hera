import logging
import glob
import os
from ... import toolkitHome
from ...utils.jsonutils import loadJSON

def stochasticLagrangian_makeDispersionFlow(arguments):
    logger = logging.getLogger("hera.bin")
    logger.execution(f"----- Start -----")
    logger.debug(f" makeDispersionFlow with arguments: {arguments}")

    if 'projectName' not in arguments:
        configurationFile = arguments.configurationFile if 'configurationFile'  in arguments else "caseConfiguration.json"

        configuration = loadJSON(configurationFile)
        projectName = configuration['projectName']
    else:
        projectName = arguments.projectName

    logger.info(f"Adding dispersion flow to project {projectName}")
    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    for flowid,flowdata in enumerate(configuration["flowFields"]["Flows"]):
        logger.debug(f'Processing flow {flowid} of {len(configuration["flowFields"]["Flows"])}')
        tk.stochasticLagrangian.createDispersionFlowField(flowData=flowdata)

    logger.execution(f"----- End -----")

def stochasticLagrangian_createDispersionCaseDirectory(arguments):
    """
        Prepares the dispersion case:

         * Copies a template to the new directory,
         * Sets up references to the root directory.

         Assumes the root directory is parallel.


    :param arguments:
             - dispersionCaseDirectory : the new case of the directory
             - dispersionFlow : the new case of the directory.
                                can be either a directory or a name of a dispersion flow field from the DB

    :return:
    """
    logger = logging.getLogger("hera.bin")
    logger.execution(f"----- Start -----")
    logger.debug(f" prepareDispersionCase with arguments: {arguments}")

    dispersionDirectory     = os.path.abspath(arguments.dispersionCaseDirectory)
    dispersionFlow           = os.path.abspath(arguments.dispersionFlowField)

    if 'projectName' not in arguments:
        configuration = loadJSON("caseConfiguration.json")
        projectName = configuration['projectName']
    else:
        projectName = arguments.projectName

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    tk.stochasticLagrangian.createDispersionCaseDirectory(dispersionDirectory,dispersionFlow=dispersionFlow)

