import logging
from ... import toolkitHome
from ...utils.jsonutils import loadJSON

def stochasticLagrangian_makeDispersionFlow(arguments):
    logger = logging.getLogger("hera.bin")
    logger.info(f"----- stochasticLagrangian_makeDispersionFlow : Start")
    logger.debug(" makeDispersionFlow with argumets: {arguments}")

    configurationFile = arguments.configurationFile if 'configurateionFile'  in arguments else "caseConfiguration.json"

    configuration = loadJSON(configurationFile)
    projectName = configuration['projectName']
    logger.info(f"Adding dispersion flow to project {projectName}")
    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    for flowid,flowdata in enumerate(configuration["flowFields"]["Flows"]):
        logger.debug(f'Processing flow {flowid} of {len(configuration["flowFields"]["Flows"])}')
        case = tk.createDispersionFlowField(flowData=flowdata)
