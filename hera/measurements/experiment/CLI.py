from ...uti

def experiments_list(arguments):
    logger = logging.getLogger("hera.bin")
    logger.execution(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")

    if 'projectName' not in arguments:
        configurationFile = arguments.configurationFile if 'configurationFile'  in arguments else "caseConfiguration.json"

        configuration = loadJSON(configurationFile)
        projectName = configuration['projectName']
    else:
        projectName = arguments.projectName
