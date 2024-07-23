import logging
import json
import os

def createRepository(arguments):
    logger = logging.getLogger("hera.CLI_riskassesment.create_repository")
    logger.debug(f" Creating a repository..")

    if not arguments.path:
        arguments.path = os.getcwd()

    repo = {}
    repo['RiskAssessment'] = {}
    repo['RiskAssessment']['DataSource'] = {}

    for file_path in os.listdir(arguments.path):
        agent_json_description = _check_if_agent(file_path)
        if agent_json_description:
            agent_name = file_path.replace(".json","").split("/")[-1]
            repo['RiskAssessment']['DataSource'][agent_name] = {"isRelativePath": "False",
                                                                       "item":{
                                                                           "dataSourceName": agent_name,
                                                                           "resource": agent_json_description,
                                                                           "dataFormat": "JSON_dict",
                                                                           "overwrite": "True"
                                                                            }
                                                                }

    with open(os.path.join(arguments.path,f"{arguments.repository_name}.json"), "w") as f:
        json.dump(repo, f, indent=4)

    logger.debug(f" Finished creating Repository")


def _check_if_agent(file_path:str):
    logger = logging.getLogger("hera.CLI_riskassesment._check_if_agent")
    logger.debug(f"Checking if {file_path} is Afgent file.")
    try:
        with open(file_path, 'r') as readFile:
            agentDescription = json.load(readFile)
    except Exception as e:
        logger.debug(f"Could not open {file_path}. Error is: {e}")
        return False

    if agentDescription.get('effectParameters'):
        return agentDescription
