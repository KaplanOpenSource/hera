import logging
import json
import os

def createRepository(arguments):
    """
    Create a risk assessment repository JSON from agent JSON files in a directory.

    Scans the given path for JSON files that contain agent descriptors
    (identified by having an ``effectParameters`` key) and generates a
    repository JSON file mapping each agent as a data source.

    Parameters
    ----------
    arguments : argparse.Namespace
        CLI arguments with ``path`` and ``repository_name`` attributes.
    """
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

    with open(os.path.join(arguments.path,f"{arguments.repository_name}"), "w") as f:
        json.dump(repo, f, indent=4)

    logger.debug(f" Finished creating Repository")


def _check_if_agent(file_path: str):
    """
    Check if a file is a valid agent JSON descriptor.

    Parameters
    ----------
    file_path : str
        Path to the file to check.

    Returns
    -------
    dict or False
        The parsed agent description if valid, False otherwise.
    """
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
