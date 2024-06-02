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
        with open(f"{os.path.join(experiment_path, 'code', arguments.experimentName)}.py", "w") as class_script:
            class_script.write(f"from hera.measurements.experiment.experiment import experimentSetupWithData")
            class_script.write("\n\n\n")
            class_script.write(f"class {arguments.experimentName}(experimentSetupWithData):")
            class_script.write("\n")
            class_script.write("\t###Implement your code here if you wish.\n")
            class_script.write("\tpass")
        logger.debug(f" finished creating an empty class for implementation..")

    def create_repository():
        logger.debug(f" Since zip file is provided, creating a repository..")
        metadata = ExperimentZipFile(arguments.zip)

        repo = {}
        perDevice = True         ##Will be defined by the updated zip format!

        repo['experiment'] = {}
        repo['experiment']['DataSource'] = {}
        repo['experiment']['DataSource'][arguments.experimentName] = {"isRelativePath": "False",
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
        for i,entity in enumerate(entities_dict_list):
            if not perDevice:
                parquet_name = entity['entityTypeName']
            else:
                parquet_name = entity['entityName']

            if parquet_name not in repo['experiment']['Measurements'].keys():
                repo['experiment']['Measurements'][parquet_name] = {"isRelativePath": "True",
                                                                      "item": {
                                                                      "type": "Experiment_rawData",
                                                                      "resource": os.path.join('data',f"{parquet_name}.parquet"),
                                                                      "dataFormat": "string",
                                                                      "desc": {
                                                                              "deviceType": entity['entityTypeName'],
                                                                              "experimentName": arguments.experimentName,
                                                                              "deviceName": entity['entityName']
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

def create_ims_experiment(arguments):
    logger = logging.getLogger("hera.bin.create_ims_experiment")
    logger.execution(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")
    if arguments.path:
        experiment_path = arguments.path
    else:
        experiment_path = os.getcwd()

    if 'token.json' not in os.listdir(experiment_path):
        raise ValueError(f"No token file found in directory {experiment_path}. For creating an IMS experiment, 'token.json' file must be in the directory.")

    f = open(f'{experiment_path}/token.json')
    token = json.load(f)

    def create_zip_for_ims():
        logger = logging.getLogger("hera.bin.create_zip_for_ims")
        logger.execution(f"----- Start -----")
        logger.debug(f" Reaching URL...")

        url = "https://api.ims.gov.il/v1/Envista/stations"
        headers = {
            "Authorization": token['Authorization']
        }
        try:
            response = requests.request("GET", url, headers=headers)
            data = json.loads(response.text.encode('utf8'))
        except ValueError as e:
            logger.error(f"Couldn't reach URL. Error message: {e}")

        logger.debug(f" Downloaded metadata from URL succecfully.")
        logger.debug(f" Will create ZIP file now..")

        zip_file = {}
        zip_file['version'] = '3.0.0'
        zip_file['name'] = 'IMS_experiment'
        zip_file['startDate'] = ''
        zip_file['endDate'] = ''
        zip_file['description'] = ''
        zip_file['trialTypes'] = []
        zip_file['imageStandalone'] = []
        zip_file['imageEmbedded'] = []

        stations = []

        for station_dict in data:
            if station_dict['active']:
                device = {}
                device['name'] = station_dict['name']
                device['description'] = ''
                device['attributes'] = [{'name': 'stationId', 'value': station_dict['stationId']},
                                        {'name': 'stationsTag', 'value': station_dict['stationsTag']},
                                        {'name': 'location', 'value': station_dict['location']},
                                        {'name': 'timebase', 'value': station_dict['timebase']},
                                        {'name': 'active', 'value': station_dict['active']},
                                        {'name': 'owner', 'value': station_dict['owner']},
                                        {'name': 'regionId', 'value': station_dict['regionId']},
                                        {'name': 'monitors', 'value': station_dict['monitors']}]

                stations.append(device)

        zip_file['deviceTypes'] = [{"name": 'STATION',
                                    'attributeTypes': [],
                                    "devices": stations}]

        logger.debug(f" Finished fetching metadata. Will save zip now.")
        json_object = json.dumps(zip_file, indent=4)
        with open(f"{experiment_path}/data.json", "w") as outfile:
            outfile.write(json_object)

        with zipfile.ZipFile(f'{experiment_path}/IMS.zip', 'w', zipfile.ZIP_DEFLATED) as zipf:
            zipf.write(f'data.json')

        logger.debug(f" Finished creating zip file.")


    def write_download_function():
        logger = logging.getLogger("hera.bin.write_download_function")
        logger.execution(f"----- Start -----")
        logger.debug(f" Will write download function to experiment class...")

        with open(f"{experiment_path}/code/IMS_experiment.py", "a") as f:
            f.write('\n\n')
            f.write("\tdef download(self,start='',end='latest', concat=True):\n")
            f.write("\t\tfrom tqdm import tqdm\n")
            f.write("\t\timport pandas as pd\n")
            f.write("\t\timport time\n")
            f.write("\t\timport os\n")
            f.write("\t\timport json\n")
            f.write("\t\timport requests\n")
            f.write("\t\tfrom datetime import datetime,timedelta\n")
            f.write("\t\turl = 'https://api.ims.gov.il/v1/Envista/stations'\n")
            f.write(f"\t\tf = open('{experiment_path}/token.json')\n")
            f.write(f"\t\ttoken = json.load(f)\n")
            f.write("\t\theaders = {'Authorization': token['Authorization']}\n")
            f.write("\t\tresponse = requests.request('GET', url, headers=headers)\n")
            f.write("\t\tdata = json.loads(response.text.encode('utf8'))\n")
            f.write("\t\tfor station_dict in tqdm(data):\n")
            f.write("\t\t\tif station_dict['active']:\n")
            f.write("\t\t\t\tstation_id = station_dict['stationId']\n")
            f.write("\t\t\t\tchances = 10\n")
            f.write(f"\t\t\t\tpath_to_data = os.path.join('{experiment_path}','data',station_dict['name'])\n")
            f.write("\t\t\t\tif concat:\n")
            f.write("\t\t\t\t\ttry:\n")
            f.write(f"\t\t\t\t\t\tstart = pd.read_csv(path_to_data)['datetime'].iloc[-1].split('T')[0].replace('-','/')\n")
            f.write("\t\t\t\t\texcept:\n")
            f.write("\t\t\t\t\t\tstart = datetime.today().strftime('%Y-%m-%d')\n")
            f.write("\t\t\t\tif end == 'latest':\n")
            f.write("\t\t\t\t\tend = (datetime.today() + timedelta(days=1)).strftime('%Y-%m-%d')\n")
            f.write("\t\t\t\twhile(chances>0):\n")
            f.write("\t\t\t\t\ttry:\n")
            f.write("\t\t\t\t\t\turl = f'https://api.ims.gov.il/v1/envista/stations/{station_id}/data?from={start}&to={end}'\n")
            f.write("\t\t\t\t\t\tresponse = requests.request('GET', url, headers=headers)\n")
            f.write("\t\t\t\t\t\tstation_data = json.loads(response.text.encode('utf8'))\n")
            f.write("\t\t\t\t\t\tcsv = pd.DataFrame()\n")
            f.write("\t\t\t\t\t\tfor i,time_stamp_dict in enumerate(station_data['data']):\n")
            f.write("\t\t\t\t\t\t\tcsv.at[i,'datetime'] = time_stamp_dict['datetime']\n")
            f.write("\t\t\t\t\t\t\tfor column in time_stamp_dict['channels']:\n")
            f.write("\t\t\t\t\t\t\t\tcsv.at[i,column['name']] = column['value']\n")
            f.write(f"\t\t\t\t\t\tif concat:\n")
            f.write(f"\t\t\t\t\t\t\tpd.concat([pd.read_csv(path_to_data),csv],axis=0).to_csv(path_to_data,index=False)\n")
            f.write(f"\t\t\t\t\t\telse:\n")
            f.write(f"\t\t\t\t\t\t\tcsv.to_csv(path_to_data,index=False)\n")
            f.write("\t\t\t\t\t\tbreak\n"),
            f.write("\t\t\t\t\texcept:\n")
            f.write("\t\t\t\t\t\tchances -= 1\n")
            f.write("\t\t\t\tif chances==0:\n")
            f.write("\t\t\t\t\tif not os.path.isfile(path_to_data):\n")
            f.write(f"\t\t\t\t\t\tpd.DataFrame().to_csv(path_to_data,index=False)\n")
            f.close()



    create_zip_for_ims()
    arguments.zip  = f"{experiment_path}/IMS.zip"
    arguments.experimentName = 'IMS_experiment'
    create_experiment(arguments)
    write_download_function()


def download(arguments):
    logger = logging.getLogger("hera.bin.download_ims")
    logger.execution(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")
    experimentToolKit = toolkitHome.getToolkit(toolkitName=toolkitHome.EXPERIMENT, projectName=arguments.project)

    experimentToolKit.getExperiment('IMS_experiment').download(start=arguments.start,end=arguments.end,concat=arguments.concat)
