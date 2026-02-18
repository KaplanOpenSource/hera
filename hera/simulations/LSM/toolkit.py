import pandas
import json
import os

from hera.utils import ConfigurationToJSON, Quantity, slurm
from hera.utils.jsonutils import JSONVariations
from hera.utils.logging.helpers import get_classMethod_logger

from .template import LSMTemplate
from itertools import product
from ... import datalayer
from ... import toolkit
from .singleSimulation import SingleSimulation
from unum.units import *
from ..utils.coordinateHandler import coordinateHandler
from pathlib import Path


class LSMToolkit(toolkit.abstractToolkit):
    """
        The LSM.old toolkit

        The datasources are the templates.

        The LSM.old template JSON structure is :


        << TO DO >>



    """
    TRUE = ".TRUE."
    FALSE = ".FALSE."
    _to_xarray = None
    _to_database = None
    _forceKeep = None
    _analysis = None

    @property
    def analysis(self):
        return self._analysis
    @property
    def to_xarray(self):
        return self._to_xarray

    @to_xarray.setter
    def to_xarray(self, value):
        if not isinstance(value,bool):
            raise ValueError("to_xarray must be boolean")
        self._to_xarray = value


    @property
    def to_database(self):
        return self._to_database

    @to_database.setter
    def to_database(self, value):
        if not isinstance(value,bool):
            raise ValueError("to_xarray must be boolean")
        self._to_database = value

    @property
    def forceKeep(self):
        return self._forceKeep

    @forceKeep.setter
    def forceKeep(self, value):
        if not isinstance(value, bool):
            raise ValueError("to_xarray must be boolean")
        self._forceKeep = value

    @property
    def singleSimulation(self):
        return SingleSimulation

    def __init__(self, projectName, filesDirectory=None, to_xarray=True, to_database=False, forceKeep=False,connectionName=None):
        """
            Initializes the LSM.old toolkit

        Parameters
        ----------

        projectName: str
            The name of the project that contains the

        filesDirectory: str
            The directory to save the simulations.old.


        to_xarray: bool
            Save the simulation results into xarray or not

        to_database: bool
            Save the simulation run in the database or not

        forceKeep: bool
            If to_xarray is true, determine wehter to keep the original files.
            if False, removes the Lagrnagian files.

        """
        super().__init__(projectName=projectName, toolkitName="LSM.old", filesDirectory=filesDirectory,connectionName=connectionName)

        self.to_xarray = to_xarray
        self.to_database = to_database
        self.forceKeep = forceKeep
        self._analysis = analysis(self)

    def getTemplates(self, **query):
        """
            Get a list of Template objects that fulfill the query

        Returns
        -------
            List of LSMTemplates

        """
        docList = self.getDataSourceDocumentsList(**query)
        return [LSMTemplate(doc,self) for doc in docList]

    def getTemplateByName(self,templateName,templateVersion=None):
        """
            Retrieve the template by its name.

        Parameters
        ----------

        templateName: str
                The name of the template.
        templateVersion: str
                The name of the template.
        Returns
        -------
            The template by the name
        """
        doc = self.getDataSourceDocument(datasourceName=templateName, version=templateVersion)
        return LSMTemplate(doc,self) if doc is not None else None

    def getTemplatesTable(self, **query):
        """
            :list the template parameters that fulfil the query
        :param query:
        :return:
        """
        docList = self.getDataSourceDocumentsList(**query)
        if len(docList) > 0:
            descList = [doc.desc.copy() for doc in docList]
            for (i, desc) in enumerate(descList):
                desc.update({'id':docList[i].id})
                desc.update({'projectName': docList[i].projectName})

            params_df_list = [pandas.DataFrame(desc.pop('params'), index=[0]) for desc in descList]
            desc_df_list = [pandas.DataFrame(desc, index=[0]) for desc in descList]
            df_list = [desc.join(params) for (desc,params) in product(desc_df_list, params_df_list)]
            ret = pandas.concat(df_list,ignore_index=True,sort=False)
        else:
            ret = pandas.DataFrame()
        return ret


    def loadData(self,fileNameOrData,saveMode=toolkit.TOOLKIT_SAVEMODE_FILEANDDB,**kwargs):
        """
            Load a template object. Possibly saves to the DB.

        Parameters
        ----------
        fileNameOrData: str or JSON str or a JSON object.

                If str , a file or an JSON str that describes a template.

        saveMode: str
                Can be either:

                    - TOOLKIT_SAVEMODE_NOSAVE   : Just create template object


                    - TOOLKIT_SAVEMODE_FILEANDDB : Creates the template and store to the DB as a source.
                                                    Raise exception if the entry exists.

                    - TOOLKIT_SAVEMODE_FILEANDDB_REPLACE: Creates the template and store to the DB as a source.
                                                    Replace the entry in the DB if it exists.

        Returns
        -------
            template.LSMTemplate object

            Return the template object.
        """

        if isinstance(fileNameOrData,str):
            if os.path.isfile(fileNameOrData):
                with open(fileNameOrData) as templateFile:
                    desc = json.load(templateFile)
            else:
                raise ValueError("fileNameOrData must be a JSON template file")
        else:
            raise ValueError("fileNameOrData must be a JSON template file")

        templateName = desc['name']
        version = kwargs.get("version",None)

        doc = None
        if saveMode in [toolkit.TOOLKIT_SAVEMODE_FILEANDDB,toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE]:

            doc = self.getDataSourceDocument(templateName)

            if doc is None:

                self.addDataSource(dataSourceName=templateName,
                                   resource=fileNameOrData,
                                   dataFormat=datalayer.datatypes.STRING,
                                   version=version,
                                   **kwargs)

            else:
                if  (saveMode == toolkit.TOOLKIT_SAVEMODE_FILEANDDB):
                    raise ValueError(f"Template {templateName} already exists in the DB")

                doc.resource = fileNameOrData
                doc.desc = kwargs
                doc.desc['version'] = version

        if doc is None:
            doc = datalayer.nonDBMetadataFrame(data=fileNameOrData,**kwargs)

        return LSMTemplate(doc,self)

    def getSimulations(self,simulationName=None,unitsTemplateVersion="v4-general", **query):
        """
        get a list of SingleSimulation objects that fulfill the query
        :param query:
        :return:
        """
        template = self.getTemplateByName(unitsTemplateVersion)

        if 'units' in template._document['desc']:
            for key in template._document['desc']["units"].keys():
                if key in query.keys():
                    query_item= query[key]
                    if isinstance(query_item, Unum):
                        query[key] = query_item.asNumber(eval(template._document['desc']["units"][key]))
                    elif isinstance(query_item, Quantity):
                        query[key] = query_item.m_as(template._document['desc']["units"][key])
                    else:
                        raise ValueError(f"query must use either pint or unum to specify units, currently type({query[key]})={type(query[key])}")
        else:
            query = ConfigurationToJSON(query)
        print(query)
        queryWithParams = {}
        for key in query.keys():
            queryWithParams[f"params__{key}"] = query[key]

        if simulationName is not None:
            queryWithParams["simulationName"] = simulationName

        docList = self.getSimulationsDocuments(type=LSMTemplate("",self).doctype_simulation, **queryWithParams)
        retList = []
        for doc in docList:
            try:
                retList.append(SingleSimulation(doc))
            except:
                print(f"Warning: could not find data with the following document: {doc.asDict()}")

        return retList

    def getSimulationsList(self, wideFormat=False, **query):
        """
            List the Simulation parameters that fulfil the query
        :param query:
        :return:
        """
        docList = self.getSimulationsDocuments(type=LSMTemplate("",self).doctype_simulation, **query)
        descList = [doc.desc.copy() for doc in docList]
        for (i, desc) in enumerate(descList):
            desc.update({'id':docList[i].id})
        params_df_list = [pandas.DataFrame(desc.pop('params'), index=[0]) for desc in descList]
        params_df_list = [df.rename(columns=dict([(x,"params__%s"%x) for x in df.columns])) for df in params_df_list]
        desc_df_list = [pandas.DataFrame(desc, index=[0]) for desc in descList]
        df_list = [desc.join(params) for (desc,params) in product(desc_df_list, params_df_list)]
        new_df_list = []
        for df in df_list:
            id = df['id'][0]
            new_df = df.copy().drop(columns=['id']).melt()
            new_df.index = [id]*len(new_df)
            new_df_list.append(new_df)
        try:
            df = pandas.concat(new_df_list)
            if wideFormat:
                return df.pivot(columns='variable', values='value')
            else:
                return df
        except ValueError:
            raise FileNotFoundError('No simulations.old found')
        

    def prepareSlurmLSMExecution(self,baseParameters,
                              jsonVariations,
                              templateName,
                              stations,
                              topography,
                              depositionRates,
                              canopy,
                              saveMode,
                              slurmExecutionFileName="submit_all.sh",
                              caseListFileName="cases.txt",
                              allocateProcessorsPerRun=None,
                              memoryInGB=None,
                              jobName="lsm_cases",
                              exclusive=False):
        """
            Creates a slurm setup to run many LSM simulations against specific variation

        Parameters
        ----------
        baseParameters : dict
                base parameters
        jsonVariations :
                Variation file (using the jsonutils variations) format to apply on the base paramaters
        slurmExecutionFileName: str
                The name of the bash file to create with the slurm batch run
        caseListFileName:
                The batchfile uses case file name, so add it.
        allocateProcessorsPerRun : int | None
                How many nodes(currently just threads) are used per job, in slurm documentation they describe "This option advises ... that job steps run ... will launch a maximum of number tasks and to provide for sufficient resources."
        memoryInGB : str | int | none
                Should memory be limited, affect depends on configuration(might kill if memory use is exceeded)
        jobName : bool
                name for slurm job
        exclusive : bool
                Should slurm run one job at a time on a GRES(Generic RESource in our case CPUs)
        Returns
        -------

        """

        logger = get_classMethod_logger(self,"prepareSlurmWorkflowExecution")
        simList = ""
        if not isinstance(baseParameters, dict):
            logger.error("Slurm preparation can only handle dictionary of parameters to run variations")

        if isinstance(jsonVariations, str):
            logger.info(f"Assuming {jsonVariations} is path to variations file")
            with open(jsonVariations, 'r') as variationsFile:
                jsonVariations = json.load(variationsFile)
        elif not isinstance(jsonVariations, dict):
            logger.error("Slurm preparation only supports json variation input as path or dict")

        SIMULATIONS_SCRIPT_DIR_NAME= "simulationsScripts"
        STATIONS_PATH = "stations.parquet"
        RUN_SIM_FILE_NAME = "run_sim.py"

        simulations_scripts_dir = Path(self.filesDirectory,SIMULATIONS_SCRIPT_DIR_NAME)
        os.makedirs(simulations_scripts_dir, exists_ok=True)
        if stations is not None:
            stations.to_parquet(simulations_scripts_dir / STATIONS_PATH)
        for i, jsonConfig in enumerate(JSONVariations(baseParameters, jsonVariations)):
            simName = f"LSM_Simulation_{i}"
            simSavePath = simulations_scripts_dir/ simName/ RUN_SIM_FILE_NAME

            os.makedirs(simulations_scripts_dir / simName, exists_ok=True)
            if isinstance(topography, str):
                topography = f'"{topography}"'
            if isinstance(canopy, str):
                canopy = f'"{canopy}"'
            if isinstance(topography, str):
                topography = f'"{topography}"'
            read_stations_line = f'stations = pandas.read_parquet("'+STATIONS_PATH+'")'

            sim_script = f"""
from hera import toolkitHome, ToolkitHome
from hera.simulations.LSM.toolkit import LSMToolkit
import pandas
old_lsm_toolkit = toolkitHome.getToolkit(toolkitName=ToolkitHome.LSM, projectName="{self.projectName}")

lsm_template = old_lsm_toolkit.getTemplates(template="{templateName}")[0]
params={jsonConfig}
{read_stations_line if stations is not None else ""}
lsm_template.run(topography={topography}, stations={None if stations is None else "stations"},canopy={canopy},depositionRates={depositionRates}, saveMode="{saveMode}",simulationName="{simName}",**params)
"""
            with open(simSavePath, "w") as f:
                f.write(sim_script)
                
            simList += f"{simName}\n"


        script =f"""
cd "$dir"
python {RUN_SIM_FILE_NAME}
        """
        caseListFilePath = simulations_scripts_dir/ caseListFileName

        logger.info(f"Writing simulation list file")
        
        with open(caseListFilePath,"w") as outputFile:
            outputFile.write(simList)
        slurm.prepareSlurmScriptExecution(script=script,
                                          slurmExecutionFilePath=simulations_scripts_dir/slurmExecutionFileName,
                                          jobDirListFilePath=caseListFilePath,
                                          allocateProcessorsPerRun=allocateProcessorsPerRun,
                                          memoryInGB=memoryInGB,
                                          jobName=jobName,
                                          quiet=False,
                                          exclusive=exclusive)



class analysis:
    """
        Analysis of the demography toolkit.
    """

    _datalayer = None
    _coordinateHandler = None

    @property
    def coordinateHandler(self):
        return self._coordinateHandler

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, dataLayer):
        self._datalayer = dataLayer
        self._coordinateHandler = coordinateHandler()
