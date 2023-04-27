import warnings
from enum import Enum, auto, unique
from typing import Union
import pandas
import shutil
import os
from ..toolkit import abstractToolkit
from ..utils import loadJSON,convertJSONtoPandas
from ..utils.query import dictToMongoQuery
from ..utils.dataframeutils import compareDataframeConfigurations
from ..datalayer import datatypes
import numpy
import pydoc
import warnings

try:
    from hermes import workflow
except ImportError:
#    raise ImportError("hermes is not installed. please install it to use the hermes workflow toolkit.")
    warnings.warn("hermes is not installed. some features will not work.")
    workflow = None


@unique
class actionModes(Enum):
    ADD = auto()
    ADDBUILD = auto()
    ADDBUILDEXECUTE = auto()

@unique
class simulationTypes(Enum):
    WORKFLOW = "hermes_workflow"
    OF_FLOWFIELD = "OF_FlowField"  # OpenFoam: calculation of the flow field.
    OF_DISPERSION = "OF_dispersion"  # OpenFoam: The dispersion itself.
    OF_FLOWDISPERSION = "OF_flowDispersion"

    @classmethod
    def value_of(cls, value):
        for k, v in cls.__members__.items():
            if value == v.value:
                return getattr(cls,k)
        else:
            raise ValueError(f"'{cls.__name__}' enum not found for '{value}'")

    @classmethod
    def isvalid(cls,value):
        for k, v in cls.__members__.items():
            if value == v.value:
                return True
        else:
            raise False


HERAMETADATA = 'heraMetaData'  # The name of the node in the workflow to which the hera meta data is added.

class workflowToolkit(abstractToolkit):
    """
        Manages the hermes worflows:

            1. Checks if they are in the DB.
            2. create a new name to them
            3. allows simple deletion
            4. allows simple comparison.
            5. retrieve by initial name.
    """
    DESC_GROUPNAME = "groupName"
    DESC_GROUPID   = "groupID"
    DESC_SIMULATIONNAME = "simulationName"
    DESC_PARAMETERS = "parameters"

    def __init__(self, projectName: str, filesDirectory: str = None,toolkitName : str="hermesWorkflowToolkit"):
        """
            Initializes the workflow toolkit.

        Parameters
        ----------
        projectName: str
            The project that the workflow will be used i.

        filesDirectory : str
            The directory to write all the workflows and the outputs. default is current directory.
        """
        super().__init__(projectName=projectName,
                         filesDirectory=filesDirectory,
                         toolkitName =toolkitName)

        ## Create the simulationType->object map

        self._simulationTypeMap = {
                        simulationTypes.WORKFLOW : "hermes.workflow",
                        simulationTypes.OF_DISPERSION : "hera.simulations.old.openFoam.datalayer.hermesWorkflow.Workflow_Dispersion",
                        simulationTypes.OF_FLOWFIELD: "hera.simulations.old.openFoam.datalayer.hermesWorkflow.Workflow_Flow",
        }


    def getHermesWorkflowFromJSON(self,workflow : Union[dict,str],simulationType : simulationTypes):
        """
            Creates a hermes workflow object from the JSON that is supplied.

            The JSON can be either file name, JSON string or a dictionary.

        Parameters
        ----------
        workflow: dict,str
            The

        simulationType : str
            The type of the workflow to create.

        Returns
        -------
            hermesWorkflow object.
        """
        workFlowJSON = loadJSON(workflow)

        if simulationType not in [x.value for x in self._simulationTypeMap]:
            raise ValueError(f"Simulation type {simulationType} not found. Must be one of : {[x for x in simulationTypes]}")
        ky = simulationTypes.value_of(simulationType)

        hermesWFObj = pydoc.locate(self._simulationTypeMap[ky])

        if hermesWFObj is None:
            raise ValueError(f"Object {self._simulationTypeMap[simulationType]} not found....")


        return hermesWFObj(workFlowJSON)


    def getHermesWorkflowFromDB(self,workflow : Union[dict,str],workflowType : simulationTypes = None):
        """
                Retrieve the hermes object from the DB.

                If the workflow is string, use it as a name. If the workflow is dict,
                use it as a filter on the paramters

        Parameters
        ----------
        workflow: str, dict
            The filtering criteria. Either name, or the parameters of the flow.

        Returns
        -------
            hermes workflow.
        """
        if isinstance(workflow,str):
            # try to find it as a name
            self.logger.debug(f"Trying to get the workflow as a name.")
            docList = self.getSimulationsDocuments(simulationName=workflow)
            if len(docList) == 0:
                self.logger.debug(f"... not found. Try to query as a json. ")
                jsn = loadJSON(workflow)
                currentQuery = dictToMongoQuery(jsn, prefix="parameters")
                docList = self.getSimulationsDocuments(**currentQuery)

        elif isinstance(workflow,dict):
            currentQuery = dictToMongoQuery(workflow, prefix="parameters")
            docList = self.getSimulationsDocuments(**currentQuery)

        if len(docList) == 0:
            self.logger.error(f"... not found. ")
            raise ValueError(f"Simulation not found.")

        if len(docList) > 0:
            self.logger.warning("Got more than 1 simulation. Using the first onle only.")
            warnings.warn("Got more than 1 simulation. Using the first onle only.")

        doc = docList[0]
        workflowType = simulationTypes.WORKFLOW.value if workflowType is None else workflowType.value
        return self.getHermesWorkflowFromJSON(doc.desc['workflow'],simulationType=workflowType)

    def getSimulationDocumentFromDB(self, nameOrWorkflowFileOrJSONOrResource : Union[dict, str]):
        """
            Returns the simulation document from the DB.

            Identify the simulation from :
             - Resource (thedirectory name)
             - Simulation name
             - Its workflow
             - workfolow dict.

            Return the first item that was found.

        Parameters
        ----------
        nameOrWorkflowFileOrJSONOrResource: str, dict

        Can be
             - Resource (thedirectory name)
             - Simulation name
             - Its workflow
             - workfolow dict.

        Returns
        -------
            A document, or None if not found. .
        """
        docList = []
        if isinstance(nameOrWorkflowFileOrJSONOrResource, str):
            # try to find it as a name
            self.logger.debug(f"Searching for {nameOrWorkflowFileOrJSONOrResource} as a name.")
            docList = self.getSimulationsDocuments(simulationName=nameOrWorkflowFileOrJSONOrResource,type=simulationTypes.WORKFLOW.value)
            if len(docList) == 0:
                self.logger.debug(f"Searching for {nameOrWorkflowFileOrJSONOrResource} as a resource.")
                docList = self.getSimulationsDocuments(resource=nameOrWorkflowFileOrJSONOrResource,type=simulationTypes.WORKFLOW.value)
                if len(docList) == 0:
                    self.logger.debug(f"... not found. Try to query as a json. ")
                    try:
                        jsn = loadJSON(nameOrWorkflowFileOrJSONOrResource)
                        wf = self.getHermesWorkflowFromJSON(jsn ,simulationTypes.WORKFLOW.value)
                        currentQuery = dictToMongoQuery(wf.parametersJSON, prefix="parameters")
                        docList = self.getSimulationsDocuments(**currentQuery)
                    except ValueError:
                        return None
            else:
                self.logger.debug(f"... Found it ")

        elif isinstance(nameOrWorkflowFileOrJSONOrResource, dict):
            currentQuery = dictToMongoQuery(nameOrWorkflowFileOrJSONOrResource, prefix="parameters")
            docList = self.getSimulationsDocuments(**currentQuery,type=simulationTypes.WORKFLOW.value)


        return docList[0] if len(docList) >0 else None


    def getSimulationsInGroup(self, simulationGroup: str, **kwargs):
        """
            Return a list of all the simulations.old with the name as a prefic, and of the requested simuationType.
            Returns the list of the documents.

            If the simuationType is None use the default simuationType (WORKFLOW).

        Parameters
        ----------

        simulationGroup : str
            The prefix name of all the runs.

        simulationType : str [optional]
            The type of the workflow.
            if None, return all.

        kwargs: additional filtering criteria.
                Use mongodb criteria.

        Returns
        -------
            list of mongo documents.

        """
        kwargs['type'] = simulationTypes.WORKFLOW.value
        return self.getSimulationsDocuments(groupName=simulationGroup, **kwargs)

    def findAvailableName(self, simulationGroup: str, **kwargs):
        """
            Finds the next availabe name of that prefix. The available name is the maximal ID + 1.

            we assume that the name of each run is:
                <prefix name>_<id>.


        Parameters
        ----------
        simulationGroup : str
            The simulation group

        simulationType : str
            The type of the workflow.
            if None, use DOCTYPE_WORKFLOW

        **kwargs : dict
                additional fileters.
                Note that these should be the full path of the parameter in the JSON.

        Returns
        -------
            int,str
            The new ID,
            The name.
        """
        simList = self.getSimulationsInGroup(simulationGroup=simulationGroup, **kwargs)
        group_ids = [int(x['desc']['groupID']) for x in simList if x['desc']['groupID'] is not None]
        if len(group_ids) == 0:
            newID = 1
        else:
            newID = int(numpy.max(group_ids)+1)

        return newID, self.getworkFlowName(simulationGroup,newID)

    @staticmethod
    def getworkFlowName(baseName,flowID):
        """
            Returns the name of the flow field from the base and the
            flow id.

            The name is <base>_<id>
            where <id> is padded.

        Parameters
        ----------
        baseName : str
                The base name
        flowID: int
                the id of the name.

        Returns
        -------

        """
        formatted_number = "{0:04d}".format(flowID)
        return f"{baseName}_{formatted_number}"


    def addToGroup(self,
                   workflowJSON: str,
                   groupName: str = None,
                   overwrite: bool = False,
                   force: bool = False,
                   assignName: bool = False,
                   action: actionModes = actionModes.ADDBUILDEXECUTE,
                   parameters: dict = dict()):
        """
            1. Adds the workflow to the database in the requested group
            2. Builds the template (.json) and python executer
            3. Runs the workflow.

            The stages are executed according to the buildMode.

            Notes:

            * If the workflow is already in the db in a different name adds to the db only if **force** is True.

            * If the simulationName already exist in the group then overwrite its definitions
              only if the **overwrite** is True.

            * If the template and python execution files exist on the disk, raise error unless overwrite is True.

            * If the group is None, parse the file name to get the group. That is, we assume that the
                file name has the structure : <groupname>_<id>. If the <id> is not an integer,
                the id in the database will be saved as None.

        Parameters
        ----------
        workflowJSON : str
            The file name that contains the workflow.

        groupName : str
            The group to assign the workflow to.
            If None, parse the name under the format
                <groupname>_<id> to get the group and the ID.

        simulationType : simulationTypes
            The type of the workflow. [optional]
            If none, use the WORKFLOW type.

        overwrite : bool
            If true, update the record if it exists.
            If false, throw an exception if the record exists.

        assignName : bool
            If true, finds the next available id and saves it in the DB.
            If groupName is None, parse the filename to get the group.

            Note that if true and group name is None, the names of the simulation will use only the string before the '_'.

            Otherwise, use the filename as the name of the simulation.

        action : buildModes enum.
            ADDBUILDEXECUTE : add the simulation to the db, builds the execution files (also saves the new tempalte) and executes it.
                          if the file already in the db -> overwrite if over
            ADDBUILD    : add to the db and build the template and execution files.
            ADD         : add to the db.

        parameters : dict

                A dictionary with the parameters to override the default parameters of the workflow.
                The structure of the dict is :

                {
                    <node name> : {
                            "parameter path 1(eg. a.b.c)" : value,
                            "parameter path 2(eg. a.b.c)" : value
                            .
                            .
                            .
                    }
                }

        Returns
        -------
            None
        """
        if workflow is None:
            raise NotImplementedError("addToGroup() requires the 'hermes' library, which is nor installed")
        self.logger.info("-- Start --")

        # 1. Getting the names of the simulation and the groups.

        #    a. Make sure that there are no extensions.
        cleanName = workflowJSON.split(".")[0]
        theType = simulationTypes.WORKFLOW.value
        self.logger.debug(f"The suggested simulation name is {cleanName} in the document {theType}")

        #   b. loading the workflow.
        self.logger.debug(f"Loading the workflow JSON {workflowJSON}")
        hermesWF = workflow(loadJSON(workflowJSON), self.FilesDirectory)
        hermesWF.updateNodes(parameters=parameters)

        #   c. Determining the simulation name, group name and group id
        groupName = groupName if groupName is not None else cleanName.split("_")[0]
        if assignName:
            self.logger.debug("Generating ID from the DB")
            groupID, simulationName = self.findAvailableName(simulationGroup=groupName, simulationType=theType)
            self.logger.debug(f" Got id : {groupID} and suggested name {simulationName}")
        else:
            simulationName = cleanName
            try:
                groupID = cleanName.split("_")[1]
                if not groupID.isdigit():
                    groupID = None
            except IndexError:
                # The name has no _ in it...
                groupID = None
            self.logger.debug(f"Use input as simulation : {simulationName} with the group {groupID}")

        self.logger.info(
            f"Simulation name is {simulationName} with type {theType} in simulation group {groupName} with id {groupID}.")

        # 2. Check if exists in the DB.

        #    a. Check if the workflow exists in the DB under different name (assume similar type)
        self.logger.debug(f"Checking if the workflow already exists in the db unde the same group and type.")
        currentQuery = dictToMongoQuery(hermesWF.parametersJSON, prefix="parameters")

        docList = self.getSimulationsInGroup(simulationGroup=groupName, simulationType=theType, **currentQuery)

        if len(docList) > 0 and (not force) and (docList[0]['desc']['simulationName'] != simulationName):
            doc = docList[0]
            wrn = f"The requested workflow {simulationName} has similar parameters to the workflow **{doc['desc']['simulationName']}** in simulation group {groupName}."
            self.logger.warning(wrn)
            raise FileExistsError(wrn)
        else:

            #   b. Check if the name of the simulation already exists in the group
            self.logger.debug(f"Check if the name of the simulation {simulationName} already exists in the group")
            docList = self.getSimulationsInGroup(simulationGroup=groupName, simulationType=theType, simulationName=simulationName)

            if len(docList) == 0:
                self.logger.info("Simulation is not in the DB, adding... ")
                doc = self.addSimulationsDocument(resource=os.path.join(self.FilesDirectory, simulationName),
                                                  dataFormat=datatypes.STRING,
                                                  type=theType,
                                                  desc=dict(
                                                      groupName=groupName,
                                                      groupID=groupID,
                                                      simulationName=simulationName,
                                                      workflow=hermesWF.json,
                                                      parameters=hermesWF.parametersJSON)
                                                  )

            elif overwrite:
                self.logger.info("Simulation in the DB, overwrite=True.  Updating... ")
                doc = docList[0]
                doc['resource'] = hermesWF.json
                doc['desc']['parameters'] = hermesWF.parametersJSON
                doc.save()
            else:
                info = f"The simulation {simulationName} with type {theType} is already in the database in group {groupName}. use the overwrite=True to update the record."
                self.logger.info(info)

        # 3.  Building and running the workflow.
        if action.value > actionModes.ADD.value:
            self.logger.info(f"Building the workflow {simulationName}")
            build = hermesWF.build(buildername=workflow.BUILDER_LUIGI)

            self.logger.info(f"Writing the workflow and the executer python {simulationName}")
            wfFileName = os.path.join(self.FilesDirectory, f"{simulationName}.json")
            # attemp to write only if the file is different than the input (that is it exists) or if overwrite (which mean it needs to be updated).
            if wfFileName != os.path.join(self.FilesDirectory,workflowJSON) or overwrite:
                hermesWF.write(wfFileName,overwrite=overwrite)
            with open(os.path.join(self.FilesDirectory, f"{simulationName}.py"), "w") as file:
                file.write(build)

        if action.value > actionModes.ADDBUILD.value:
            self.logger.info(f"Run the workflow {simulationName}")

            if overwrite:
                # delete the run files if exist.
                executionfileDir = os.path.join(self.FilesDirectory, f"{simulationName}_targetFiles")
                shutil.rmtree(executionfileDir, ignore_errors=True)

            pythonPath = os.path.join(self.FilesDirectory, f"{simulationName}")
            executionStr = f"python3 -m luigi --module {os.path.basename(pythonPath)} finalnode_xx_0 --local-scheduler"
            self.logger.debug(executionStr)
            os.system(executionStr)

    def compare(self, workFlows: Union[list, str], docType: str = None, nodes: Union[list, str] = None,
                allParameters: bool = True, JSON: bool = False) -> Union[dict, pandas.DataFrame]:
        """
            Compares two or more simulations.old.


        Parameters
        ----------
        workFlows : str,list
                A single input uses it as a group name,
                a list is the list of cases to compare.

        docType : str
                The type of simulations.old. If not supplied use DOCTYPE_WORKFLOW

        nodes  : str,list
                The node name or a list of nodes to dislay

        allParameters: bool
                If true, then display all the parameters and not only those that are differ between two simulations.old.

        JSON: bool
                If true, return the results as a JSON and not pandas.DataFrame.

        Returns
        -------
            pandas.DataFrame, json (depends on the input flags).
            Return the differences between the parametrs of the requested cases.
        """
        if workflow is None:
            raise NotImplementedError("compare() requires the 'hermes' library, which is nor installed")
        self.logger.info("--- Start ---")

        # 1. Get all the simulations.old
        if isinstance(workFlows, Iterable):
            self.logger.debug("Workflow is iterable. Trying to retrieve the parameters for each item individually. ")

            workflowList = []

            for simulation in workFlows:
                if os.path.exists(simulation):
                    workflowJSON = loadJSON(simulation)
                else:
                    theType = self.DOCTYPE_WORKFLOW if docType is None else docType
                    docList = self.getSimulationDocuments(simulationName=simulation, type=theType)
                    if len(docList) == 0:
                        raise ValueError(f"Simulation {simulation} was not found on disk or in the project. ")
                    workflowJSON = docList[0]['desc']['workflow']

                workflowList.append(workflow(workflowJSON, WD_path=self.FilesDirectory))

        else:
            self.logger.debug("Workflow is a groupName. Get the simulations.old from the group")
            simulationList = self.getSimulationsInGroup(simulationGroup=workFlows, simulationType=docType)
            workflowList = [workflow(x['desc']['workflow'], WD_path=self.FilesDirectory) for x in simulationList]

        return workflowList

    def listSimulations(self,
                        simulationGroup:str,
                        parametersOfNodes:list = None,
                        allParams:bool = False,
                        longFormat:bool=False,
                        jsonFormat:bool = False
                        ) -> Union[pandas.DataFrame,dict]:
        """
            Lists all the simulations.old in the simulation group (of this project).

            Allows additional filters using the simulationType.

            If parameters is not None, return the list of parameters.
            return the parameters of all the nodes if the paraleters is an empty List, or the requested parameters.
            The default behaviour is to return only the parameters that are different from each other, unless allParams
            is True.

            The output is either pandas.DataFrame (if jsonFormat is False) or a JSON (if JSON is True).

        Parameters
        ----------
        simulationGroup : str
            The name of the group
        simulationType : str, optional
            Additional filter according to the simulation type.

        parametersOfNodes  : list[str]
            If None, just return the names of the simulations.old. Otherwise add the parameters from the requested nodes.

        allParams: bool
            If true, list all the parameters and not just the parameters that were different between the simulations.old.
        jsonFormat: bool
            If true, return JSON and not a normalized pandas.DataFrame.

        Returns
        -------
            pandas.DataFrame or dict
            A list of the simulations.old and their values.

        """
        simulationList = self.getSimulationsInGroup(simulationGroup=simulationGroup)
        if not simulationList:
            return {} if jsonFormat else pandas.DataFrame()
        if parametersOfNodes is not None:
            if workflow is None:
                raise NotImplementedError(
                    "listSimulations() with parametersOfNodes requires the 'hermes' library, which is nor installed"
                )
            qry = None
            if len(parametersOfNodes) > 0:
                qry = "nodeName in @parametersOfNodes"

            simParamList = []
            for simulationDoc in simulationList:
                wf = workflow(simulationDoc['desc']['workflow'])
                simulationParameters = convertJSONtoPandas(wf.parametersJSON).assign(simulationName=simulationDoc.desc[self.DESC_SIMULATIONNAME])
                if qry is not None:
                    simulationParameters = simulationParameters.query(qry)

                simParamList.append(simulationParameters)

            res = pandas.concat(simParamList)
            res = res.assign(nodeName=res.apply(lambda x: x.parameterName.split(".")[0],axis=1))
            if not allParams:
                ret =  compareDataframeConfigurations(res,
                                                      datasetName="simulationName",
                                                      parameterName="parameterName",
                                                      indexList="nodeName",
                                                      longFormat=longFormat)
            else:
                ret = res.pivot(index="simulationName",columns=["nodeName","parameterName"],values="value")
        else:
            ret = pandas.DataFrame([x.desc[self.DESC_SIMULATIONNAME] for x in simulationList],columns=[self.DESC_SIMULATIONNAME])

        if jsonFormat:
            ret = ret.to_json()

        return ret


