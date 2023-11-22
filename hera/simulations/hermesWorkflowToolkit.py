from enum import Enum, auto, unique
from typing import Union
import pandas
import shutil
import os
from ..toolkit import abstractToolkit
from ..utils import loadJSON,compareJSONS
from ..utils.query import dictToMongoQuery
from ..utils.dataframeutils import compareDataframeConfigurations
from ..datalayer import datatypes
import numpy
import pydoc
import warnings
from ..utils.logging import with_logger,get_classMethod_logger

try:
    from hermes import workflow
    from hermes.utils.workflowAssembly import handler_build,handler_buildExecute,handler_expand,handler_execute
except ImportError:
#    raise ImportError("hermes is not installed. please install it to use the hermes workflow toolkit.")
    warnings.warn("hermes is not installed. some features will not work.")
    workflow = None


@unique
class actionModes(Enum):
    ADD = auto()
    ADDBUILD = auto()
    ADDBUILDEXECUTE = auto()


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
    DESC_WORKFLOWNAME = "workflowName"
    DESC_PARAMETERS = "parameters"

    DOCTYPE_WORKFLOW = "hermesWorkflow"

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

        # ## Create the simulationType->object map
        # self._simulationTypeMap = {
        #                 workflowsTypes.WORKFLOW.value : "hermes.workflow",
        #                  workflowsTypes.OF_DISPERSION.value  : "hera.simulations.openFoam.datalayer.hermesWorkflow.Workflow_Dispersion",
        #                  workflowsTypes.OF_FLOWFIELD.value  : "hera.simulations.openFoam.datalayer.hermesWorkflow.Workflow_Flow"
        # }


    def getHemresWorkflowFromDocument(self,documentList,returnFirst=True):
        """
            Return a hermes-workflow (or a list of hermes workflows) to the user.

        Parameters
        ----------
        documentList : list, document
            A hera.datalayer document or a list of documents.

        returnFirst : bool
            If true, return obj of only the first iterm in the list (if it is a list).

        Returns
        -------
            hermes workflow object (or one of its derivatives).
        """

        docList = numpy.atleast_1d(documentList)

        if returnFirst:
            doc = docList[0]
            ret =  self.getHermesWorkflowFromJSON(doc.desc['workflow'],name=doc.desc['workflowName'])
        else:
            ret = [self.getHermesWorkflowFromJSON(doc.desc['workflow'],name=doc.desc['workflowName']) for doc in docList]

        return ret


    def getHermesWorkflowFromJSON(self,workflow : Union[dict,str],name=None):
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
        ky = workFlowJSON['workflow'].get('workflowType',None)

        if ky is None:
            hermesWFObj = pydoc.locate("hermes.workflow")
        else:
            hermesWFObj = pydoc.locate(f"hera.simulations.openFoam.OFWorkflow.workflow_{ky}")

        if hermesWFObj is None:
            err = f"The workflow type {ky} not found"
            self.logger.error(err)
            raise ValueError(err)

        return hermesWFObj(workFlowJSON,name=name)


    def getHermesWorkflowFromDB(self,nameOrWorkflowFileOrJSONOrResource : Union[dict, str,list,workflow],returnFirst=True,**query):
        """
                Retrieve workflows from the DB as hermes.workflow objects (or its derivatives).

                If the workflow is string, use it as a name. If the workflow is dict,
                use it as a filter on the paramters

                If returnFirst is False, return a list with all the results of the query.
                else, returns a single hermesworkflow.

        Parameters
        ----------
        workflow: str, dict
            The filtering criteria. Either name, or the parameters of the flow.

        returnFirst : bool
            If true, return only the first object (if found several results in the DB)

        query: arguments
            Additional query criteria.

        Returns
        -------
            list (returnFirst is False)
            hermes workflow.
        """

        docList = self.getWorkflowDocumentFromDB(nameOrWorkflowFileOrJSONOrResource,**query)

        if len(docList) == 0:
            self.logger.error(f"... not found. ")
            ret = None
        else:
            ret = self.getHemresWorkflowFromDocument(documentList=docList,returnFirst=returnFirst)
        return ret



    def getItemsFromDB(self,nameOrWorkflowFileOrJSONOrResource,doctype=None,**query):
        """
            Tries to find item as name, workflow directory , groupname or through the resource.
            Additional queries are also applicable.

        Parameters
        ----------
        nameOrWorkflowFileOrJSONOrResource : string or dict
                The name/dict that defines the item
        type  : string
            document type.

        query : dict
            Additional criteria.
        Returns
        -------
            doc or empty list if not found.
        """
        doctype = self.DOCTYPE_WORKFLOW if doctype is None else doctype
        # try to find it as a name
        mongo_crit = dictToMongoQuery(query)
        if isinstance(nameOrWorkflowFileOrJSONOrResource, str):
            self.logger.debug(f"Searching for {nameOrWorkflowFileOrJSONOrResource} as a name.")
            docList = self.getSimulationsDocuments(workflowName=nameOrWorkflowFileOrJSONOrResource, type=doctype,**mongo_crit)
            if len(docList) == 0:
                self.logger.debug(f"Searching for {nameOrWorkflowFileOrJSONOrResource} as a resource.")
                docList = self.getSimulationsDocuments(resource=nameOrWorkflowFileOrJSONOrResource, type=doctype,**mongo_crit)
                if len(docList) == 0:
                    self.logger.debug(f"Searching for {nameOrWorkflowFileOrJSONOrResource} as a workflow group.")
                    docList = self.getSimulationsDocuments(groupName=nameOrWorkflowFileOrJSONOrResource,type=doctype,**mongo_crit)
                    if len(docList) == 0:
                        self.logger.debug(f"... not found. Try to query as a json. ")
                        try:
                            jsn = loadJSON(nameOrWorkflowFileOrJSONOrResource)
                            wf = self.getHermesWorkflowFromJSON(jsn)
                            currentQuery = dictToMongoQuery(wf.parametersJSON, prefix="parameters")
                            currentQuery.update(mongo_crit)
                            docList = self.getSimulationsDocuments(type=self.DOCTYPE_WORKFLOW, **currentQuery)
                        except ValueError:
                            # self.logger.debug(f"Searching for {nameOrWorkflowFileOrJSONOrResource} as a file.")
                            # if os.path.isfile(nameOrWorkflowFileOrJSONOrResource):
                            #     from ..datalayer.document import nonDBMetadataFrame
                            #     workflowName = os.path.basename(nameOrWorkflowFileOrJSONOrResource).split(".")[0]
                            #     grpTuple = workflowName.split("_")
                            #     groupName = grpTuple[0]
                            #     groupID = grpTuple[1] if len(grpTuple) > 1 else 0
                            #     hermesWF = self.getHermesWorkflowFromJSON(nameOrWorkflowFileOrJSONOrResource)
                            #     res = nonDBMetadataFrame(data=None,
                            #                              projecName=self.projectName,
                            #                              resource=os.path.join(self.FilesDirectory,
                            #                                                    nameOrWorkflowFileOrJSONOrResource),
                            #                              dataFormat=datatypes.STRING,
                            #                              type=self.DOCTYPE_WORKFLOW,
                            #                              groupName=groupName,
                            #                              groupID=groupID,
                            #                              workflowName=workflowName,
                            #                              workflowType=hermesWF.workflowType,
                            #                              workflow=hermesWF.json,
                            #                              parameters=hermesWF.parametersJSON
                            #                              )
                            #     docList = [res]
                            # else:
                            self.logger.debug(f"not found")
                            docList = []
                    else:
                        self.logger.debug(f"... Found it as workflow group ")
                else:
                    self.logger.debug(f"... Found it as resource ")
            else:
                self.logger.debug(f"... Found it as name")

        elif isinstance(nameOrWorkflowFileOrJSONOrResource, dict) or isinstance(nameOrWorkflowFileOrJSONOrResource, workflow):

            qryDict = nameOrWorkflowFileOrJSONOrResource.parametersJSON if isinstance(nameOrWorkflowFileOrJSONOrResource, workflow) else nameOrWorkflowFileOrJSONOrResource

            self.logger.debug(f"Searching for {qryDict} using parameters")
            currentQuery = dictToMongoQuery(qryDict, prefix="parameters")
            currentQuery.update(mongo_crit)
            docList = self.getSimulationsDocuments(**currentQuery, type=self.DOCTYPE_WORKFLOW)
        else:
            docList = []

        return docList


    def getWorkflowDocumentFromDB(self, nameOrWorkflowFileOrJSONOrResource : Union[dict, str, list,workflow],**query):
        """
            Returns the simulation document from the DB.
            The nameOrWorkflowFileOrJSONOrResource can be either group name

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

        query : dict
            Additional query cireteria to the DB.

        Returns
        -------
            A document, or None if not found. .
        """

        if isinstance(nameOrWorkflowFileOrJSONOrResource,list):
            docList = []
            for simulationItem in nameOrWorkflowFileOrJSONOrResource:
                docList += self.getItemsFromDB(simulationItem)
        else:
            docList = self.getItemsFromDB(nameOrWorkflowFileOrJSONOrResource)

        return docList

    def getSimulationsInGroup(self, groupName: str, **kwargs):
        """
            Return a list of all the simulations.old with the name as a prefic, and of the requested simuationType.
            Returns the list of the documents.

            If the simuationType is None use the default simuationType (WORKFLOW).

        Parameters
        ----------

        groupName : str
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
        return self.getSimulationsDocuments(groupName=groupName, type=self.DOCTYPE_WORKFLOW, **kwargs)

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

        **kwargs : dict
                additional fileters.
                Note that these should be the full path of the parameter in the JSON.

        Returns
        -------
            int,str
            The new ID,
            The name.
        """
        simList = self.getSimulationsInGroup(groupName=simulationGroup, **kwargs)
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
                   execute: bool = False,
                   parameters: dict = dict()):
        """
            1. Adds the workflow to the database in the requested group
            2. Builds the template (.json) and python executer
            3. Runs the workflow.

            The stages are executed according to the buildMode.

            Notes:

            * If the workflow is already in the db in a different name adds to the db only if **force** is True.

            * If the workflowName already exist in the group then overwrite its definitions
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

        overwrite : bool
            If true, update the record if it exists.
            If false, throw an exception if the record exists.

        force : bool
            Allow duplicate workflows in the project.

        assignName : bool
            If true, finds the next available id and saves it in the DB.
            If groupName is None, parse the filename to get the group.

            Note that if true and group name is None, the names of the simulation will use only the string before the '_'.

            Otherwise, use the filename as the name of the simulation.

        execute : bool
            If true, execute the workflow

        buildModes enum.
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
        logger = get_classMethod_logger(self,"addToGroup")
        if workflow is None:
            raise NotImplementedError("addToGroup() requires the 'hermes' library, which is nor installed")
        logger.info("-- Start --")

        # 1. Getting the names of the simulation and the groups.

        #    a. Make sure that there are no extensions.
        cleanName = workflowJSON.split(".")[0]

        #   b. loading the workflow.
        self.logger.debug(f"Loading the workflow JSON {workflowJSON}")
        hermesWF = workflow(loadJSON(workflowJSON), self.FilesDirectory)
        hermesWF.updateNodes(parameters=parameters)
        theType = hermesWF.workflowType

        logger.debug(f"The suggested simulation name is {cleanName} as a workflow type {theType} (in file {workflowJSON})")

        #   c. Determining the simulation name, group name and group id
        groupName = groupName if groupName is not None else cleanName.split("_")[0]
        if assignName:
            self.logger.debug("Generating ID from the DB")
            groupID, workflowName = self.findAvailableName(simulationGroup=groupName, workflowType=theType)
            self.logger.debug(f" Got id : {groupID} and suggested name {workflowName}")
        else:
            workflowName = cleanName
            try:
                groupID = cleanName.split("_")[1]
                if not groupID.isdigit():
                    groupID = None
            except IndexError:
                # The name has no _ in it...
                groupID = None
            self.logger.debug(f"Use input as simulation : {workflowName} with the group {groupID}")

        self.logger.info(f"Simulation name is {workflowName} with type {theType} in simulation group {groupName} with id {groupID}.")

        # 2. Check if exists in the DB.

        #    a. Check if the workflow exists in the DB under different name (assume similar type)
        logger.debug(f"Checking if the workflow already exists in the db unde the same group and type.")
        currentQuery = dictToMongoQuery(hermesWF.parametersJSON, prefix="parameters")

        docList = self.getSimulationsInGroup(groupName=groupName, workflowType=theType, **currentQuery)

        if len(docList) > 0 and (not force) and (docList[0]['desc']['workflowName'] != workflowName):
            doc = docList[0]
            wrn = f"The requested workflow {workflowName} has similar parameters to the workflow **{doc['desc']['workflowName']}** in simulation group {groupName}."
            self.logger.warning(wrn)
            raise FileExistsError(wrn)
        else:

            #   b. Check if the name of the simulation already exists in the group
            self.logger.debug(f"Check if the name of the simulation {workflowName} already exists in the group")
            docList = self.getSimulationsInGroup(groupName=groupName, workflowType=theType, workflowName=workflowName)

            if len(docList) == 0:
                self.logger.info("Simulation is not in the DB, adding... ")
                doc = self.addSimulationsDocument(resource=os.path.join(self.FilesDirectory, workflowName),
                                                  dataFormat=datatypes.STRING,
                                                  type=self.DOCTYPE_WORKFLOW,
                                                  desc=dict(
                                                      groupName=groupName,
                                                      groupID=groupID,
                                                      workflowName=workflowName,
                                                      workflowType=theType,
                                                      workflow=hermesWF.json,
                                                      parameters=hermesWF.parametersJSON)
                                                  )

            elif overwrite:
                self.logger.info("Simulation in the DB, overwrite=True.  Updating... ")
                doc = docList[0]
                doc['desc']['workflow'] = hermesWF.json
                doc['desc']['parameters'] = hermesWF.parametersJSON
                doc.save()
            else:
                info = f"The simulation {workflowName} with type {theType} is already in the database in group {groupName}. use the overwrite=True to update the record."
                self.logger.info(info)

        # 3.  Building and running the workflow.
        if execute:
            logger.info(f"Building and executing the workflow {workflowName}")
            build = hermesWF.build(buildername=workflow.BUILDER_LUIGI)

            logger.info(f"Writing the workflow and the executer python {workflowName}")
            wfFileName = os.path.join(self.FilesDirectory, f"{workflowName}.json")
            # attemp to write only if the file is different than the input (that is it exists) or if overwrite (which mean it needs to be updated).
            if wfFileName != os.path.join(self.FilesDirectory,workflowJSON) or overwrite:
                hermesWF.write(wfFileName,overwrite=overwrite)
            with open(os.path.join(self.FilesDirectory, f"{workflowName}.py"), "w") as file:
                file.write(build)

            # delete the run files if exist.
            logger.debug(f"Removing the targetfiles and execute")
            executionfileDir = os.path.join(self.FilesDirectory, f"{workflowName}_targetFiles")
            shutil.rmtree(executionfileDir, ignore_errors=True)

            pythonPath = os.path.join(self.FilesDirectory, f"{workflowName}")
            executionStr = f"python3 -m luigi --module {os.path.basename(pythonPath)} finalnode_xx_0 --local-scheduler"
            self.logger.debug(executionStr)
            os.system(executionStr)



    def compareWorkflowsObj(self,
                         workflowList,
                         longFormat : bool = False):
        """
            Compares the parameters of the workflows to each other.


        Parameters
        ----------
        workflowList : list of hermes workflow objects
                The list of workflows to compare.
        Returns
        -------

        """
        return compareJSONS(**dict([(wf.name,wf.parametersJSON) for wf in workflowList]),
                            longFormat=longFormat)

    def workflow_compare(self,
                         workFlows: Union[list, str],
                         longFormat : bool = False,
                         transpose  : bool = False) -> Union[dict, pandas.DataFrame]:
        """
            Compares two or more hermes workflows.

        Parameters
        ----------
        workFlows : str,list
                A single input uses it as a group name,
                a list is the list of case names to compare.

        diffParams: bool
                If true display only differences.

        JSON: bool
                If true, return the results as a JSON and not pandas.DataFrame.

        Returns
        -------
            pandas.DataFrame, json (depends on the input flags).
            Return the differences between the parametrs of the requested cases.
        """
        if workFlows is None:
            raise NotImplementedError("compare() requires the 'hermes' library, which is nor installed")
        self.logger.info("--- Start ---")

        workflowList = []
        for workflowName in list(workFlows):
            if os.path.exists(workflowName):
                workflowList.append(self.getHermesWorkflowFromJSON(workflowName,name=workflowName))
            else:
                simulationList = self.getWorkflowDocumentFromDB(workflowName)
                groupworkflowList = [workflow(simulationDoc['desc']['workflow'],WD_path=self.FilesDirectory,name=simulationDoc.desc[self.DESC_WORKFLOWNAME]) for simulationDoc in simulationList]
                workflowList+=groupworkflowList

        res = self.compareWorkflowsObj(workflowList,longFormat=longFormat)

        return res.T if transpose else res

    def workflow_compareInGroup(self,workflowGroup,longFormat=False,transpose=False) :
        """
            Compares all the workflows in the group name.

        Parameters
        ----------
        workflowGroup : str
            The group name.

        Returns
        -------
            Pandas with the difference in the parameter names.
        """
        simulationList = self.getSimulationsInGroup(groupName=workflowGroup)
        workflowList = [workflow(simulationDoc['desc']['workflow'],WD_path=self.FilesDirectory,name=simulationDoc.desc[self.DESC_WORKFLOWNAME]) for simulationDoc in simulationList]
        res = self.compareWorkflowsObj(workflowList,longFormat=longFormat)
        return res.T if transpose else res

    def workflow_list(self,
                      workflowGroup:str,
                      listNodes : bool = False,
                      listParameters : bool = False) -> Union[pandas.DataFrame,dict]:
        """
            Lists all the simulations in the simulation group (of this project).

            Allows additional filters using the simulationType.

            If parameters is not None, return the list of parameters.
            return the parameters of all the nodes if the paraleters is an empty List, or the requested parameters.
            The default behaviour is to return only the parameters that are different from each other, unless allParams
            is True.

            The output is either pandas.DataFrame (if jsonFormat is False) or a JSON (if JSON is True).

        Parameters
        ----------
        workflowGroup : str
            The name of the group

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
        simulationList = self.getSimulationsInGroup(groupName=workflowGroup)
        ret = []
        for simdoc in simulationList:
            val = dict(workflowName=simdoc['desc']['workflowName'])

            if listNodes:
                val['nodes'] = simdoc['desc']['workflow']['workflow']['nodeList']

            if listParameters:
                val['parameters'] = simdoc['desc']['parameters']

            ret.append(val)

        return ret

    def workflow_listInGroups(self, workflowType=None, workflowName=True):
        """
            Lists all the simulation groups of the current project.

        Parameters
        ----------

        workflowType : str
                    The type of workflow to list.
                    If None, print all of them.

        workflowName : bool
                    if true, also lists all the simulations in that group.

        Returns
        -------

        """
        qry = dict(type=self.DOCTYPE_WORKFLOW)
        if workflowType is not None:
            qry['workflowType'] = workflowType

        docLists = self.getSimulationsDocuments(**qry)
        if len(docLists)==0:
            print(f"There are no workflow-groups in project {self.projectName}")
        else:
            data = pandas.DataFrame([dict(type=doc['desc']['workflowType'],workflowName=doc['desc']['workflowName'],groupName=doc['desc']['groupName']) for doc in docLists])

            for (groupType,groupName),grpdata in data.groupby(["type","groupName"]):
                ttl = f"{groupType}"
                print(ttl)
                print("-"*(len(ttl)))
                print(f"\t* {groupName}")
                if workflowName:
                    for simName in grpdata.workflowName.unique():
                        print(f"\t\t + {simName}")




















