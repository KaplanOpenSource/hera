import json
from enum import Enum, auto, unique
from typing import Union
import pandas
import shutil
import os

from hera.toolkit import abstractToolkit
from hera.utils import loadJSON, compareJSONS
from hera.utils.query import dictToMongoQuery
from hera.datalayer import datatypes
import numpy
import pydoc
import warnings
from ..utils.logging import with_logger, get_classMethod_logger
from collections.abc import Iterable

try:
    from hermes import workflow
    from hermes.utils.workflowAssembly import handler_build, handler_buildExecute, handler_expand, handler_execute
except ImportError:
    #    raise ImportError("hermes is not installed. please install it to use the hermes workflow toolkit.")
    warnings.warn("hermes is not installed. some features will not work.")
    workflow = None


@unique
class actionModes(Enum):
    ADD = auto()
    ADDBUILD = auto()
    ADDBUILDEXECUTE = auto()


class hermesWorkflowToolkit(abstractToolkit):
    """
        Manages the hermes worflows:

            1. Checks if they are in the DB.
            2. create a new name to them
            3. allows simple deletion
            4. allows simple comparison.
            5. retrieve by initial name.
    """
    DESC_GROUPNAME = "groupName"
    DESC_GROUPID = "groupID"
    DESC_WORKFLOWNAME = "workflowName"
    DESC_PARAMETERS = "parameters"

    DOCTYPE_WORKFLOW = "hermesWorkflow"

    DOCKIND_CACHE = "Cache"
    DOCKIND_SIMULATIONS = "Simulations"

    def __init__(self, projectName: str, filesDirectory: str = None, toolkitName: str = "hermesWorkflowToolkit"):
        """
            Initializes the workflow toolkit.

        Parameters
        ----------
        projectName: str
            The project that the workflow will be used i.

        filesDirectory : str
            The directory to write all the Workflow and the outputs. default is current directory.
        """
        super().__init__(projectName=projectName,
                         filesDirectory=filesDirectory,
                         toolkitName=toolkitName)

        # ## Create the simulationType->object map
        # self._simulationTypeMap = {
        #                 WorkflowTypes.WORKFLOW.value : "hermes.workflow",
        #                  WorkflowTypes.OF_DISPERSION.value  : "hera.simulations.openFoam.datalayer.hermesWorkflow.Workflow_Dispersion",
        #                  WorkflowTypes.OF_FLOWFIELD.value  : "hera.simulations.openFoam.datalayer.hermesWorkflow.Workflow_Flow"
        # }

    def listHermesSolverTemplates(self, solverName):
        """
            Returns a list of all the templates that were loaded for that specific solver.

        Parameters
        ----------
        solverName : str


        Returns
        -------

        """
        retList = []
        for doc in self.getDataSourceDocumentsList(solver=solverName):
            data = dict(doc.desc)  # ['desc'])
            retList.append(data)

        if len(retList) > 0:
            return pandas.DataFrame(retList).set_index("name")
        else:
            return pandas.DataFrame()

    def getHermesFlowTemplate(self, hermesFlowName):
        """
            Get a hermes flow template
        Parameters
        ----------
        solverName

        Returns
        -------

        """
        return self.getDataSourceData(hermesFlowName, desc__component="Flow")

    def listHermesNodesTemplates(self):
        """
            Returns a list of all the templates that were loaded for that specific solver.

        Parameters
        ----------
        solverName : str


        Returns
        -------

        """
        retList = []
        for doc in self.getDataSourceDocumentsList(desc__component="Node"):
            data = dict(doc.desc['desc'])
            data['nodeName'] = doc.desc['datasourceName']
            retList.append(data)

        if len(retList) > 0:
            return pandas.DataFrame(retList).set_index("nodeName")
        else:
            return pandas.DataFrame()

    def getHermesNodeTemplate(self, hermesNodeName):
        """
            Get a hermes flow template
        Parameters
        ----------
        solverName

        Returns
        -------

        """
        return self.getDataSourceData(hermesNodeName, desc__component="Node")

    def getHemresWorkflowFromDocument(self, documentList, returnFirst=True):
        """
            Return a hermes-workflow (or a list of hermes Workflow) to the user.

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

        docList = documentList if isinstance(documentList, Iterable) else [documentList]

        if returnFirst:
            doc = docList[0]
            ret = self.getHermesWorkflowFromJSON(doc.desc['workflow'], name=doc.desc['workflowName'])
        else:
            ret = [self.getHermesWorkflowFromJSON(doc.desc['workflow'], name=doc.desc['workflowName']) for doc in
                   docList]

        return ret

    def getHermesWorkflowFromJSON(self, workflow: Union[dict, str], name=None):
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
        logger = get_classMethod_logger(self, "getHermesWorkflowFromJSON")
        workFlowJSON = loadJSON(workflow)

        ky = workFlowJSON['workflow'].get('solver', None)

        if ky is None:
            hermesWFObj = pydoc.locate("hermes.workflow")
        else:
            hermesWFObj = pydoc.locate(f"hera.simulations.openFoam.OFWorkflow.workflow_{ky}")

        if hermesWFObj is None:
            err = f"The workflow type {ky} not found"
            logger.error(err)
            raise ValueError(err)

        return hermesWFObj(workFlowJSON, name=name)

    def getHermesWorkflowFromDB(self, nameOrWorkflowFileOrJSONOrResource: Union[dict, str, list, workflow],
                                returnFirst=True, **query):
        """
                Retrieve Workflow from the DB as hermes.workflow objects (or its derivatives).

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
        logger = get_classMethod_logger(self, "getHermesWorkflowFromDB")
        docList = self.getWorkflowListDocumentFromDB(nameOrWorkflowFileOrJSONOrResource, **query)

        if len(docList) == 0:
            logger.error(f"... not found. ")
            ret = None
        else:
            ret = self.getHemresWorkflowFromDocument(documentList=docList, returnFirst=returnFirst)
        return ret

    def getWorkflowDocumentByName(self,name,doctype=None, dockind="Simulations",**query):
        """
            Retrieve the simulation using only the name.

            For a more sophisticated retrieval (that tries to retrieve by group name and by the content of the workflow)  use getWorkflowDocumentFromDB
        Parameters
        ----------
        name : str
            The name of the workflow.

        doctype  : string
            document type.

        dockind  : string
            Whether the document is cachaed or Simulation.

        query : dict
            Additional criteria.

        Returns
        -------

        """
        logger = get_classMethod_logger(self, "getWorkflowDocumentFromDB")
        doctype = self.DOCTYPE_WORKFLOW if doctype is None else doctype
        # try to find it as a name
        mongo_crit = dictToMongoQuery(query)

        retrieve_func = getattr(self, f"get{dockind}Documents")

        logger.info(f"Searching for {name} as a name of kind {dockind}")
        docList = retrieve_func(workflowName=name, type=doctype, **mongo_crit)
        return None if len(docList) == 0 else docList[0]


    def getWorkflowDocumentFromDB(self, nameOrWorkflowFileOrJSONOrResource, doctype=None, dockind="Simulations",
                                  **query):
        """
            Tries to find item as name, workflow directory , groupname or through the resource.
            Additional queries are also applicable.

        Parameters
        ----------
        nameOrWorkflowFileOrJSONOrResource : string or dict
                The name/dict that defines the item
        doctype  : string
            document type.

        dockind  : string
            Whether the document is cachaed or Simulation.

        query : dict
            Additional criteria.
        Returns
        -------
            doc or empty list if not found.
        """
        logger = get_classMethod_logger(self, "getWorkflowDocumentFromDB")
        doctype = self.DOCTYPE_WORKFLOW if doctype is None else doctype
        # try to find it as a name
        mongo_crit = dictToMongoQuery(query)

        retrieve_func = getattr(self, f"get{dockind}Documents")

        if isinstance(nameOrWorkflowFileOrJSONOrResource, str):
            logger.info(f"Searching for {nameOrWorkflowFileOrJSONOrResource} as a name of kind {dockind}")
            docList = retrieve_func(workflowName=nameOrWorkflowFileOrJSONOrResource, type=doctype, **mongo_crit)
            if len(docList) == 0:
                logger.info(f"Searching for {nameOrWorkflowFileOrJSONOrResource} as a resource of kind {dockind}.")
                docList = retrieve_func(resource=nameOrWorkflowFileOrJSONOrResource, type=doctype, **mongo_crit)
                if len(docList) == 0:
                    logger.info(
                        f"Searching for {nameOrWorkflowFileOrJSONOrResource} as a workflow group of kind {dockind}.")
                    docList = retrieve_func(groupName=nameOrWorkflowFileOrJSONOrResource, type=doctype, **mongo_crit)
                    if len(docList) == 0:
                        logger.info(f"... not found. Try to query as a json. ")
                        try:
                            jsn = loadJSON(nameOrWorkflowFileOrJSONOrResource)
                            wf = self.getHermesWorkflowFromJSON(jsn)
                            currentQuery = dictToMongoQuery(wf.parametersJSON, prefix="parameters")
                            currentQuery.update(mongo_crit)
                            docList = retrieve_func(type=self.DOCTYPE_WORKFLOW, **currentQuery)
                        except ValueError as e:

                            # logger.debug(f"Searching for {nameOrWorkflowFileOrJSONOrResource} as a file.")
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
                            err = f"Error {e} when trying to load as JSON."
                            logger.error(err)
                            docList = []
                        except IsADirectoryError:
                            logger.debug(f"not found")
                            docList = []
                    else:
                        logger.info(f"... Found it as workflow group ")
                else:
                    logger.info(f"... Found it as resource ")
            else:
                logger.info(f"... Found it as name")

        elif isinstance(nameOrWorkflowFileOrJSONOrResource, dict) or isinstance(nameOrWorkflowFileOrJSONOrResource,
                                                                                workflow):
            qryDict = nameOrWorkflowFileOrJSONOrResource.parametersJSON if isinstance(
                nameOrWorkflowFileOrJSONOrResource, workflow) else nameOrWorkflowFileOrJSONOrResource
            logger.debug(f"Searching for {qryDict} using parameters")
            currentQuery = dictToMongoQuery(qryDict, prefix="parameters")
            currentQuery.update(mongo_crit)
            docList = retrieve_func(**currentQuery, type=self.DOCTYPE_WORKFLOW)
        else:
            docList = []

        return docList

    def getWorkflowListDocumentFromDB(self, nameOrWorkflowFileOrJSONOrResource: Union[dict, str, list, workflow],
                                      **query):
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

        if isinstance(nameOrWorkflowFileOrJSONOrResource, list):
            docList = []
            for simulationItem in nameOrWorkflowFileOrJSONOrResource:
                docList += self.getWorkflowDocumentFromDB(simulationItem)
        else:
            docList = self.getWorkflowDocumentFromDB(nameOrWorkflowFileOrJSONOrResource)

        return docList

    def getWorkflowListOfSolvers(self, solverName: str, **query):
        """
            Returns all the documents of the requested solver.

        Parameters
        ----------
        solverName : str
            The name of the solver
        query : param-list
            Additional query

        Returns
        -------
            list of documnts.
        """
        return self.getSimulationsDocuments(solver=solverName, type=self.DOCTYPE_WORKFLOW, **query)

    def getWorkflowDocumentsInGroup(self, groupName: str, **kwargs):
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
        newID = self.getCounterAndAdd(f"simulations_{simulationGroup}")
        return newID, self.getworkFlowName(simulationGroup, newID)

    @staticmethod
    def getworkFlowName(baseName, flowID):
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

    def addUpdateWorkflowFileInGroup(self,workflowFileName):
        """
            Updates the content of the filename in the database.

            Use the file name as the simulation name.

        Parameters
        ----------
        workflowFileName

        addToGroup: bool
            Add the workflow to the group if False.

        Returns
        -------

        """
        logger = get_classMethod_logger(self, "updateWorkflowFileInGroup")
        if workflow is None:
            raise NotImplementedError("addUpdateWorkflowFileInGroup() requires the 'hermes' library, which is nor installed")

        workflowName = workflowFileName.split(".")[0]
        doc = self.getWorkflowDocumentByName(workflowName)
        if doc is None:
            workflowJSON = loadJSON(workflowFileName)
            hermesWF = workflow(workflowJSON, self.FilesDirectory)

            doc = self.addSimulationsDocument(resource=os.path.join(self.FilesDirectory, workflowName),
                                              dataFormat=datatypes.STRING,
                                              type=self.DOCTYPE_WORKFLOW,
                                              desc=dict(
                                                  groupName=workflowName.split("_")[0],
                                                  groupID=workflowName.split("_")[-1],
                                                  workflowName=workflowName,
                                                  solver=hermesWF.solver,
                                                  workflow=hermesWF.json,
                                                  parameters=hermesWF.parametersJSON)
                                              )
        return doc


    def addWorkflowToGroup(self, workflowJSON: str, groupName: str, writeWorkflowToFile:bool=False):
        """
            Adds the workflow to the database, or updates an existing document.

            If the fullName is true:
                The groupOrFullName is the name of the simulation.
                Hence, try to get the name from the DB. If the document exists update it with the
                new data if the overwrite is True.
                If the record is not in the DB, aadd it with that name.

            If the fullName is False:
                Check if the workflow is in the DB. If it is not,
                generate a new name with the counter and add it.


        Parameters
        ----------
        workflowJSON : str
            The file name that contains the workflow.

        groupName : str
            The group to assign the workflow to.

        fullName: str
            `Treat the groupOrFullName as a full name if True.

        overwrite : bool
            If true, then update the json workflow in the document if it exists.

        writeWorkflowToFile : bool
            If true, then write the JSON file to the disk.

        Returns
        -------
            The document of the new workflow.
        """
        logger = get_classMethod_logger(self, "addWorkflowToGroup")
        if workflow is None:
            raise NotImplementedError("addWorkflowToGroup() requires the 'hermes' library, which is nor installed")

        logger.debug(f"The name is a groupName, check if the workflow is in the DB and if not, generate a name and add it")
        docList = self.getWorkflowDocumentFromDB(loadJSON(workflowJSON))
        if len(docList) > 0:
            logger.info(f"...Found. Returning the document.")
            doc = docList[0]
        else:
            logger.info("...Not Found, adding the input to the DB")
            groupID = self.getCounterAndAdd(groupName)
            workflowData = loadJSON(workflowJSON)
            workflowName = self.getworkFlowName(groupName, groupID)
            hermesWF = workflow(workflowData, self.FilesDirectory)
            doc = self.addSimulationsDocument(resource=os.path.join(self.FilesDirectory, workflowName),
                                              dataFormat=datatypes.STRING,
                                              type=self.DOCTYPE_WORKFLOW,
                                              desc=dict(
                                                  groupName=groupName,
                                                  groupID=groupID,
                                                  workflowName=workflowName,
                                                  solver=hermesWF.solver,
                                                  workflow=hermesWF.json,
                                                  parameters=hermesWF.parametersJSON)
                                              )
            if writeWorkflowToFile:
                with open(os.path.join(self.FilesDirectory, f"{workflowName}.json"), "w") as outFile:
                    json.dump(hermesWF.json, outFile, indent=4)


        return doc


    def executeWorkflowFromDB(self, nameOrWorkflowFileOrJSONOrResource):
        """
            Building and Executing the workflow.

            Note that it is only executing the workflow.
            For OpenFOAM simulations, use the runOFSimulation method that will build the case
            and then run it.

            Note the procedure removes the [name]_targetFiles

        Parameters
        ----------
                nameOrWorkflowFileOrJSONOrResource: str, dict

        Can be
             - Resource (thedirectory name)
             - Simulation name
             - Its workflow
             - workfolow dict.

        build : bool [default = True]
            If true, also builds the workflow.

        Returns
        -------
            None
        """
        logger = get_classMethod_logger(self, "executeWorkflow")
        docList = self.getWorkflowListDocumentFromDB(nameOrWorkflowFileOrJSONOrResource)

        for doc in docList:
            workflowJSON = doc.desc['workflow']
            workflowName = doc.desc['workflowName']
            logger.info(f"Processing {workflowName}")

            hermesWF = self.getHermesWorkflowFromJSON(workflowJSON)

            logger.info(f"Building and executing the workflow {workflowName}")
            build = hermesWF.build(buildername=workflow.BUILDER_LUIGI)

            logger.info(f"Writing the workflow and the executer python {workflowName}")
            wfFileName = os.path.join(self.FilesDirectory, f"{workflowName}.json")
            hermesWF.write(wfFileName)

            pythonFileName = os.path.join(self.FilesDirectory, f"{workflowName}.py")
            with open(pythonFileName, "w") as outFile:
                outFile.write(build)

            # delete the run files if exist.
            logger.debug(f"Removing the targetfiles and execute")
            executionfileDir = os.path.join(self.FilesDirectory, f"{workflowName}_targetFiles")
            shutil.rmtree(executionfileDir, ignore_errors=True)

            pythonPath = os.path.join(self.FilesDirectory, f"{workflowName}")
            executionStr = f"python3 -m luigi --module {os.path.basename(pythonPath)} finalnode_xx_0 --local-scheduler"
            logger.debug(executionStr)
            os.system(executionStr)
            os.remove(pythonFileName)


    def compareWorkflowObj(self,
                           workflowList,
                           longFormat: bool = False):
        """
            Compares the parameters of the Workflow to each other.


        Parameters
        ----------
        workflowList : list of hermes workflow objects
                The list of Workflow to compare.
        Returns
        -------

        """
        return compareJSONS(**dict([(wf.name, wf.parametersJSON) for wf in workflowList]),
                            longFormat=longFormat,changeDotToUnderscore=True)


    def compareWorkflows(self,
                         Workflow: Union[list, str],
                         longFormat: bool = False,
                         transpose: bool = False) -> Union[dict, pandas.DataFrame]:
        """
            Compares two or more hermes Workflow.

        Parameters
        ----------
        Workflow : str,list
                A single input uses it as a group name,
                a list is the list of Workflow names to compare.

        diffParams: bool
                If true display only differences.

        JSON: bool
                If true, return the results as a JSON and not pandas.DataFrame.

        Returns
        -------
            pandas.DataFrame, json (depends on the input flags).
            Return the differences between the parametrs of the requested Workflow.
        """
        logger = get_classMethod_logger(self, "compareWorkflow")
        if Workflow is None:
            raise NotImplementedError("compare() requires the 'hermes' library, which is nor installed")
        logger.info("--- Start ---")

        workflowList = []
        for workflowName in numpy.atleast_1d(Workflow):
            simulationList = self.getWorkflowListDocumentFromDB(workflowName)
            if len(simulationList) == 0:
                if os.path.isfile(workflowName):
                    workflowList.append(self.getHermesWorkflowFromJSON(workflowName, name=workflowName))
                else:
                    err = f"Cannog find simulations for {workflowName}"
                    logger.error(err)
            else:
                groupworkflowList = [workflow(simulationDoc['desc']['workflow'], WD_path=self.FilesDirectory,
                                              name=simulationDoc.desc[self.DESC_WORKFLOWNAME]) for simulationDoc in
                                     simulationList]
                workflowList += groupworkflowList

        res = self.compareWorkflowObj(workflowList, longFormat=longFormat)

        return res.T if transpose else res


    def compareWorkflowInGroup(self, workflowGroup, longFormat=False, transpose=False):
        """
            Compares all the Workflow in the group name

            Each parameter that has different value across the workgroup is in the row, each simulation
            is in the column.

        Parameters
        ----------
        workflowGroup : str
            The group name.
        longFormat : bool
            If True, return the results in long format rather than in a wide table.
        transpose : bool
            If True return the simulation names as rows

        Returns
        -------
            Pandas with the difference in the parameter names.
        """
        simulationList = self.getWorkflowDocumentsInGroup(groupName=workflowGroup)
        workflowList = [workflow(simulationDoc['desc']['workflow'], WD_path=self.FilesDirectory,
                                 name=simulationDoc.desc[self.DESC_WORKFLOWNAME]) for simulationDoc in simulationList]
        if len(workflowList) == 0:
            ret = None
        else:
            res = self.compareWorkflowObj(workflowList, longFormat=longFormat).T
            ret =  res.T if transpose else res

        return ret

    def deleteWorkflowInGroup(self,workflowGroup,deepDelete=False,resetCounter=True):
        """
            Deletes all the workflows in the group.
        Parameters
        ----------
        workflowGroup : str
            The name of the workgroup.

        deepDelete: bool [default = False]
            If true, delete the resources.

        resetCounter : bool [default = True]
            Reset the counter of the group to 1.

        Returns
        -------

        """
        simulationList = self.getWorkflowDocumentsInGroup(groupName=workflowGroup)
        for doc in simulationList:
            if os.path.exists(doc.resource) and deepDelete:
                shutil.rmtree(doc.resource)

            doc.delete()

        if resetCounter:
            self.setCounter(counterName=workflowGroup)




    def listWorkflows(self,
                      workflowGroup: str,
                      listNodes: bool = False,
                      listParameters: bool = False) -> Union[pandas.DataFrame, dict]:
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
        simulationList = self.getWorkflowDocumentsInGroup(groupName=workflowGroup)
        ret = []
        for simdoc in simulationList:
            val = dict(workflowName=simdoc['desc']['workflowName'])

            if listNodes:
                val['nodes'] = simdoc['desc']['workflow']['workflow']['nodeList']

            if listParameters:
                val['parameters'] = simdoc['desc']['parameters']

            ret.append(val)

        return ret


    def listGroups(self, solver=None, workflowName=True):
        """
            Lists all the simulation groups of the current project.

        Parameters
        ----------

        solver : str
                    The name of the solver of this workflow.
                    If None, print all of them.

        workflowName : bool
                    if true, also lists all the simulations in that group.

        Returns
        -------

        """
        logger = get_classMethod_logger(self, "listGroups")
        qry = dict(type=self.DOCTYPE_WORKFLOW)
        logger.info("Getting the groups in the project.")
        if solver is not None:
            logger.info(f"....Using solver {solver}")
            qry['solver'] = solver

        docLists = self.getSimulationsDocuments(**qry)
        if len(docLists) == 0:
            logger.info(f"There are no workflow-groups in project {self.projectName}")
        else:
            data = pandas.DataFrame([dict(solver=doc['desc']['solver'], workflowName=doc['desc']['workflowName'],
                                          groupName=doc['desc']['groupName']) for doc in docLists])

            for (solverType, groupName), grpdata in data.groupby(["solver", "groupName"]):
                ttl = f"{solverType}"
                print(ttl)
                print("-" * (len(ttl)))
                print(f"\t* {groupName}")
                if workflowName:
                    for simName in grpdata.workflowName.unique():
                        print(f"\t\t + {simName}")


    def workflowTable(self, workflowGroup, longFormat=False, transpose=False):
        """
            Compares all the Workflow in the group name

            Each parameter that has different value across the workgroup is in the row, each simulation
            is in the column.

            Identical to the method compareWorkflowInGroup.

        Parameters
        ----------
        workflowGroup : str
            The group name.
        longFormat : bool
            If True, return the results in long format rather than in a wide table.
        transpose : bool
            If True return the simulation names as rows

        Returns
        -------
            Pandas with the difference in the parameter names.
        """

        return self.compareWorkflowInGroup(workflowGroup, longFormat, transpose)
