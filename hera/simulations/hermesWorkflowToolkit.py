from enum import Enum, auto, unique
from typing import Union
from collections.abc import Iterable
import pandas
import shutil
import os
from ..toolkit import abstractToolkit
from hermes import workflow
from ..utils import loadJSON
from ..utils.query import dictToMongoQuery
from ..datalayer import datatypes
import numpy

try:
    from hermes import workflow
except ImportError:
    raise ImportError("hermes is not installed. please install it to use the hermes workflow toolkit.")


@unique
class actionModes(Enum):
    ADD = auto()
    ADDBUILD = auto()
    ADDBUILDRUN = auto()

DOCTYPE_WORKFLOW = "hermes_workflow"

class workflowToolkit(abstractToolkit):
    """
        Manages the hermes worflows:

            1. Checks if they are in the DB.
            2. create a new name to them
            3. allows simple deletion
            4. allows simple comparison.
            5. retrieve by initial name.
    """



    def __init__(self, projectName: str, filesDirectory: str = None):
        """
            Initializes the workflow toolkit.

        Parameters
        ----------
        projectName: str
            The project that the workflow will be used i.

        filesDirectory : str
            The directory to write all the workflows and the outputs. default is current directory.
        """
        super().__init__(projectName=projectName, filesDirectory=filesDirectory, toolkitName="hermesWorkflowToolkit")

    def getSimulationsInGroup(self, groupName: str, docType: str = None, **kwargs):
        """
            Return a list of all the simulations with the name as a prefic, and of the requested docTYPE.
            if the docType is None use the default doctype.

        Parameters
        ----------

        groupName : str
            The prefix name of all the runs.

        docType : str
            The type of the workflow.
            if None, use DOCTYPE_WORKFLOW

        Returns
        -------
            list
                A list of JSON files.
        """
        theType = self.DOCTYPE_WORKFLOW if docType is None else docType
        return self.getSimulationsDocuments(groupName=groupName, type=theType,**kwargs)

    def _findAvailableName(self, prefixName: str, docType: str = None, **kwargs):
        """
            Finds the next availabe name of that prefix. The available name is the maximal ID + 1.

            we assume that the name of each run is:
                <prefix name>_<id>.


        Parameters
        ----------
        prefixName : str
            The prefix name of all the runs.

        docType : str
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
        simList = self.getSimulationsInGroup(name=prefixName, docType=docType, **kwargs)
        group_ids = [x['desc']['groupID'] for x in simList if x['desc']['groupID'] is not None]
        if len(group_ids) == 0:
            newID = 1
        else:
            newID = numpy.max(group_ids)

        return newID, f"{prefixName}_{newID}"

    def addToGroup(self,
                   workflowJSON: str,
                   groupName: str = None,
                   docType: str = None,
                   overwrite: bool = False,
                   force: bool = False,
                   assignName: bool = False,
                   action: actionModes = actionModes.ADDBUILDRUN,
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

        docType : str
            The type of the workflow. [optional]
            If none, use the DOCTYPE_WORKFLOW type.

        overwrite : bool
            If true, update the record if it exists.
            If false, throw an exception if the record exists.

        assignName : bool
            If true, finds the next available id and saves it in the DB.
            If groupName is None, parse the filename to get the group.

            Note that if true and group name is None, the names of the simulation will use only the string before the '_'.

            Otherwise, use the filename as the name of the simulation.

        action : buildModes enum.
            ADDBUILDRUN : add the simulation to the db, builds the execution files (also saves the new tempalte) and executes it.
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
        self.logger.info("-- Start --")

        # 1. Getting the names of the simulation and the groups.

        #    a. Make sure that there are no extensions.
        cleanName = workflowJSON.split(".")[0]
        theType = self.DOCTYPE_WORKFLOW if docType is None else docType
        self.logger.debug(f"The suggested simulation name is {cleanName} in the document {theType}")

        #   b. loading the workflow.
        self.logger.debug(f"Loading the workflow JSON {workflowJSON}")
        hermesWF = workflow(loadJSON(workflowJSON), self.FilesDirectory)
        hermesWF.updateNodes(parameters=parameters)

        #   c. Determining the simulation anem, group name and group id
        groupName = groupName if groupName is not None else cleanName.split("_")[0]
        if assignName:
            self.logger.debug("Generating ID from the DB")
            groupID, simulationName = self._findAvailableName(prefixName=groupName, docType=theType)
            self.logger.debug(f" Got id : {groupID} and suggested name {newName}")
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
            f"Simulation name is {simulationName} with type {theType} in group {groupName} with id {groupID}.")

        # 2. Check if exists in the DB.

        #    a. Check if the workflow exists in the DB under different name (assume similar type)
        self.logger.debug(f"Checking if the workflow already exists in the db unde the same group and type.")
        currentQuery = dictToMongoQuery(hermesWF.parametersJSON, prefix="parameters")

        docList = self.getSimulationsInGroup(groupName=groupName, docType=theType, **currentQuery)



        if len(docList) > 0 and (not force) and (docList[0]['desc']['name'] != simulationName):
            doc = docList[0]
            wrn = f"The requested workflow {simulationName} has similar parameters to the workflow **{doc['desc']['name']}** in group {groupName}."
            self.logger.warning(wrn)
            raise FileExistsError(wrn)
        else:

            #   b. Check if the name of the simulation already exists in the group
            self.logger.debug(f"Check if the name of the simulation {simulationName} already exists in the group")
            docList = self.getSimulationsInGroup(groupName=groupName, docType=theType, name=simulationName)

            if len(docList) == 0:
                self.logger.info("Simulation is not in the DB, adding... ")
                doc = self.addSimulationsDocument(resource=os.path.join(self.FilesDirectory, simulationName),
                                                  dataFormat=datatypes.STRING,
                                                  type=theType,
                                                  desc=dict(
                                                      groupName=groupName,
                                                      groupID=groupID,
                                                      name=simulationName,
                                                      workflow=hermesWF.workflowJSON,
                                                      parameters=hermesWF.parametersJSON)
                                                  )


            elif overwrite:
                self.logger.info("Simulation in the DB, overwrite=True.  Updating... ")
                doc = docList[0]
                doc['desc']['workflow'] = workflowJSON
                doc['desc']['parameters'] = hermesWF.parametersJSON
                doc.save()
            else:
                err = f"The simulation {simulationName} with type {theType} is already in the database in group {groupName}. use the overwrite=True to update the record."
                self.logger.error(err)
                raise FileExistsError(err)



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
            Compares two or more simulations.


        Parameters
        ----------
        workFlows : str,list
                A single input uses it as a group name,
                a list is the list of cases to compare.

        docType : str
                The type of simulations. If not supplied use DOCTYPE_WORKFLOW

        nodes  : str,list
                The node name or a list of nodes to dislay

        allParameters: bool
                If true, then display all the parameters and not only those that are differ between two simulations.

        JSON: bool
                If true, return the results as a JSON and not pandas.DataFrame.

        Returns
        -------
            pandas.DataFrame, json (depends on the input flags).
            Return the differences between the parametrs of the requested cases.
        """
        self.logger.infor("--- Start ---")

        # 1. Get all the simulations
        if isinstance(workFlows, Iterable):
            self.logger.debug("Workflow is iterable. Trying to retrieve the parameters for each item individually. ")

            workflowList = []

            for simulation in workFlows:
                if os.path.exists(simulation):
                    workflowJSON = loadJSON(simulation)
                else:
                    theType = self.DOCTYPE_WORKFLOW if docType is None else docType
                    docList = self.getSimulationDocuments(name=simulation, type=theType)
                    if len(docList) == 0:
                        raise ValueError(f"Simulation {simulation} was not found on disk or in the project. ")
                    workflowJSON = docList[0]['desc']['workflow']

                workflowList.append(workflow(workflowJSON, WD_path=self.FilesDirectory))





        else:
            self.logger.debug("Workflow is a groupName. Get the simulations from the group")
            simulationList = self.getSimulationsInGroup(groupName=workFlows, docType=docType)
            workflowList = [workflow(x['desc']['workflow'], WD_path=self.FilesDirectory) for x in simulationList]

        return workflowList
