from typing import Union
from collections.abc import Iterable
import pandas
import os
from ..toolkit import abstractToolkit
from hermes import workflow
from ..utils import loadJSON
from ..utils.query import dictToMongoQuery
from ..datalayer import datatypes
import numpy


class workflowToolkit(abstractToolkit):
    """
        Manages the hermes worflows:

            1. Checks if they are in the DB.
            2. create a new name to them
            3. allows simple deletion
            4. allows simple comparison.
            5. retrieve by initial name.
    """

    DOCTYPE_WORKFLOW = "hermes_workflow"

    def __init__(self, projectName: str, FilesDirectory: str = None):
        """
            Initializes the workflow toolkit.

        Parameters
        ----------
        projectName: str
            The project that the workflow will be used i.

        FilesDirectory : str
            The directory to write all the workflows and the outputs. default is current directory.
        """
        super().__init__(projectName=projectName, FilesDirectory=FilesDirectory)

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
        return self.getSimulationDocuments(groupName=groupName, type=theType)

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
        if len(group_ids) ==0:
            newID = 1
        else:
            newID = numpy.max(group_ids)

        return newID, f"{prefixName}_{newID}"

    def addToGroup(self, workflowJSON: str, groupName: str, docType: str = None, parameters=dict()):
        """
            Adds the workflow to the requested group, after the chane in the parameters
            if it does not exist in the DB.

            Assigns the new simlation a name from te group and then saves the
            expaned workflow and python execute to the disk with that name.


        Parameters
        ----------
        workflowJSON  : str
                JSON file, or json configuration

        groupName : str
                The base name of the simulation. Used to find the available names.

        docType  : str
                The type of the document. if None use general docWORKFLOW.

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
            str
            The name ofthe new
        """
        self.logger.info("--- Start ---")
        # 1. expand workflow:
        hermesWF = workflow(loadJSON(workflowJSON), self.FilesDirectory)
        hermesWF.updateNodes(parameters=parameters)

        # 2. Check if that simulation is already in the database.
        currentQuery = dictToMongoQuery(newTemplate.workflowJSON, prefix="workflow")
        theType = self.DOCTYPE_WORKFLOW if docType is None else docType

        docList = self.getSimulationsInGroup(groupName=groupName, docType=theType, **currentQuery)
        self.logger.info("Check if the current workflow exists in the database.")
        if len(docList) == 0:
            self.logger.info("Does not exist in the db. find a new name.")
            newID, newName = self._findAvailableName(prefixName=groupName, docType=theType)

            self.logger.info(f"Adding the new simulation to the database using the name {newName}")

            doc = self.addSimulationDocument(resource=os.path.join(self.FilesDirectory, newName),
                                             format=datatypes.STRING,
                                             type=theType,
                                             desc=dict(groupName=groupName,
                                                       groupID=newID,
                                                       name=newName,
                                                       workflow=newTemplate.workflowJSON)
                                             )

            build = hermesWF.build(buildername=workflow.BUILDER_LUIGI)

            self.logger.info(f"Writing the workflow and the executer python {newName}")
            with open(os.path.join(self.FilesDirectory, f"{newName}.json"), "w") as file:
                file.write(newTemplate.workflowJSON)
            with open(os.path.join(self.FilesDirectory, f"{newName}.py"), "w") as file:
                file.write(build)

            ret = newName
        else:
            ret = docList[0]['desc']['name']

    def addToGroupFromFile(self, fileName: str,groupName : str =None, docType: str = None, overwrite=False,force=False):
        """
            Loads the workflow in the file to the database to the requested group.

            If the group is None, use the file name as the group. We assume that the
            file name has the structure : <groupname>_<id>. If the <id> is not an integer,
            the id in the database will be saved as None.

            If the simulationName already exist in the group then overwrite its definitions only if
            the overwrite is True.

            Check if the workflow is in the group. if it does, add it with the new requested name (thus, create duplicity)
            only if force is True.

        Parameters
        ----------
        fileName : str
            The file name that contains the workflow.

        docType : str
            The type of the workflow. [optional]
            If none, use the DOCTYPE_WORKFLOW type.

        overwrite : bool
            If true, update the record if it exists.
            If false, throw an exception if the record exists.

        groupName : str
            The group to assign the workflow to.
            If None, use the file name as the group and the ID.


        Returns
        -------
            None
        """
        self.logger.info("-- Start --")

        cleanName = fileName.split(".")[0]
        theType = self.DOCTYPE_WORKFLOW if docType is None else docType
        workflow = loadJSON(fileName)

        groupName = groupName if groupName is not None else cleanName.split("_")[0]
        groupID   = cleanName.split("_")[1]
        if not groupID.isdigit():
            groupID = None

        self.logger.info(f"Loading {fileName} with type {theType} to simulation group {groupName} with group ID {groupID}")

        # 1. Check if the workflow itself is already in the DB.
        currentQuery = dictToMongoQuery(workflow, prefix="workflow")
        docList = self.getSimulationsInGroup(groupName=groupName, docType=theType, **currentQuery)

        if len(docList) > 0 and not force:
            doc = docList[0]
            raise FileExistsError(f"The requested workflow already exists in the group {groupName} with the name {doc['desc']['name']}")

        # 2. Check if the name of the simulation already exists
        docList = self.getSimulationsInGroup(groupName=groupName, docType=theType, name=cleanName)

        if len(docList) ==0:
            # add to the db
            doc = self.addSimulationDocument(resource=os.path.join(self.FilesDirectory, cleanName),
                                             format=datatypes.STRING,
                                             type=theType,
                                             desc=dict(
                                                 groupName=groupName,
                                                 groupID=groupID,
                                                 name=cleanName,
                                                 workflow=workflow)
                                             )

        elif overwrite:
            # update the record.
            doc = docList[0]
            doc['desc']['workflow'] = workflow
            doc.save()
        else:
            # print error and show the differences between the two cases.
            raise FileExistsError(
                f"The simulation {fileName} is already in the database. use the overwrite=True to update the record.")


    def compare(self,workFlows : Union[list,str],docType: str=None,nodes : Union[list,str]=None,allParameters:bool=True,JSON:bool=False)-> Union[dict,pandas.DataFrame]:
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
        if isinstance(workFlows,Iterable):
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

                workflowList.append(workflow(workflowJSON,WD_path=self.FilesDirectory))





        else:
            self.logger.debug("Workflow is a groupName. Get the simulations from the group")
            simulationList = self.getSimulationsInGroup(groupName=workFlows,docType=docType)
            workflowList = [workflow(x['desc']['workflow'],WD_path=self.FilesDirectory) for x in simulationList]

        return workflowList