import warnings

import numpy
import os
import glob
import shutil
from distutils.dir_util import copy_tree
from ....utils import loadJSON
from ....utils import loggedObject
from ...hermesWorkflowToolkit import workflowToolkit,simulationTypes
from ...openFoam import ofObjectHome
from hera.datalayer import datatypes

class OFworkflowToolkit(workflowToolkit):
    """
        Adds some functionality to handle workflow with specialized funtions for OF.

        For example, add the code to create and handle the baseFlows for the dispersion workflow.
    """

    def __init__(self,projectName,filesDirectory):
        super().__init__(projectName=projectName,
                         filesDirectory=filesDirectory,
                         toolkitName="OFworkflowToolkit")




    def prepareFlowFieldForDispersion(self, flowData, simulationGroup:str):
        """
            Prepares the case directory of the flow for the dispersion, and assigns a local name to it.
            Currently, assumes the case is parallel.

            Steps:
            ------

            1. First checks in the DB if the simulation was already defined.
            2. If it is uses is id.
               If not:
                2.1. add it to the DB and get the id.
                2.2 Perpare the case:
                    1. Copies the base directory to the simulationsDirectory.
                    2. Creates the symbolic link for the mesh.
                    3. Creates the Hmix and ustar.
                    4. build the change directory and run it.
                    5. run the create cell heights.
                2.3 return the id

        Parameters
        ----------
        flowData : dict
            The configuration of the base flow.

            Has the structure:
                    "baseFlow" : {
                        type: "directory|simulationName|workflowFile|ID"
                        parameters:  {
                                .. depend on the type
                        },
                        "useTime" : [float], optional, default: last timestep
                        fromTimestep : [float], optional,  default: 0
                        maximalDispersionTime : [float], required.
                        copyMesh : [bool], default false

                    },
                    "dispersionField" : {

                    }

            base flow parameters:
            =====================

                * useTime : float, the time step to copy from the simulation.
                          If None (or not specified) use the last time step.

                * fromTimestep: float, the first time step in the simulation.
                * maximalDispersionTime: float, required, the maximal timestep of the dispersion simulation (of the flow field).
                * copyMesh: bool, if true, then copy the mesh. otherwise, make symbolic links.

                base flow type
                --------------
                        - directory: The name of the directory to use.
                            Parameters:
                                    + name
                        - simulationName : query the hera db by the simulation name,
                            Parameters
                                    + name : the base name of the run.

                        - simulationGroup: query the simulation group by the filters.
                                    simulationGroup : the group name.
                                    + filters : { nodename : { jpath : value}}
                                            where nodename is the name of the node in the workflow.

                        - workflowFile: query the hera db using the workflow.
                            Parameters:
                                    + filters : { nodename : { jpath : value}}
                                            where nodename is the name of the node in the workflow.
                        - ID : the id of the document in the heradb.

                        Implemented:
                            - directory
                            - simulationName (without filters).

                        To implement the rst, we might need to extend the scripts to use multiple baseFlow simulations.

            simulationGroup : str
                    The name of the group.

        Returns
        -------
            A list of str

            The ids of the documents that will be used as the directories.

        """
        self.logger.info("-------- Start ---------")
        baseFlow = flowData['baseFlow']

        # 1. Get the case of the run
        self.logger.debug(f"Getting the base flow directory from handler  {baseFlow['type']}")
        handler = getattr(self,f"getBaseFlow_{baseFlow['type']}")
        orig = handler(baseFlow['parameters'])
        self.logger.info(f"Found the base flow directory: {orig}")

        TS = [float(os.path.basename(ts)) for ts in glob.glob(os.path.join(orig, "processor0", "*")) if
              os.path.basename(ts).replace(".", "").isdigit()]
        TS.sort()

        useTime = baseFlow.get("useTime",None)
        fromTimestep = baseFlow.get("fromTimestep","0")
        maximalDispersionTime = baseFlow["maximalDispersionTime"]
        copyMesh = baseFlow.get("copyMesh",False)

        if useTime is None:
            # find maximal TS, assume it is parallel:
            uts = str(TS[-1])
        else:
            # find the closes TS.
            request = float(useTime)
            uts = TS[min(range(len(TS)), key=lambda i: abs(float(TS[i]) - request))]

        fromTime = fromTimestep

        self.logger.debug(f"Using Time step {uts} for Steady state")

        ## Now we should find out if there is a similar run.
        # 1. same original run.
        # 2. same Hmix/ustar/... other fields.
        # 3. the requested start and end are within the existing start and end.

        querydict = dict(
            groupName=simulationGroup,
            flowParameters=dict(
                baseFlowDirectory = orig,
                flowFields=flowData['dispersionFields'],
                usingTimeStep = uts,
                baseFlow = baseFlow
            )
        )

        self.logger.debug(f"Test if the requested flow field already exists in the project")
        docList = self.getSimulationsDocuments(type=simulationTypes.OF_FLOWDISPERSION.value,**querydict)

        if len(docList) == 0:
            self.logger.info(f"Flow field not found, creating new and adding to the DB. ")
            ofhome = ofObjectHome()

            groupID,simulationName = self.findAvailableName(simulationGroup,simuationType=simulationTypes.OF_FLOWDISPERSION.value)
            self.logger.info(f"Creating new name for the simulation.... Got ID {groupID} with simulation name {simulationName}")
            dest = os.path.abspath(os.path.join(self.FilesDirectory,simulationName))

            try:
                self.logger.info(f"Creating the directory {dest}")
                os.makedirs(dest,exist_ok=True)
            except FileExistsError:
                raise FileExistsError("The case already exists, use --overwrite ")

            self.logger.info("Copying the configuration directories")
            # copy constant, 0 and system.
            for general in ["constant", "system", "0"]:
                self.logger.debug(f"\tCopying {general} in {orig} directory --> {dest}")
                orig_general = os.path.join(orig, general)
                dest_general = os.path.join(dest, general)
                if os.path.exists(dest_general):
                    copy_tree(orig_general, dest_general)
                else:
                    shutil.copytree(orig_general, dest_general)

            self.logger.info(f"Copy the parallel directories")
            for proc in glob.glob(os.path.join(orig,"processor*")):

                orig_proc = os.path.join(proc, str(uts))
                # The directories in round times are int and not float.
                # that is 10 and not 10.0. Therefore, we must check and correct if does not exist.
                if not os.path.exists(orig_proc):
                    self.logger.debug("Source does not exist, probably de to float/int issues. Recreating the source with int format")
                    orig_proc = os.path.join(proc, str(int(uts)))
                    if not os.path.exists(orig_proc):
                        self.logger.error(f"Source {orig_proc} does not exist... make sure the simulation is OK.")
                        raise ValueError(f"Source {orig_proc} does not exist... make sure the simulation is OK.")

                for dest_time in [fromTime, maximalDispersionTime]:
                    self.logger.debug(f"Handling {proc} - time {dest_time}")
                    dest_proc = os.path.join(dest,os.path.basename(proc),str(dest_time))
                    shutil.copytree(orig_proc,dest_proc)

                self.logger.info(f"Copying the mesh")
                if copyMesh:
                    orig_constant = os.path.join(proc,"constant")
                    dest_constant = os.path.join(dest,os.path.basename(proc),"constant")
                    shutil.copytree(orig_constant, dest_constant)
                else:
                    fullpath = os.path.abspath(os.path.join(proc, "constant", "polyMesh"))
                    destination = os.path.join(dest, os.path.basename(proc), "constant", "polyMesh")
                    os.makedirs(os.path.dirname(destination), exist_ok=True)
                    os.system(f"ln -s {fullpath} {destination}")

            self.logger.info(f"Creating the flow specific fields in the flow needed for the dispersion")
            dispersionFields = flowData['dispersionFields']
            for dispersionFieldName in dispersionFields.keys():
                fieldDimensions = dispersionFields[dispersionFieldName].get("dimensions",None)
                fieldComponents = dispersionFields[dispersionFieldName].get("components", None)
                self.logger.debug(f"Creating the flow specific field: {dispersionFieldName}. ")
                field = ofhome.getField(fieldName=dispersionFieldName,
                                              fieldGroup=ofObjectHome.GROUP_DISPERSION,
                                              dimensions=fieldDimensions,
                                              componentNames=fieldComponents)

                for proc in glob.glob(os.path.join(orig, "processor*")):
                    procName = os.path.split(proc)[-1]
                    self.logger.debug(f"Writing the flow specific field {dispersionFieldName} to processor {procName} . ")
                    for dest_time in [fromTime, maximalDispersionTime]:
                        field.emptyParallelField(caseDirectory=dest,
                                                 timeName=str(dest_time),
                                                 processor=procName,
                                                 boundaryField=dispersionFields[dispersionFieldName]["boundaryField"],
                                                 data=dispersionFields[dispersionFieldName].get("internalField"))

                # We should look into it more closly, why parallel case doesn't recognize the time steps of the
                # processors. For now, just create these directories in the main root as well.
                for dest_time in [fromTime, maximalDispersionTime]:
                    os.makedirs(os.path.join(dest,str(dest_time)),exist_ok=True)


            self.logger.info("Finished creating the flow field for the dispersion. Adding to the database. ")
            self.logger.debug("Updating the metadata of the record with the new group ID and simulation name")
            querydict.update(dict(
                groupID=groupID,
                name=simulationName,
            ))

            self.logger.debug("Adding record to the database")
            self.addSimulationsDocument(resource=dest,
                                       type=simulationTypes.OF_FLOWDISPERSION.value,
                                       dataFormat=datatypes.STRING,
                                       desc=querydict)

            ret = dest
        else:
            self.logger.info(f"Found the requested flow in the flowFields of the project. Returning {docList[0].resource}")
            ret = docList[0].resource

        return ret

    ########################################## baseTemplateHandler
    #
    #  Retrieves the base workflows.
    #
    def getBaseFlow_directory(self, parameters : dict):
        """
            Return the name of the requested directory.

        Parameters
        ----------
        parameters: dict
            The json of the directory:
                {
                    name : ...
                }

        projectName : str
            The name of the project (not used in this procedure).

        Returns
        -------

        """
        return parameters['name']

    def getBaseFlow_simulationName(self, parameters):
        """
            Check if the directory is already the db.
            If it is, then return the id.

        Parameters
        ----------
        parameters: dict
            The json of the directory:
                {
                    name : ...
                }
        projectName : str
            The name of the project.

        Returns
        -------

        """
        docList = self.getSimulationsDocuments(name=parameters['name'])
        if len(docList) >0:
            if len(docList) > 1:
                warnings.warn(f"Found more than 1 simulation with the name {parameters['name']}. Return the first one")
            return docList[0].resource
        else:
            raise ValueError(f"Cannot find flows with the name {parameters['name']}. Use hera-OF-flows list to see the names of existing simulations.")




