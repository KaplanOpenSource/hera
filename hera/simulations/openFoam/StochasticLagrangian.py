import os
import glob
import shutil
import json
from ...datalayer import datatypes
from . import DECOMPOSED_CASE,RECONSTRUCTED_CASE,TYPE_VTK_FILTER
from .datalayer.OFObjects import OFMeshBoundary
from ..hermesWorkflowToolkit import workflowToolkit,workflowsTypes
from ...utils.query import dictToMongoQuery
from ...utils.jsonutils import loadJSON

class stochasticLagrangianDataLayer:
    """
        A base class that handles all the stochastic lagrangian things (datalauer
    """
    toolkit = None
    logger  = None # take the logger from the OF toolkit.

    analysis = None # link to the analysis of the StochasticLagrnagian
    presentation = None # link to the presentation of the stochasticLagrangian

    def __init_(self,toolkit):
        self.toolkit = toolkit
        self.logger = toolkit.logger

    def createDispersionFlowField(self, flowData, overwrite:bool=False, useDBSupport: bool = True):
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
                    "originalFlow" : {
                        "source" : <name>,
                        "time" : {
                            type : "steadyState|dynamic",
                            "timestep" : <time>
                        },
                        linkMeshSymbolically : True
                    },
                    dispersionDuration : <duration>
                    dispersionFields : {

                    }

            base flow parameters:
            =====================

                * dynamicaltype: steady-state or dynamic.
                                 if dynamic, use the timesteps of the simulation (from timeStep and on). The time in the simulation will be
                                 [time - firsttime]

                                 if static, use the  timeStep and copy it to dispersionDuration.


                * timeStep: float, the first time step in the simulation. In case of steady-state, the only time step.
                                If
                * dispersionDuration: float, required, the maximal timestep of the dispersion simulation (of the flow field).
                                         required only in steady-state simulations.
                * copyMesh: bool, if true, then copy the mesh. otherwise, make symbolic links.

                base flow type
                --------------
                        - directory: The name of the directory to use.
                            Parameters:
                                    + name
                        - dispersionFlowFieldName : query the hera db by the simulation name,
                            Parameters
                                    + name : the base name of the run.

                        - workflowGroup: query the simulation group by the filters.
                                    workflowGroup : the group name.
                                    + filters : { nodename : { jpath : value}}
                                            where nodename is the name of the node in the workflow.

                        - workflowFile: query the hera db using the workflow.
                            Parameters:
                                    + filters : { nodename : { jpath : value}}
                                            where nodename is the name of the node in the workflow.
                        - ID : the id of the document in the heradb.

                        Implemented:
                            - directory
                            - dispersionFlowFieldName (without filters).

                        To implement the rst, we might need to extend the scripts to use multiple originalFlow simulations.old.

            useDBSupport : bool
                    If False, does not access the DB to check if the flow exists.
                    and does not add it. Used in the hermes module, when each dispersion case has its own flow (as it is soft link, it does not take space).

        Returns
        -------
            A list of str

            The ids of the documents that will be used as the directories.

        """
        def toTimeFormat(basedir,timeDirectory,logger):
            orig_proc = os.path.join(basedir, str(timeDirectory))
            # The directories in round times are int and not float.
            # that is 10 and not 10.0. Therefore, we must check and correct if does not exist.
            if not os.path.exists(orig_proc):
                logger.debug("Source does not exist, probably due to float/int issues. Recreating the source with int format")
                orig_proc = os.path.join(basedir, str(int(timeDirectory)))
                if not os.path.exists(orig_proc):
                    logger.error(f"Source {orig_proc} does not exist... make sure the simulation is OK.")
                    raise ValueError(f"Source {orig_proc} does not exist... make sure the simulation is OK.")

        self.logger.info(f"Creating dispersion flow field : {flowData['originalFlow']['source']} ")
        originalFlow = flowData['originalFlow']

        # 1. Get the case directory of the original flow.
        self.logger.debug(f"tying to find the flow {originalFlow['source']} in the DB")
        docList = self.toolkit.getWorkflowDocumentFromDB(originalFlow['source'])

        if len(docList) == 0:
            self.logger.deug(f"that flow is not in the DB. trying to interpret as a directory ")
            if not os.path.isdir(originalFlow['source']):
                err = f"The original flow {originalFlow['source']} is not a directory and does not exist in the DB. Use [hera-workflows list workflows] to list all the workflows"
                self.logger.critical(err)
                raise FileNotFoundError(err)
            else:
                if os.path.isdir(os.path.join(originalFlow['source'],'system')) and os.path.isdir(os.path.join(originalFlow['source'],'constant')):
                    originalFlowCaseDir = originalFlow['source']
                    workflowGroup = os.path.basename(originalFlowCaseDir)
                else:
                    err = f"The directory {originalFlow['source']} is not a case directory (doesn't have system or constant subdirs)"
                    self.logger.critical(err)
                    raise ValueError(err)
        elif len(docList) > 1:
            err = f"The name {originalFlow['source']} has more than one simulations to it. Did you supply a workflow group?"
            self.logger.error(err)
            raise ValueError(err)
        else:
            originalFlowCaseDir = docList[0].resource
            workflowGroup       = docList[0].desc['workflowName']

        workflowGroup = f"{workflowGroup}_Dispersion"

        self.logger.info(f"Found the original flow directory: {originalFlowCaseDir}. Using {workflowGroup} as the workflow group for the disperison flow")

        self.logger.debug("Getting the time in the original flow. Determine whether the simulation is parallel or not.")
        if os.path.exists(os.path.join(originalFlowCaseDir, "processor0")):
            self.logger.debug("Found directory 'processor0' assuming parallel")
            ptPath = ["processor0", "*"]
            parallelOriginal = True
        else:
            self.logger.debug("Directory 'processor0' not found!.  assuming single processor")
            ptPath = ["*"]
            parallelOriginal = False

        TS = [float(os.path.basename(ts)) for ts in glob.glob(os.path.join(originalFlowCaseDir, *ptPath)) if
              os.path.basename(ts).replace(".", "").isdigit()]
        TS.sort()

        self.logger.info(f"Found timesteps : {TS} in original flow")

        dynamicType = originalFlow['time']['type']

        timeStep = originalFlow.get("timeStep",None)
        dispersionDuration = flowData["dispersionDuration"]

        copyMesh = not originalFlow.get("linkMeshSymbolically", True)

        if timeStep is None:
            # find maximal TS, assume it is parallel:
            uts = TS[-1]
        else:
            uts = TS[min(range(len(TS)), key=lambda i: abs(TS[i] - timeStep))]

        if dynamicType==self.toolkit.TIME_STEADYSTATE:
            self.logger.debug(f"Using Time step {uts} for Steady state")
            # steady state, only 2 time steps
            timeList = [uts, str(float(dispersionDuration) + float(uts))]
        else:
            self.logger.debug(f"Using Time step {uts} as first time step for dynamic simulation")
            timeList = [x for x in TS if x > uts]

        self.logger.info(f"The simulation type is: {dynamicType}: Using time steps: {timeList}")

        ## Now we should find out if there is a similar run.
        # 1. same original run.
        # 2. same Hmix/ustar/... other fields.
        # 3. the requested start and end are within the existing start and end.

        if useDBSupport:
            querydict = dict(
                groupName=workflowGroup,
                flowParameters=dict(
                    baseFlowDirectory = originalFlowCaseDir,
                    dynamicType = dynamicType,
                    flowFields=flowData['dispersionFields'],
                    timeStep = uts,
                    originalFlowParameters = originalFlow
                )
            )
            self.logger.debug(f"Trying to find the simulation in the database. The run is: \n {json.dumps(querydict,indent=4)}")
            docList = self.toolkit.getSimulationsDocuments(type=workflowsTypes.OF_FLOWDISPERSION.value, **dictToMongoQuery(querydict))
            self.logger.debug(f"Found {len(docList)} in the database")
        else:
            self.logger.debug(f"Running without DB support, so does not query the db ")
            docList = []

        if len(docList) == 0 or overwrite:
            self.logger.info(f"Flow field not found, creating new and adding to the DB. ")
            ofhome = self.toolkit.OFObjectHome()
            if len(docList) == 0:
                self.logger.debug("Find the first name that is free")
                groupID = 1
                dispersionFlowFieldName = workflowGroup
                while not os.path.exists(os.path.join(self.toolkit.FilesDirectory,f"{dispersionFlowFieldName}_{groupID}")):
                    groupID += 1
                groupID = str(groupID)

            else: # overwrite existing
                self.logger.info(f"Found simulation {docList[0]['name']} on the disk , overwriting")
                resource = docList[0].resource
                dispersionFlowFieldName = docList[0].desc['name']
                groupID= docList[0].desc['groupID']
                self.logger.debug(f"The flow exists with the name {dispersionFlowFieldName}, so deleting {resource}, and writing over it")
                shutil.rmtree(resource)

            dispersionFlowFieldName = f"{dispersionFlowFieldName}_{groupID}"
            dispersionFlowFieldDirectory = os.path.abspath(os.path.join(self.toolkit.FilesDirectory,dispersionFlowFieldName))
            self.logger.info(f"Creating Dispersion flow simulation {dispersionFlowFieldName} in {os.path.abspath(self.toolkit.FilesDirectory)}")
            os.makedirs(dispersionFlowFieldDirectory,exist_ok=True)

            self.logger.info(f"Creating the flow specific fields in the flow needed for the dispersion")
            dispersionFields = flowData['dispersionFields']

            meshBoundary = OFMeshBoundary(originalFlowCaseDir).getBoundary()
            dispersionFieldList = []
            for dispersionFieldName, dispersionFieldData in dispersionFields.item():
                self.logger.debug(f"Creating the flow specific field: {dispersionFieldName}. ")
                field = ofhome.getFieldFromJSON(fieldName=dispersionFieldName,configuration=dispersionFieldData,meshBoundary=meshBoundary)
                dispersionFieldList.append( field )

            self.logger.info("Copying the configuration directories from the original to the new configuration")
            # copy constant, 0 and system.
            for general in ["constant", "system", "0"]:
                self.logger.debug(f"\tCopying {general} in {originalFlowCaseDir} directory --> {dispersionFlowFieldDirectory}")
                orig_general = os.path.join(originalFlowCaseDir, general)
                dest_general = os.path.join(dispersionFlowFieldDirectory, general)
                if os.path.exists(dest_general):
                    self.logger.debug(f"path {dest_general} exists... removing before copy")
                    shutil.rmtree(dest_general)
                shutil.copytree(orig_general, dest_general)

            self.logger.info(f"Copy directories. The run is {'parallel' if parallelOriginal else 'single-core'}")
            if parallelOriginal:
                origDirsList =  glob.glob(os.path.join(originalFlowCaseDir, "processor*"))
            else:
                origDirsList = [originalFlowCaseDir]

            for orig_time in timeList:
                dest_time = orig_time - timeList[0]
                if (dynamicType == self.toolkit.TIME_STEADYSTATE) and (orig_time == timeList[-1]):
                    orig_time = timeList[0]

                # We should look into it more closly, why parallel case doesn't recognize the time steps of the
                # processors. For now, just create these directories in the main root as well.
                os.makedirs(os.path.join(dispersionFlowFieldDirectory, str(dest_time)), exist_ok=True)

                self.logger.info(f"Mapping Dispersion time {dest_time} to {orig_time}")

                for origDir in origDirsList:
                    # The directories in round times are int and not float.
                    # that is 10 and not 10.0. Therefore, we must check and correct if does not exist.
                    orig_proc_timestep = toTimeFormat(basedir=origDir, timeDirectory=str(orig_time), logger=self.logger)

                    parallelOrSinglePathTime = [os.path.basename(origDir),str(dest_time)] if parallelOriginal else [str(dest_time)]

                    dest_proc_timestep = os.path.join(dispersionFlowFieldDirectory,*parallelOrSinglePathTime)
                    self.logger.info(f"\tMapping {orig_proc_timestep} --> {dest_proc_timestep}")
                    if os.path.exists(dest_proc_timestep):
                        self.logger.debug(f"path {dest_proc_timestep} exists... removing before copy")
                        shutil.rmtree(dest_proc_timestep)
                    shutil.copytree(orig_proc_timestep,dest_proc_timestep)

                    if copyMesh:
                        orig_constant = os.path.join(origDir,"constant")
                        parallelOrSinglePathConstant = [os.path.basename(origDir), "constant"] if parallelOriginal else ["constant"]
                        dest_constant = os.path.join(dispersionFlowFieldDirectory,*parallelOrSinglePathConstant)
                        self.logger.info(f"Copying the mesh in {orig_constant} to  {dest_constant}")
                        shutil.copytree(orig_constant, dest_constant)
                    else:
                        orig_constant_polymesh = os.path.abspath(os.path.join(origDir, "constant", "polyMesh"))
                        parallelOrSinglePathConstant = [os.path.basename(origDir), "constant","polyMesh"] if parallelOriginal else ["constant","polyMesh"]
                        destination_constant_polymesh = os.path.join(dispersionFlowFieldDirectory, os.path.basename(origDir), *parallelOrSinglePathConstant)
                        self.logger.info(f"Linking mesh in {orig_constant_polymesh} to {destination_constant_polymesh}")
                        if not os.path.exists(destination_constant_polymesh):
                            self.logger.debug(f"Linking {orig_constant_polymesh} -> {destination_constant_polymesh}")
                            os.makedirs(os.path.dirname(destination_constant_polymesh), exist_ok=True)
                            os.system(f"ln -s {orig_constant_polymesh} {destination_constant_polymesh}")

                    for field in dispersionFieldList:
                        self.logger.info(f"Writing field {field.name} to {dispersionFlowFieldDirectory} in time step {str(dest_time)}")
                        field.write(caseDirectory=dispersionFlowFieldDirectory,
                                    location=str(dest_time),
                                    parallel=parallelOriginal)


            self.logger.info("Finished creating the flow field for the dispersion. ")
            if useDBSupport:
                self.logger.info("Adding to the database.")
                if len(docList) ==0:
                    self.logger.debug("Updating the metadata of the record with the new group ID and simulation name")
                    querydict.update(dict(
                        groupID=groupID,
                        name=dispersionFlowFieldName,
                    ))
                    self.logger.debug("Adding record to the database")
                    self.toolkit.addSimulationsDocument(resource=dispersionFlowFieldDirectory, type=workflowsTypes.OF_FLOWDISPERSION.value, dataFormat=datatypes.STRING, desc=querydict)

            ret = dispersionFlowFieldDirectory
        else:
            self.logger.info(f"Found the requested flow in the flowFields of the project. Returning {docList[0].resource}")
            ret = docList[0].resource

        self.logger.info("... Done")
        return ret

    def createDispersionCaseDirectory(self,dispersionCaseDirectory,dispersionWorkflow,updateDB=True,exportFromDB=False,allowDuplicate=False,rewrite=False):

        self.logger.execution(f"----- Start -----")

        if (updateDB and exportFromDB):
            err = "Cannot use both --updateDB and --exportFromDB"
            self.logger.error(err)
            raise ValueError(err)

        try:
            dispersionWorkFlow = loadJSON(dispersionWorkflow)
        except ValueError:
            dispersionWorkFlow = None

        if dispersionWorkFlow is None:
            err = f"The dispersionCaseDirectory must be a workflow file. Got {dispersionCaseDirectory}"
            self.logger.error(err)
            raise ValueError(err)

        hermes_dispersionWorkflow = self.toolkit.getHermesWorkflowFromJSON(dispersionWorkFlow)


        # 1. Check if the dispersionFlowField is a workflow name or a directory
        dispersionFlowFieldDocumentList = self.getWorkflowDocumentFromDB(dispersionFlowField)
        if len(dispersionFlowFieldDocumentList) == 0:
            self.logger.debug(
                f"Not found in the DB, trying to use the argument {dispersionFlowField} as a directory")
            if os.path.isdir(os.path.abspath(dispersionFlowField)):
                dispersionFlowFieldName = os.path.abspath(dispersionFlowField)
                dispersionFlowFieldDirectory = dispersionFlowFieldName
                self.logger.info(
                    f"Using the dispersion flow {dispersionFlowFieldName} as name and directory. Parameters will not be available since it is not in DB")
            else:
                err = f"Cannot find the dispersion flow field {dispersionFlowField}"
                self.logger.error(err)
                raise ValueError(err)
        else:
            dispersionFlowFieldName = dispersionFlowFieldDocumentList[0].desc['workflowName']
            dispersionFlowFieldDirectory = dispersionFlowFieldDocumentList[0].resource
            self.logger.debug(
                f"Found the dispersion flow {dispersionFlowFieldName} in the DB. Using direcotry {dispersionFlowFieldDirectory}")

        currentDispersionName = os.path.basename(dispersionWorkflow).split(".")[0]
        self.logger.info(
            f"The dispersion workflow {currentDispersionName} and the dispersionFlowField Name {dispersionFlowFieldName} (direccotry {dispersionFlowFieldDirectory})")

        updateDispersionFlowField = False
        updateWorkflow = False
        needRewrite = False
        needDBAdd = False

        # 2. Check if the dispersion workflow name is in the DB.
        DBDispersionWorkflowDocument = self.toolkit.getWorkflowDocumentFromDB(currentDispersionName)
        if len(DBDispersionWorkflowDocument) > 0:
            self.logger.execution(
                f" the workflow {currentDispersionName} in the DB. First see if the dispersion flow field is similar")

            # 2.1 Check if the dispersion flow fields are similar
            if (DBDispersionWorkflowDocument[0]['desc']['dispersionFlowFieldName'] != dispersionFlowFieldName):
                self.logger.debug(
                    f"The dispersion with the name {currentDispersionName} exists in DB and is associated to the disperesion flow {DBDispersionWorkflowDocument[0]['desc']['dispersionFlowFieldName']}.")
                if updateDB:
                    self.logger.debug(
                        "updateDB and rewrite flags exist, so mark for removing the dispersion case and updating the record")
                    updateDispersionFlowField = True
                    needRewrite = True
                else:
                    err = f"The dispersion with the name {currentDispersionName} exists in DB and is associated to the disperesion flow {DBDispersionWorkflowDocument[0]['desc']['dispersionFlowFieldName']}. In order to change it to dispersion {dispersionFlowFieldName}, either change the dispersion name or use the --updateDB"
                    self.logger.error(err)
                    raise ValueError(err)
            else:
                self.logger.debug("dispersion not found in DB, add it")
                needDBAdd = True

            self.logger.execution(f"Check if the disk workflow is identical to the DB workflow")
            diskDispersionWorkflow = self.toolkit.getHermesWorkflowFromJSON(dispersionWorkFlow)
            DBDispersionWorkflow = self.toolkit.getHermesWorkflowFromJSON(DBDispersionWorkflowDocument[0]['desc']['workflow'])

            # 2.1 Check if the workflows are similar
            res = self.toolkit.compareWorkflowsObj([DBDispersionWorkflow, diskDispersionWorkflow])
            if len(res.columns) == 1:
                self.logger.info(f"The workflow in the DB is identical to the disk.")
                DispersionWorkflow = diskDispersionWorkflow
            else:
                self.logger.info(
                    f"The workflow in the DB is different than the disk. Must specify whether to use the disk version  (--updateDB), or the DB version (--exportFromDB). ")
                if exportFromDB:
                    self.logger.execution(f"Exporting the workflow to file {dispersionWorkflow}")
                    with open(dispersionWorkflow, 'w') as JSONOut:
                        json.dump(DBDispersionWorkflowDocument[0]['desc']['workflow'], JSONOut, indent=4)

                    DispersionWorkflow = DBDispersionWorkflow
                elif updateDB:
                    self.logger.execution("Updating the DB with the workflow from disk")
                    DispersionWorkflow = diskDispersionWorkflow
                    updateWorkflow = True
                    needRewrite = True
                else:
                    err = f"The workflow in {dispersionWorkflow} file and in the DB are different. Use the --updateDB to update the DB or --exportFromDB to rewrite the file"
                    self.logger.error(err)
                    raise ValueError(err)

        # 3. Check if workflow already in DB under a different name
        self.logger.execution("Check if workflow already in DB under a different name")
        #    Stop unless allowDuplicate flag exist
        DBDispersionWorkflow = self.toolkit.getHermesWorkflowFromDB(dispersionWorkFlow,
                                                          dispersionFlowFieldName=dispersionFlowFieldName)

        if DBDispersionWorkflow is not None:
            self.logger.debug(
                f"Found the workflow in the DB under the name {DBDispersionWorkflow[0]['desc']['workflowName']}")
            if currentDispersionName != DBDispersionWorkflow[0]['desc'][
                'workflowName'] and not allowDuplicate:
                info = f"Current simulation name is {currentDispersionName}. The workflow already in DB under the name {DBDispersionWorkflow[0]['desc']['workflowName']}. use the --allowDuplicate to continue"
                self.logger.error(info)
                raise ValueError(info)

        dispersionName = dispersionWorkflow.split(".")[0]
        dispersionDirectoryName = os.path.abspath(dispersionName)
        self.logger.debug(
            f"Getting the dispersion directory name from {dispersionWorkflow}: using {dispersionDirectoryName}")

        # 4.  check if the dispersion directory exists.
        #    if it exist, remove if the --rewrite flag exists.
        self.logger.execution(
            "Check if dispersion already exists on the disk. If it is remove for fresh start only in rewrite flag exists")
        if os.path.isdir(dispersionDirectoryName):
            if needRewrite:
                self.logger.execution(
                    f"The directory {dispersionDirectoryName} exists and differ than DB. Rewriting is needed as disk version was selected with --updateDB flag")
            if rewrite:
                self.logger.info(f"rewrite flag exists: removing the directory {dispersionDirectoryName}")
                shutil.rmtree(dispersionDirectoryName)
            else:
                err = f"The dispersion directory {dispersionDirectoryName} already exists. use --rewrite to force removing and recreating"
                self.logger.error(err)
                raise ValueError(err)

        # 3. Create the dispersion case and link to the workflow.
        self.logger.info(
            f"Creating dispersion case {dispersionDirectoryName}  and linking to {dispersionFlowFieldDirectory}")
        self.toolkit.stochasticLagrangian.createDispersionCaseDirectory(dispersionDirectoryName,
                                                              dispersionFlowDirectory=dispersionFlowFieldDirectory)


        # 5. Add/update to the DB.
        self.logger.execution("add to DB if not exist. If exists, check what needs to be updated (dispersionFlowField or the entire code)")
        if needDBAdd:
            groupName = dispersionName.split("_")[0]
            groupID = dispersionName.split("_")[1]
            self.logger.debug(f"Adding new document with group {groupName} and {groupID}")

            parameters = DispersionWorkflow.parametersJSON

            if len(DBDispersionWorkflowDocument) > 0:
                parameters['flowParameters'] = DBDispersionWorkflowDocument[0].desc['flowParameters']

            doc = self.toolkit.addSimulationsDocument(resource=dispersionFlowFieldDirectory,
                                            dataFormat=datatypes.STRING,
                                            type=workflowsTypes.OF_DISPERSION.value,
                                            desc=dict(
                                                dispersionFlowFieldName=dispersionFlowFieldName,
                                                groupName=groupName,
                                                groupID=groupID,
                                                workflowName=dispersionDirectoryName,
                                                workflowType=DispersionWorkflow.workflowType,
                                                workflow=DispersionWorkflow.json,
                                                parameters=parameters)
                                            )

        elif (updateDispersionFlowField or updateWorkflow):
            doc = dispersionFlowFieldDocumentList[0]
            if updateDispersionFlowField:
                self.logger.debug("Updating dispersion flow field")
                doc.desc['dispersionFlowFieldName'] = dispersionFlowFieldName

                if len(DBDispersionWorkflowDocument) > 0:
                    doc.desc['parameters']['flowParameters'] = DBDispersionWorkflowDocument[0].desc['flowParameters']

            if updateWorkflow:
                self.logger.debug("Updating workflow")
                doc.desc['workflow'] = DispersionWorkflow.json
                doc.desc['parameters'] = DispersionWorkflow.parametersJSON

                if len(DBDispersionWorkflowDocument) > 0:
                    doc.desc['parameters']['flowParameters'] = DBDispersionWorkflowDocument[0].desc['flowParameters']

            self.logger.debug("Save to DB")
            doc.save()

    # def createDispersionCaseDirectory(self, dispersionDirectory, dispersionFlowDirectory):
    #     dispersionDirectory = os.path.abspath(dispersionDirectory)
    #     dispersionFlowDirectory = os.path.abspath(dispersionFlowDirectory)
    #     self.logger.info(f"Creating dispersion at {dispersionDirectory} with base flow {dispersionFlowDirectory}")
    #
    #     if not os.path.isfile(dispersionFlowDirectory):
    #         err = f"The {dispersionFlowDirectory} is not not a directory."
    #         self.logger.error(err)
    #         raise ValueError(err)
    #
    #     constantDir = os.path.join(dispersionDirectory, "constant")
    #     systemDir = os.path.join(dispersionDirectory, "system")
    #     os.makedirs(constantDir, exist_ok=True)
    #     os.makedirs(systemDir, exist_ok=True)
    #
    #     self.logger.debug(f"Copying constant and system from the base flow")
    #     shutil.copytree(os.path.join(dispersionFlowDirectory, "constant"), os.path.join(dispersionDirectory, "constant"))
    #     shutil.copytree(os.path.join(dispersionFlowDirectory, "system"), os.path.join(dispersionDirectory, "system"))
    #
    #     for proc in glob.glob(os.path.join(dispersionFlowDirectory, "processor*")):
    #         self.logger.execution(f"Working on processor {proc}")
    #         fullpath = os.path.abspath(os.path.join(proc, "constant", "polyMesh"))
    #         destination = os.path.join(dispersionDirectory, os.path.basename(proc), "constant", "polyMesh")
    #         os.makedirs(os.path.dirname(destination), exist_ok=True)
    #         self.logger.debug(f"\t Linking: ln -s {fullpath} {destination}")
    #         self.logger.debug(f"\t Linking root case : ln -s {os.path.abspath(proc)} {os.path.join(dispersionDirectory, os.path.basename(proc))}/rootCase")
    #         os.system(f"ln -s {fullpath} {destination}")
    #         os.system(f"ln -s {os.path.abspath(proc)} {os.path.join(dispersionDirectory, os.path.basename(proc))}/rootCase")
    #
    #         # create the 0 directory in all processors.
    #         os.makedirs(os.path.join(dispersionDirectory, os.path.basename(proc), '0'), exist_ok=True)
    #         self.logger.debug(f"Creating the 0 time in {os.path.join(dispersionDirectory, os.path.basename(proc), '0')} ")
    #
    #     # linking the decpomposePar dict from the root.
    #     # root_decomposePar = os.path.abspath(os.path.join(dispersionFlow,"system","decomposeParDict"))
    #     # decompose_dest    = os.path.abspath(os.path.join(dispersionDirectory,"system"))
    #     # os.system(f"ln -s {root_decomposePar} {decompose_dest}")
    #
    #     # linking the rootCase in the root directory of the dispersion dispersionDirectory.
    #     self.logger.debug(f"Linking the root case: ln -s {dispersionFlowDirectory} {os.path.join(dispersionDirectory, 'rootCase')}")
    #     os.system(f"ln -s {dispersionFlowDirectory} {os.path.join(dispersionDirectory, 'rootCase')}")
    #
    #     # create the 0 directory in the root.
    #     self.logger.debug(f"Making the 0 in {os.path.join(dispersionDirectory, '0')}")
    #     os.makedirs(os.path.join(dispersionDirectory, '0'), exist_ok=True)
    #
    #     self.logger.info("... Done")