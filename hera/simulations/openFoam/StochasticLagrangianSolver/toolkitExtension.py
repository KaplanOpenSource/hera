import numpy
import pandas
import os
import glob
import shutil
import json
from ....datalayer import datatypes
from ..OFObjects import OFMeshBoundary
from ....utils.query import dictToMongoQuery
from ....utils.logging import get_classMethod_logger

class toolkitExtension_LagrangianSolver:
    """
        A base class that handles all the stochastic lagrangian things (datalauer
    """
    toolkit = None
    logger  = None # take the logger from the OF toolkit.

    analysis = None # link to the analysis of the StochasticLagrnagian
    presentation = None # link to the presentation of the stochasticLagrangian

    def __init__(self,toolkit):
        self.toolkit = toolkit
        self.logger = toolkit.logger

    def createDispersionFlowField(self,flowName, flowData, overwrite:bool=False, useDBSupport: bool = True):
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

        flowName : str
            Suffix of the string.

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

        logger = get_classMethod_logger(self,"createDispersionFlowField")
        logger.info(f"Creating dispersion flow field : {flowData['originalFlow']['source']} ")
        originalFlow = flowData['originalFlow']

        # 1. Get the case directory of the original flow.
        logger.debug(f"tying to find the flow {originalFlow['source']} in the DB")
        docList = self.toolkit.getWorkflowDocumentFromDB(originalFlow['source'])

        if len(docList) == 0:
            logger.deug(f"that flow is not in the DB. trying to interpret as a directory ")
            if not os.path.isdir(originalFlow['source']):
                err = f"The original flow {originalFlow['source']} is not a directory and does not exist in the DB. Use [hera-workflows list workflows] to list all the workflows"
                logger.critical(err)
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
            err = f"The name {originalFlow['source']} has more than one simulations to it. "
            logger.error(err)
            raise ValueError(err)
        else:
            originalFlowCaseDir = docList[0].resource
            workflowGroup       = docList[0].desc['workflowName']

        workflowGroup = f"{workflowGroup}_DFF"
        dispersionFlowFieldName   = f"{workflowGroup}_{flowName}"

        logger.info(f"Found the original flow directory: {originalFlowCaseDir}. Using {workflowGroup} as the workflow group for the disperison flow, and {dispersionFlowFieldName} as its new name")

        logger.debug("Getting the time in the original flow. Determine whether the simulation is parallel or not.")
        if os.path.exists(os.path.join(originalFlowCaseDir, "processor0")):
            logger.debug("Found directory 'processor0' assuming parallel")
            ptPath = ["processor0", "*"]
            parallelOriginal = True
        else:
            logger.debug("Directory 'processor0' not found!.  assuming single processor")
            ptPath = ["*"]
            parallelOriginal = False

        TS = [float(os.path.basename(ts)) for ts in glob.glob(os.path.join(originalFlowCaseDir, *ptPath)) if
              os.path.basename(ts).replace(".", "").isdigit()]
        TS.sort()

        logger.info(f"Found timesteps : {TS} in original flow")
        dynamicType = originalFlow['time']['type']

        timeStep = originalFlow.get("timeStep",None)
        dispersionDuration = flowData["dispersionDuration"]

        copyMesh = not originalFlow.get("linkMeshSymbolically", True)

        if timeStep is None:
            logger.debug("timeStep is None: find maximal TS and assume it is parallel")
            uts = TS[-1]
        else:
            uts = TS[min(range(len(TS)), key=lambda i: abs(TS[i] - timeStep))]

        if dynamicType==self.toolkit.TIME_STEADYSTATE:
            logger.debug(f"Using Time step {uts} for Steady state")
            # steady state, only 2 time steps
            timeList = [uts, str(float(dispersionDuration) + float(uts))]
        else:
            logger.debug(f"Using Time step {uts} as first time step for dynamic simulation")
            timeList = [x for x in TS if x > uts]

        logger.info(f"The simulation type is: {dynamicType}: Using time steps: {timeList}")

        ## Now we should find out if there is a similar run.
        # 1. same original run.
        # 2. same Hmix/ustar/... other fields.
        # 3. the requested start and end are within the existing start and end.

        if useDBSupport:
            logger.debug(f"Check if {dispersionFlowFieldName} is in the DB with the required parameters")
            querydict = dict(
                groupName=workflowGroup,
                flowParameters=dict(
                    baseFlowDirectory = originalFlowCaseDir,
                    dynamicType = dynamicType,
                    flowFields=flowData['dispersionFields'],
                    timeStep = uts,
                    originalFlow = flowData['originalFlow']
                )
            )
            logger.debug(f"Trying to find the dispersion workflow  in the database. The run is: \n {json.dumps(querydict,indent=4)}")
            docList = self.toolkit.getSimulationsDocuments(type=self.toolkit.OF_FLOWDISPERSION, **dictToMongoQuery(querydict),workflowName=dispersionFlowFieldName)
            logger.debug(f"Found {len(docList)} in the database")

            if len(docList) > 1:
                err = f"Found more than one {dispersionFlowFieldName} with the same parameters set. Please fix it manually (using Project and deleteSimulationsDocuments)"
                logger.error(err)
                raise ValueError(err)

            if len(docList) == 0:
                logger.debug(f"Did not find dispersion field {dispersionFlowFieldName} with the requested paraeters. Checking to see if that name exists with other parameters. ")
                docList = self.toolkit.getSimulationsDocuments(type=self.toolkit.OF_FLOWDISPERSION,workflowName=dispersionFlowFieldName)
                logger.debug(f"Found {len(docList)} in the database")

                if len(docList) == 0:
                    logger.debug("Not found. Creating a new workflow")
                    doc = None
                else:
                    logger.debug("Found the name but with different paramters. overwrite if overwrite=True")
                    doc = docList[0]
            else:
                logger.debug("Found the name with the same paramters. overwrite if overwrite=True")
                doc = docList[0]
        else:
            logger.debug(f"Running without DB support, so does not query the db ")
            doc = None

        ofhome = self.toolkit.OFObjectHome()
        if doc is not None and overwrite:
            logger.info(f"Starting to overwrite existing dispersion field. Remove existing workflow field.")

            logger.info(f"Found Dispersion flow field {doc.desc['name']} on the disk , overwriting the same directory")
            resource = docList[0].resource
            logger.debug(f"Deleting {resource}, and writing over it")
            shutil.rmtree(resource)

        dispersionFlowFieldDirectory = os.path.abspath(os.path.join(self.toolkit.FilesDirectory,dispersionFlowFieldName))
        logger.info(f"Creating Dispersion flow simulation {dispersionFlowFieldName} in {dispersionFlowFieldDirectory}")
        os.makedirs(dispersionFlowFieldDirectory,exist_ok=True)

        logger.info(f"Creating the flow specific fields in the flow needed for the dispersion")
        dispersionFields = flowData['dispersionFields']

        meshBoundary = OFMeshBoundary(originalFlowCaseDir).getBoundary()
        dispersionFieldList = []
        for dispersionFieldName, dispersionFieldData in dispersionFields.item():
            logger.debug(f"Creating the flow specific field: {dispersionFieldName}. ")
            field = ofhome.getFieldFromJSON(fieldName=dispersionFieldName,configuration=dispersionFieldData,meshBoundary=meshBoundary)
            dispersionFieldList.append( field )

        logger.info("Copying the configuration directories from the original to the new configuration")
        # copy constant, 0 and system.
        for general in ["constant", "system", "0"]:
            logger.debug(f"\tCopying {general} in {originalFlowCaseDir} directory --> {dispersionFlowFieldDirectory}")
            orig_general = os.path.join(originalFlowCaseDir, general)
            dest_general = os.path.join(dispersionFlowFieldDirectory, general)
            if os.path.exists(dest_general):
                logger.debug(f"path {dest_general} exists... removing before copy")
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


        logger.info("Finished creating the flow field for the dispersion. ")
        if useDBSupport:
            logger.info("Adding to the database.")
            querydict['workflowName'] = dispersionFlowFieldName
            if doc is None:
                logger.debug("Updating the metadata of the record with the new group ID and simulation name")

                self.logger.debug("Adding record to the database")
                self.toolkit.addSimulationsDocument(resource=dispersionFlowFieldDirectory, type=self.toolkit.OF_FLOWDISPERSION, dataFormat=datatypes.STRING, desc=querydict)

            ret = dispersionFlowFieldDirectory
        else:
            logger.info(f"Found the requested flow in the flowFields of the project. Updating the description. Returning {docList[0].resource}")
            doc.desc =querydict
            doc.save()
            ret = doc.resource

        logger.info("... Done")
        return ret

    def createDispersionCaseDirectory(self, hermes_dispersionWorkflow, updateDB=True,exportFromDB=False, allowDuplicate=False, rewrite=False):
        """
            Creates a dispersion case directory and linking to the dispersion workflow.

            This function also updates the database and ensures that the hermes_dispersionWorkflow and the DB record are consistent.
            It also searches for a different workflow with the same name.

            If the DB and the input are not consistent then update the DB if the
            updateDB flag is true. Otherwise, update the input from the DB.

            If the directory already exists on the disk, there are inconsistencies and the updateDB flag is True (i.e update db),
            then we assume that the current direcotry represent the record that was in the DB and it will be
            and rebuild only if the rewrite flag is true.

            If another dispersion with the same parameters exists in the DB under a different name,
            then the current workflow will be added only if the allowDuplicate flag is True.

            If the workflow exists in the DB under a different name, then raise a FileExistsError.

        Parameters
        ----------
        hermes_dispersionWorkflow : hera.simulations.openFOAM.Workflow_StochasticLagrangianSolver
            An instance of the workflow to create.

        updateDB : bool
            If true, use the input there is inconsistencies between DB and  the input workflow.

            If the case directory already exists, update DB and relink it only if  --rewrite flag also exists.

        exportFromDB :  bool
            If true, then use the DB record if there is inconsistencies with the input workflow

        allowDuplicate : bool
            Add to DB even if the same workflow exists under a different name

        rewrite : bool
            Rebuild the directory if it exists, found inconsistencies
        Returns
        -------
            The workflow (synchronized with the DB).
            The disk will also be definitly synchronized only if rewrite flag is true.
        """
        if (updateDB and exportFromDB):
            err = "Cannot use both --updateDB and --exportFromDB"
            self.logger.error(err)
            raise ValueError(err)

        self.logger.execution(f"----- Start -----")

        if hermes_dispersionWorkflow.name is None:
            self.logger.error("Must set the name property in the  dispersionWorkflow object")
            raise ValueError("Must set the name property in the  dispersionWorkflow object")

        dispersionDirectoryName = os.path.abspath(os.path.join(self.toolkit.filesDirectory,hermes_dispersionWorkflow.name))

        self.logger.info(f"The dispersion workflow {hermes_dispersionWorkflow.name} and the dispersionFlowField Name {hermes_dispersionWorkflow.dispersionFlowField}.")

        self.logger.debug(f"Trying to find the dispersionFlowField in the DB to determine the directory")
        dispersionFlowFieldDocumentList = self.toolkit.getWorkflowDocumentFromDB(hermes_dispersionWorkflow.dispersionFlowField)
        if len(dispersionFlowFieldDocumentList) == 0:
            self.logger.error(f"Could not find dispersion workflow {hermes_dispersionWorkflow.dispersionFlowField} in DB.")
            self.logger.debug(f"Trying again to look for the directory of the dispersion flow in the DB")
            dispersionFlowFieldDocumentList = self.toolkit.getWorkflowDocumentFromDB(os.path.abspath(hermes_dispersionWorkflow.dispersionFlowField))
            if len(dispersionFlowFieldDocumentList) == 0:
                err = f"The {hermes_dispersionWorkflow.dispersionFlowField} is not in the DB. Add it before you can continue"
                self.logger.critical(err)
                raise ValueError(err)

        self.logger.execution(f"Found dispersion flow field in DB under the name {dispersionFlowFieldDocumentList[0].desc['workflowName']}.")
        if hermes_dispersionWorkflow.dispersionFlowField != dispersionFlowFieldDocumentList[0].desc['workflowName']:
            self.logger.debug(f"The current dispersion name is {hermes_dispersionWorkflow.dispersionFlowField} and probably a directory. Updating to {dispersionFlowFieldDocumentList[0].desc['workflowName']}")
            hermes_dispersionWorkflow.dispersionFlowField = dispersionFlowFieldDocumentList[0].desc['workflowName']

        dispersionFlowFieldName = dispersionFlowFieldDocumentList[0].desc['workflowName']
        dispersionFlowFieldDirectory = dispersionFlowFieldDocumentList[0].resource
        self.logger.info(f"Using dispersion flow field {dispersionFlowFieldName} with directory {dispersionFlowFieldDirectory}")

        updateWorkflow = False
        foundInconsistency = None
        needDBAdd = False

        # 1. Check if the dispersion workflow name is in the DB under a different name.
        self.logger.execution("Querying DB to see if the a workflow exists with a similar name")
        doc_nameQueryList = self.toolkit.getWorkflowDocumentFromDB(hermes_dispersionWorkflow.name)
        if len(doc_nameQueryList) > 0:
            self.logger.execution("Found a workflow with the same name. Check for inconsistency between DB and input")

            doc_similarWorkflow = self.toolkit.getHemresWorkflowFromDocument(doc_nameQueryList[0])
            res = self.toolkit.compareWorkflowsObj([doc_similarWorkflow, hermes_dispersionWorkflow])
            if len(res.columns) == 1:
                self.logger.info(f"The workflow in the DB is identical to the disk.")
                hermes_dispersionWorkflow = hermes_dispersionWorkflow
            else:
                self.logger.execution(f"The workflow in the DB is different than the disk: \n\n {res}")
                if exportFromDB:
                    hermes_dispersionWorkflow = doc_similarWorkflow
                elif updateDB:
                    self.logger.execution("Updating the DB with the workflow from disk")
                    hermes_dispersionWorkflow = hermes_dispersionWorkflow
                    updateWorkflow = True
                    foundInconsistency = res
                else:
                    err = f"The workflow in {hermes_dispersionWorkflow} file and in the DB are different. Use the --updateDB to update the DB or --exportFromDB to rewrite the file"
                    self.logger.error(err)
                    raise ValueError(err)
        else:
            self.logger.execution("Record not found. Querying DB to see if that workflow exists under a different name")

            otherWorkflows = self.toolkit.getHermesWorkflowFromDB(hermes_dispersionWorkflow)
            self.logger.debug(f"Found {len(otherWorkflows)} records that match workflow.")
            if len(otherWorkflows) > 0:
                foundNames = [x.desc['workflowName'] for x in otherWorkflows]
                self.logger.info(f"The workflows {','.join(foundNames)} match the workflow that was given as input, but have different names.")
                if not allowDuplicate:
                    err = f"Found the input workflow under the names {','.join(foundNames)}. The supplied name is {hermes_dispersionWorkflow.name}. " \
                          f"Force the addition of a duplicate case by using the --allowDuplicate flag"
                    self.logger.error(err)
                    raise FileExistsError(err)
                else:
                    self.logger.execution("--allowDuplicate was supplied, so allowing duplication addition to DB")
                    needDBAdd = True

        self.logger.execution(f"Check if dispersion {dispersionDirectoryName} already exists. If it is remove for fresh start only if it needs rewrite and the  --rewrite was supplied")
        if os.path.isdir(dispersionDirectoryName):
            self.logger.debug("Directory already exists")
            if foundInconsistency is not None:
                self.logger.execution(
                    f"The directory {dispersionDirectoryName} exists and differ than DB. Rewriting is needed as disk version was selected with --updateDB flag")
                if rewrite:
                    self.logger.info(f"rewrite flag exists: removing the directory {dispersionDirectoryName}")
                    shutil.rmtree(dispersionDirectoryName)
                else:
                    err = f"The dispersion directory {dispersionDirectoryName} already exists, and the input workflow does not match DB: \n{res}\n\n Use --rewrite to force removing and recreating"
                    self.logger.error(err)
                    raise ValueError(err)

        # 3. Create the dispersion case and link to the workflow.
        self.logger.info(f"Creating dispersion case {dispersionDirectoryName}  and linking to {dispersionFlowFieldDirectory}")
        self.toolkit.stochasticLagrangian._createAndLinkDispersionCaseDirectory(dispersionDirectoryName,dispersionFlowDirectory=dispersionFlowFieldDirectory)

        # 5. Add/update to the DB.
        self.logger.execution("add to DB if not exist. If exists, check what needs to be updated (dispersionFlowField or the entire code)")
        if needDBAdd:
            groupName = hermes_dispersionWorkflow.name.split("_")[0]
            groupID = hermes_dispersionWorkflow.name.split("_")[1]

            self.logger.execute(f"Adding new document with group {groupName} and {groupID}")
            self.toolkit.addSimulationsDocument(resource=dispersionDirectoryName,\
                                                dataFormat=datatypes.STRING,\
                                                type=self.toolkit.DOCTYPE_WORKFLOW,\
                                                desc=dict(\
                                                    groupName=groupName,\
                                                    groupID=groupID,\
                                                    workflowName=dispersionFlowFieldName,\
                                                    workflowType=hermes_dispersionWorkflow.workflowType,\
                                                    workflow=hermes_dispersionWorkflow.json,\
                                                    parameters=hermes_dispersionWorkflow.parametersJSON\
                                                    )\
                                                )
        elif updateWorkflow:
                self.logger.execute("Updating workflow")
                doc = doc_nameQueryList[0]
                doc.desc['workflow'] = hermes_dispersionWorkflow.json
                doc.desc['parameters'] = hermes_dispersionWorkflow.parametersJSON

                jsonstr = json.dumps(doc.desc,indent=4)
                self.logger.debug(f"Saving desc \n{jsonstr}\n to DB")
                doc.save()

    def _createAndLinkDispersionCaseDirectory(self, dispersionDirectory, dispersionFlowDirectory):
        dispersionDirectory = os.path.abspath(dispersionDirectory)
        dispersionFlowDirectory = os.path.abspath(dispersionFlowDirectory)
        self.logger.info(f"Creating dispersion at {dispersionDirectory} with base flow {dispersionFlowDirectory}")

        if not os.path.isfile(dispersionFlowDirectory):
            err = f"The {dispersionFlowDirectory} is not not a directory."
            self.logger.error(err)
            raise ValueError(err)

        constantDir = os.path.join(dispersionDirectory, "constant")
        systemDir = os.path.join(dispersionDirectory, "system")
        os.makedirs(constantDir, exist_ok=True)
        os.makedirs(systemDir, exist_ok=True)

        self.logger.debug(f"Copying constant and system from the base flow")
        shutil.copytree(os.path.join(dispersionFlowDirectory, "constant"), os.path.join(dispersionDirectory, "constant"))
        shutil.copytree(os.path.join(dispersionFlowDirectory, "system"), os.path.join(dispersionDirectory, "system"))

        for proc in glob.glob(os.path.join(dispersionFlowDirectory, "processor*")):
            self.logger.execution(f"Working on processor {proc}")
            fullpath = os.path.abspath(os.path.join(proc, "constant", "polyMesh"))
            destination = os.path.join(dispersionDirectory, os.path.basename(proc), "constant", "polyMesh")
            os.makedirs(os.path.dirname(destination), exist_ok=True)
            self.logger.debug(f"\t Linking: ln -s {fullpath} {destination}")
            self.logger.debug(f"\t Linking root case : ln -s {os.path.abspath(proc)} {os.path.join(dispersionDirectory, os.path.basename(proc))}/rootCase")
            os.system(f"ln -s {fullpath} {destination}")
            os.system(f"ln -s {os.path.abspath(proc)} {os.path.join(dispersionDirectory, os.path.basename(proc))}/rootCase")

            # create the 0 directory in all processors.
            os.makedirs(os.path.join(dispersionDirectory, os.path.basename(proc), '0'), exist_ok=True)
            self.logger.debug(f"Creating the 0 time in {os.path.join(dispersionDirectory, os.path.basename(proc), '0')} ")

        # linking the decpomposePar dict from the root.
        # root_decomposePar = os.path.abspath(os.path.join(dispersionFlow,"system","decomposeParDict"))
        # decompose_dest    = os.path.abspath(os.path.join(dispersionDirectory,"system"))
        # os.system(f"ln -s {root_decomposePar} {decompose_dest}")

        # linking the rootCase in the root directory of the dispersion dispersionDirectory.
        self.logger.debug(f"Linking the root case: ln -s {dispersionFlowDirectory} {os.path.join(dispersionDirectory, 'rootCase')}")
        os.system(f"ln -s {dispersionFlowDirectory} {os.path.join(dispersionDirectory, 'rootCase')}")

        # create the 0 directory in the root.
        self.logger.debug(f"Making the 0 in {os.path.join(dispersionDirectory, '0')}")
        os.makedirs(os.path.join(dispersionDirectory, '0'), exist_ok=True)

        self.logger.info("... Done")


    @property
    def sourcesTypeList(self):
        return [x.split("_")[1] for x in dir(self) if "makeSource" in x and "_" in x]


    def makeSource(self, x, y, z, nParticles, type="Point", **kwargs):
        slist = self.sourcesTypeList
        if type not in slist:
            raise ValueError(f"The type must be [{','.join(slist)}]. Got {type} instead")

        return getattr(self, f"makeSource_{type}")(x, y, z, nParticles, **kwargs)


    def makeSource_Point(self, x, y, z, nParticles, **kwargs):
        return pandas.DataFrame(
            {"x": [x for i in range(nParticles)], "y": [y for i in range(nParticles)], "z": [z for i in range(nParticles)]})


    def makeSource_Circle(self, x, y, z, nParticles, radius, distribution="uniform", **kwargs):
        Rs = getattr(numpy.random, distribution)(0, radius, nParticles)
        thetas = numpy.random.uniform(0, 2 * numpy.pi, nParticles)
        xs = []
        ys = []
        for R, theta in zip(Rs, thetas):
            xs.append(x + R * numpy.cos(theta))
            ys.append(y + R * numpy.sin(theta))
        return pandas.DataFrame({"x": xs, "y": ys, "z": [z for i in range(nParticles)]})


    def makeSource_Sphere(self, x, y, z, nParticles, radius, distribution="uniform", **kwargs):
        Rs = getattr(numpy.random, distribution)(0, radius, nParticles)
        thetas = numpy.random.uniform(0, 2 * numpy.pi, nParticles)
        phis = numpy.random.uniform(0, 2 * numpy.pi, nParticles)
        xs = []
        ys = []
        zs = []
        for R, theta, phi in zip(Rs, thetas, phis):
            xs.append(x + R * numpy.sin(theta) * numpy.cos(phi))
            ys.append(y + R * numpy.sin(theta) * numpy.sin(phi))
            zs.append(z + R * numpy.cos(theta))
        return pandas.DataFrame({"x": xs, "y": ys, "z": zs})


    def makeSource_Cylinder(self, x, y, z, nParticles, radius, height, horizontalDistribution="uniform",
                            verticalDistribution="uniform", **kwargs):
        Rs = getattr(numpy.random, horizontalDistribution)(0, radius, nParticles)
        Hs = getattr(numpy.random, verticalDistribution)(-height / 2, height / 2, nParticles)
        thetas = numpy.random.uniform(0, 2 * numpy.pi, nParticles)
        xs = []
        ys = []
        zs = []
        for R, theta, H in zip(Rs, thetas, Hs):
            xs.append(x + R * numpy.cos(theta))
            ys.append(y + R * numpy.sin(theta))
            zs.append(z + H)
        return pandas.DataFrame({"x": xs, "y": ys, "z": zs})


    def makeSource_Rectangle(self, x, y, z, nParticles, lengthX, lengthY, rotateAngle=0, **kwargs):
        xdist = numpy.random.uniform(-lengthX / 2, lengthX / 2, nParticles)
        ydist = numpy.random.uniform(-lengthY / 2, lengthY / 2, nParticles)
        xs = []
        ys = []
        for i in range(nParticles):
            xs.append(x + xdist[i] * numpy.cos(rotateAngle) + ydist[i] * numpy.sin(rotateAngle))
            ys.append(y - xdist[i] * numpy.sin(rotateAngle) + ydist[i] * numpy.cos(rotateAngle))
        return pandas.DataFrame({"x": xs, "y": ys, "z": [z for i in range(nParticles)]})


    def makeSource_Cube(self, x, y, z, nParticles, lengthX, lengthY, lengthZ, rotateAngle=0, **kwargs):
        xdist = numpy.random.uniform(-lengthX / 2, lengthX / 2, nParticles)
        ydist = numpy.random.uniform(-lengthY / 2, lengthY / 2, nParticles)
        zdist = numpy.random.uniform(-lengthZ / 2, lengthZ / 2, nParticles)
        xs = []
        ys = []
        zs = []
        for i in range(nParticles):
            xs.append(x + xdist[i] * numpy.cos(rotateAngle) + ydist[i] * numpy.sin(rotateAngle))
            ys.append(y - xdist[i] * numpy.sin(rotateAngle) + ydist[i] * numpy.cos(rotateAngle))
            zs.append(z + zdist[i])
        return pandas.DataFrame({"x": xs, "y": ys, "z": zs})
