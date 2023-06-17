import json
import os
import glob
from .datalayer.OFObjects import OFField,OFMeshBoundary
from .datalayer.hermesWorkflow import Workflow_Flow
from .datalayer.OFObjects import ofObjectHome
from . import DECOMPOSED_CASE,RECONSTRUCTED_CASE,TYPE_VTK_FILTER
from ..hermesWorkflowToolkit import workflowToolkit,workflowsTypes
from ...utils.query import dictToMongoQuery
from ..hermesWorkflowToolkit import workflowToolkit
from ...datalayer import datatypes
from .analysis.VTKPipeline import VTKpipeline
import shutil


class OFToolkit(workflowToolkit):
    """
        The goal of this toolkit is to provide the functions that are required to run workflows.
        and to mange the workflows in the DB.

        This toolkit might relay on the hermes project in order to manipulate the nodes
        of the workflow. (TBD).

    """
    TIME_STEADYSTATE = "steadyState"
    TIME_DYNAMIC = "dynamic"

    _ofObjectHome = None

    @property
    def ofObjectHome(self):
        return self._ofObjectHome

    def __init__(self, projectName, filesDirectory=None):
        super().__init__(projectName=projectName,
                         filesDirectory=filesDirectory,
                         toolkitName="OFworkflowToolkit")

        self._ofObjectHome = ofObjectHome()
        self._analysis = analysis(self)

    def processorList(self,caseDirectory):
        """
            Returns the list of processors directories in the case
        Parameters
        ----------
        caseDirectory : str
            Path to the directory.

        Returns
        -------

        """
        return [os.path.basename(proc) for proc in glob.glob(os.path.join(caseDirectory, "processor*"))]

    def getHermesWorkflow_Flow(self,workflowfile):
        """
            Returns the workflow of the requested JSON file.
        Parameters
        ----------
        workflowfile

        Returns
        -------

        """
        return Workflow_Flow(workflowfile)

    def getMesh(self,caseDirectory,parallel=True,time=0):
        """
            Reads the mesh from the mesh directory.

            Reads the decomposed case if it exists and parallel is true,
            otherwise, reads just the single case.

            Unfortunately, we have to use the OF postProcess utility in order to interpolate the
            mesh points to their centers.

        Parameters
        ----------
            caseDirectory: str
                    The path to the case. Should be absolute in order to determine whether we need to add the -case tot he postProcess.

            parallel: bool
                    If parallel case exists, read it .
        Returns
        -------
            pandas dataframe with the points in the columns             x,y,z
            the index column (don't mix up with the index of pandas)  is the sequential number of the point.

            If the case is decomposed, return processorNumber and index columns.
            The index is the internal order in the processor.
        """

        # 1. Run the postProcess utility to set the cell centers

        self.logger.info(f"Start. case {caseDirectory}. Current directory is : {os.getcwd()}.")

        casePointer = "" if caseDirectory == os.getcwd() else f"-case {caseDirectory}"

        useParallel= False
        if parallel:
            self.logger.debug(f"Attempt to load parallel case")
            # Check if the case is decomposed, if it is, run it.
            proc0dir = os.path.join(caseDirectory,"processor0")

            if os.path.exists(proc0dir):
                self.logger.debug(f"Found parallel case, using decomposed case")
                useParallel = True
            else:
                self.logger.debug(f"parallel case NOT found. Using composed case")

        # Calculating the cell centers
        checkPath = os.path.join(caseDirectory,"processor0",str(time),"C") if useParallel else os.path.join(caseDirectory,str(time),"C")
        parallelExec = "-parallel" if useParallel else ""
        caseType = "decomposed" if useParallel else "composed"
        if not os.path.exists(checkPath):
            self.logger.debug(f"Cell centers does not exist in {caseType} case. Calculating...")
            os.system(f"foamJob {parallelExec} -wait postProcess -func writeCellCentres {casePointer}")
            self.logger.debug(f"done: foamJob {parallelExec} -wait postProcess -func writeCellCentres {casePointer}")
            if not os.path.exists(checkPath):
                self.logger.error("Error running the writeCellCentres. Check mesh")
                raise RuntimeError("Error running the writeCellCentres. Check mesh")
        else:
            self.logger.debug(f"Cell centers exist in {caseType} case.")

        cellCenters = OFField(name="C",dimensions="",componentNames=['x','y','z'])
        self.logger.debug(f"Loading the cell centers in time {time}. Usint {caseType}")
        ret =  cellCenters.load(caseDirectory,times=time,parallelCase=useParallel)

        self.logger.info(f"--- End ---")
        return ret

    def runWorkflow(self,workflowNameOrJSON,saveMode):
        """
            Runs a workflow. Checks if it is in the DB first for the appropriate saveMode.


        :param workflowNameOrJSON:
        :return:
        """
        pass

    def getCaseBoundary(self,caseDirectory:str,filterProcessor:bool = True):
        """
            Returns a list of the boundary fields.
        Parameters
        ----------
        caseDirectory : str
                the directory of the flow.

        filterProcessor : bool
                If true, filter the processor* boundary out.

        Returns
        -------
            A list of field names.
        """

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
                        type: "directory|dispersionFlowFieldName|workflowFile|ID"
                        parameters:  {
                                .. depend on the type
                        },
                        "useTime" : [float], optional, default: last timestep
                        timeStep : [float], optional,  default: 0
                        dispersionDuration : [float], required.
                        copyMesh : [bool], default false

                    },
                    "dispersionField" : {

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
        def findTimeDirectory(basedir,timeDirectory,logger):
            orig_proc = os.path.join(basedir, str(timeDirectory))
            # The directories in round times are int and not float.
            # that is 10 and not 10.0. Therefore, we must check and correct if does not exist.
            if not os.path.exists(orig_proc):
                logger.debug("Source does not exist, probably due to float/int issues. Recreating the source with int format")
                orig_proc = os.path.join(basedir, str(int(timeDirectory)))
                if not os.path.exists(orig_proc):
                    logger.error(f"Source {orig_proc} does not exist... make sure the simulation is OK.")
                    raise ValueError(f"Source {orig_proc} does not exist... make sure the simulation is OK.")

        self.logger.info("-------- createDispersionFlowField: Start ---------")
        originalFlow = flowData['originalFlow']

        # 1. Get the case directory of the original flow.
        self.logger.debug(f"tying to find the flow {originalFlow['source']} in the DB")
        docList = self.getSimulationDocumentFromDB(originalFlow['source'])

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

        if dynamicType==self.TIME_STEADYSTATE:
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
            self.logger.debug(f"Trying to find the simulation in the database. The run is: \n {json.dumps(querydict,indent=4)")
            docList = self.getSimulationsDocuments(type=workflowsTypes.OF_FLOWDISPERSION.value, **dictToMongoQuery(querydict))
            self.logger.debug(f"Found {len(docList)} in the database")
        else:
            self.logger.debug(f"Running without DB support, so does not query the db ")
            docList = []

        if len(docList) == 0 or overwrite:
            self.logger.info(f"Flow field not found, creating new and adding to the DB. ")
            ofhome = ofObjectHome()
            if len(docList) == 0:
                self.logger.debug("Find the first name that is free")
                groupID = 1
                dispersionFlowFieldName = workflowGroup
                while not os.path.exists(os.path.join(self.FilesDirectory,f"{dispersionFlowFieldName}_{groupID}")):
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
            dispersionFlowFieldDirectory = os.path.abspath(os.path.join(self.FilesDirectory,dispersionFlowFieldName))
            self.logger.info(f"Creating Dispersion flow simulation {dispersionFlowFieldName} in {os.path.abspath(self.FilesDirectory)}")
            os.makedirs(dispersionFlowFieldDirectory,exist_ok=True)

            self.logger.info(f"Creating the flow specific fields in the flow needed for the dispersion")
            dispersionFields = flowData['dispersionFields']

            meshBoundary = OFMeshBoundary(originalFlowCaseDir).getBoundary()
            newFieldList = []
            for dispersionFieldName in dispersionFields.keys():
                fieldDimensions = dispersionFields[dispersionFieldName].get("dimensions",None)
                fieldComponents = dispersionFields[dispersionFieldName].get("components", None)
                fieldBoundary   = dispersionFields[dispersionFieldName].get("boundaryField", dict())
                self.logger.debug(f"Creating the flow specific field: {dispersionFieldName}. ")
                field = ofhome.getField(fieldName=dispersionFieldName,
                                        simulationType=ofObjectHome.GROUP_DISPERSION,
                                        dimensions=fieldDimensions,
                                        componentNames=fieldComponents)

                self.logger.debug(f"Adding the boundaries condition 'zerGradient' to all the boundaries that were not specified")
                for bnd in meshBoundary:
                    fieldBoundary.setdefault(bnd,dict(type='zeroGradient'))

                newFieldList.append( (field,fieldBoundary) )

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
                if (dynamicType == self.TIME_STEADYSTATE) and (orig_time == timeList[-1]):
                    orig_time = timeList[0]

                self.logger.info(f"Mapping Dispersion time {dest_time} to {orig_time}")

                for origDir in origDirsList:
                    # The directories in round times are int and not float.
                    # that is 10 and not 10.0. Therefore, we must check and correct if does not exist.
                    orig_proc_timestep = findTimeDirectory(basedir=origDir, timeDirectory=str(orig_time), logger=self.logger)

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


                for origDir in origDirsList:
                    procName = os.path.split(origDir)[-1]
                    self.logger.debug(f"Writing the flow specific field {dispersionFieldName} to processor {procName} . ")
                    for dest_time in timeList:
                        field.emptyParallelField(caseDirectory=dispersionFlowFieldDirectory,
                                                 timeName=str(dest_time),
                                                 processor=procName,
                                                 boundaryField=fieldBoundary,
                                                 data=dispersionFields[dispersionFieldName].get("internalField"))

                # We should look into it more closly, why parallel case doesn't recognize the time steps of the
                # processors. For now, just create these directories in the main root as well.
                if parallelOriginal:
                    for dest_time in timeList:
                        os.makedirs(os.path.join(dispersionFlowFieldDirectory,str(dest_time)),exist_ok=True)


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
                    self.addSimulationsDocument(resource=dispersionFlowFieldDirectory, type=workflowsTypes.OF_FLOWDISPERSION.value, dataFormat=datatypes.STRING, desc=querydict)

            ret = dispersionFlowFieldDirectory
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
        return parameters['name'],os.path.basename(parameters['name'])

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
        docList = self.getSimulationsDocuments(workflowName=parameters['name'],type=workflowToolkit.DOC_TYPE)
        if len(docList) >0:
            if len(docList) > 1:
                import warnings
                warnings.warn(f"Found more than 1 simulation with the name {parameters['name']}. Return the first one")
            return docList[0].resource,parameters['name']
        else:
            raise ValueError(f"Cannot find flows with the name {parameters['name']}. Use hera-OF-flows list to see the names of existing simulations.old.")



class analysis:
    """
        The analysis of the OpenFOAM.

        Handles both eulerian (flow) and lagrangian fields.
    Returns
    -------

    """

    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self,datalayer):
        """
            Initializes the analysis layer.
        Parameters
        ----------
        datalayer : OFToolkit
                The datalayer to use.

        """
        self._datalayer = datalayer


    def makeVTKPipeline(self, nameOrWorkflowFileOrJSONOrResource, vtkPipeline, caseType=DECOMPOSED_CASE, servername=None, fieldNames=None):
        """
            Creates a new VTK pipeline from the simulation.

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

        vtkPipeline : str,dict

            The JSON that describes the VTK pipeline.

        Returns
        -------
            analysis.VTKPipeline.VTKpipeline
        """
        return VTKpipeline(datalayer=self.datalayer,
                           pipelineJSON=vtkPipeline,
                           nameOrWorkflowFileOrJSONOrResource=nameOrWorkflowFileOrJSONOrResource,
                           caseType = caseType,
                           serverName = servername,
                           fieldNames = fieldNames)

    def getFiltersDocuments(self,nameOrWorkflowFileOrJSONOrResource,filterName=None):
        """
                Returns the cache documents of the filters.

                nameOrWorkflowFileOrJSONOrResource identifies the simulation and can be:
                    - The path to the directory,
                    - The name of the simulation
                    - The workflow file.
                    - The workflow dict.

                If filterName is None, will return all the filters of the case.
                If nameOrWorkflowFileOrJSONOrResource is None, will return all the

        Parameters
        ----------
        filterName : str
                The name of the filter to retrieve
        nameOrWorkflowFileOrJSONOrResource : str
                The identifier of the simulation.


        Returns
        -------
            list of DB documents.
        """
        wrkflow = self.datalayer.getSimulationDocumentFromDB(nameOrWorkflowFileOrJSONOrResource)
        if wrkflow is None:
            raise ValueError(f"The case {nameOrWorkflowFileOrJSONOrResource} was not found in the project {self.datalayer.projectName}")
        simulationName = wrkflow['desc']['simulationName']

        qry = dict(simulationName=simulationName)
        if filterName is not None:
            qry['filterName'] = filterName

        return self.datalayer.getCacheDocuments(type=TYPE_VTK_FILTER,**qry)