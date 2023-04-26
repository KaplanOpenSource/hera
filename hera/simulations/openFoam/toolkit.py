from ...toolkit import abstractToolkit
import os
import glob
from .datalayer.OFObjects import OFField,OFMeshBoundary
from .datalayer.hermesWorkflow import Workflow_Flow
from .datalayer.OFObjects import ofObjectHome
from . import DECOMPOSED_CASE,RECONSTRUCTED_CASE,TYPE_VTK_FILTER
from ..hermesWorkflowToolkit import workflowToolkit,simulationTypes
from ...utils import loadJSON
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

    def prepareFlowFieldForDispersion(self, flowData, simulationGroup:str,overwrite:bool=False, useDBSupport: bool = True):
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

                        To implement the rst, we might need to extend the scripts to use multiple baseFlow simulations.old.

            simulationGroup : str
                    The name of the group.

            useDBSupport : bool
                    If False, does not access the DB to check if the flow exists.
                    and does not add it. Used in the hermes module, when each dispersion case has its own flow (as it is soft link, it does not take space).

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
            # find the closest TS.
            request = float(useTime)
            uts = TS[min(range(len(TS)), key=lambda i: abs(float(TS[i]) - request))]

        fromTime = fromTimestep

        self.logger.debug(f"Using Time step {uts} for Steady state")

        ## Now we should find out if there is a similar run.
        # 1. same original run.
        # 2. same Hmix/ustar/... other fields.
        # 3. the requested start and end are within the existing start and end.

        if useDBSupport:

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
        else:
            self.logger.debug(f"Running without DB support, so does not query the db ")
            docList = []

        if len(docList) == 0 or overwrite:
            self.logger.info(f"Flow field not found, creating new and adding to the DB. ")
            ofhome = ofObjectHome()

            if useDBSupport:
                if len(docList) == 0:
                    self.logger.debug("Getting a new name")
                    groupID,simulationName = self.findAvailableName(simulationGroup, simulationType=simulationTypes.OF_FLOWDISPERSION.value)
                else:
                    resource = docList[0].resource
                    simulationName = docList[0].desc['name']
                    groupID= docList[0].desc['groupID']
                    self.logger.debug(f"The flow exists with the name {simulationName}, so deleting {resource}, and writing over it")
                    shutil.rmtree(resource)

                self.logger.info(
                    f"Creating new name for the flow simulation.... Got ID {groupID} with simulation name {simulationName}")
            else:
                groupID = "1"
                simulationName = simulationGroup
                self.logger.info(f"Creating flow simulation with {simulationName}")

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
                    self.logger.debug(f"path {dest_general} exists... removing before copy")
                    shutil.rmtree(dest_general)
                shutil.copytree(orig_general, dest_general)

            self.logger.info(f"Copy the parallel directories")
            for proc in glob.glob(os.path.join(orig,"processor*")):

                orig_proc = os.path.join(proc, str(uts))
                # The directories in round times are int and not float.
                # that is 10 and not 10.0. Therefore, we must check and correct if does not exist.
                if not os.path.exists(orig_proc):
                    self.logger.debug("Source does not exist, probably due to float/int issues. Recreating the source with int format")
                    orig_proc = os.path.join(proc, str(int(uts)))
                    if not os.path.exists(orig_proc):
                        self.logger.error(f"Source {orig_proc} does not exist... make sure the simulation is OK.")
                        raise ValueError(f"Source {orig_proc} does not exist... make sure the simulation is OK.")

                for dest_time in [fromTime, maximalDispersionTime]:
                    self.logger.debug(f"Handling {proc} - time {dest_time}")
                    dest_proc = os.path.join(dest,os.path.basename(proc),str(dest_time))
                    if os.path.exists(dest_proc):
                        self.logger.debug(f"path {dest_proc} exists... removing before copy")
                        shutil.rmtree(dest_proc)
                    shutil.copytree(orig_proc,dest_proc)

                self.logger.info(f"Copying the mesh")
                if copyMesh:
                    orig_constant = os.path.join(proc,"constant")
                    dest_constant = os.path.join(dest,os.path.basename(proc),"constant")
                    shutil.copytree(orig_constant, dest_constant)
                else:
                    fullpath = os.path.abspath(os.path.join(proc, "constant", "polyMesh"))
                    destination = os.path.join(dest, os.path.basename(proc), "constant", "polyMesh")
                    if not os.path.exists(destination):
                        self.logger.debug(f"Linking {fullpath} -> {destination}")
                        os.makedirs(os.path.dirname(destination), exist_ok=True)
                        os.system(f"ln -s {fullpath} {destination}")

            self.logger.info(f"Creating the flow specific fields in the flow needed for the dispersion")
            dispersionFields = flowData['dispersionFields']

            meshBoundary = OFMeshBoundary(orig).getBoundary()

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


                for proc in glob.glob(os.path.join(orig, "processor*")):
                    procName = os.path.split(proc)[-1]
                    self.logger.debug(f"Writing the flow specific field {dispersionFieldName} to processor {procName} . ")
                    for dest_time in [fromTime, maximalDispersionTime]:
                        field.emptyParallelField(caseDirectory=dest,
                                                 timeName=str(dest_time),
                                                 processor=procName,
                                                 boundaryField=fieldBoundary,
                                                 data=dispersionFields[dispersionFieldName].get("internalField"))

                # We should look into it more closly, why parallel case doesn't recognize the time steps of the
                # processors. For now, just create these directories in the main root as well.
                for dest_time in [fromTime, maximalDispersionTime]:
                    os.makedirs(os.path.join(dest,str(dest_time)),exist_ok=True)


            self.logger.info("Finished creating the flow field for the dispersion. ")
            if useDBSupport:
                self.logger.info("Adding to the database.")
                if len(docList) ==0:
                    self.logger.debug("Updating the metadata of the record with the new group ID and simulation name")
                    querydict.update(dict(
                        groupID=groupID,
                        name=simulationName,
                    ))
                    self.logger.debug("Adding record to the database")
                    self.addSimulationsDocument(resource=dest,type=simulationTypes.OF_FLOWDISPERSION.value,dataFormat=datatypes.STRING,desc=querydict)

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
                import warnings
                warnings.warn(f"Found more than 1 simulation with the name {parameters['name']}. Return the first one")
            return docList[0].resource
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