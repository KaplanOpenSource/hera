import numpy
import os
import glob
import shutil
from distutils.dir_util import copy_tree
from hera import toolkitHome
from hermes.workflow import workflow,hermesNode
from ....utils import loadJSON
from ....utils import loggedObject
from ....datalayer import Project
from itertools import product
import json
from shutil import copyfile
from ...openFoam import ofObjectHome

try:
    from hermes import workflow
except ImportError:
    raise ImportError("Cannot use this module without hermes... Install it. ")


class abstractWorkflow(workflow):
    """
            An abstract specialization of the hermes workflow to the
            openfoam workflow.

    """

    DOCTYPE_DISPERSION_FLOWFIELD = "dispersionFlowField" # The flow field of the dispersion.
    DOCTYPE_FLOWFIELD = "FlowField" # calculation of the flow field.
    DOCTYPE_DISPERSION = "dispersion" # The dispersion itself.

    @property
    def parameters(self):
        return self['Parameters']

    def __init__(self,workflowJSON):
        super().__init__(workflowJSON=workflowJSON)
        self.logger = loggedObject(loggerName=None).logger

    def buildCaseExecutionScript(self,caseDirectory,configuration):
        """
            Writes the allRun file that executes the workflow and the allClean file
            that cleans the case to the case directory.

        Parameters
        ----------
        caseDirectory: str
            The directory to write the files to.

        configuration: dict
            A configuration file that is used to build the execution of the node.
            Specifically, includes a node 'caseExecution' with the structure:

             {
                "parallelCase": true|false             # Write the execution for parallel execution.
                "runFile": [
                            --------------- A list of nodes.
                  {
                    "name": "blockMesh",               # The name of the program to execute.
                    "couldRunInParallel": false,       # Write as parallel (only if parallel case is True).
                    "parameters": null                 # Parameters for each run.
                  }
                ]
                            --------------- A list of nodes.
             }

        """
        self.logger.info("---- Start ----")
        try:
            execConfiguration = configuration['caseExecution']
        except KeyError:
            err = "'caseExecution' node not found in configuration"
            self.logger.error(err)
            raise ValueError(err)


        isParallel = execConfiguration['parallelCase']
        self.logger.info(f"Running this case as parallel? {isParallel}")

        execLine = ""

        for execNode in execConfiguration['runFile']:
            self.logger.execution(f"Processing Node {execNode['name']}")

            parallelFlag = "-parallel" if (isParallel and execNode['couldRunInParallel']) else ""
            progName = execNode['name']
            parameters = execNode.get('parameters',None)

            if parameters is not None:
                params   = " ".join(numpy.atleast_1d(execNode['parameters']))
            else:
                params = ""

            foamJob = execNode.get("foamJob",True)

            if foamJob:
                execLine += f"foamJob {parallelFlag} -screen -wait {progName} {params}\n"
            else:
                execLine += f"{progName} {params}\n"

        allrunFile = os.path.join(caseDirectory,"Allrun")
        with open(allrunFile,'w') as execFile:
            execFile.write(execLine)
        os.chmod(allrunFile, 0o777)

        # Now write the allClean file.
        allCleanContent = """
#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

# Remove surface and features
rm -f constant/triSurface/motorBike.obj.gz > /dev/null 2>&1
rm -rf constant/extendedFeatureEdgeMesh > /dev/null 2>&1
rm -f constant/triSurface/motorBike.eMesh > /dev/null 2>&1

cleanCase
        """
        allcleanFile = os.path.join(caseDirectory,"Allclean")
        with open(allcleanFile,'w') as allclean:
            allclean.write(allCleanContent)

        os.chmod(allcleanFile, 0o777)

    @classmethod
    def baseTemplateHandler_directory(cls,jsonConfiguration):
        """
            Check if the directory is already the db.
            If it is, then return the id.

        Parameters
        ----------
        json: str
            The json of the directory:
                {
                    name : ...
                }


        Returns
        -------

        """
        return jsonConfiguration['name']

class Workflow_Flow(abstractWorkflow):
    """
        This class manages the hermes workflow of openFOAM that is designated to
        calculate flow fields.


        It contains additional functions that are specific to openFOAM workflows
        such as controlDict, fvOptions, fvSchemes, and fvSolutions that return the relevant
        nodes.

        Also, it contains a function change the from composed to decomposed run.

        Each node has a list of default flags that could be modified when the case execution script is built.

        The case is considered as running in parallel if there is an execution that is required to run
        in parallel (for example the decompose par).

    """
    def __init__(self,workflowJSON,parallelNodes=None):
        """
            Initializes a openFOAM hermes workflow.

            workflowExecutionJSON is a JSON that describes how to execute the workflow.

            Its structure is an ordered list of dict.
            Each dict represents a command line that has to be executed to run the case.
                [
                    {
                        ... node description ...
                    },

                ]

            The node description is:
            {
                name:       The name of the program to execute.
                parallel:   Is it necessary for parallel run.
                parameters: a list of string to append as parameters.
            }


        Parameters
        ----------
        workflow : str, dict, file type.
            Can be the file name, a string of the JSON, or a file type.


        parallelNodes: str, list
            A name of node, or a list of nodes in the workflow that are required to build the case in parallel.
            Will be removed if the workflow is executed as a unified case.
        """
        super().__init__(workflowJSON=workflowJSON)
        self._parallelNodes = [] if parallelNodes is None else numpy.atleast_1d(parallelNodes)

        # examine here that all the nodes exist, if not - it is not a flow
        for node in ['controlDict','fvSolution','fvSchemes','blockMesh','fileWriter','defineNewBoundaryConditions']:
            if node not in self.workflowJSON['nodes']:
                raise ValueError(f"The node {node} does not exist in the flow. Not a flow workflow.")



    @property
    def controlDict(self):
        return self['controlDict']

    @property
    def fvSolution(self):
        return self['fvSolution']

    @property
    def fvScheme(self):
        return self['fvScheme']

    @property
    def blockMesh(self):
        return self['blockMesh']


    @property
    def snappyHexMesh(self):
        return self['snappyHexMesh'] if 'snappyHexMesh' in self.nodes else None

    @property
    def fileWriter(self):
        return self['fileWriter']


    @property
    def defineNewBoundaryConditions(self):
        return self['defineNewBoundaryConditions']

    def changeWorkflowToRunAsComposed(self):
        """
            Removes the decompose node from the workflow.

        Parameters
        ----------

        Returns
        -------
            None
        """

        # 1. Remove from the nodelist.
        for decmdNode in self._parallelNodes:
            if decmdNode in self.nodeList:
                self.nodeList.remove(decmdNode)
                del self.nodes[decmdNode]

    def prepareGeometry(self, workingDirectory, configurationFile,overwrite=False):
        """
            Adapts the geometry of the workflow.

        Parameters
        ----------
        caseDirectory : str
            The directory to save the STL in

        configurationFile: str,file like, dict
            Reads the JSON  from the directory.
            Could be a file like object, a path to the file, a JSON string or a dict.

        overwrite: bool
            If true, then generate the STL even it exists.
            if False, skip the genereation.

        Returns
        -------
            dict,
            The name of the object and the names that were written.

        """
        self.logger.info("-- Start --")
        configuration = loadJSON(configurationFile)

        geom_handlers = [x.split("_")[1] for x in dir(self) if x.startswith("geometryHandler")]

        for name,geometryObject in configuration['geometry'].items():
            self.logger.debug(f"Processing: \n {json.dumps(geometryObject,indent=4)}")

            # 1. Handle the gemetry-type specifics. (create STL and ect.)
            if geometryObject['source']['type'] not in geom_handlers:
                err = f"{geometryObject['source']['type']} must be one of {','.join(geom_handlers)}"
                self.logger.error(err)
                raise ValueError(err)

            handler = getattr(self, f"geometryHandler_{geometryObject['source']['type']}")
            handler(name, geometryObject, workingDirectory, configuration,overwrite)

            # 2. Change the snappyHexMesh node.
            # snappyHexNode = self.snappyHexMesh
            #
            # if snappyHexNode is None:
            #     self.logger.warning("snappyHexMesh node does not exist. Geometry might not be used!")
            # else:
            #     snappyHexNode["geometry"]['objects'][name] = geometryObject["meshing"]

        self.logger.info("-- End --")

    def prepareMesh(self, configurationFile):
        """
            Reads the configuration file and:

            Adapt the following workflow nodes:

                1. blockMesh
                2. Snappy
                        Setup the locatin In Mesh property in the snappy hex mesh (if exists).

        Parameters
        ----------
        configurationFile

        Returns
        -------

        """
        self.logger.info("-- Start --")
        configuration = loadJSON(configurationFile)


        blockMeshNode = configuration['mesh'].get('blockMesh',None)

        if blockMeshNode is not None:
            self.logger.execution("Found blockMesh Node")

            ## 1. Block mesh handler.
            blockMesh_handlers = [x.split("_")[1] for x in dir(self) if x.startswith("blockMeshHandler")]
            if blockMeshNode['type'] not in blockMesh_handlers:
                err = f"{blockMeshNode['type']} must be one of {','.join(blockMesh_handlers)}"
                self.logger.error(err)
                raise ValueError(err)

            handler = getattr(self, f"blockMeshHandler_{blockMeshNode['type']}")
            self.logger.debug(f"Calling the blockMeshHandler_{blockMeshNode['type']} handler")
            handler(configuration)

            locationNode =  configuration['mesh'].get('locationInMesh',None)

            if locationNode is not None:
                self.logger.execution("Found location Node")
                ## 2. Location in mesh handler.
                locationInMesh_handlers = [x.split("_")[1] for x in dir(self) if x.startswith("locationInMeshHandler")]
                if locationNode['type'] not in locationInMesh_handlers:
                    err = f"{blockMeshNode['locationInMesh']} must be one of {','.join(locationInMesh_handlers)}"
                    self.logger.error(err)
                    raise ValueError(err)

                handler = getattr(self, f"locationInMeshHandler_{locationNode['type']}")
                self.logger.execution(f"Calling the locationInMeshHandler_{locationNode['type']} handler")
                handler(configuration)
        else:
            self.logger.execution("blockMesh node was not Found. ")

        self.logger.info("-- End -- ")

    def prepareICandBC(self,configurationFile):
        """
            Handles the initial conditions of the file
            calls the iniital conditions handler.


        Parameters
        ----------
        configurationFile:
                The customization json file.

        Returns
        -------

        """

        for icbc in configurationFile['ICandBC']:
            icbc_type = icbc['type']
            func = getattr(self,f"ICandBCHandler_{icbc_type}")(icbc)

    ##########################################
    ##
    ##                  Handlers
    ##
        ##################


    ############## geometryHandler
    # Defines the geometry objects and sets up the snappy hex mesh.

    def geometryHandler_topography(self,regionName, geometryJSON, workingDirectory, configuration,overwrite=False):
        """
            Build the topography STL using the topography toolkit.


        Parameters
        ----------
        geometryJSON: dict
                The geometry object.
                Has the structure:

                "source" :    {
                        "type" : "topography",  <-- the type
                        "datasource" : "BNTL",  <-- The name of the datasource to use.
                        "region" : "canopy"     <-- The name of the region in the configuration file to use.
                        "dxdy" : 50             <-- The resolution of the conversion.
                    }

        workingDirectory: str
                The path to save the STL to.

        regions: dict
                A dict with the regions as polygons.

        Returns
        -------

        """
        self.logger.info(" -- Start -- ")
        tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_TOPOGRAPHY,
                                    projectName=configuration['projectName'])

        stlFileName = f"{geometryJSON['meshing']['objectName']}.stl"
        fullFileName = os.path.join(workingDirectory, stlFileName)

        region = geometryJSON['source']['region']
        regionList = configuration.get('regions', {})
        regionCoords = regionList.get(region, None)

        if regionCoords is None:
            raise ValueError(f"The region: {region} is not found. Found the regions: {','.join(regionList.keys())}")

        bx = [regionCoords['parameters']['xmin'], regionCoords['parameters']['ymin'], regionCoords['parameters']['xmax'], regionCoords['parameters']['ymax']]
        self.logger.debug(f"Found the region {region} with coords {bx}")

        if os.path.exists(fullFileName) and not overwrite:
            self.logger.info(f"File {stlFileName} in {regionName} exists. and overwrite=False")
        else:
            self.logger.info(f"Generating {stlFileName} for geometry {regionName}")
            stlcontent = tk.regionToSTL(shapeDataOrName=bx, dxdy=geometryJSON['source']['dxdy'],
                                        datasourceName=geometryJSON['source']['datasourceName'])

            self.logger.debug(f"Converting to STL complete. Now writing to the file {fullFileName}")
            with open(fullFileName, 'w') as stloutputfile:
                stloutputfile.write(stlcontent)

        ## Here we should update the snappyHexMesh node

        self.logger.debug(f"Adding the region {regionName} to the workflow")
        self.parameters["domains"] = {regionName : dict(shapeDataOrName=bx,datasourceName=geometryJSON['source']['datasourceName'])}
        self.logger.info(" -- End -- ")

    def geometryHandler_buildings(self,regionName, geometryJSON, workingDirectory, configuration):
        """
            Creates the STL of the buildings in the region using the building toolkit.

            Sets the parameter buildingsBBOX to the bounding box of the buildings area.

            We assume that the region to crop is bbox. (and given in the configuration as dict with the xmax,xmin,ymin,ymax keys).


        Parameters
        ----------
        regionName: str
            The name of the region.

            Not used here, only for completeness of the interface.

        geometryJSON: dict
                The geometry object is a dict with the following keys:

                * source - parameters for building the STL of the buildings.
                           used by the buildings toolkit.

                    has the strucutre:
                        "source" :    {
                                "type" : "buildings",  <-- the type
                                "datasource" : "BNTL",  <-- The name of the datasource to use.
                                "region" : "canopy"     <-- The name of the region in the configuration file to use.
                            }
                * meshing - parameters used by the blockMesh and the snappyHexMesh to build the code.
                    has the structure of the snappyHexMesh node ["geometry"]['objects'][name]

        configuration: dict
                The overall configuration of the case.

        Returns
        -------

        """
        tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_BUILDINGS, projectName=configuration['projectName'])

        region = geometryJSON['source']['region']
        regionList = configuration.get('regions', {})
        regionCoords = regionList.get(region, None)

        if regionCoords is None:
            raise ValueError(f"The region{region} is not found. Found the regions: {','.join(regionList.keys())}")

        bx = [regionCoords['parameters']['xmin'], regionCoords['parameters']['ymin'], regionCoords['parameters']['xmax'], regionCoords['parameters']['ymax']]

        stlFileName = f"{geometryJSON['meshing']['objectName']}.stl"
        fullFileName = os.path.join(workingDirectory, stlFileName)

        tk.regionToSTL(regionNameOrData=bx, outputFileName=fullFileName,flat=None,saveMode=toolkitHome.TOOLKIT_SAVEMODE_NOSAVE)

        self.parameters['buildingsBBOX'] = bx
        ## Here we should update the snappyHexMesh node

    def geometryHandler_indoorBuilding(self,regionName, geometryJSON, workingDirectory, configuration,overwrite=None):
        """
            Updating the snappyHexMesh node.

            Assumes that the name of the region is the name of the object in the
            workflow node. .

        Parameters
        ----------
        geometryJSON: dict
                The geometry object that defines the object.
                Has the structure:

                "source" :    {
                        "type": "indoorBuilding",
                        "fileName": "StationOpenDoors",
                        "fileType" : "obj"
                    }
                "mesh" : {
                    ... snappy hex mesh data
                }

        workingDirectory: str
                The path to copy the obj file to.

        overwrite:
            Not used

        Returns
        -------

        """
        self.logger.info(" -- Start -- ")

        objectName = geometryJSON['source']['fileName']
        objectType = geometryJSON['source']['fileType']
        # Updating the snappyHexMesh object:
        snappyNode = self.snappyHexMesh
        if snappyNode is None:
            self.logger.error("There is no snappyHexMesh Node in this workflow")
            raise ValueError("There is no snappyHexMesh Node in this workflow")

        # # 1. change the object name to the regionName
        # try:
        #     objNode =  snappyNode['geometry']['objects'][regionName]
        # except KeyError:
        #     err = f"The object {regionName} was not found in snappyHexMesh node. Found : {','.join(snappyNode['geometry'].keys())}"
        #     self.logger.error(err)
        #     raise ValueError(err)
        #
        # # 2. change all the keys from the casefile:
        # for keyName,keyValue in geometryJSON['meshing'].items():
        #     objNode[keyName] = keyValue
        self.parameters['objectFile'] = f"{objectName}.{objectType}"
        snappyNode['geometry']['objects'][regionName].update(geometryJSON['meshing'])

    ############## blockMeshHandler
    ## Defines the block mesh node .

    def blockMeshHandler_region(self, configuration,overwrite=False):
        """
            Creates the vertices in the the blockMesh, and updates the workflow.

            The location of the region was taken from the region that was defined in the configuration.
            Currently, we only support box-like blockMesh domains.

            We assume that the region to crop is bbox. (and given in the configuration as dict with the xmax,xmin,ymin,ymax keys).

            Does not change the boundaries or makes the mesh cyclic.

        Parameters
        ----------
        configuration : dict
            The configuration data.


        Returns
        -------
            None
        """
        regionName = configuration['mesh']['blockMesh']['region']
        height = configuration['mesh']['blockMesh']['height']

        corners = configuration['regions'][regionName]['parameters']

        Xlist = [corners['xmin'],corners['xmax'],corners['xmax'],corners['xmin']]
        Ylist = [corners['ymin'], corners['ymin'], corners['ymax'], corners['ymax']]

        VerticesList = []
        for h,xy in product(height,zip(Xlist,Ylist)):
            VerticesList.append([xy[0],xy[1],h])


        blockMesh = self.blockMesh
        blockMesh['vertices'] = VerticesList

    def blockMeshHandler_objFile(self, configuration,overwrite=False):
        """
            Creates the vertices in the the blockMesh, and updates the workflow using the
            bounds of the obj file.
            The location of the region was taken from the region that was defined in the configuration.
            Currently, we only support box-like blockMesh domains.

            Does not change the boundaries or makes the mesh cyclic.

            Note that the obj file is given in millimeters and it is converted to meters.

        Parameters
        ----------
        configuration : dict
            The configuration data.


        Returns
        -------
            None, updates the blockMesh node.
        """
        try:
            from freecad import app as FreeCAD
            import Mesh
        except ImportError:
            self.logger.error("freecad module  is not installed in the environment")
            raise ImportError("freecad is not installed. please install it before trying again.")

        # Load the file
        fileName = configuration['mesh']['blockMesh']['fileName']
        Mesh.open(fileName)
        objFile = FreeCAD.getDocument("Unnamed")

        bboxes = [x.Mesh.BoundBox for x in objFile.findObjects()]

        maxPropList = ['XMax','YMax','ZMax']
        corners = dict()
        for propName in maxPropList:
             corners[propName] = numpy.max([getattr(x,propName) for x in bboxes])/1000

        minPropList = ['XMin', 'YMin', 'ZMin']
        for propName in minPropList:
             corners[propName] = numpy.min([getattr(x,propName) for x in bboxes])/1000

        Xlist = [corners['XMin'],corners['XMax'],corners['XMax'],corners['XMin']]
        Ylist = [corners['YMin'], corners['YMin'], corners['YMax'], corners['YMax']]
        Zlist = [corners['ZMin'], corners['ZMax']]

        VerticesList = []
        for h,xy in product(Zlist,zip(Xlist,Ylist)):
            VerticesList.append([xy[0],xy[1],h])

        blockMesh = self.blockMesh
        blockMesh['vertices'] = VerticesList

        # copy the obj file to the template directory.
        #copyfile( fileName,    os.path.join(configuration["specializedTemplateDirectory"],fileName))


    ############## location in the mesh
    ## Handles the location mesh.

    def locationInMeshHandler_relative(self,configuration,overwrite=False):
        """
            Calculates the location in mesh point according to the fractions
            in each dimension.

        Parameters
        ----------
        configuration : dict
                The configurtation of the change.

        Returns
        -------

        """
        regionName = configuration['mesh']['blockMesh']['region']
        height = configuration['mesh']['blockMesh']['height']

        corners = configuration['regions'][regionName]['parameters']

        x_frac = configuration['mesh']['locationInMesh']['x']
        y_frac = configuration['mesh']['locationInMesh']['y']
        z_frac = configuration['mesh']['locationInMesh']['z']


        x = corners['xmin'] + x_frac *(corners['xmax']-corners['xmin'])
        y = corners['ymin'] + y_frac *(corners['ymax']-corners['ymin'])
        z = height[0] + z_frac * (height[1] -height[0])

        snappyNode = self.snappyHexMesh

        if snappyNode is not None:
            snappyNode["castellatedMeshControls"]["locationInMesh"] = [x,y,z]


    def locationInMeshHandler_absolute(self,configuration,overwrite=False):
        """
            Calculates the location in mesh point according to the fractions
            in each dimension.



        Parameters
        ----------
        configuration : dict
                The configurtation of the change.

        Returns
        -------

        """
        x = configuration['mesh']['locationInMesh']['x']
        y = configuration['mesh']['locationInMesh']['y']
        z = configuration['mesh']['locationInMesh']['z']

        snappyNode = self.snappyHexMesh

        if snappyNode is not None:
            snappyNode["castellatedMeshControls"]["locationInMesh"] = [x,y,z]


    ############## initial and boundry conditions.
    ## creates the initial and the boundary conditions

    def ICandBCHandler_node(self,icnode):
        """
            write the defineNewBoundaryConditions node.

        Returns
        -------

        """
        self.defineNewBoundaryConditions['fields'] = icnode['data']

class Workflow_Dispersion(abstractWorkflow):

    def __init__(self ,workflowJSON):
        super().__init__(workflowJSON=workflowJSON)

        # examine here that all the nodes exist, if not - it is not a flow

    @staticmethod
    def getFlowFieldName(baseName,flowID):
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

        formatted_number = "{0:04d}".format(flowID+1)
        return f"{baseName}_{formatted_number}"

    @classmethod
    def prepareFlowField(cls,projectName,flowData,suggsetedName):
        """
            Prepares the case directory of the flow for the dispersion.
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
        configurationJSON : dict
            The configuration of the base flow.

            Has the structure:
                    "baseFlow" : {
                        type: "directory|simulationName|workflowFile|ID"
                        parameters:  {
                                .. depend on the type
                        },
                        "useTimestep" : [float], optional, default: last timestep
                        fromTimestep : [float], optional,  default: 0
                        maximalDispersionTime : [float], required.
                        copyMesh : [bool], default false

                    },
                    "dispersionField" : {

                    }

            base flow parameters:
            =====================

                * useTimestep : float, the time step to copy from the simulation.
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
                                    + filters : { nodename : { jpath : value}}
                                            where nodename is the name of the node in the workflow.

                        - workflowFile: query the hera db using the workflow.
                            Parameters:
                                    + filters : { nodename : { jpath : value}}
                                            where nodename is the name of the node in the workflow.
                        - ID : the id of the document in the heradb.

                        ONLY directory is implemented now.

        Returns
        -------
            A list of str

            The ids of the documents that will be used as the directories.

        """

        logger = loggedObject(loggerName="simulations.openFoam.datalayer.hermesWorkflow").logger
        db = Project(projectName=projectName)


        logger.info("-------- Start ---------")
        baseFlow = flowData['baseFlow']

        # 1. Get the case of the run
        logger.debug(f"Getting the base flow directory from handler  {baseFlow['type']}")
        handler = getattr(cls,f"baseTemplateHandler_{baseFlow['type']}")
        orig = handler(baseFlow['parameters'])
        logger.debug(f"Getting the base flow directory: {orig}")

        TS = [os.path.basename(ts) for ts in glob.glob(os.path.join(orig, "processor0", "*")) if
              os.path.basename(ts).replace(".", "").isdigit()]
        TS.sort()

        useTimestep = baseFlow.get("useTimestep",None)
        fromTimestep = baseFlow.get("fromTimestep","0")
        maximalDispersionTime = baseFlow["maximalDispersionTime"]
        copyMesh = baseFlow.get("copyMesh",False)

        if useTimestep is None:
            # find maximal TS, assume it is parallel:
            uts = TS[-1]
        else:
            # find the closes TS.
            request = float(useTimestep)
            uts = TS[min(range(len(TS)), key=lambda i: abs(float(TS[i]) - request))]

        fromTime = fromTimestep

        logger.debug(f"Using Time step {uts} for Steady state")

        ## Now we should find out if there is a similar run.
        # 1. same original run.
        # 2. same Hmix/ustar/... other fields.
        # 3. the requested start and end are within the existing start and end.

        querydict = dict(
            flowParameters=dict(
                baseFlowDirectory=orig,
                flowFields=flowData['dispersionFields'],
                usingTimeStep = uts
            )
        )

        logger.debug(f"Test if the requested flow field already exists in the project")
        docList = db.getCacheDocuments(type=cls.DOCTYPE_DISPERSION_FLOWFIELD,**querydict)

        if len(docList) == 0:
            logger.info(f"Flow field not found, creating new and adding to the DB")
            ofhome = ofObjectHome()
            ## Add new document to the db.


            # Copy

            ### Tempoarary
            DispersionCase = suggsetedName # this will be replaced with the document id from the db
            ########
            dest = os.path.abspath(DispersionCase)
            try:
                logger.info(f"Creating the directory {dest}")
                os.makedirs(dest,exist_ok=True)
            except FileExistsError:
                raise FileExistsError("The case already exists, use --overwrite ")

            # copy constant, 0 and system.
            for general in ["constant", "system", "0"]:
                logger.info(f"Copying {general} in {orig} directory --> {dest}")
                orig_general = os.path.join(orig, general)
                dest_general = os.path.join(dest, general)
                if os.path.exists(dest_general):
                    copy_tree(orig_general, dest_general)
                else:
                    shutil.copytree(orig_general, dest_general)

            logger.debug(f"Copy the parallel directories")
            for proc in glob.glob(os.path.join(orig,"processor*")):
                for dest_time in [fromTime, maximalDispersionTime]:
                    logger.info(f"Handling time {dest_time}/Process {proc}")
                    orig_proc = os.path.join(proc,str(uts))
                    dest_proc = os.path.join(dest,os.path.basename(proc),str(dest_time))
                    shutil.copytree(orig_proc,dest_proc)

                logger.debug(f"Copying the mesh")
                if copyMesh:
                    orig_constant = os.path.join(proc,"constant")
                    dest_constant = os.path.join(dest,os.path.basename(proc),"constant")
                    shutil.copytree(orig_constant, dest_constant)
                else:
                    fullpath = os.path.abspath(os.path.join(proc, "constant", "polyMesh"))
                    destination = os.path.join(dest, os.path.basename(proc), "constant", "polyMesh")
                    os.makedirs(os.path.dirname(destination), exist_ok=True)
                    os.system(f"ln -s {fullpath} {destination}")

            logger.debug(f"Creating the flow specific fields in the flow needed for the dispersion")

            dispersionFields = flowData['dispersionFields']
            for dispersionFieldName in dispersionFields.keys():
                fieldDimensions = dispersionFields[dispersionFieldName].get("dimensions",None)
                fieldComponents = dispersionFields[dispersionFieldName].get("components", None)

                logger.debug(f"Creating the flow specific field: {dispersionFieldName}. ")
                field = ofhome.getField(fieldName=dispersionFieldName,
                                              fieldGroup=ofObjectHome.GROUP_DISPERSION,
                                              dimensions=fieldDimensions,
                                              componentNames=fieldComponents)

                for proc in glob.glob(os.path.join(orig, "processor*")):
                    procName = os.path.split(proc)[-1]
                    logger.debug(f"Writing the flow specific field {dispersionFieldName} to processor {procName} . ")
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
            ret = dest
        else:
            ret = docList[0].resource

        return ret

