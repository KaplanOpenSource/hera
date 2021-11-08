from ..hermesWorkflow import hermesWorkflow,hermesNode
import numpy
import os
from ... import toolkitHome
from ...utils.jsonutils import loadJSON
from itertools import product
import json


class hermesOpenFOAMWorkflow(hermesWorkflow):
    """
        This class manages the hermes workflow of openFOAM.

        It contains additional functions that are specific to openFOAM workflows
        such as controlDict, fvOptions, fvSchemes, and fvSolutions that return the relevant
        nodes.

        Also, it contains a function change the from composed to decomposed run.

        Each node has a list of default flags that could be modified when the case execution script is built.

        The case is considered as running in parallel if there is an execution that is required to run
        in parallel (for example the decompose par).

    """

    _workflowExecutionJSON = None # Holds the datastructure that allows the workflow to build the execution
                              # of the workflow.

    def __init__(self,workflow,workflowExecutionJSON=None,parallelNodes=None):
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

        """
        super().__init__(workflow)

        self._workflowExecutionJSON = workflowExecutionJSON
        self._parallelNodes = [] if parallelNodes is None else numpy.atleast_1d(parallelNodes)



    @property
    def parameters(self):
        return hermesNode('controlDict',self.nodes['Parameters'])

    @property
    def controlDict(self):
        return hermesNode('controlDict',self.nodes['controlDict'])

    @property
    def fvSolution(self):
        return hermesNode('fvSolution',self.nodes['fvSolution'])

    @property
    def fvScheme(self):
        return hermesNode('fvScheme',self.nodes['fvScheme'])

    @property
    def blockMesh(self):
        return hermesNode('blockMesh',self.nodes['blockMesh'])


    @property
    def snappyHexmesh(self):
        return hermesNode('snappyHexMesh',self.nodes['snappyHexMesh']) if 'snappyHexMesh' in self.nodes else None


    @property
    def fileWriter(self):
        return hermesNode('fileWriter',self.nodes['fileWriter'])


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
            execLine += f"foamJob {parallelFlag} -screen -wait {progName} {params}\n"


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
            snappyHexNode = self.snappyHexmesh

            if snappyHexNode is None:
                self.logger.warning("snappyHexMesh node does not exist. Geometry might not be used!")
            else:
                snappyHexNode["geometry"]['objects'][name] = geometryObject["meshing"]

        self.logger.info("-- End --")

    def prepareMesh(self, configurationFile, workingDirectory):
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
                if blockMeshNode['locationInMesh'] not in locationInMesh_handlers:
                    err = f"{blockMeshNode['locationInMesh']} must be one of {','.join(locationInMesh_handlers)}"
                    self.logger.error(err)
                    raise ValueError(err)

                handler = getattr(self, f"locationInMeshHandler_{locationNode['type']}")
                self.logger.execution(f"Calling the locationInMeshHandler_{locationNode['type']} handler")
                handler(configuration)
        else:
            self.logger.execution("blockMesh node was not Found. ")

        self.logger.info("-- End -- ")


    ##########################################
    ##
    ##                  Handlers
    ##
    ##################
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
            self.logger.info(f"Geometry {regionName} exists. and overwrite=False")
        else:
            self.logger.info(f"Generating geometry {regionName}")
            stlcontent = tk.regionToSTL(shapeDataOrName=bx, dxdy=geometryJSON['source']['dxdy'],
                                        dataSourceName=geometryJSON['source']['datasource'])

            self.logger.debug(f"Converting to STL complete. Now writing to the file {fullFileName}")
            with open(fullFileName, 'w') as stloutputfile:
                stloutputfile.write(stlcontent)

        self.parameters["domains"] = dict(regionName=dict(shapeDataOrName=bx,dataSourceName=geometryJSON['source']['datasource']))
        self.logger.info(" -- End -- ")

    def geometryHandler_buildings(self,regionName, geometryJSON, workingDirectory, configuration,overwrite=False):
        """
            Build the buildings STL using the buildings toolkit.

            We assume that the region to crop is bbox. (and given in the configuration as dict with the xmax,xmin,ymin,ymax keys).

        Parameters
        ----------
        geometryJSON: dict
                The geometry object.
                Has the structure:

                "source" :    {
                        "type" : "buildings",  <-- the type
                        "datasource" : "BNTL",  <-- The name of the datasource to use.
                        "region" : "canopy"     <-- The name of the region in the configuration file to use.
                    }

        workingDirectory: str
                The path to save the STL to.

        regions: dict
                A dict with the regions as polygons.

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
        Ylist = [corners['ymin'], corners['ymin'], corners['ymax'], corners['ymin']]

        VerticesList = []
        for xy,h in product(zip(Xlist,Ylist),height):
            VerticesList.append([xy[0],xy[1],h])


        blockMesh = self.blockMesh
        blockMesh['vertices'] = VerticesList

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

        snappyNode = self.snappyHexmesh

        if snappyNode is not None:
            snappyNode["castellatedMeshControls"]["locationInMesh"] = [x,y,z]






