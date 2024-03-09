import os
import json
from ... import toolkitHome
from ...utils import loadJSON
from ...utils.logging import helpers as hera_logging
from ...utils.freeCAD import getObjFileBoundaries
from itertools import product

try:
    import hermes
except ImportError:
    raise ImportError("Cannot use this module without hermes... Install it. ")


class abstractWorkflow(hermes.workflow):
    """
            An abstract specialization of the hermes workflow to the
            openfoam workflow.

            Add 2 capabilities to the hermes.workflow:

            1.  Access to the group and the project names.
            1.  Adds the specific functionalities for an OpenfFOAM workflow (determining initial conditions, and ect).
    """

    workflowDescirption = None

    _requiredNodeList = None

    def __init__(self, workflowJSON, workflowHeraDocument=None,name=None):
        """
            Initializes the abstract workflow.

        Parameters
        ----------
        workflowJSON : json
                The JSON that describes the workflow.

        workflowDoc : hera document.
                Holds the data of the database (optional).
        """
        super().__init__(workflowJSON=workflowJSON,name=name)
        self.logger = hera_logging.get_logger(self)  # was: loggedObject(name=None).logger, ignores actual class
        self.workflowHeraDocument = workflowHeraDocument

        self._requiredNodeList = ['controlDict', 'fvSolution', 'fvSchemes', 'fileWriter', 'defineNewBoundaryConditions','Parameters']

    @property
    def workflowGroup(self):
        ret = None
        if self.workflowHeraDocument:
            return self.workflowHeraDocument['desc']['groupName']

    @workflowGroup.setter
    def workflowGroup(self, value):
        ret = None
        if not isinstance(value,str):
            raise ValueError("Group name must be a string")
        if self.workflowHeraDocument:
            self.workflowHeraDocument['desc']['groupName'] = value
            self.workflowHeraDocument.save()
        return ret

    @property
    def workflowGroupID(self):
        ret = None
        if self.workflowHeraDocument:
            ret =  self.workflowHeraDocument['desc']['groupID']
        return ret

    @workflowGroup.setter
    def workflowGroupID(self, value):
        ret = None
        if not isinstance(value,str):
            raise ValueError("Group name must be a string")
        if self.workflowHeraDocument:
            self.workflowHeraDocument['desc']['groupID'] = value
            self.workflowHeraDocument.save()
        return ret

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
    def fileWriter(self):
        return self['fileWriter']

    @property
    def parameters(self):
        return self['parameters']

    @property
    def defineNewBoundaryConditions(self):
        return self['defineNewBoundaryConditions']

##############################################################################
##                          Workflow Eulerian/Lagrangian
##############################################################################

#### Eulerian
class workflow_Eulerian(abstractWorkflow):
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
    def __init__(self,workflowJSON,workflowHeraDocument=None,name=None):
        """
            Initializes a hermes workflow for openflow flow simulations.


        Parameters
        ----------
        workflow : str, dict, file type.
            Can be the file name, a string of the JSON, or a file type.


        parallelNodes: str, list
            A name of node, or a list of nodes in the workflow that are required to build the case in parallel.
            Will be removed if the workflow is executed as a unified case.
        """
        super().__init__(workflowJSON=workflowJSON, workflowHeraDocument=workflowHeraDocument,name=name)
        self._requiredNodeList.append("blockMesh")

        # examine here that all the nodes exist, if not - it is not a flow
        for node in self._requiredNodeList:
            if node not in self.workflowJSON['nodes']:
                raise ValueError(f"The node {node} does not exist in the flow. Not a flow workflow.")

    @property
    def blockMesh(self):
        return self['blockMesh']

    @property
    def snappyHexMesh(self):
        return self['snappyHexMesh'] if 'snappyHexMesh' in self.nodes else None

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
    ############## initial and boundry conditions.
    ## creates the initial and the boundary conditions

    def ICandBCHandler_node(self,icnode):
        """
            write the defineNewBoundaryConditions node.

        Returns
        -------

        """
        self.defineNewBoundaryConditions['fields'] = icnode['data']


#### Lagrangian

class workflow_Lagrangian(abstractWorkflow):

    def __init__(self ,workflowJSON,workflowHeraDocument=None,name=None):
        super().__init__(workflowJSON=workflowJSON, workflowHeraDocument=workflowHeraDocument,name=name)


##############################################################################
##                          Workflow_StochasticLagrangianSolver
##############################################################################
class workflow_StochasticLagrangianSolver(workflow_Lagrangian):

    def __init__(self ,workflowJSON,workflowHeraDocument=None,name=None):
        super().__init__(workflowJSON=workflowJSON, workflowHeraDocument=workflowHeraDocument,name=name)

        # Make sure that the
        # dispersionFlowField exists
        if 'dispersionFlowField' not in self.parameters.parameters:
            err = "The StochasticLagrangianSolver must have a dispersionFlowField specification in the parameters node"
            self.logger.error(err)
            raise ValueError(err)

    @property
    def dispersionFlowField(self):
        return

    @property
    def dispersionFlowField(self):
        return self.parameters.parameters['dispersionFlowField']

    @dispersionFlowField.setter
    def dispersionFlowField(self, value):
        self.parameters.parameters['dispersionFlowField'] = str(value)


##########################################################
##                          Workflow_simpleFoam
##########################################################
class workflow_simpleFoam(workflow_Eulerian):

    def __init__(self ,workflowJSON,workflowHeraDocument=None,name=None):
        super().__init__(workflowJSON=workflowJSON, workflowHeraDocument=workflowHeraDocument,name=name)

##########################################################
##                          Workflow_indoorFOAMBoussinesq
##########################################################
class workflow_indoorFOAMBoussinesq(workflow_Eulerian):

    def __init__(self ,workflowJSON,workflowHeraDocument=None,name=None):
        super().__init__(workflowJSON=workflowJSON, workflowHeraDocument=workflowHeraDocument,name=name)





#
# def prepareGeometry(self, workingDirectory, configurationFile,overwrite=False):
#     """
#         Adapts the geometry of the workflow.
#
#     Parameters
#     ----------
#     caseDirectory : str
#         The directory to save the STL in
#
#     configurationFile: str,file like, dict
#         Reads the JSON  from the directory.
#         Could be a file like object, a path to the file, a JSON string or a dict.
#
#     overwrite: bool
#         If true, then generate the STL even it exists.
#         if False, skip the genereation.
#
#     Returns
#     -------
#         dict,
#         The name of the object and the names that were written.
#
#     """
#     self.logger.info("-- Start --")
#     configuration = loadJSON(configurationFile)
#
#     geom_handlers = [x.split("_")[1] for x in dir(self) if x.startswith("geometryHandler")]
#
#     for name,geometryObject in configuration['geometry'].items():
#         self.logger.debug(f"Processing: \n {json.dumps(geometryObject,indent=4)}")
#
#         # 1. Handle the gemetry-type specifics. (create STL and ect.)
#         if geometryObject['source']['type'] not in geom_handlers:
#             err = f"{geometryObject['source']['type']} must be one of {','.join(geom_handlers)}"
#             self.logger.error(err)
#             raise ValueError(err)
#
#         handler = getattr(self, f"geometryHandler_{geometryObject['source']['type']}")
#         handler(name, geometryObject, workingDirectory, configuration,overwrite)
#
#         # 2. Change the snappyHexMesh node.
#         # snappyHexNode = self.snappyHexMesh
#         #
#         # if snappyHexNode is None:
#         #     self.logger.warning("snappyHexMesh node does not exist. Geometry might not be used!")
#         # else:
#         #     snappyHexNode["geometry"]['objects'][name] = geometryObject["meshing"]
#
#     self.logger.info("-- End --")
#
# def prepareMesh(self, configurationFile):
#     """
#         Reads the configuration file and:
#
#         Adapt the following workflow nodes:
#
#             1. blockMesh
#             2. Snappy
#                     Setup the locatin In Mesh property in the snappy hex mesh (if exists).
#
#     Parameters
#     ----------
#     configurationFile
#
#     Returns
#     -------
#
#     """
#     self.logger.info("-- Start --")
#     configuration = loadJSON(configurationFile)
#     blockMeshNode = configuration['mesh'].get('blockMesh',None)
#
#     if blockMeshNode is not None:
#         self.logger.execution("Found blockMesh Node")
#
#         ## 1. Block mesh handler.
#         blockMesh_handlers = [x.split("_")[1] for x in dir(self) if x.startswith("blockMeshHandler")]
#         if blockMeshNode['type'] not in blockMesh_handlers:
#             err = f"{blockMeshNode['type']} must be one of {','.join(blockMesh_handlers)}"
#             self.logger.error(err)
#             raise ValueError(err)
#
#         handler = getattr(self, f"blockMeshHandler_{blockMeshNode['type']}")
#         self.logger.debug(f"Calling the blockMeshHandler_{blockMeshNode['type']} handler")
#         handler(configuration)
#
#         locationNode =  configuration['mesh'].get('locationInMesh',None)
#
#         if locationNode is not None:
#             self.logger.execution("Found location Node")
#             ## 2. Location in mesh handler.
#             locationInMesh_handlers = [x.split("_")[1] for x in dir(self) if x.startswith("locationInMeshHandler")]
#             if locationNode['type'] not in locationInMesh_handlers:
#                 err = f"{blockMeshNode['locationInMesh']} must be one of {','.join(locationInMesh_handlers)}"
#                 self.logger.error(err)
#                 raise ValueError(err)
#
#             handler = getattr(self, f"locationInMeshHandler_{locationNode['type']}")
#             self.logger.execution(f"Calling the locationInMeshHandler_{locationNode['type']} handler")
#             handler(configuration)
#     else:
#         self.logger.execution("blockMesh node was not Found. ")
#
#     self.logger.info("-- End -- ")
#
# ##########################################
# ##
# ##                  Handlers
# ##
#     ##################
#
#
# ############## geometryHandler
# # Defines the geometry objects and sets up the snappy hex mesh.
# def geometryHandler_topography(self,regionName, geometryJSON, workingDirectory, configuration,overwrite=False):
#     """
#         Build the topography STL using the topography toolkit.
#
#
#     Parameters
#     ----------
#     geometryJSON: dict
#             The geometry object.
#             Has the structure:
#
#             "source" :    {
#                     "type" : "topography",  <-- the type
#                     "datasource" : "BNTL",  <-- The name of the datasource to use.
#                     "region" : "canopy"     <-- The name of the region in the configuration file to use.
#                     "dxdy" : 50             <-- The resolution of the conversion.
#                 }
#
#     workingDirectory: str
#             The path to save the STL to.
#
#     regions: dict
#             A dict with the regions as polygons.
#
#     Returns
#     -------
#
#     """
#     self.logger.info(" -- Start -- ")
#     tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_TOPOGRAPHY,
#                                 projectName=configuration['projectName'])
#
#     stlFileName = f"{geometryJSON['meshing']['objectName']}.stl"
#     fullFileName = os.path.join(workingDirectory, stlFileName)
#
#     region = geometryJSON['source']['region']
#     regionList = configuration.get('regions', {})
#     regionCoords = regionList.get(region, None)
#
#     if regionCoords is None:
#         raise ValueError(f"The region: {region} is not found. Found the regions: {','.join(regionList.keys())}")
#
#     bx = [regionCoords['parameters']['xmin'], regionCoords['parameters']['ymin'], regionCoords['parameters']['xmax'], regionCoords['parameters']['ymax']]
#     self.logger.debug(f"Found the region {region} with coords {bx}")
#
#     if os.path.exists(fullFileName) and not overwrite:
#         self.logger.info(f"File {stlFileName} in {regionName} exists. and overwrite=False")
#     else:
#         self.logger.info(f"Generating {stlFileName} for geometry {regionName}")
#         stlcontent = tk.regionToSTL(shapeDataOrName=bx, dxdy=geometryJSON['source']['dxdy'],
#                                     datasourceName=geometryJSON['source']['datasourceName'])
#
#         self.logger.debug(f"Converting to STL complete. Now writing to the file {fullFileName}")
#         with open(fullFileName, 'w') as stloutputfile:
#             stloutputfile.write(stlcontent)
#
#     ## Here we should update the snappyHexMesh node
#
#     self.logger.debug(f"Adding the region {regionName} to the workflow")
#     self.parameters["domains"] = {regionName : dict(shapeDataOrName=bx,datasourceName=geometryJSON['source']['datasourceName'])}
#     self.logger.info(" -- End -- ")
#
# def geometryHandler_buildings(self,regionName, geometryJSON, workingDirectory, configuration):
#     """
#         Creates the STL of the buildings in the region using the building toolkit.
#
#         Sets the parameter buildingsBBOX to the bounding box of the buildings area.
#
#         We assume that the region to crop is bbox. (and given in the configuration as dict with the xmax,xmin,ymin,ymax keys).
#
#
#     Parameters
#     ----------
#     regionName: str
#         The name of the region.
#
#         Not used here, only for completeness of the interface.
#
#     geometryJSON: dict
#             The geometry object is a dict with the following keys:
#
#             * source - parameters for building the STL of the buildings.
#                        used by the buildings toolkit.
#
#                 has the strucutre:
#                     "source" :    {
#                             "type" : "buildings",  <-- the type
#                             "datasource" : "BNTL",  <-- The name of the datasource to use.
#                             "region" : "canopy"     <-- The name of the region in the configuration file to use.
#                         }
#             * meshing - parameters used by the blockMesh and the snappyHexMesh to build the code.
#                 has the structure of the snappyHexMesh node ["geometry"]['objects'][name]
#
#     configuration: dict
#             The overall configuration of the case.
#
#     Returns
#     -------
#
#     """
#     tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_BUILDINGS, projectName=configuration['projectName'])
#
#     region = geometryJSON['source']['region']
#     regionList = configuration.get('regions', {})
#     regionCoords = regionList.get(region, None)
#
#     if regionCoords is None:
#         raise ValueError(f"The region{region} is not found. Found the regions: {','.join(regionList.keys())}")
#
#     bx = [regionCoords['parameters']['xmin']*0.95, regionCoords['parameters']['ymin']*0.95, regionCoords['parameters']['xmax']*1.05, regionCoords['parameters']['ymax']*1.05]
#
#     stlFileName = f"{geometryJSON['meshing']['objectName']}.stl"
#     fullFileName = os.path.join(workingDirectory, stlFileName)
#
#     tk.regionToSTL(regionNameOrData=bx, outputFileName=fullFileName,flat=None,saveMode=toolkitHome.TOOLKIT_SAVEMODE_NOSAVE)
#
#     self.parameters['buildingsBBOX'] = bx
#     ## Here we should update the snappyHexMesh node
#
# def geometryHandler_indoorBuilding(self,regionName, geometryJSON, workingDirectory, configuration,overwrite=None):
#     """
#         Updating the snappyHexMesh node.
#
#         Assumes that the name of the region is the name of the object in the
#         workflow node. .
#
#     Parameters
#     ----------
#     geometryJSON: dict
#             The geometry object that defines the object.
#             Has the structure:
#
#             "source" :    {
#                     "type": "indoorBuilding",
#                     "fileName": "StationOpenDoors",
#                     "fileType" : "obj"
#                 }
#             "mesh" : {
#                 ... snappy hex mesh data
#             }
#
#     workingDirectory: str
#             The path to copy the obj file to.
#
#     overwrite:
#         Not used
#
#     Returns
#     -------
#
#     """
#     self.logger.info(" -- Start -- ")
#     objectName = geometryJSON['source']['fileName']
#     objectType = geometryJSON['source']['fileType']
#     # Updating the snappyHexMesh object:
#     snappyNode = self.snappyHexMesh
#     if snappyNode is None:
#         self.logger.error("There is no snappyHexMesh Node in this workflow")
#         raise ValueError("There is no snappyHexMesh Node in this workflow")
#
#     # # 1. change the object name to the regionName
#     # try:
#     #     objNode =  snappyNode['geometry']['objects'][regionName]
#     # except KeyError:
#     #     err = f"The object {regionName} was not found in snappyHexMesh node. Found : {','.join(snappyNode['geometry'].keys())}"
#     #     self.logger.error(err)
#     #     raise ValueError(err)
#     #
#     # # 2. change all the keys from the casefile:
#     # for keyName,keyValue in geometryJSON['meshing'].items():
#     #     objNode[keyName] = keyValue
#     self.parameters['objectFile'] = f"{objectName}.{objectType}"
#     snappyNode['geometry']['objects'][regionName].update(geometryJSON['meshing'])
#
# ############## blockMeshHandler
# ## Defines the block mesh node .
#
# def blockMeshHandler_region(self, configuration):
#     """
#         Creates the vertices in the blockMesh, and updates the workflow.
#
#         The location of the region was taken from the region that was defined in the configuration.
#         Currently, we only support box-like blockMesh domains.
#
#         We assume that the region to crop is bbox. (and given in the configuration as dict with the xmax,xmin,ymin,ymax keys).
#
#         Does not change the boundaries or makes the mesh cyclic.
#
#     Parameters
#     ----------
#     configuration : dict
#         The configuration data.
#
#     Returns
#     -------
#         None
#     """
#     regionName = configuration['mesh']['blockMesh']['region']
#     height = configuration['mesh']['blockMesh']['height']
#
#     corners = configuration['regions'][regionName]['parameters']
#
#     Xlist = [corners['xmin'],corners['xmax'],corners['xmax'],corners['xmin']]
#     Ylist = [corners['ymin'], corners['ymin'], corners['ymax'], corners['ymax']]
#
#     VerticesList = []
#     for h,xy in product(height,zip(Xlist,Ylist)):
#         VerticesList.append([xy[0],xy[1],h])
#
#
#     blockMesh = self.blockMesh
#     blockMesh['vertices'] = VerticesList
#
# def blockMeshHandler_objFile(self, configuration,overwrite=False):
#     """
#         Creates the vertices in the the blockMesh, and updates the workflow using the
#         bounds of the obj file. However, if the 'domainbounds' key exists, then replace the
#         domain value by the input.
#
#         The location of the region was taken from the region that was defined in the configuration.
#         Currently, we only support box-like blockMesh domains.
#
#         Does not change the boundaries or makes the mesh cyclic.
#
#         Note that the obj file is assumed to be given in millimeters and it is converted to meters.
#         Therefore, we divide by 1000.
#
#
#
#     Parameters
#     ----------
#     configuration : dict
#         The configuration data.
#
#         A dict with the strucutre:
#         'mesh': {
#             <other nodes>
#             'blockmesh' : {
#                 "type": "objFile",
#                 "fileName": "buildings.obj",
#                 "absoluteDomainBounds" : {    <- optional
#                         [coordinate ] : 10000
#                 }
#             }
#         }
#             where [coordinate] is 'XMin','XMax','YMin','YMax','ZMin','ZMax'
#
#
#     Returns
#     -------
#         None, updates the blockMesh node.
#     """
#     fileName = configuration['mesh']['blockMesh']['fileName']
#     corners = getObjFileBoundaries(fileName)
#
#     # check if the user supplied corners explicitly
#     determinedCorners =  configuration['mesh']['blockMesh'].get("absoluteDomainBounds",{})
#     corners.update(determinedCorners)
#
#     Xlist = [corners['XMin'],corners['XMax'],corners['XMax'],corners['XMin']]
#     Ylist = [corners['YMin'], corners['YMin'], corners['YMax'], corners['YMax']]
#     Zlist = [corners['ZMin'], corners['ZMax']]
#
#     VerticesList = []
#     for h,xy in product(Zlist,zip(Xlist,Ylist)):
#         VerticesList.append([xy[0],xy[1],h])
#
#     blockMesh = self.blockMesh
#     blockMesh['vertices'] = VerticesList
#
#     # copy the obj file to the template directory.
#     #copyfile( fileName,    os.path.join(configuration["specializedTemplateDirectory"],fileName))
#
#
# ############## location in the mesh
# ## Handles the location mesh.
#
# def locationInMeshHandler_relative(self,configuration,overwrite=False):
#     """
#         Calculates the location in mesh point according to the fractions
#         in each dimension.
#
#     Parameters
#     ----------
#     configuration : dict
#             The configurtation of the change.
#
#     Returns
#     -------
#
#     """
#     meshType = configuration['mesh']['blockMesh']['type']
#
#     func = getattr(self,f"locationInMeshHandlerRelative_{meshType}")
#     func(configuration)
#
# def locationInMeshHandlerRelative_objFile(self,configuration):
#     """
#         Finds the relative point in the objFile, and update the snappyHexMesh node.
#     Returns
#     -------
#
#     """
#     self.logger.info("Calculating relative location for objFile")
#
#     fileName = configuration['mesh']['blockMesh']['fileName']
#     corners = getObjFileBoundaries(fileName)
#     self.logger.debug("Got the corners {corners}")
#     x_frac = configuration['mesh']['locationInMesh']['x']
#     y_frac = configuration['mesh']['locationInMesh']['y']
#     z_frac = configuration['mesh']['locationInMesh']['z']
#
#     self.logger.debug(f"Got the fractions x: {x_frac}, y: {y_frac}, z: {z_frac}")
#
#     x = corners['XMin'] + x_frac *(corners['XMax']-corners['XMin'])
#     y = corners['YMin'] + y_frac *(corners['YMax']-corners['YMin'])
#     z = corners['ZMin'] + z_frac * (corners['ZMax'] - corners['ZMin'])
#
#     self.logger.debug(f"Got the ceners x: {x}, y: {y}, z: {z}")
#
#     snappyNode = self.snappyHexMesh
#
#     if snappyNode is not None:
#         snappyNode["castellatedMeshControls"]["locationInMesh"] = [x,y,z]
#
# def locationInMeshHandlerRelative_region(self,configuration):
#     """
#         Finds the relative point in the region, and update the snappyHexMesh node.
#     Returns
#     -------
#
#     """
#
#     regionName = configuration['mesh']['blockMesh']['region']
#     height = configuration['mesh']['blockMesh']['height']
#
#     corners = configuration['regions'][regionName]['parameters']
#     x_frac = configuration['mesh']['locationInMesh']['x']
#     y_frac = configuration['mesh']['locationInMesh']['y']
#     z_frac = configuration['mesh']['locationInMesh']['z']
#
#     x = corners['XMin'] + x_frac *(corners['XMax']-corners['XMin'])
#     y = corners['YMin'] + y_frac *(corners['YMax']-corners['YMin'])
#     z = height[0] + z_frac * (height[1] -height[0])
#
#     snappyNode = self.snappyHexMesh
#
#     if snappyNode is not None:
#         snappyNode["castellatedMeshControls"]["locationInMesh"] = [x,y,z]
#
# def locationInMeshHandler_absolute(self,configuration,overwrite=False):
#     """
#         Calculates the location in mesh point according to the fractions
#         in each dimension.
#
#
#
#     Parameters
#     ----------
#     configuration : dict
#             The configurtation of the change.
#
#     Returns
#     -------
#
#     """
#     x = configuration['mesh']['locationInMesh']['x']
#     y = configuration['mesh']['locationInMesh']['y']
#     z = configuration['mesh']['locationInMesh']['z']
#
#     snappyNode = self.snappyHexMesh
#
#     if snappyNode is not None:
#         snappyNode["castellatedMeshControls"]["locationInMesh"] = [x,y,z]
#
