import os
import json
import numpy

from ... import toolkitHome
from ...utils import loadJSON
from ...utils.logging import get_classMethod_logger
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

    def __init__(self, workflowJSON, workflowHeraDocument=None,name=None, **kwargs):
        """
            Initializes the abstract workflow.

        Parameters
        ----------
        workflowJSON : json
                The JSON that describes the workflow.

        workflowDoc : hera document.
                Holds the data of the database (optional).
        """
        super().__init__(workflowJSON=workflowJSON,name=name, **kwargs)
        self.workflowHeraDocument = workflowHeraDocument

        self._requiredNodeList = ['controlDict', 'fvSolution', 'fvSchemes', 'fileWriter','Parameters']

    @property
    def workflowGroup(self):
        """Return the workflow group name from the hera document."""
        ret = None
        if self.workflowHeraDocument:
            return self.workflowHeraDocument['desc']['groupName']

    @workflowGroup.setter
    def workflowGroup(self, value):
        """Set the workflow group name and save to the hera document."""
        ret = None
        if not isinstance(value,str):
            raise ValueError("Group name must be a string")
        if self.workflowHeraDocument:
            self.workflowHeraDocument['desc']['groupName'] = value
            self.workflowHeraDocument.save()
        return ret

    @property
    def workflowGroupID(self):
        """Return the workflow group ID from the hera document."""
        ret = None
        if self.workflowHeraDocument:
            ret =  self.workflowHeraDocument['desc']['groupID']
        return ret

    @workflowGroupID.setter
    def workflowGroupID(self, value):
        """Set the workflow group ID and save to the hera document."""
        ret = None
        if not isinstance(value,str):
            raise ValueError("Group name must be a string")
        if self.workflowHeraDocument:
            self.workflowHeraDocument['desc']['groupID'] = value
            self.workflowHeraDocument.save()
        return ret

    @property
    def controlDict(self):
        """Return the ControlDict workflow node."""
        return self['ControlDict']

    @property
    def fvSolution(self):
        """Return the fvSolution workflow node."""
        return self['fvSolution']


    @property
    def fvScheme(self):
        """Return the fvScheme workflow node."""
        return self['fvScheme']

    @property
    def fileWriter(self):
        """Return the fileWriter workflow node."""
        return self['fileWriter']

    @property
    def parameters(self):
        """Return the Parameters workflow node."""
        return self['Parameters']

    @property
    def buildAllRun(self):
        """Return the buildAllRun workflow node."""
        return self['buildAllRun']

    @property
    def defineNewBoundaryConditions(self):
        """Return the defineNewBoundaryConditions workflow node."""
        return self['defineNewBoundaryConditions']


    def __setitem__(self, key, value):
        """
            Adds the node to the workflow.

            This procedure also updates the nodeList and the buildAllRun node.
            It will add it after the 'nodeListPosition' object.
            The value is a dict with the fields:

            - buildAllRun : dict
                An object that specifies the all run entry:

                - "name": The command that activates this node.
                - "parameters": additional parameters
                - "couldRunInParallel": bool , Can this command run in parallel.
                - "foamJob": bool , is it a foamJob commandor just a regular script.

            - fileWriter


            - nodeListPosition : string or int, [optional]
                    If exists, write the new node after the node name (if string) or its position (if int).
                    if does not exist, add at the end.
            - node             :  dict
                    The data of the node
        Parameters
        ----------
        key
        value

        Returns
        -------

        """
        self.buildAllRun[key] = value['buildAllRun']
        self.fileWriter[key] = value['fileWriter']

        super().__setitem(key=key,value=value)

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
    def __init__(self,workflowJSON,workflowHeraDocument=None,name=None, **kwargs):
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
        super().__init__(workflowJSON=workflowJSON, workflowHeraDocument=workflowHeraDocument,name=name, **kwargs)
        #self._requiredNodeList.append("blockMesh")

        # examine here that all the nodes exist, if not - it is not a flow
        for node in self._requiredNodeList:
            if node not in self.workflowJSON['nodes']:
                raise ValueError(f"The node {node} does not exist in the flow. Not a flow workflow.")


    @property
    def blockMesh(self):
        """Return the blockMesh workflow node."""
        return self['blockMesh']

    @property
    def snappyHexMesh(self):
        """Return the snappyHexMesh node, or None if it does not exist."""
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

    def set_blockMesh_blockBoundaries(self,minx,maxx,miny,maxy,minz,maxz,dx,dy,dz):
        """
            Sets the blockmesh dictionary
        Parameters
        ----------
        minx : float

        maxx : float
        miny : float
        maxy : float
        minz : float
        maxz : float
        dx : float
            The x spacing in m.
        dy : float
            The y spacing in m.
        dz: float
            The z = spacing in m
        Returns
        -------
            None
        """
        logger = get_classMethod_logger(self,"set_blockMesh_boundaries")
        logger.info(f"Setting the blockMesh with block that spans from  [{minx},{miny},{minz}] to [{maxx},{maxy},{maxz}]")
        # Define the 4 XY corner points of the hexahedron base in CCW order
        # (OpenFOAM blockMesh convention: vertices 0-3 bottom, 4-7 top).
        Xlist = [minx,maxx,maxx,minx]
        Ylist = [miny, miny, maxy, maxy]
        Zlist = [minz, maxz]

        # Build 8 vertices via Cartesian product: for each Z height,
        # pair with each (X,Y) corner. Result: bottom 4 + top 4 vertices.
        VerticesList = []
        for h,xy in product(Zlist,zip(Xlist,Ylist)):
            VerticesList.append([xy[0],xy[1],h])

        blockMesh = self.blockMesh
        blockMesh['vertices'] = VerticesList
        # Hex vertex indices: 0..7 for a single block (standard OF ordering).
        blockMesh["hex"] = [x for x in range(8)]
        # Cell counts: ceiling division ensures at least one cell per direction.
        blockMesh["cellCount"] = [int(numpy.ceil((maxx-minx)/dx)),
                                  int(numpy.ceil((maxy - miny) / dy)),
                                  int(numpy.ceil((maxz - minz) / dz))]



    def set_blockMesh_blockHeight(self,Z,dz):
        """
            Sets the blockmesh dictionary
        Parameters
        ----------
        minx : float

        maxx : float
        miny : float
        maxy : float
        minz : float
        maxz : float
        dx : float
            The x spacing in m.
        dy : float
            The y spacing in m.
        dz: float
            The z = spacing in m
        Returns
        -------
            None
        """
        blockMesh = self.blockMesh
        verticsList=  blockMesh['vertices']
        for i in range(len(verticsList[4:])):
            verticsList[i][2] = Z

        minZ = verticsList[0][2]
        blockMesh["cellCount"][2] = (Z-minZ)/dz


    def foam_snappyhexmesh_addobject(self,dataOrFile,objectFile):
        """
            Adds a
        Parameters
        ----------
        dataOrFile
        objectFile

        Returns
        -------

        """

#### Lagrangian

class workflow_Lagrangian(abstractWorkflow):
    """Hermes workflow specialization for Lagrangian particle simulations."""

    def __init__(self ,workflowJSON,workflowHeraDocument=None,name=None, **kwargs):
        """Initialize the Lagrangian workflow."""
        super().__init__(workflowJSON=workflowJSON, workflowHeraDocument=workflowHeraDocument,name=name, **kwargs)


##############################################################################
##                          Workflow_StochasticLagrangianSolver
##############################################################################
class workflow_StochasticLagrangianSolver(workflow_Lagrangian):
    """Hermes workflow for stochastic Lagrangian dispersion simulations."""

    def __init__(self ,workflowJSON,workflowHeraDocument=None,name=None, **kwargs):
        """Initialize the stochastic Lagrangian solver workflow and validate parameters."""
        super().__init__(workflowJSON=workflowJSON, workflowHeraDocument=workflowHeraDocument,name=name, **kwargs)
        logger = get_classMethod_logger(self,"init")
        # Make sure that the
        # dispersionFlowField exists
        if 'dispersionFlowField' not in self.parameters.parameters:
            err = "The StochasticLagrangianSolver must have a dispersionFlowField specification in the parameters node"
            logger.error(err)
            raise ValueError(err)

    @property
    def dispersionName(self):
        """Return the name of the dispersion workflow."""
        return self.name


    @property
    def originalFlowFieldName(self):
        """Return the name of the original flow field."""
        return self.parameters.parameters['originalFlowField']

    @property
    def dispersionFlowFieldName(self):
        """Return the composite dispersion flow field name."""
        return f"{self.parameters.parameters['originalFlowField']}_DFF_{self.parameters.parameters['dispersionFlowField']}"

    @dispersionFlowFieldName.setter
    def dispersionFlowFieldName(self, value):
        """Set the dispersion flow field name parameter."""
        self.parameters.parameters['dispersionFlowField'] = str(value)

    @property
    def dispersionDuration(self):
        """Return the dispersion duration from the parameters."""
        return self.parameters.parameters['dispersionDuration']

##########################################################
##                          Workflow_simpleFoam
##########################################################
class workflow_simpleFoam(workflow_Eulerian):
    """Hermes workflow for the simpleFoam steady-state solver."""

    def __init__(self ,workflowJSON,workflowHeraDocument=None,name=None, **kwargs):
        """Initialize the simpleFoam workflow."""
        super().__init__(workflowJSON=workflowJSON, workflowHeraDocument=workflowHeraDocument,name=name, **kwargs)

##########################################################
##                          workflow_buoyantReactingFoam
##########################################################
class workflow_buoyantReactingFoam(workflow_Eulerian):
    """Hermes workflow for the buoyantReactingFoam compressible solver."""

    def __init__(self ,workflowJSON,workflowHeraDocument=None,name=None, **kwargs):
        """Initialize the buoyantReactingFoam workflow."""
        super().__init__(workflowJSON=workflowJSON, workflowHeraDocument=workflowHeraDocument,name=name, **kwargs)


class workflow_scalarTransportFoam(workflow_Eulerian):
    """Hermes workflow for the scalarTransportFoam solver."""

    def __init__(self ,workflowJSON,workflowHeraDocument=None,name=None, **kwargs):
        """Initialize the scalarTransportFoam workflow."""
        super().__init__(workflowJSON=workflowJSON, workflowHeraDocument=workflowHeraDocument,name=name, **kwargs)


##########################################################
##                          Workflow_indoorFOAMBoussinesq
##########################################################
class workflow_indoorFOAMBoussinesq(workflow_Eulerian):
    """Hermes workflow for the indoor Boussinesq approximation solver."""

    def __init__(self ,workflowJSON,workflowHeraDocument=None,name=None, **kwargs):
        """Initialize the indoorFOAMBoussinesq workflow."""
        super().__init__(workflowJSON=workflowJSON, workflowHeraDocument=workflowHeraDocument,name=name, **kwargs)

##########################################################
##                          Workflow_indoorFOAMBoussinesq
##########################################################
class workflow_homogenousWindLogProfile(workflow_Eulerian):
    """Hermes workflow for the homogeneous wind logarithmic profile solver."""

    def __init__(self ,workflowJSON,workflowHeraDocument=None,name=None, **kwargs):
        """Initialize the homogenousWindLogProfile workflow."""
        super().__init__(workflowJSON=workflowJSON, workflowHeraDocument=workflowHeraDocument,name=name, **kwargs)

