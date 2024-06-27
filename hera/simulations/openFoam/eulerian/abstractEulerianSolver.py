from ..OFWorkflow import workflow_Eulerian
from ....utils.freeCAD import getObjFileBoundaries
from ....utils.logging import get_classMethod_logger
class absractEulerianSolver_toolkitExtension:

    toolkit = None

    analysis = None # link to the analysis of the StochasticLagrnagian
    presentation = None # link to the presentation of the stochasticLagrangian

    DOCTYPE_EULERIAN_CACHE = "eulerianCacheDocType"

    solverName = None
    incompressible = None

    def __init__(self,toolkit,solverName,incompressible):
        self.toolkit = toolkit
        self.solverName = solverName
        self.incompressible = incompressible
        #self.analysis = analysis(self)

    @property
    def flowType(self):
        return self.toolkit.SIMULATIONTYPE_INCOMPRESSIBLE if self.incompressible else SIMULATIONTYPE_COMPRESSIBLE


    def blockMesh_setBoundFromBounds(self, eulerianWF, minx,maxx,miny,maxy,minz,maxz,dx,dy,dz):
        """
            Set the boundaries of the blockMesh from the boundaries of the obj.

            Assuming a single block.
        Parameters
        ----------
        eulerianWorkFlow : eulerain workflow
                The workflow to change.

        minx : float
            The min x coordinate
        maxx : float
            The max x coordinate
        miny : float
            The min x coordinate
        maxy : float
            The max y coordinate
        minz : float
            The min z coordinate
        maxz : float
            The max z coordinate
        dx : float
            The x spacing in m.
        dy : float
            The y spacing in m.
        dz: float
            The z = spacing in m
        Returns
        -------

        """
        logger = get_classMethod_logger(self,blockMesh_setBoundFromBounds)
        if not isinstance(eulerianWF,workflow_Eulerian):
            err = f"The workflow is not eulerian (does not inherit the workflow_Eulerian class, cannot change the mesh. Got {type(eulerianWF)}"
            logger.error(err)
            raise ValueError(err)

        eulerianWF.set_blockMesh_boundaries(minx=minx,maxx=maxx,miny=miny,maxy=maxy,minz=minz,maxz=maxz,dx=dx,dy=dy,dz=dz)


    def blockMesh_setBoundFromFile(self,eulerianWorkFlow,fileName,dx,dy,dz):
        """
            Setting the bounds of the mesh according to the obj.

            Assuming a single block.
        Parameters
        ----------
        eulerianWorkFlow

        Returns
        -------

        """
        logger = get_classMethod_logger(self, "blockMesh_setBoundFromBounds")
        if not isinstance(eulerianWF, workflow_Eulerian):
            err = f"The workflow is not eulerian (does not inherit the workflow_Eulerian class, cannot change the mesh. Got {type(eulerianWF)}"
            logger.error(err)
            raise ValueError(err)

        corners = getObjFileBoundaries(fileName)
        eulerianWF.set_blockMesh_boundaries(**corners, dx=dx,dy=dy, dz=dz)


    def blockMesh_setDomainHeight(self,eulerianWorkFlow,Z,dz):
        """
            Sets the domain height and updates the number of

            Assuming a single block.
        Parameters
        ----------
        eulerianWorkFlow

        Returns
        -------

        """
        logger = get_classMethod_logger(self, "blockMesh_setDomainHeight")
        if not isinstance(eulerianWF, workflow_Eulerian):
            err = f"The workflow is not eulerian (does not inherit the workflow_Eulerian class, cannot change the mesh. Got {type(eulerianWF)}"
            logger.error(err)
            raise ValueError(err)

        corners = getObjFileBoundaries(fileName)
        eulerianWF.set_blockMesh_boundaries(**corners, dx=dx,dy=dy, dz=dz)

    def writeFieldInCase(self, fieldName,caseDirectory,components,xarrayData, boundaryConditions=None,data=None):
        """
            Writes
        Parameters
        ----------
        fieldName
        caseDirectory
        components
        xarrayData
        boundaryConditions
        data

        Returns
        -------

        """

        """
            Returns an OF Object.

        Parameters
        ----------
        fieldName  : string
            The name of the field.
            The field must exist in the field descriptions.

        boundaryConditions
        data : pandas, list
            The value of the field (ordered with the cell centers).

        Returns
        -------

        """

        return self.toolkit.OFObjectHome.getField(fieldName=fieldName,flowType=self.flowType,boundaryConditions=boundaryConditions,data=data)