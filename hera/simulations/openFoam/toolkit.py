from ...toolkit import abstractToolkit
import os
import glob
from .datalayer.OFObjects import OFField
from .datalayer.hermesWorkflow import Workflow_Flow
from .datalayer.OFObjects import ofObjectHome

class OFToolkit(abstractToolkit):
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

    def __init__(self, projectName, filesDirectory=None, toolkitName="OFToolkit"):
        super().__init__(toolkitName=toolkitName, projectName=projectName, filesDirectory=filesDirectory)
        self._ofObjectHome = ofObjectHome()

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
