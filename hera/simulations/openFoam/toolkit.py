import json
import os
import glob
from .datalayer.OFObjects import OFField
from .datalayer.hermesWorkflow import Workflow_Flow
from .datalayer.OFObjects import OFObjectHome
from . import DECOMPOSED_CASE,RECONSTRUCTED_CASE,TYPE_VTK_FILTER
from ..hermesWorkflowToolkit import workflowToolkit
from ...datalayer import datatypes
from .analysis.VTKPipeline import VTKpipeline
import shutil
from . import StochasticLagrangian

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

    stochasticLagrangian = None

    @property
    def ofObjectHome(self):
        return self._ofObjectHome

    def __init__(self, projectName, filesDirectory=None):
        super().__init__(projectName=projectName,
                         filesDirectory=filesDirectory,
                         toolkitName="OFworkflowToolkit")

        self._ofObjectHome = OFObjectHome()
        self._analysis = analysis(self)

        self.stochasticLagrangian = StochasticLagrangian.stochasticLagrangianDataLayer(self)
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