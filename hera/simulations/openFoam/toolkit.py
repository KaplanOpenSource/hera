import numpy
import os
import glob
from .OFObjects import OFField
from .OFWorkflow import workflow_Eulerian
from .OFObjects import OFObjectHome
from . import DECOMPOSED_CASE, TYPE_VTK_FILTER
from ..hermesWorkflowToolkit import workflowToolkit
from .VTKPipeline import VTKpipeline
from .lagrangian.StochasticLagrangianSolver import StochasticLagrangianSolver_toolkitExtension
from ...utils.jsonutils import loadJSON,compareJSONS
from ...utils.logging import get_classMethod_logger
from evtk import hl as evtk_hl
import dask.dataframe as dd
from itertools import chain

class OFToolkit(workflowToolkit):
    """
        The goal of this toolkit is to provide the functions that are required to run workflows.
        and to mange the workflows in the DB.

        This toolkit might relay on the hermes project in order to manipulate the nodes
        of the workflow. (TBD).

    """
    TIME_STEADYSTATE = "steadyState"
    TIME_DYNAMIC = "dynamic"

    SIMULATIONTYPE_COMPRESSIBLE = "compressible"
    SIMULATIONTYPE_INCOMPRESSIBLE = "incompressible"
    SIMULATIONTYPE_DISPERSION = "dispersion"

    OF_FLOWDISPERSION = "flowDispersion"

    stochasticLagrangian = None

    def __init__(self, projectName, filesDirectory=None):
        super().__init__(projectName=projectName,
                         filesDirectory=filesDirectory,
                         toolkitName="OFworkflowToolkit")

        self.OFObjectHome = OFObjectHome()
        self._analysis = Analysis(self)
        self._presentation = Presentation(self,self.analysis)
        self.stochasticLagrangian = StochasticLagrangianSolver_toolkitExtension(self)

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
        return workflow_Eulerian(workflowfile)

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
        logger = get_classMethod_logger(self,"getMesh")
        logger.info(f"Start. case {caseDirectory}. Current directory is : {os.getcwd()}.")

        casePointer = "" if caseDirectory == os.getcwd() else f"-case {caseDirectory}"

        useParallel= False
        if parallel:
            logger.debug(f"Attempt to load parallel case")
            # Check if the case is decomposed, if it is, run it.
            proc0dir = os.path.join(caseDirectory,"processor0")

            if os.path.exists(proc0dir):
                logger.debug(f"Found parallel case, using decomposed case")
                useParallel = True
            else:
                logger.debug(f"parallel case NOT found. Using composed case")

        # Calculating the cell centers
        checkPath = os.path.join(caseDirectory,"processor0",str(time),"C") if useParallel else os.path.join(caseDirectory,str(time),"C")
        parallelExec = "-parallel" if useParallel else ""
        caseType = "decomposed" if useParallel else "composed"
        if not os.path.exists(checkPath):
            logger.debug(f"Cell centers does not exist in {caseType} case. Calculating...")
            os.system(f"foamJob {parallelExec} -wait postProcess -func writeCellCentres {casePointer}")
            logger.debug(f"done: foamJob {parallelExec} -wait postProcess -func writeCellCentres {casePointer}")
            if not os.path.exists(checkPath):
                logger.error("Error running the writeCellCentres. Check mesh")
                raise RuntimeError("Error running the writeCellCentres. Check mesh")
        else:
            logger.debug(f"Cell centers exist in {caseType} case.")

        cellCenters = OFField(name="C",dimensions="",componentNames=['x','y','z'])
        logger.debug(f"Loading the cell centers in time {time}. Usint {caseType}")
        ret =  cellCenters.load(caseDirectory,times=time,parallelCase=useParallel)

        logger.info(f"--- End ---")
        return ret

    def createEmptyCase(self, caseDirectory :str, fieldList : list, simulationType:str, additionalFieldsDescription  =dict()):
        """
            Creates an empty case directory for the simulation.
            fields is a list of fields to create in the case directory.

            The simulation type (copressible, incompressible, dispersion) is needed to get the dimensions and components
            of the fields. If the fields are not in the standard list then their description can be supplied in the
            additionalFieldsDescription parameters

        Parameters
        ----------
        caseDirectory : str
            The case directory to create

        fieldList : list
            The list of field names to create

        simulationType : str
            compressible, incompressible or dispersion.
            The dimension of the fields is determined by the type of simulation

        additionalFieldsDescription : dict | str
            Definition of additional fields:
            has the structure :
            {
                dimensions : {kg : .., m : ..},
                componentNames : None|list
            )
            the keys for the dimensions are kg,m,s,K,mol,A,cd

            Can also be a JSON file name.

        Returns
        -------

        """
        print(f"Making case {caseDirectory} with fields {','.join(fieldList)}")

        # Make the case :
        if os.path.isfile(caseDirectory):
            raise ValueError(
                f"The file {caseDirectory} exists as a file. Cannot create a directory. Please remove/rename it and rerun. ")

        os.makedirs(os.path.join(caseDirectory, "constant"), exist_ok=True)
        os.makedirs(os.path.join(caseDirectory, "system"), exist_ok=True)
        os.makedirs(os.path.join(caseDirectory, "constant", "triSurface"), exist_ok=True)
        os.makedirs(os.path.join(caseDirectory, "0"), exist_ok=True)
        os.makedirs(os.path.join(caseDirectory, "0.orig"), exist_ok=True)
        os.makedirs(os.path.join(caseDirectory, "0.parallel"), exist_ok=True)

        fileaddition = dict()
        if additionalFieldsDescription is not None:
                fileaddition = loadJSON(additionalFieldsDescription)

        # Makes the empty fields
        for fieldName in fieldList:
            field = self.OFObjectHome.getField(fieldName, flowType=simulationType, additionalFieldsDescription=fileaddition)
            field.write(caseDirectory=caseDirectory, location=0)
            field.write(caseDirectory=caseDirectory, location="0.orig")
            field.write(caseDirectory=caseDirectory, location="0.parallel",parallel=False,parallelBoundary=True)


    def template_add(self,name,objFile,workflowObj=None):
        """
            Adds a templates to the toolkit.

            Templates can be
                - Flow : Holds Hermes flow templates.
                - Node : Holds a hermes node objects
                - Field : Holds a field templates.
                        This can be
                            * xarray
                            * pandas/dask
                            * constant

        Parameters
        ----------
        name
        objFile
        workflowObj

        Returns
        -------

        """
        pass

class Analysis:
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
        wrkflow = self.datalayer.getCaseListDocumentFromDB(nameOrWorkflowFileOrJSONOrResource)
        if wrkflow is None:
            raise ValueError(f"The case {nameOrWorkflowFileOrJSONOrResource} was not found in the project {self.datalayer.projectName}")
        simulationName = wrkflow['desc']['simulationName']

        qry = dict(simulationName=simulationName)
        if filterName is not None:
            qry['filterName'] = filterName

        return self.datalayer.getCacheDocuments(type=TYPE_VTK_FILTER,**qry)


class Presentation:

    def __init__(self,datalayer,analysis):
        self.datalayer = datalayer
        self.analysis = analysis

    def to_paraview_CSV(self, data, outputdirectory, filename, timeFactor=1):
        """
            Writes the globalPositions (globalX,globalY,globalZ) as  CSV for visualization in paraview.
            In paraview, each timestep is a different file.

        Parameters
        -----------
        data: dask.dataframe or pandas.dataframe
            The data to present

        outputdirectory: str
            The directory to write the files in

        timeFactor : int
            Multiply the time by a factro to make the time step round (so that paraview will recognize it).

        filename: str
            The filename to write.

        Returns
        -------
            None
        """
        for times, timedata in data.groupby("time"):
            with open(os.path.join(outputdirectory, f"{filename}_{str(int(timeFactor * times)).replace('.', '_')}.csv"),
                      "w") as outputfile:
                outputfile.writelines(timedata[['globalX', 'globalY', 'globalZ']].to_csv(index=False))

    def toUnstructuredVTK(self, data, outputdirectory, filename, timeNameOutput=True,fieldList=None,xcoord="x",ycoord="y",zcoord="z",timecoord="time"):
        """
            Writes the data as a VTK vtu file.



        :param data: panas.Dataframe.
                The data in a dataframe of pandas.
        :param outputdirectory:
        :param filename: str
                    The filename
        :param timeNameOutput: bool
                    If true, use the time as the suffix to the filename. Else, use running number.
        :return:
        """
        logger = get_classMethod_logger(self, "toUnstructuredVTK")
        namePath = os.path.join(outputdirectory, filename)
        os.makedirs(outputdirectory, exist_ok=True)
        logger.info(f"Writing unstructured VTK file to base filename: {namePath}")

        logger.debug("Getting the field names")
        if fieldList is None:
            fieldList = [x for x in data.columns if x not in [xcoord,ycoord,zcoord]]

        if isinstance(data,dd.DataFrame):
            logger.debug("The input type is dask dataframe. Getting timesteps and then using partition")
            timeList = [z for z in chain(*[x.index.unique().compute() for x in data.partitions])]
            for indx,time in enumerate(timeList):
                timeData = data.loc[time].compute()
                outdata = dict((x,timeData[x].values) for x in fieldList)
                finalfile = f"{namePath}_{indx}"
                logger.info(f"Writing time step {time} to {finalfile}")

                x = timeData[xcoord].values
                y = timeData[ycoord].values
                z = timeData[zcoord].values
                evtk_hl.pointsToVTK(finalfile, x, y, z, outdata)
        else:
            # all the data is loaded:
            for indx, (timeName, timeData) in enumerate(data.groupby("time")):
                outdata = dict((x,timeData[x].values) for x in fieldList)
                finalfile = f"{namePath}_{indx}"
                logger.info(f"Writing time step {time} to {finalfile}")
                
                x = timeData[xcoord].values
                y = timeData[ycoord].values
                z = timeData[zcoord].values
                evtk_hl.pointsToVTK(finalfile, x, y, z, outdata)

    def toStructuredVTK(self, data, outputdirectory, filename, extents, dxdydz, timeNameOutput=True):
        """
            Converts the data to structured grid, and calcualtes the cocnentration

            ** for now, the ppm is of SF6.

        :param data:
        :param outputdirectory:
        :param filename:
        :param timeNameOutput: bool
                    Use the time or sequence for the output file name.
        :param extents: dict
                    with keys: xmin,xmax,ymin,ymax,zmin,zmax of the entire domain.
        :param dxdydz: float
                    The mesh steps.

        :return:
        """
        namePath = os.path.join(outputdirectory, filename)
        os.makedirs(outputdirectory, exist_ok=True)

        for indx, (timeName, timeData) in enumerate(data.groupby("time")):
            # assign the timestep into the large mesh of that time step.
            fulldata = self.analysis.calcConcentrationTimeStepFullMesh(timeData, extents, dxdydz=dxdydz)

            finalFile = f"{namePath}_{int(timeName)}" if timeNameOutput else f"{namePath}_{indx}"

            X, Y, Z = numpy.meshgrid(fulldata.xI, fulldata.yI, fulldata.zI)

            C = numpy.ascontiguousarray(fulldata.transpose("yI", "xI", "zI").values)

            data = dict(C_kg_m3=C)  # 1kg/m**3=160000ppm
            evtk_hl.structuredToVTK(finalFile, X, Y, Z, pointData=data)




    # def compareWorkflows_AllGroups(self,workflowsTypes):
    #     """
    #         Lists all the simulations of the type (flow, dispersion, flowdispersion)
    #         and return the differences between each group
    #     Parameters
    #     ----------
    #     workflowsTypes
    #
    #     Returns
    #     -------
    #
    #     """
    #     groupsList = pandas.DataFrame(dict(groups=[doc.desc['groupName'] for doc\
    #                                                in self.getSimulationDocuments(type=workflowsTypes)])).drop_duplicates()
    #
    #     resList = []
    #     for group in groupsList:
    #         wfDict = dict([(doc.desc['name'], doc.desc['name']['flowParameters']) for doc in
    #                        self.getSimulationDocuments(type=workflowsTypes,
    #                                                    groupName=group)])
    #         res = compareJSONS(wfDict)
    #         resList.append(res.assign(groupName=group))
    #
    #     return pandas.concat(resList)

    ########################################## baseTemplateHandler
    #
    #  Retrieves the base workflows.
    #
    # def getBaseFlow_directory(self, parameters : dict):
    #     """
    #         Return the name of the requested directory.
    #
    #     Parameters
    #     ----------
    #     parameters: dict
    #         The json of the directory:
    #             {
    #                 name : ...
    #             }
    #
    #     projectName : str
    #         The name of the project (not used in this procedure).
    #
    #     Returns
    #     -------
    #
    #     """
    #     return parameters['name'],os.path.basename(parameters['name'])

    # def getBaseFlow_simulationName(self, parameters):
    #     """
    #         Check if the directory is already the db.
    #         If it is, then return the id.
    #
    #     Parameters
    #     ----------
    #     parameters: dict
    #         The json of the directory:
    #             {
    #                 name : ...
    #             }
    #     projectName : str
    #         The name of the project.
    #
    #     Returns
    #     -------
    #
    #     """
    #     docList = self.getSimulationsDocuments(workflowName=parameters['name'],type=workflowToolkit.DOC_TYPE)
    #     if len(docList) >0:
    #         if len(docList) > 1:
    #             import warnings
    #             warnings.warn(f"Found more than 1 simulation with the name {parameters['name']}. Return the first one")
    #         return docList[0].resource,parameters['name']
    #     else:
    #         raise ValueError(f"Cannot find flows with the name {parameters['name']}. Use hera-OF-flows list to see the names of existing simulations.old.")
    #
