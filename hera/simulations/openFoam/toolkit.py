import numpy
import os
import glob
from .OFWorkflow import workflow_Eulerian
from .preprocessOFObjects import OFObjectHome
from . import TYPE_VTK_FILTER
from ..hermesWorkflowToolkit import hermesWorkflowToolkit
from .postProcess.VTKPipeline import VTKPipeLine
from .lagrangian.StochasticLagrangianSolver import StochasticLagrangianSolver_toolkitExtension
from ...utils.jsonutils import loadJSON
from ...utils.logging import get_classMethod_logger
from evtk import hl as evtk_hl
import dask.dataframe as dask_dataframe
from itertools import chain
from itertools import product
from collections.abc import Iterable
from . import FLOWTYPE_DISPERSION, FIELDTYPE_SCALAR, FIELDTYPE_TENSOR, \
    FIELDTYPE_VECTOR, CASETYPE_DECOMPOSED, CASETYPE_RECONSTRUCTED,FLOWTYPE_COMPRESSIBLE, FLOWTYPE_INCOMPRESSIBLE
from .eulerian.buoyantReactingFoam import buoyantReactingFoam_toolkitExtension
import pandas
from dask.delayed import delayed
import dask

class OFToolkit(hermesWorkflowToolkit):
    """
        The goal of this toolkit is to provide the functions that are required to run workflows.
        and to mange the workflows in the DB.

        This toolkit might relay on the hermes project in order to manipulate the nodes
        of the workflow. (TBD).

    """
    TIME_STEADYSTATE = "steadyState"
    TIME_DYNAMIC = "dynamic"

    FLOWTYPE_COMPRESSIBLE = FLOWTYPE_COMPRESSIBLE
    FLOWTYPE_INCOMPRESSIBLE = FLOWTYPE_INCOMPRESSIBLE
    FLOWTYPE_DISPERSION = FLOWTYPE_DISPERSION

    FIELDTYPE_SCALAR = FIELDTYPE_SCALAR
    FIELDTYPE_VECTOR = FIELDTYPE_VECTOR
    FIELDTYPE_TENSOR = FIELDTYPE_TENSOR

    CASETYPE_DECOMPOSED = CASETYPE_DECOMPOSED
    CASETYPE_RECONSTRUCTED = CASETYPE_RECONSTRUCTED

    OF_FLOWDISPERSION = "flowDispersion"

    stochasticLagrangian = None

    buoyantReactingFoam = None

    def __init__(self, projectName, filesDirectory=None):
        super().__init__(projectName=projectName,
                         filesDirectory=filesDirectory,
                         toolkitName="OFworkflowToolkit")

        self.OFObjectHome = OFObjectHome()
        self._analysis = Analysis(self)
        self._presentation = Presentation(self, self.analysis)
        self.stochasticLagrangian = StochasticLagrangianSolver_toolkitExtension(self)
        self.buoyantReactingFoam  = buoyantReactingFoam_toolkitExtension(self)

    def processorList(self, caseDirectory):
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

    def getHermesWorkflow_Flow(self, workflowfile):
        """
            Returns the workflow of the requested JSON file.
        Parameters
        ----------
        workflowfile

        Returns
        -------

        """
        return workflow_Eulerian(workflowfile)

    def getMesh(self, caseDirectory, readParallel=True, time=0):
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

            readParallel: bool
                    If parallel case exists, read it .
        Returns
        -------
            pandas dataframe with the points in the columns             x,y,z
            the index column (don't mix up with the index of pandas)  is the sequential number of the point.

            If the case is decomposed, return processorNumber and index columns.
            The index is the internal order in the processor.
        """

        # 1. Run the postProcess utility to set the cell centers
        logger = get_classMethod_logger(self, "getMesh")
        logger.info(f"Start. case {caseDirectory}. Current directory is : {os.getcwd()}.")

        casePointer = "" if caseDirectory == os.getcwd() else f"-case {caseDirectory}"

        useParallel = False
        if readParallel:
            logger.debug(f"Attempt to load parallel case")
            # Check if the case is decomposed, if it is, run it.
            proc0dir = os.path.join(caseDirectory, "processor0")

            if os.path.exists(proc0dir):
                logger.debug(f"Found parallel case, using decomposed case")
                useParallel = True
            else:
                logger.debug(f"parallel case NOT found. Using composed case")

        # Calculating the cell centers
        checkPath = os.path.join(caseDirectory, "processor0", str(time), "C") if useParallel else os.path.join(
            caseDirectory, str(time), "C")
        parallelExec = "-parallel" if useParallel else ""
        caseType = "decomposed" if useParallel else "composed"
        if not os.path.exists(checkPath):
            logger.debug(f"Cell centers does not exist in {caseType} case. Calculating...")
            os.system(f"foamJob {parallelExec} {casePointer} -wait postProcess -func writeCellCentres ")
            logger.debug(f"done: foamJob {parallelExec} -wait postProcess -func writeCellCentres {casePointer}")
            if not os.path.exists(checkPath):
                logger.error("Error running the writeCellCentres. Check mesh")
                raise RuntimeError("Error running the writeCellCentres. Check mesh")
        else:
            logger.debug(f"Cell centers exist in {caseType} case.")

        logger.debug(f"Loading the cell centers in time {time}. Usint {caseType}")
        cellCenters = self.OFObjectHome.readFieldFromCase(fieldName="cellCenters",
                                                         flowType=FLOWTYPE_INCOMPRESSIBLE,
                                                         caseDirectory=caseDirectory,
                                                         timeStep=time,
                                                         readParallel=readParallel)
        return cellCenters

    def createEmptyCase(self, caseDirectory: str, fieldList: list, flowType: str, additionalFieldsDescription=dict()):
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

        flowType : str
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
        logger = get_classMethod_logger(self,"createEmptyCase")
        logger.info(f"Making case {caseDirectory} with fields {','.join(fieldList)}")

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

        for fieldName, fieldDefs in fileaddition.items():
            logger.info(f"Adding temporary field {fieldName} to the field directory")
            self.OFObjectHome.addFieldDefinitions(fieldName=fieldName, **fieldDefs)

        # Makes the empty fields
        for fieldName in fieldList:
            logger.info(f"Creating field {fieldName}")
            self.writeEmptyField(fieldName=fieldName,flowType=flowType,caseDirectory=caseDirectory,timeOrLocation=0)
            self.writeEmptyField(fieldName=fieldName,flowType=flowType,caseDirectory=caseDirectory,timeOrLocation="0.orig")
            self.writeEmptyField(fieldName=fieldName,flowType=flowType,caseDirectory=caseDirectory,timeOrLocation="0.parallel",writeProcBoundary=True)

    def writeEmptyField(self,fieldName,flowType,caseDirectory,timeOrLocation=0,readBoundaryFromCase=False,writeProcBoundary=False):
        """
            Writes an empty field in the case.

            If the readBoundaryField is True, then the field is written with the relevant boundaries (that are red from the case).

        Parameters
        ----------
        fieldName : str
            The name of the field
        flowType : str
            The flow type (compressible/incompressible)
        timeOrLocation : float  [default : 0]
            The time to write the new field in
        caseDirectory : str
            The name of te new case directory
        readParallel : bool
            If true, then attempt to write as parallel fields.
        readBoundaryFromCase : bool
            If True, tries to read the boundary names from the case.
            Otherwise write the general ".*" boundary field.

        writeProcBoundary : bool
           If true writes the proc boundaries.
           If readBoundaryFromCase is True, then write the specific proc for each processor (when it is parallel)
           Otherwise, write the "proc.*" boundary field.

        Returns
        -------

        """
        logger = get_classMethod_logger(self,"writeEmptyField")
        logger.info(f"Creating the field: {fieldName}. ")

        field = self.OFObjectHome.getEmptyField(fieldName=fieldName, flowType=flowType)

        if readBoundaryFromCase:
            field.readBoundariesFromCase(caseDirectory,readParallel=True)

        if writeProcBoundary:
            field.addProcBoundary()
        field.writeToCase(caseDirectory=caseDirectory, timeOrLocation=timeOrLocation)

    #############################################################

    def template_add(self, name, objFile, workflowObj=None):
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

    def xarrayToSetFieldsDictDomain(self, xarrayData, xColumnName="x", yColumnName="y", zColumnName="z", time=None,
                                    timeColumn="time", **kwargs):
        """
            Converts the xarray to the setFields dict of the internal domain.

            Not debugged.

            Use

        Parameters
        ----------
        xarrayData : xarray dataset
                The data set as xarray.
        xColumnName : string
            The coordinate name of the x-axis
            If None, ignore that coordinate (the mesh should be in the right dimension).
        yColumnName: string
            The coordinate name of the y-axis
            If None, ignore that coordinate (the mesh should be  in the right dimension).
        zColumnName: string
            The coordinate name of the z-axis
            If None, ignore that coordinate (the mesh should be  in the right dimension).
        time : string
            The time is not None, select the closest (nearest neighbours).

            If None ignore the index.
        kwargs : dict
            A mapping of OpenFoamField -> dataset Field
            If the dataset Field is a tuple, then write it as vector.
            Each component can be either a field, float (to keep fixed value).
            If you want to map a function, just create its value as a Field name.

            so
                U = ("u","v",0)
                will map the feld 'u' as Ux, 'v' as Uy and Uz=0


        Returns
        -------
            string (setField dict).

        """
        logger = get_classMethod_logger(self, "xarrayToSetFieldsDict")
        logger.debug(f"------------ Start  : {logger.name} ")

        coordList = []
        coordNames = []
        if xColumnName is not None:
            coordList.append(xarrayData.coords[xColumnName].shape[0] - 1)
            coordNames.append(xColumnName)
        if yColumnName is not None:
            coordList.append(xarrayData.coords[yColumnName].shape[0] - 1)
            coordNames.append(yColumnName)
        if zColumnName is not None:
            coordList.append(xarrayData.coords[zColumnName].shape[0] - 1)
            coordNames.append(zColumnName)

        ret = []
        arryToOFvector = lambda arry: f"({' '.join(arry)} )"
        for X in product(*coordList):
            lowerLeft = [xarrayData.coords[coordNames[coordIndx]][ind].item() for coordIndx, ind in enumerate(X)]
            upperRight = [xarrayData.coords[coordNames[coordIndx]][ind + 1].item() for coordIndx, ind in enumerate(X)]

            fielaValueStr = ""
            for OFField, fieldName in kwargs.items():

                accessDict = dict([(coordname, indx) for (coordname, indx) in zip(coordNames, X)])
                if time is not None:
                    timeDict = {timeColumn: time}
                    valArray = xarrayData.sel(**timeDict, method='nearest')

                else:
                    valArray = xarrayData

                if isinstance(fieldName, Iterable):
                    valArray = []
                    for mappedField in fieldName:
                        if isinstance(mappedField, str):
                            valArray.append(valArray[mappedField].isel(**accessDict).item())
                        elif isinstance(mappedField, float) or isinstance(mappedField, int):
                            valArray.append(mappedField)
                        else:
                            err = f"The mapping {OFField}->{fieldName} contains a value that is not a string or number. "
                            raise ValueError(err)
                    if len(fieldName) == 3:
                        fielaValueStr += f"volVectorFieldValue {OFField} {arryToOFvector(valArray)}\n"
                    elif len(fieldName) == 9:
                        fielaValueStr += f"volTensorFieldValue {OFField} {arryToOFvector(valArray)}\n"
                    else:
                        err = f"The number of components in the  mapping {OFField}->{fieldName} must be 1,3 or 9. got {len(fieldName)}."
                        raise ValueError(err)
                else:
                    value = valArray[fieldName].isel(**accessDict).item()

                    fielaValueStr += f"volScalarFieldValue {OFField} {value}\n"

            BoxRecord = f"""
                boxToCell 
                {{
                    box {arryToOFvector(lowerLeft)} {arryToOFvector(upperRight)};
                    fieldValues 
                    (
                        {fielaValueStr}
                    
                    );                  
                }}
            """
            ret.append(BoxRecord)
        return "\n".join(ret)


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

    def __init__(self, datalayer):
        """
            Initializes the analysis layer.
        Parameters
        ----------
        datalayer : OFToolkit
                The datalayer to use.

        """
        self._datalayer = datalayer

    def getVTKPipeline(self):
        """
            Creates a new VTK pipeline from the simulation.

        Parameters
        ----------
        """
        return VTKPipeLine()

    def getFiltersDocuments(self, nameOrWorkflowFileOrJSONOrResource, filterName=None):
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
        wrkflow = self.datalayer.getWorkflowDocumentFromDB(nameOrWorkflowFileOrJSONOrResource)
        if wrkflow is None:
            raise ValueError(
                f"The case {nameOrWorkflowFileOrJSONOrResource} was not found in the project {self.datalayer.projectName}")
        simulationName = wrkflow['desc']['simulationName']

        qry = dict(simulationName=simulationName)
        if filterName is not None:
            qry['filterName'] = filterName

        return self.datalayer.getCacheDocuments(type=TYPE_VTK_FILTER, **qry)

    def executeAndLoad(self,vtkPipeline,sourceOrName=None,timeList=None, tsBlockNum=50, overwrite=False,append=False,overwriteMetadata=False):
        """

        Parameters
        ----------
        vtkPipeline : VTKPipeLine

        sourceOrName : Non or string
            The name of the filter/reader to start with. None will be the default reader.


        timeList : list
                The list of timesteps to read.

        tsBlockNum: int
                The block number
        overwrite : bool  [default False]
               If true, overwrite on the parquet/netcdf file.

        append : bool [default False]
               If true, append to existing parquet

        overwriteMetadata : bool [default False]
                If True, overwrite the meatadata of the filters.

        Returns
        -------

        """
        logger = get_classMethod_logger(self,"executeAndLoad")
        logger.info(f"Executing the VTK pipeline")
        vtkPipeline.execute(sourceOrName=sourceOrName,timeList=timeList, tsBlockNum=tsBlockNum, overwrite=overwrite,append=append)
        logger.info(f"Loading the VTK pipeline to project {self.datalayer.projectName}")
        vtkPipeline.loadToProject(datalayer =self.datalayer, overwrite=overwriteMetadata)

class Presentation:

    def __init__(self, datalayer, analysis):
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

    def toUnstructuredVTK(self, data, outputdirectory, filename, fieldList=None, xcoord="globalX",ycoord="globalY",zcoord="globalZ", timecoord="time",overwrite=False):
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
            fieldList = [x for x in data.columns if x not in [xcoord, ycoord, zcoord]]

        if isinstance(data, dask_dataframe.DataFrame):
            logger.debug("The input type is dask dataframe. Getting timesteps and then using partition")
            timeList = [z for z in chain(*[x.index.unique().compute() for x in data.partitions])]
            for indx, time in enumerate(timeList):
                timeData = data.loc[time].compute()
                outdata = dict((x, timeData[x].values) for x in fieldList)
                finalfile = f"{namePath}_{indx}"
                logger.info(f"Writing time step {time} to {finalfile}")

                x = timeData[xcoord].values
                y = timeData[ycoord].values
                z = timeData[zcoord].values

                if os.path.exists(finalfile):
                    logger.debug(f"The file exists, remove it if overwrite is True")
                    if overwrite:
                        os.remove(finalfile)
                    else:
                        err = f"File {finalfile} exists. Use --overwrite to overwrite the files with newer version"
                        logger.error(err)
                        raise FileExistsError(err)

                evtk_hl.pointsToVTK(finalfile, x, y, z, outdata)
        else:
            # all the data is loaded:
            for indx, (timeName, timeData) in enumerate(data.groupby(timecoord)):
                outdata = dict((x, timeData[x].values) for x in fieldList)
                finalfile = f"{namePath}_{indx}"
                logger.info(f"Writing time step {timeName} to {finalfile}")

                x = timeData[xcoord].values
                y = timeData[ycoord].values
                z = timeData[zcoord].values

                if os.path.exists(finalfile):
                    logger.debug(f"The file exists, remove it if overwrite is True")
                    if overwrite:
                        os.remove(finalfile)
                    else:
                        err = f"File {finalfile} exists. Use --overwrite to overwrite the files with newer version"
                        logger.error(err)
                        raise FileExistsError(err)

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

    def loadLagrangianDataParallel(self,
                                   casePath,
                                   times=None,
                                   parallelCase=True,
                                   withVelocity=False,
                                   withReleaseTimes=False,
                                   withMass=False,
                                   cloudName="kinematicCloud"):
        """
            Extracts results of an OpenFOAM cloud particles.

            TODO : need to be debugged as valid only in the old version.

        Parameters
        -----------

            times: list
                 a list of time steps to extract.
                 If None, it extracts all time steps in the casePath.
            file: str
                The name for a file, in which the data is saved. Default is "<cloud name>_data.parquet", in the current working directory.

            withVelocity: bool
                    True: extract the particles' velocities in addition to their positions.
                    Default is False.

            **kwargs:
                Any additional parameters to add to the description in the DB.
        Returns
        --------
            document of the data.
        """

        finalCasePath = os.path.abspath(casePath)
        loader = lambda timeName: readLagrangianRecord(timeName,
                                                       casePath=finalCasePath,
                                                       withVelocity=withVelocity,
                                                       withReleaseTimes=withReleaseTimes,
                                                       withMass=withMass,
                                                       cloudName=cloudName)

        if parallelCase:
            processorList = [os.path.basename(proc) for proc in glob.glob(os.path.join(finalCasePath, "processor*"))]
            if len(processorList) == 0:
                raise ValueError(f"There are no processor* directories in the case {finalCasePath}. Is it parallel?")

            if times is None:
                timeList = sorted([x for x in os.listdir(os.path.join(finalCasePath, processorList[0])) if (
                        os.path.isdir(os.path.join(finalCasePath, processorList[0], x)) and
                        x.isdigit() and
                        (not x.startswith("processor") and x not in ["constant", "system", "rootCase", 'VTK']))],
                                  key=lambda x: int(x))
            else:
                timeList = numpy.atleast_1d(times)

            data = dask.dataframe.from_delayed(
                [delayed(loader)(os.path.join(processorName, timeName)) for processorName, timeName in
                 product(processorList, timeList[1:])])
        else:

            if times is None:
                timeList = sorted([x for x in os.listdir(finalCasePath) if (
                        os.path.isdir(x) and
                        x.isdigit() and
                        (not x.startswith("processor") and x not in ["constant", "system", "rootCase", 'VTK']))],
                                  key=lambda x: int(x))
            else:
                timeList = numpy.atleast_1d(times)

            data = dask.dataframe.from_delayed([delayed(loader)(timeName) for timeName in timeList])

        return data


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
