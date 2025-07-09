
from itertools import product
import pandas
import numpy
import dask.dataframe as dd
import glob
import shutil
import os
import xarray
import tqdm
from hera.simulations.openFoam import CASETYPE_DECOMPOSED,CASETYPE_RECONSTRUCTED
from hera import get_classMethod_logger
from deprecated import deprecated
from dask.diagnostics import ProgressBar

#### import the simple module from the paraview
try:
    import paraview.vtk.numpy_interface.dataset_adapter as dsa
    import paraview.simple as pvsimple
    from paraview import servermanager
    #### disable automatic camera reset on 'Show'
    pvsimple._DisableFirstRenderCameraReset()
except ImportError:
    print("paraview module is not Found!. VTK pipeline wont work")

from hera.utils.logging import helpers as hera_logging

class paraviewOpenFOAM:
    """
        A class to extract openFOAM file format
        using VTK filters and write as parquet or netcdf files.
    """

    _componentsNames = None  # names of components for reading.

    def __init__(self, casePath,caseType=CASETYPE_DECOMPOSED, servername=None,name="mainreader"):
        """
            Initializes the paraviewOpenFOAM class.

            Supports single case or decomposed case and
            works with paraview server if initializes.

        Parameters
        -----------

        casePath: str
                    A full path to the case directory.

        outputPath : str
                    The path to save the results in

        CaseType:  str
                Either 'Decomposed Case' for parallel cases or 'Reconstructed Case'
                for single processor cases.

        servername: str
                if None, work locally.
                connection string to the paraview server.

                The connection string is printed when the server is initialized.

        """
        logger = get_classMethod_logger(self,"__init__")
        if servername is not None:
            pvsimple.Connect(servername)

        self._componentsNames = {}

        self.casePath = casePath
        self.caseType = caseType

        # Array shape length 1 - scalar.
        #					 2 - vector.
        #					 3 - tensor.
        # the dict holds their names.
        self._componentsNames = {(): "",
                                 (0,): "_x",
                                 (1,): "_y",
                                 (2,): "_z",
                                 (0, 0): "_xx",
                                 (0, 1): "_xy",
                                 (0, 2): "_xz",
                                 (1, 0): "_yx",
                                 (1, 1): "_yy",
                                 (1, 2): "_yz",
                                 (2, 0): "_zx",
                                 (2, 1): "_zy",
                                 (2, 2): "_zz"}

    def initializeReader(self, readerName="reader"):
        """
            Constructs a reader and register it in the vtk pipeline.

            Handles either parallel or single format.

        Parameters
        -----------

        readerName: str
                The name of the reader.  (of the pipline).
                When using server, then you can have different pipelines with different names.

        casePath: str
                a full path to the case directory.
        CaseType: str
                Either 'Decomposed Case' for parallel cases or 'Reconstructed Case'
                for single processor cases.
        fieldnames: list of str
                List of field names to load.
                if None, read all the fields.

        servername: str
                The address of pvserver. If None, use the local single threaded case.
        :return:
                the reader
        """
        reader = pvsimple.OpenFOAMReader(FileName="%s/tmp.foam" % self.casePath, CaseType=self.caseType,guiName=readerName)
        reader.MeshRegions.SelectAll()
        possibleRegions = list(reader.MeshRegions)
        reader.MeshRegions = ['internalMesh']
        reader.UpdatePipeline()

        # setting the local variable.
        self.reader = reader
        self.readerName  = readerName

        return reader

    def readTimeSteps(self, datasourcenamedict, timelist=None, fieldnames=None, regularMesh=False):
        """
            reads a list of datasource lists to a dictionary

        Parameters
        ----------

        datasourcenamedict: dict
                filtername -> output path .

        timelist: list
                The list of times to read.
        fieldnames:
                The list of fields to write.
        regtularMesh
                convert pandas results to xarray (works only for regular grids).

        Return
        ------

        For each time step.
                    A map datasourcename -> pandas
        """
        logger = get_classMethod_logger(self, "readTimeSteps")
        for timeslice in timelist:
            # read the timestep.
            logger.info("\r Reading time step %s" % timeslice)
            ret = {}
            for filterName,outputFile in datasourcenamedict.items():
                datasource = pvsimple.FindSource(filterName)
                logger.debug(f"Reading source {filterName}")
                rt = self._readTimeStep(datasource, timeslice, fieldnames, regularMesh)
                if rt is not None:
                    ret[filterName] = rt
            yield ret

    def _readTimeStep(self, datasource, timeslice, fieldnames=None, regularMesh=False):
        # read the timestep.
        logger = get_classMethod_logger(self, "_readTimeStep")

        datasource.UpdatePipeline(timeslice)
        rawData = servermanager.Fetch(datasource)
        data = dsa.WrapDataObject(rawData)

        if isinstance(data.Points, dsa.VTKNoneArray):
            logger.debug("No data exists for filter... return with None")
            return None
        elif isinstance(data.Points, dsa.VTKArray):
            points = numpy.array(data.Points).squeeze()
        else:
            points = numpy.concatenate([numpy.array(x) for x in data.Points.GetArrays()]).squeeze()

        logger.debug(f"Filter has {points.shape[0]} points. Building basic dataFrame. ")
        curstep = pandas.DataFrame()

        # create index
        curstep['x'] = points[:, 0]
        curstep['y'] = points[:, 1]
        curstep['z'] = points[:, 2]
        curstep['time'] = timeslice

        fieldlist = data.PointData.keys() if fieldnames is None else fieldnames
        for field in fieldlist:
            if isinstance(data.PointData[field], dsa.VTKNoneArray):
                continue
            elif isinstance(data.PointData[field], dsa.VTKArray):
                arry = numpy.array(data.PointData[field]).squeeze()
            else:
                arry = numpy.concatenate([numpy.array(x) for x in data.PointData[field].GetArrays() if not isinstance(x,dsa.VTKNoneArray)]).squeeze()

            # Array shape length 1 - scalar.
            #					 2 - vector.
            #					 3 - tensor.
            # the dict holds their names.
            TypeIndex = len(arry.shape) - 1
            for indxiter in product(*([range(3)] * TypeIndex)):
                L = tuple([slice(None, None, None)] + list(indxiter))
                try:
                    curstep["%s%s" % (field, self._componentsNames[indxiter])] = arry[L]
                except ValueError:
                    logger.warning("Field %s is problematic... ommiting" % field)

        curstep = curstep.assign(x=curstep.x.round(7), y=curstep.y.round(7), z=curstep.z.round(7), time=curstep.time.round(7))
        curstep = curstep.set_index(['time', 'x', 'y', 'z']).to_xarray() if regularMesh else curstep

        return curstep

    def writeCase(self, filtersDict, regularMesh, timeList=None, fieldnames=None, tsBlockNum=50, overwrite=False):
        """
                Writes the requested fileters as parquet files.

                if overwrite=False then  find the time that is saved. Assume that all the filters have the same timelist.
        Parameters
        ----------
        filtersDict : dict : filterName -> output file name
                List of VTK filters to write as parquets and the corresponding file name. 

        timeList : None, list or float.
                * None :  Use all time series. of the solver
                * List : get the listed timesteps.
                * float : Get only timesteps larger than this number
        fieldnames : list
                List of the fields to write in the filters (i.e the variables).

        FileBaseName : str
                The basic name of the file. Used to distinguish several filters with the same name.
        tsBlockNum : int
                number of block

        overwrite : bool


        append : bool
                The

        Returns
        -------

        """
        logger = get_classMethod_logger(self,"writeNonRegularCase")
        logger.info(f"Starting writing to parquet filters {','.join(filtersDict.keys())}")

        slice_filext = "zarr" if regularMesh else "parquet"

        maxTime = -1  # take all
        if overwrite:
            logger.info(f"Removing the old results")
            for filterName,outputPath in filtersDict.items():
                logger.debug(f"The data for {filterName} : {outputPath}")
                if os.path.isfile(outputPath):
                    logger.debug(f"\tParquet file {outputPath} is a file. Removing it")
                    os.remove(outputPath)
                elif os.path.isdir(outputPath):
                    logger.debug(f"\tParquet file {outputPath} is a directory. Removing the tree")
                    shutil.rmtree(outputPath)

        logger.info(f"Making sure that the output directories exist")
        for filterName, outputFile in filtersDict.items():
            outputPath = os.path.dirname(outputFile)
            logger.debug(f"{filterName} for directory {outputPath}")
            if not os.path.isdir(outputPath):
                logger.debug(f"\t Does not exist. Creating {outputPath}")
                os.makedirs(outputPath)


        readTimesList = self.reader.TimestepValues if timeList is None else timeList
        logger.info(f"Getting timelist {readTimesList}")

        blockID = 0
        tempList = []

        append = False if overwrite else True
        for filtersData in tqdm.tqdm(self.readTimeSteps(datasourcenamedict=filtersDict,
                                              timelist=readTimesList,
                                              fieldnames=fieldnames,
                                              regularMesh=regularMesh)):

            tempList.append(filtersData)
            logger.debug(f"Current dataFrames in memory  {len(tempList)}")
            if len(tempList) == tsBlockNum:
                self.writeList(tempList,blockID,filtersDict,regularMesh,slice_filext)
                tempList=[]
                blockID += 1

        if len(tempList) > 0:
            self.writeList(tempList,blockID,filtersDict,regularMesh,slice_filext)

        logger.info("Repartitioning to 100MB per partition")
        for filterName,outputFile in filtersDict.items():
            logger.debug(f"Working on {filterName}")

            outputFilterName = filterName
            outputPath = os.path.dirname(outputFile)
            outputFileList = [x for x in glob.glob(os.path.join(outputPath, f"tmp_{outputFilterName}_*.{slice_filext}"))]

            logger.debug(f"Saving all data to {outputFile}: {outputFileList}")
            with (ProgressBar()):
                if regularMesh:
                        # input_data is the directory path here
                        lazy_ds = xarray.open_mfdataset(outputFileList, chunks='auto', engine="zarr")
                        if append and os.path.exists(outputFile):
                            old_data = xarray.open_mfdataset(outputFile, chunks='auto', engine="zarr")
                            lazy_ds = xarray.concat([lazy_ds,old_data],dim="time").sortby("time")

                        try:
                            lazy_ds.to_zarr(f"{outputFile}.final", mode='w')
                        except NotImplementedError:
                            # somethimes this works and sometimes the other. not clear yet when...
                            lazy_ds.chunk("auto").to_zarr(f"{outputFile}.final", mode='w')

                else:
                    newDataList = [dd.read_parquet(fileName) for fileName in outputFileList]
                    if append and os.path.exists(outputFile):
                        newDataList.append(dd.read_parquet(outputFile))

                    allData = dd.concat(newDataList).repartition(partition_size="100MB")\
                        .sort_values("time")\
                        .set_index("time")\
                        .to_parquet(f"{outputFile}.final")

            if os.path.exists(outputFile):
                if os.path.isfile(outputFile):
                    os.remove(outputFile)
                else:
                    shutil.rmtree(outputFile)

            os.rename(f"{outputFile}.final",outputFile)

            logger.debug(f"Removing the old tmp files. ")
            for fileTodelete in outputFileList:
                if os.path.isfile(fileTodelete):
                    os.remove(fileTodelete)
                else:
                    shutil.rmtree(fileTodelete)

    def writeList(self,theList,blockID,filtersDict,regularMesh,fileExt):
        logger = get_classMethod_logger(self,"writeList")
        filterList = [x for x in theList[0].keys()]
        for filterName in filterList:
            outputFilterName = filterName.replace(".","-")
            outputPath = os.path.dirname(filtersDict[filterName])
            outputFile = os.path.join(outputPath,f"tmp_{outputFilterName}_{blockID:06d}.{fileExt}")

            logger.info(f"\tWriting filter {filterName} in temporary file {outputFile} ")

            if regularMesh:
                ds_slice = xarray.concat([item[filterName] for item in theList], dim='time')
                ds_slice.to_zarr(outputFile,mode='w')
            else:
                block_data = pandas.concat([item[filterName] for item in theList], ignore_index=True,sort=True)
                data = dd.from_pandas(block_data, npartitions=1)
                data.sort_values("time").to_parquet(outputFile)


    ############################################################################################
    ####        Depracated
    ############################################################################################

    @deprecated(reason="Use writeRegularCase instead")
    def write_netcdf(self, datasourcenamelist, timeList=None, fieldnames=None, tsBlockNum=50, overwrite=False,
                     append=False):
        self.writeRegularCase(datasourcenamelist, timeList, fieldnames, tsBlockNum, overwrite,append)


    @deprecated(reason="Old Name, use readTimeSteps with regularMesh=False")
    def to_pandas(self, datasourcenamelist, timelist=None, fieldnames=None):
        return self.readTimeSteps(datasourcenamelist, timelist, fieldnames, regtularMesh=False)

    @deprecated(reason="Old Name, use readTimeSteps with regularMesh=True")
    def to_xarray(self, datasourcenamelist, timelist=None, fieldnames=None):
        return self.readTimeSteps(datasourcenamelist, timelist, fieldnames, regtularMesh=True)

    @deprecated(reason="Old Name, use readTimeSteps with regularMesh=False")
    def to_dataFrame(self, datasourcenamelist, timelist=None, fieldnames=None):
        return self.readTimeSteps(datasourcenamelist, timelist, fieldnames, regtularMesh=False)

    @deprecated(reason="Old Name, use readTimeSteps with regularMesh=True")
    def to_dataArray(self, datasourcenamelist, timelist=None, fieldnames=None):
        return self.readTimeSteps(datasourcenamelist, timelist, fieldnames, regtularMesh=True)
    @deprecated(reason="Use writeNonRegularCase instead")
    def write_parquet(self, datasourcenamelist, timeList=None, fieldnames=None, tsBlockNum=50, overwrite=False,
                      append=False, filterList=None):
        writeNonRegularCase(datasourcenamelist, timeList, fieldnames, tsBlockNum, overwrite, append, filterList)

    @deprecated(reason="Use writeNonRegularCase instead")
    def write_parquet(self, datasourcenamelist, timeList=None, fieldnames=None, tsBlockNum=50, overwrite=False,
                      append=False, filterList=None):
        writeNonRegularCase(datasourcenamelist, timeList, fieldnames, tsBlockNum, overwrite, append, filterList)


    #
    # def writeRegularCase(self, filtersDict, timeList=None, fieldnames=None, tsBlockNum=50, overwrite=False,append=False):
    #     """
    #         Writes a list of datasources (vtk filters) to netcdf (with xarray).
    #         The grid data **must** be regular!!!.
    #         Todo: add a an option for regularization function.
    #
    #     Parameters
    #     ----------
    #
    #     readername: str
    #             The name of the reader to use.
    #     filtersDict:
    #             The name of the datasources to write.,
    #     outfile: str
    #             the directory to write the files.
    #     timeList: list
    #             the times to write
    #     fieldnames: list
    #             the fields to write
    #     tsBlockNum: int
    #             the number of
    #
    #     Returns
    #     -------
    #
    #     None
    #
    #     """
    #
    #     def checkIfExist(self,dataChunk,blockID,fileDirectory):
    #         filterList = [k for k in dataChunk.keys()]
    #         blockfrmt = ('{:0%dd}' % blockDig).format(blockID)
    #         for filtername in filterList:
    #             curfilename = os.path.join(fileDirectory, "%s_%s.nc" % (filtername, blockfrmt))
    #             if os.path.exists(curfilename):
    #                 if not overwrite:
    #                     raise Exception('NOTE: "%s" is alredy exists and will be not overwitten' % curfilename)
    #
    #     timeList = self.reader.TimestepValues if timeList is None else numpy.atleast_1d(timeList)
    #     os.makedirs(self.netcdfdir,exist_ok=True)
    #
    #     blockDig = max(5, numpy.ceil(numpy.log10(len(timeList))) + 1)
    #     blockID = 0
    #
    #     L = []
    #     checkExist=True
    #
    #     for xray in self.readTimeSteps(datasourcenamelist=filtersDict, timelist=timeList, fieldnames=fieldnames,regularMesh=True):
    #
    #         if checkExist:
    #             checkExist =False
    #             checkIfExist(xray, blockID, self.netcdfdir)
    #
    #         L.append(xray)
    #         if len(L) == tsBlockNum:
    #             if isinstance(L[0],dict):
    #                 filterList = [k for k in L[0].keys()]
    #                 for filtername in filterList:
    #                     writeList([item[filtername] for item in L], blockID, blockDig, overwrite, self.netcdfdir)
    #
    #             else:
    #                 writeList(L, blockID, blockDig, overwrite, self.netcdfdir)
    #             L = []
    #             blockID += 1
    #             checkExist = False
    #
    #     if len(L)>0:
    #         checkIfExist(xray, blockID, self.netcdfdir)
    #         if isinstance(L[0],dict):
    #             filterList = [k for k in L[0].keys()]
    #             for filtername in filterList:
    #                 writeList([item[filtername] for item in L], blockID, blockDig, overwrite, self.netcdfdir)
    #         else:
    #             writeList(L,blockID,blockDig,self.netcdfdir)
    #
