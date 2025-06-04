from itertools import product
import pandas
import numpy
import dask.dataframe as dd

import os
import xarray

from hera.simulations.openFoam import CASETYPE_DECOMPOSED,CASETYPE_RECONSTRUCTED
from hera import get_classMethod_logger
from deprecated import deprecated

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
    datalayer = None

    def __init__(self, casePath,datalayer,caseType=CASETYPE_DECOMPOSED, servername=None,name="mainreader"):
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

        self.datalayer = datalayer

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
                rt = self._readTimeStep(datasource, timeslice, fieldnames, regtularMesh)
                if rt is not None:
                    ret[filterName] = rt
            yield ret

    def _readTimeStep(self, datasource, timeslice, fieldnames=None, xarray=False):
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

        curstep = curstep.set_index(['time', 'x', 'y', 'z']).to_xarray() if xarray else curstep

        return curstep

    def writeRegularCase(self, datasourcenamelist, timeList=None, fieldnames=None, tsBlockNum=50, overwrite=False,append=False):
        """
            Writes a list of datasources (vtk filters) to netcdf (with xarray).
            The grid data **must** be regular!!!.
            Todo: add a an option for regularization function.

        Parameters
        ----------

        readername: str
                The name of the reader to use.
        datasourcenamelist:
                The name of the datasources to write.,
        outfile: str
                the directory to write the files.
        timeList: list
                the times to write
        fieldnames: list
                the fields to write
        tsBlockNum: int
                the number of

        Returns
        -------

        None

        """
        def writeList(theList, blockID, blockDig, overwrite, fileDirectory):
            data = xarray.concat(theList, dim="time")
            blockfrmt = ('{:0%dd}' % blockDig).format(blockID)
            curfilename = os.path.join(fileDirectory, "%s_%s.nc" % (filtername, blockfrmt))
            if os.path.exists(curfilename):
                if not overwrite:
                    raise FileExistsError(f'NOTE: {curfilename} exists and will be not overwitten')

            data.to_netcdf(curfilename)
            blockID += 1

        def checkIfExist(self,dataChunk,blockID,fileDirectory):
            filterList = [k for k in dataChunk.keys()]
            blockfrmt = ('{:0%dd}' % blockDig).format(blockID)
            for filtername in filterList:
                curfilename = os.path.join(fileDirectory, "%s_%s.nc" % (filtername, blockfrmt))
                if os.path.exists(curfilename):
                    if not overwrite:
                        raise Exception('NOTE: "%s" is alredy exists and will be not overwitten' % curfilename)

        timeList = self.reader.TimestepValues if timeList is None else numpy.atleast_1d(timeList)
        os.makedirs(self.netcdfdir,exist_ok=True)

        blockDig = max(5, numpy.ceil(numpy.log10(len(timeList))) + 1)
        blockID = 0

        L = []
        checkExist=True

        for xray in self.readTimeSteps(datasourcenamelist=datasourcenamelist, timelist=timeList, fieldnames=fieldnames,regularMesh=True):

            if checkExist:
                checkExist =False
                checkIfExist(xray, blockID, self.netcdfdir)

            L.append(xray)
            if len(L) == tsBlockNum:
                if isinstance(L[0],dict):
                    filterList = [k for k in L[0].keys()]
                    for filtername in filterList:
                        writeList([item[filtername] for item in L], blockID, blockDig, overwrite, self.netcdfdir)

                else:
                    writeList(L, blockID, blockDig, overwrite, self.netcdfdir)
                L = []
                blockID += 1
                checkExist = False

        if len(L)>0:
            checkIfExist(xray, blockID, self.netcdfdir)
            if isinstance(L[0],dict):
                filterList = [k for k in L[0].keys()]
                for filtername in filterList:
                    writeList([item[filtername] for item in L], blockID, blockDig, overwrite, self.netcdfdir)
            else:
                writeList(L,blockID,blockDig,self.netcdfdir)


    def writeNonRegularCase(self, filtersDict, timeList=None, fieldnames=None, tsBlockNum=50, overwrite=False):
        """
                Writes the requested fileters as parquet files.

                if overwrite=False then  find the time that is saved. Assume that all the filters have the same timelist.
        Parameters
        ----------
        filtersDict : dict : filterName -> output file name
                List of VTK filters to write as parquets and the corresponding file name. 

        timeList : list
                List of timesteps to write. optional, if None, the use all time series.

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
        logger.info(f"Starting writing to parquet filters {','.join(numpy.atleast_1d(datasourcenamelist))}")

        def writeList(theList,blockID,filtersDict):
            filterList = [x for x in theList[0].keys()]
            for filterName in filterList:
                outputFilterName = filterName.replace(".","-")
                outputPath = os.path.dirname(filtersDict[filterName])
                outputFile = os.path.join(outputPath,f"tmp_{outputFilterName}_{blockID}.parquet")

                logger.info("\tWriting filter %s in file %s" % (filterName, outputFile))
                block_data = pandas.concat([pandas.DataFrame(item[filterName]['data']) for item in theList], ignore_index=True,sort=True)
                data = dd.from_pandas(block_data, npartitions=1)
                data.sort_values("Time").to_parquet(outputFile)

        if not os.path.isdir(outputPath):
            logger.debug(f"Creating output directory {outputPath}")
            os.makedirs(outputPath)

        maxTime = -1  # take all
        if overwrite:
            logger.info(f"Removing the old results")
            for filterName,outputPath in filtersDict:
                if os.path.isfile(outputPath):
                    logger.debug(f"Parquet file {outputPath} is a file. Removing it")
                    os.remove(outputPath)
                elif os.path.isdir(outputPath):
                    logger.debug(f"Parquet file {outputPath} is a directory. Removing the tree")
                    shutil.rmtree(outputPath)

        logger.debug(f"Getting the time list.")
        readTimesList = self.reader.TimestepValues
        if isinstance(timeList,list):
            logger.debug(f"Got a list. Then getting it exactly: {timeList}")
            readTimesList = timeList
        elif isinstance(timeList,float) or isinstance(timeList,int):
            logger.debug(f"Got a number. Then getting only timesteps greater than it: {timeList}")
            readTimesList = [x for x in readTimesList if x > timeList]
            logger.debug(f"Getting {readTimesList}")

        import pdb
        pdb.set_trace()

        blockID = 0
        tempList = []
        append = False if overwrite else True
        for filtersData in self.readTimeSteps(datasourcenamelist=filtersDict,
                                              timelist=timeList,
                                              fieldnames=fieldnames,
                                              outputPath=outputPath,
                                              regularMesh=False):

            tempList.append(filtersData)
            logger.debug(f"Current dataFrames in memory  {len(L)}")
            if len(tempList) == tsBlockNum:
                writeList(tempList,blockID,filtersDict)
                tempList=[]
                blockID += 1

        if len(tempList) > 0:
            writeList(tempList,blockID,filtersDict)

        logger.info("Repartitioning to 100MB per partition")
        for filterItem in datasourcenamelist:
            logger.debug(f"Working on {filterItem['filter']}")

            outputFilterName = filterName.replace(".", "-")
            outputPath = os.path.dirname(filtersDict[filterName])
            outputFile = os.path.join(outputPath, f"tmp_{outputFilterName}_*.parquet")

            logger.debug(f"Saving  partitions from the files {outputFile}")
            dd.read_parquet(outputFile).repartition(partition_size="100MB")\
                .sort_values("time")\
                .set_index("time")\
                .to_parquet(filtersDict[filterName])

            logger.debug(f"Removing the old tmp files. ")
            for fileTodelete in glob.glob(outputFile):
                os.remove(fileTodelete)





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
