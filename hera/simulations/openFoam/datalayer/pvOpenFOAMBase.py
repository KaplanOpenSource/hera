from ....utils import loggedObject
from itertools import product
import pandas
import numpy
import dask.dataframe as dd

import os
import vtk.numpy_interface.dataset_adapter as dsa
import xarray

#### import the simple module from the paraview
import paraview.simple as pvsimple
from paraview import servermanager

#### disable automatic camera reset on 'Show'
pvsimple._DisableFirstRenderCameraReset()


class paraviewOpenFOAM(loggedObject):
    """
        A class to extract openFOAM file format
        using VTK filters and write as parquet or netcdf files.
    """

    _componentsNames = None  # names of components for reading.

    _hdfdir        = None    # path to save the hdf.
    _netcdfdir     = None    # path to save the netcdf.
    _parquetdir    = None


    _readerName = None
    _father = None

    @property
    def reader(self):
        return self._reader

    @property
    def readerName(self):
        return self._readerName


    @property
    def outputPath(self):
        return self._outputPath

    @property
    def hdfdir(self):
        return self._hdfdir

    @hdfdir.setter
    def hdfdir(self, hdfdir):
        self._hdfdir = hdfdir

    @property
    def netcdfdir(self):
        return self._netcdfdir
    @property
    def possibleRegions(self):
        return self._possibleRegions

    @netcdfdir.setter
    def netcdfdir(self, netcdfdir):
        self._netcdfdir = netcdfdir

    @property
    def parquetdir(self):
        return self._parquetdir

    @parquetdir.setter
    def parquetdir(self, parquetdir):
        self._parquetdir = parquetdir

    def __init__(self, casePath,caseType='Decomposed Case', servername=None,fieldNames =None,name="mainreader"):
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
        super().__init__()
        if servername is not None:
            pvsimple.Connect(servername)

        self._componentsNames = {}

        self.netcdfdir = "netcdf"
        self.hdfdir = "hdf"
        self.parquetdir = "parquet"


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


        self._ReadCase(readerName=name, casePath=casePath, CaseType=caseType, fieldNames=fieldNames)

    def _ReadCase(self, readerName, casePath, CaseType='Decomposed Case', fieldNames=None):
        """
            Constructs a reader and register it in the vtk pipeline.

            Handles either parallel or single format.

        Parameters
        -----------

        readerName:
                the name of the reader.
        casePath:
                a full path to the case directory.
        CaseType: str
                Either 'Decomposed Case' for parallel cases or 'Reconstructed Case'
                for single processor cases.
        fieldnames: list of str
                List of field names to load.
                if None, read all the fields.
        :return:
                the reader
        """
        self._readerName  = readerName
        self._reader = pvsimple.OpenFOAMReader(FileName="%s/tmp.foam" % casePath, CaseType=CaseType, guiName=readerName)
        self.reader.MeshRegions.SelectAll()
        self._possibleRegions = list(self._reader.MeshRegions)
        self._reader.MeshRegions = ['internalMesh']
        if fieldNames is not None:
            self._reader.CellArrays = fieldNames

        self._reader.UpdatePipeline()

    def to_pandas(self, datasourcenamelist, timelist=None, fieldnames=None):
        return self.readTimeSteps(datasourcenamelist, timelist, fieldnames, xarray=False)

    def to_xarray(self, datasourcenamelist, timelist=None, fieldnames=None):
        return self.readTimeSteps(datasourcenamelist, timelist, fieldnames, xarray=True)

    def to_dataFrame(self, datasourcenamelist, timelist=None, fieldnames=None):
        return self.readTimeSteps(datasourcenamelist, timelist, fieldnames, xarray=False)

    def to_dataArray(self, datasourcenamelist, timelist=None, fieldnames=None):
        return self.readTimeSteps(datasourcenamelist, timelist, fieldnames, xarray=True)


    def readTimeSteps(self, datasourcenamelist, timelist=None, fieldnames=None, xarray=False):
        """
            reads a list of datasource lists to a dictionary

        Parameters
        ----------

        readername: VTK filter, str
                The reader filter (or its name)

        datasourcenamelist: list
                A list of names of filters to get.

        timelist: list
                The list of times to read.
        fieldnames:
                The list of fields to write.
        xarray
                convert pandas results to xarray (works only for regular grids).

        Return
        ------

        For each time step.
                    A map datasourcename -> pandas
        """
        datasourcenamelist = numpy.atleast_1d(datasourcenamelist)

        timelist = self.reader.TimestepValues if timelist is None else numpy.atleast_1d(timelist)
        for timeslice in timelist:
            # read the timestep.
            self.logger.info("\r Reading time step %s" % timeslice)

            ret = {}
            for datasourcename in datasourcenamelist:
                datasource = pvsimple.FindSource(datasourcename)
                self.logger.execution(f"Reading source {datasourcename}")
                rt = self._readTimeStep(datasource, timeslice, fieldnames, xarray)
                if rt is not None:
                    ret[datasourcename] = rt
            yield ret

    def _readTimeStep(self, datasource, timeslice, fieldnames=None, xarray=False):
        # read the timestep.
        datasource.UpdatePipeline(timeslice)
        rawData = servermanager.Fetch(datasource)
        data = dsa.WrapDataObject(rawData)

        if isinstance(data.Points, dsa.VTKNoneArray):
            self.logger.execution("No data exists for filter... return with None")
            return None
        elif isinstance(data.Points, dsa.VTKArray):
            points = numpy.array(data.Points).squeeze()
        else:
            points = numpy.concatenate([numpy.array(x) for x in data.Points.GetArrays()]).squeeze()

        self.logger.execution(f"Filter has {points.shape[0]} points. Building basic dataFrame. ")
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
                    self.logger.warning("Field %s is problematic... ommiting" % field)

        curstep = curstep.set_index(['time', 'x', 'y', 'z']).to_xarray() if xarray else curstep

        return curstep

    def write_netcdf(self, datasourcenamelist, timeList=None, fieldnames=None, tsBlockNum=50, overwrite=False,append=False):
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

        def writeList(theList, blockID, blockDig,overWrite,fileDirectory):
            data = xarray.concat(theList, dim="time")
            blockfrmt = ('{:0%dd}' % blockDig).format(blockID)
            curfilename = os.path.join(fileDirectory, "%s_%s.nc" % (filtername, blockfrmt))
            if os.path.exists(curfilename):
                if not overWrite:
                    raise Exception ('NOTE: "%s" is alredy exists and will be not overwitten' % curfilename)

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

        for xray in self.to_xarray(datasourcenamelist=datasourcenamelist, timelist=timeList, fieldnames=fieldnames):

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

    def write_parquet(self, datasourcenamelist, timeList=None, fieldnames=None, tsBlockNum=50, overwrite=False, append=False, filterList=None):
        """
                Writes the requested fileters as parquet files.

                if overwrite=False then  find the time that is saved. Assume that all the filters have the same timelist.
        Parameters
        ----------
        datasourcenamelist : list
                List of VTK filters to write as parquets

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
        self.logger.info(f"Starting writing to parquet filters {','.join(numpy.atleast_1d(datasourcenamelist))}")

        def writeList(theList,filePath,append,overwrite):
            filterList = [x for x in theList[0].keys()]
            for filtername in filterList:
                outfile = os.path.join(filePath,f"{filtername}.parquet")
                self.logger.info("\tWriting filter %s in file %s" % (filtername, outfile))

                if os.path.exists(outfile) and not append and not overwrite:
                    self.logger.error(f"The output for {filtername} ({outfile}) exists. use --append to append or --overwrite to overwrite. Ignoring")
                    continue
                block_data = pandas.concat([pandas.DataFrame(item[filtername]) for item in theList], ignore_index=True,sort=True)
                data = dd.from_pandas(block_data, npartitions=1)
                data.set_index("time").to_parquet(outfile,append=append,overwrite=overwrite)


        if not os.path.isdir(self.parquetdir):
            self.logger.debug(f"Creating output directory {self.parquetdir}")
            os.makedirs(self.parquetdir)

        maxTime = -1  # take all time steps.
        if not overwrite:
            # find the time that is saved. Assume that all the filters have the same timelist.
            # find the first.
            self.logger.execution("Appending to existing data, find the maximal time. ")
            for filtername in datasourcenamelist:
                prqtFile = os.path.join(self.parquetdir, f"{filtername}.parquet")
                if not os.path.exists(prqtFile):
                    continue
                maxTime = dd.read_parquet(prqtFile).index.max().compute()
                break
            self.logger.debug(f"The maximal time found is {maxTime}. Skipping all the timesteps beofre that.")
        else:
            self.logger.execution("Overwriting existing data... Deleting current parquets, if they exist")
            import shutil
            for filtername in datasourcenamelist:
                prqtFile = os.path.join(self.parquetdir, f"{filtername}.parquet")
                if os.path.isfile(prqtFile):
                    self.logger.debug(f"Parquet file {prqtFile} is a file. Removing it")
                    os.remove(prqtFile)
                elif os.path.isdir(prqtFile):
                    self.logger.debug(f"Parquet file {prqtFile} is a directory. Removing the tree")
                    shutil.rmtree(prqtFile)

        timeList = self.reader.TimestepValues if timeList is None else numpy.atleast_1d(timeList)
        self.logger.debug(f"Filtering out the timesteps. Original list includes: {timeList}")
        timeList = [x for x in timeList if x > maxTime]
        self.logger.debug(f"Filtering out the timesteps. After filteration: {timeList}")
        blockID = 0
        L = []
        for pnds in self.to_pandas(datasourcenamelist=datasourcenamelist, timelist=timeList, fieldnames=fieldnames):

            L.append(pnds)
            self.logger.debug(f"Current dataFrames in memory  {len(L)}")
            if len(L) == tsBlockNum:
                writeList(L, self.parquetdir,append=append,overwrite=overwrite)

                # From the second iteration, we must append to the newly created file.
                overwrite = False
                append = True
                L=[]
                blockID += 1
        if len(L) > 0:
            writeList(L, self.parquetdir,append=append,overwrite=overwrite)

        filterList = [x for x in L[0].keys()]
        for filtername in filterList:
            outfile = os.path.join(self.parquetdir,f"{filtername}.parquet")
            if os.path.exists(outfile):
                self.logger.execution("Repartitioning to 100MB per partition")
                dd.read_parquet(outfile).repartition(partition_size="100MB").reset_index().sort_values("time").set_index(
                    "time").to_parquet(outfile)
