
from itertools import product
import json
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
    import vtkmodules.numpy_interface.dataset_adapter as dsa
    from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet
    import paraview.simple as pvsimple
    from paraview import servermanager
    #### disable automatic camera reset on 'Show'
    from paraview.servermanager import Proxy
    from paraview.servermanager import ProxyProperty

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
        reader.MeshRegions = possibleRegions
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
        def debug_proxy_data(proxy: Proxy, nested=True, depth = 0):
            """Use to make sure all filter data is set correctly(right before execution)

            Parameters
            ----------
                proxy (Proxy): every paraview object with data is a proxy and so is the paraview filter
                nested (bool, optional): should the debug log go into each sub object(hence creating a tree). Defaults to True.
                depth (int, optional): the starting depth in the print(mostly internal use). Defaults to 0.
            """
            for prop_name in proxy.ListProperties():
                prop = proxy.GetProperty(prop_name)
                tree_sign = '├' if prop_name != proxy.ListProperties()[-1] else '╰'
                print_prefix = f"{'.│'*depth} {f'{tree_sign}─' if depth != 0 else ''}"
                if isinstance(prop,ProxyProperty):
                    prop_data = prop.GetData()
                    if isinstance(prop_data,Proxy):
                        logger.debug(f"{print_prefix}{prop_name} -> Proxy of type {prop} with {len(proxy.ListProperties())} properties")
                        if nested:
                            depth += 1
                            debug_proxy_data(prop_data, depth=depth)
                            depth -= 1
                    else:
                        logger.debug(f"{print_prefix}{prop_name} -> {prop_data}")
                else:
                    logger.debug(f"{print_prefix}{prop_name} -> {prop}")
        
        logger = get_classMethod_logger(self, "readTimeSteps")
        for timeslice in timelist:
            # read the timestep.
            logger.info("\r Reading time step %s" % timeslice)
            ret = {}
            for filterName,_ in datasourcenamedict.items():
                datasource = pvsimple.FindSource(filterName)
                logger.debug(f"Reading source {filterName}")
                assert isinstance(datasource, Proxy)
                logger.debug("Calculating filter with the following parameters:")
                debug_proxy_data(datasource)
                rt = self._readTimeStep(datasource, timeslice, fieldnames, regularMesh)
                if rt is not None:
                    ret[filterName] = rt
            yield ret

    def _parsePointSet(self, pointSet, timeslice, fieldnames, regularMesh):
        """Parse a VTK point set into a pandas DataFrame or xarray Dataset."""
        logger = get_classMethod_logger(self, "_parsePointSet")
        if isinstance(pointSet.Points, dsa.VTKNoneArray):
            logger.debug("No data exists for filter... return with None")
            return None
        elif isinstance(pointSet.Points, dsa.VTKArray):
            points = numpy.array(pointSet.Points).squeeze()
        else:
            points = numpy.concatenate([numpy.array(x) for x in pointSet.Points.GetArrays()]).squeeze()

        logger.debug(f"Filter has {points.shape[0]} points. Building basic dataFrame. ")

        # For filters like integrate, the len points.shape==1, and therefore [:,0] will break the code.
        # However, it means that there are no points, and therefore we will skip it.
        if len(points.shape)==2:
            # create index
            curstep = pandas.DataFrame()
            curstep['x'] = points[:, 0]
            curstep['y'] = points[:, 1]
            curstep['z'] = points[:, 2]
            curstep['time'] = timeslice
            curstep = curstep.assign(x=curstep.x.round(7), y=curstep.y.round(7), z=curstep.z.round(7),
                                     time=curstep.time.round(7))
        else:
            curstep = pandas.DataFrame([dict(time=timeslice)])

        fieldlist = pointSet.PointData.keys() if fieldnames is None else fieldnames
        for field in fieldlist:
            if isinstance(pointSet.PointData[field], dsa.VTKNoneArray):
                continue
            elif isinstance(pointSet.PointData[field], dsa.VTKArray):
                arry = numpy.array(pointSet.PointData[field]).squeeze()
            else:
                arry = numpy.concatenate([numpy.array(x) for x in pointSet.PointData[field].GetArrays() if not isinstance(x,dsa.VTKNoneArray)]).squeeze()

            # Array shape length 0 - Integration - a scalar that is not a field.
            #                    1 - scalar.
            #					 2 - vector.
            #					 3 - tensor.
            # the dict holds their names.
            TypeIndex = len(arry.shape) - 1
            if TypeIndex >=0:
                for indxiter in product(*([range(3)] * TypeIndex)):
                    L = tuple([slice(None, None, None)] + list(indxiter))
                    try:
                        curstep["%s%s" % (field, self._componentsNames[indxiter])] = arry[L]
                    except ValueError:
                        logger.warning("Field %s is problematic... ommiting" % field)
            else:
                # For integrate variables filter, there is only scalar.
                val = numpy.atleast_1d(arry)[0] # 0 dim did some troubles.
                curstep[field] = val

        if 'x' in curstep:
            curstep = curstep.set_index(['time', 'x', 'y', 'z']).to_xarray() if regularMesh else curstep
        else:
            curstep = curstep.set_index(['time']).to_xarray() if regularMesh else curstep
        
        return curstep

    def _parseMultiBlockDataSet(self, multiBlock, timeslice, fieldnames, regularMesh):
        """Parse a VTK multi-block dataset, preserving block names."""
        num_blocks = multiBlock.GetNumberOfBlocks()
        final_data = []
        if num_blocks == 0:
            return final_data
        
        # name doesn't matter and list doesn't have to be used
        elif num_blocks == 1:
            return self._parseVTKData(multiBlock.GetBlock(0), timeslice, fieldnames, regularMesh)

        for i in range(num_blocks):
            # Check if block exists and get its name
            block = multiBlock.GetBlock(i)
            if block:
                blockName = self._getBlockName(multiBlock, i)
                parsed_block = self._parseVTKData(block, timeslice, fieldnames, regularMesh)
                if isinstance(parsed_block, list):
                    # we assumee that list elements already have block names assigned to them or they are irrelevent
                    final_data.extend(parsed_block)
                else:
                    parsed_block = self._assignBlockName(regularMesh, blockName, parsed_block)
                    final_data.append(parsed_block)
            else:
                raise RuntimeError(f"Failed to get block number {i}, raw VTK object:\n{multiBlock}")
        return final_data

    def _assignBlockName(self, regularMesh, blockName, parsed_block):
        """Assign a block name column/variable to a parsed data block."""
        if regularMesh:
            assert isinstance(parsed_block, xarray.Dataset) or isinstance(parsed_block, xarray.DataArray), f"while processing regular mesh got produced {type(parsed_block)} instead of xarray"
            blockNameMap=numpy.full(tuple(parsed_block.sizes.values()), fill_value=blockName)
            if isinstance(parsed_block, xarray.DataArray):
                parsed_block = parsed_block.to_dataset()
            parsed_block['blockName'] = (tuple(parsed_block.sizes.keys()), blockNameMap)
        else:
            assert isinstance(parsed_block, pandas.DataFrame), f"while processing non regular mesh got {type(parsed_block)} instead of pandas.DataFrame"
            parsed_block['blockName']=blockName
        return parsed_block

    def _getBlockName(self, multiBlock, i):
        """Retrieve the metadata name of a block by its index."""
        blockName = numpy.nan
        if multiBlock.HasMetaData(i):
            blockMetaData = multiBlock.GetMetaData(i)
            blockMetaDataNameKey = multiBlock.NAME()
            if blockMetaData.Has(blockMetaDataNameKey):
                blockName = blockMetaData.Get(blockMetaDataNameKey)
        return blockName

    def _parseVTKData(self, data, timeslice, fieldnames, regularMesh):
        """Recursively going over the vtk data and generating a final dataframe trying to keep the block name

        :param data: vtk data,
        :type data: Any
        :param timeslice: timestep at which the data is calculated for
        :type timeslice: int
        :param fieldnames: names of the fields of interest
        :type fieldnames: List[str]
        :param regularMesh: decides whether to save the data with regular format(xarray) or generic table(parquet with pandas), notice all data gets processed as pandas initially
        :type regularMesh: bool
        :raises NotImplementedError: Some vtk data's parsing might not be implemented in which case raises
        :return: final list/sole dataframe or xarray containing VTK data
        :rtype: List[pandas.Dataframe] | List[xarray] | Dataframe | xarray
        """
        logger = get_classMethod_logger(self, "_parseVTKWrappedData")
        logger.debug(f"Parsing VTK raw data of type {type(data)}")
        
        # composite data needs to be parsed as ktk objects to retain the block names
        data = dsa.WrapDataObject(data) if not isinstance(data, vtkMultiBlockDataSet) else data
        # make sure we bring back composite data to their vtkMultiBlockDataSet structure 
        data = data.VTKObject if isinstance(data, dsa.CompositeDataSet) else data
        logger.debug(f"Normalized raw data to type {type(data)}")
        
        # instead of using the wrapped data we iterate over it to get the name of the blocks
        # required when block data has similar structure and the name is required to associate the values
        if isinstance(data, vtkMultiBlockDataSet):
            logger.debug("Parsing composite data, combining all blocks with their names...")
            return self._parseMultiBlockDataSet(data, timeslice, fieldnames, regularMesh)
        elif isinstance(data, dsa.PointSet) or isinstance(data, dsa.PolyData) or isinstance(data, dsa.UnstructuredGrid):
            logger.info("Parsing point data")
            #  the cases we encountered with dsa.PolyData just needed the PointSet data
            # hence JUST FOR NOW they all get funneled to _pointSetParser
            return self._parsePointSet(data, timeslice, fieldnames, regularMesh)
        elif isinstance(data, dsa.Table):
            logger.info("Parsing table data")
            return self._parseTable(data, timeslice, regularMesh)
        else:
            raise NotImplementedError(f"Currently reading timestamp doesn't implement parsing of {type(data)}")

    def _parseTable(self, table, timeslice, regularMesh):
        """Parse a VTK table into a pandas DataFrame or xarray Dataset."""
        df = pandas.DataFrame()
        for columnName in table.RowData.keys():
            df[columnName]=table.RowData[columnName]
        df['time'] = timeslice
        ret =  df.to_xarray() if regularMesh else df
        return ret

    def _readTimeStep(self, datasource, timeslice, fieldnames=None, regularMesh=False):
        """Fetch and parse data from a VTK source at a given time slice."""
        # read the timestep.
        logger = get_classMethod_logger(self, "_readTimeStep")

        datasource.UpdatePipeline(timeslice)
        # parview works as a server always so we must request the data from the datasource using servermanager.Fetch
        # pre and post processing algorithms may be applied by secifying arg1(pre) and arg2(post) and client port with idx
        # it should return a vtk object(rawData)
        rawData = servermanager.Fetch(datasource)
        
        # dsa is a utility that "provides classes that allow Numpy-type access to VTK datasets and arrays"
        # https://docs.vtk.org/en/latest/api/python/vtkmodules/vtkmodules.numpy_interface.dataset_adapter.html
        # dsa.WraoDataObject returns "a Numpy friendly wrapper of a vtkDataObject" which is dsa.DataObject
        # 
        # There are 8 classes that are polymorphic forms of dsa.DataObject and then multiple polymorphic forms of them(10 total) listed here:
        #    dsa.DataObject(FieldData), dsa.Table(RowData), dsa.HyperGridTree(CellData(AMR)),
        #    dsa.DataSet(base form for dsa.PointSet(Points, base for either dsa.PolyData(Polygons) or dsa.UnstructuredGrid(Cells, CellsType, CellsLocations)),
        #    dsa.CompositeDataSet(composite here means seperated dsa.DataSet's that need to be iterated over),
        #    dsa.Graph(vertices nd edeges), dsa.Molecule(Atoms and Bonds)

        return self._parseVTKData(rawData, timeslice, fieldnames, regularMesh)

    def writeCase(self, filtersDict, regularMesh, timeList=None, fieldnames=None, tsBlockNum=50, overwrite=False, latestTimestamp=False):
        """
        Write VTK filter results to parquet (unstructured) or zarr (regular) files.

        Algorithm overview:
          1. Prepare the filesystem (clean old outputs, create directories).
          2. Resolve which timesteps to process.
          3. Stream timesteps in fixed-size blocks to temporary files (controls memory).
          4. Merge all temporary files into the final output per filter.

        Parameters
        ----------
        filtersDict : dict
            Mapping of filter name -> output file path.
        regularMesh : bool
            If True, write zarr (xarray). If False, write parquet (dask DataFrame).
        timeList : list or None
            Timesteps to process. None = all available from the reader.
        fieldnames : list, optional
            VTK field names to extract.
        tsBlockNum : int
            Number of timesteps to accumulate before flushing a temporary block.
        overwrite : bool
            If True, remove existing output files before writing.
        latestTimestamp : bool
            If True, only process the latest available timestep.
        """
        logger = get_classMethod_logger(self, "writeNonRegularCase")
        logger.info(f"Starting writing to parquet filters {','.join(filtersDict.keys())}")

        # Choose file extension based on mesh type: zarr for regular, parquet otherwise.
        slice_filext = "zarr" if regularMesh else "parquet"
        # In non-overwrite mode, new data is appended to any existing output.
        append = not overwrite

        # Step 1: Prepare the filesystem for writing.
        if overwrite:
            self._removeOldOutputs(filtersDict)
        self._ensureOutputDirs(filtersDict)

        # Step 2: Determine which timesteps to read from the simulation.
        readTimesList = self._resolveTimeList(timeList, latestTimestamp)

        # Step 3: Stream timesteps in blocks, writing each full block to a temp file.
        self._writeTimeStepBlocks(filtersDict, readTimesList, fieldnames,
                                  regularMesh, slice_filext, tsBlockNum)

        # Step 4: For each filter, merge its temp files into a single final output.
        logger.info("Repartitioning to 100MB per partition")
        for filterName, outputFile in filtersDict.items():
            tmpFiles = self._collectTmpFiles(filterName, outputFile, slice_filext)
            self._mergeToFinalOutput(outputFile, tmpFiles, regularMesh, append)
            self._atomicReplace(outputFile)
            self._cleanupTmpFiles(tmpFiles)

    # ------------------------------------------------------------------
    # Private helpers for writeCase
    # ------------------------------------------------------------------

    def _removeOldOutputs(self, filtersDict):
        """Remove existing output files or directories before overwriting."""
        logger = get_classMethod_logger(self, "_removeOldOutputs")
        logger.info("Removing the old results")
        for filterName, outputPath in filtersDict.items():
            logger.debug(f"The data for {filterName} : {outputPath}")
            if os.path.isfile(outputPath):
                logger.debug(f"\tParquet file {outputPath} is a file. Removing it")
                os.remove(outputPath)
            elif os.path.isdir(outputPath):
                logger.debug(f"\tParquet file {outputPath} is a directory. Removing the tree")
                shutil.rmtree(outputPath)

    def _ensureOutputDirs(self, filtersDict):
        """Create output directories for each filter if they do not already exist."""
        logger = get_classMethod_logger(self, "_ensureOutputDirs")
        logger.info("Making sure that the output directories exist")
        for filterName, outputFile in filtersDict.items():
            outputPath = os.path.dirname(outputFile)
            logger.debug(f"{filterName} for directory {outputPath}")
            if not os.path.isdir(outputPath):
                logger.debug(f"\t Does not exist. Creating {outputPath}")
                os.makedirs(outputPath)

    def _resolveTimeList(self, timeList, latestTimestamp):
        """Determine which timesteps to process from the reader or caller input.

        Returns the full reader timestep list when timeList is None,
        or trims to only the latest entry when latestTimestamp is True.
        """
        logger = get_classMethod_logger(self, "_resolveTimeList")
        readTimesList = self.reader.TimestepValues if timeList is None else timeList
        logger.info(f"Getting timelist {readTimesList}")
        if latestTimestamp and len(readTimesList) != 0:
            readTimesList = [readTimesList[-1]]
        return readTimesList

    def _writeTimeStepBlocks(self, filtersDict, readTimesList, fieldnames,
                             regularMesh, slice_filext, tsBlockNum):
        """Stream timesteps from the VTK pipeline, flushing to temp files in blocks.

        Accumulates up to tsBlockNum timesteps in memory, then writes them
        to a numbered temporary file via writeList. Any leftover timesteps
        that do not fill a complete block are flushed at the end.
        """
        logger = get_classMethod_logger(self, "_writeTimeStepBlocks")
        blockID = 0
        tempList = []

        for filtersData in tqdm.tqdm(self.readTimeSteps(
                datasourcenamedict=filtersDict, timelist=readTimesList,
                fieldnames=fieldnames, regularMesh=regularMesh)):
            tempList.append(filtersData)
            logger.debug(f"Current dataFrames in memory  {len(tempList)}")
            # Flush the block to disk once it reaches the target size.
            if len(tempList) == tsBlockNum:
                self.writeList(tempList, blockID, filtersDict, regularMesh, slice_filext)
                tempList = []
                blockID += 1

        # Flush any remaining timesteps that did not fill a complete block.
        if len(tempList) > 0:
            self.writeList(tempList, blockID, filtersDict, regularMesh, slice_filext)

    def _collectTmpFiles(self, filterName, outputFile, slice_filext):
        """Glob all numbered temporary block files produced for a given filter."""
        outputPath = os.path.dirname(outputFile)
        tmpPattern = f"tmp_{filterName.replace('.', '-')}_*.{slice_filext}"
        return glob.glob(os.path.join(outputPath, tmpPattern))

    def _mergeToFinalOutput(self, outputFile, tmpFiles, regularMesh, append):
        """Merge temporary block files into a single '.final' staging file.

        For regular meshes (zarr): opens all blocks as a lazy multi-file dataset,
        optionally concatenates with previously saved data, and writes to zarr.
        For unstructured meshes (parquet): concatenates all blocks with dask,
        repartitions to ~100 MB chunks indexed by time, and writes to parquet.
        """
        logger = get_classMethod_logger(self, "_mergeToFinalOutput")
        logger.debug(f"Saving all data to {outputFile}: {tmpFiles}")

        with ProgressBar():
            if regularMesh:
                self._mergeZarr(outputFile, tmpFiles, append)
            else:
                self._mergeParquet(outputFile, tmpFiles, append)

    def _mergeZarr(self, outputFile, tmpFiles, append):
        """Concatenate temporary zarr blocks, optionally appending old data."""
        lazy_ds = xarray.open_mfdataset(tmpFiles, chunks='auto', engine="zarr")
        # If appending, include previously saved data so nothing is lost.
        if append and os.path.exists(outputFile):
            old_data = xarray.open_mfdataset(outputFile, chunks='auto', engine="zarr")
            lazy_ds = xarray.concat([lazy_ds, old_data], dim="time").sortby("time")
        try:
            lazy_ds.to_zarr(f"{outputFile}.final", mode='w')
        except NotImplementedError:
            # somethimes this works and sometimes the other. not clear yet when...
            lazy_ds.chunk("auto").to_zarr(f"{outputFile}.final", mode='w')

    def _mergeParquet(self, outputFile, tmpFiles, append):
        """Concatenate temporary parquet blocks, repartition to ~100 MB, index by time."""
        newDataList = [dd.read_parquet(fileName) for fileName in tmpFiles]
        # If appending, include the previously saved parquet data.
        if append and os.path.exists(outputFile):
            newDataList.append(dd.read_parquet(outputFile))
        dd.concat(newDataList).repartition(partition_size="100MB") \
            .reset_index() \
            .set_index("time") \
            .to_parquet(f"{outputFile}.final")

    def _atomicReplace(self, outputFile):
        """Atomically swap the '.final' staging file into the target output path.

        Removes the old output (file or directory) first, then renames.
        """
        if os.path.exists(outputFile):
            if os.path.isfile(outputFile):
                os.remove(outputFile)
            else:
                shutil.rmtree(outputFile)
        os.rename(f"{outputFile}.final", outputFile)

    def _cleanupTmpFiles(self, tmpFiles):
        """Remove all temporary block files after a successful merge."""
        logger = get_classMethod_logger(self, "_cleanupTmpFiles")
        logger.debug("Removing the old tmp files. ")
        for fileTodelete in tmpFiles:
            if os.path.isfile(fileTodelete):
                os.remove(fileTodelete)
            else:
                shutil.rmtree(fileTodelete)

    def writeList(self,theList,blockID,filtersDict,regularMesh,fileExt):
        """Write a list of time step data blocks to temporary files."""
        logger = get_classMethod_logger(self,"writeList")
        filterList = [x for x in theList[0].keys()]
        for filterName in filterList:
            outputFilterName = filterName.replace(".","-")
            outputPath = os.path.dirname(filtersDict[filterName])
            outputFile = os.path.join(outputPath,f"tmp_{outputFilterName.replace('.','-')}_{blockID:06d}.{fileExt}")

            logger.info(f"\tWriting filter {filterName} in temporary file {outputFile} ")

            concat_list = []
            for item in theList:
                if isinstance(item[filterName], list):
                    concat_list.extend(item[filterName])
                else:
                    concat_list.append(item[filterName])
            if regularMesh:
                ds_slice = xarray.concat(concat_list, dim='time')
                ds_slice.to_zarr(outputFile,mode='w')
            else:
                block_data = pandas.concat(concat_list, ignore_index=True,sort=True)
                data = dd.from_pandas(block_data, npartitions=1)
                data.sort_values("time").set_index("time").to_parquet(outputFile)


    ############################################################################################
    ####        Depracated
    ############################################################################################

    @deprecated(reason="Use writeRegularCase instead")
    def write_netcdf(self, datasourcenamelist, timeList=None, fieldnames=None, tsBlockNum=50, overwrite=False,
                     append=False):
        """Write datasources as netcdf files (deprecated)."""
        self.writeRegularCase(datasourcenamelist, timeList, fieldnames, tsBlockNum, overwrite,append)


    @deprecated(reason="Old Name, use readTimeSteps with regularMesh=False")
    def to_pandas(self, datasourcenamelist, timelist=None, fieldnames=None):
        """Read time steps as pandas DataFrames (deprecated)."""
        return self.readTimeSteps(datasourcenamelist, timelist, fieldnames, regtularMesh=False)

    @deprecated(reason="Old Name, use readTimeSteps with regularMesh=True")
    def to_xarray(self, datasourcenamelist, timelist=None, fieldnames=None):
        """Read time steps as xarray Datasets (deprecated)."""
        return self.readTimeSteps(datasourcenamelist, timelist, fieldnames, regtularMesh=True)

    @deprecated(reason="Old Name, use readTimeSteps with regularMesh=False")
    def to_dataFrame(self, datasourcenamelist, timelist=None, fieldnames=None):
        """Read time steps as pandas DataFrames (deprecated)."""
        return self.readTimeSteps(datasourcenamelist, timelist, fieldnames, regtularMesh=False)

    @deprecated(reason="Old Name, use readTimeSteps with regularMesh=True")
    def to_dataArray(self, datasourcenamelist, timelist=None, fieldnames=None):
        """Read time steps as xarray DataArrays (deprecated)."""
        return self.readTimeSteps(datasourcenamelist, timelist, fieldnames, regtularMesh=True)
    @deprecated(reason="Use writeNonRegularCase instead")
    def write_parquet(self, datasourcenamelist, timeList=None, fieldnames=None, tsBlockNum=50, overwrite=False,
                      append=False, filterList=None):
        """Write datasources as parquet files (deprecated)."""
        writeNonRegularCase(datasourcenamelist, timeList, fieldnames, tsBlockNum, overwrite, append, filterList)

    @deprecated(reason="Use writeNonRegularCase instead")
    def write_parquet(self, datasourcenamelist, timeList=None, fieldnames=None, tsBlockNum=50, overwrite=False,
                      append=False, filterList=None):
        """Write datasources as parquet files (deprecated)."""
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
