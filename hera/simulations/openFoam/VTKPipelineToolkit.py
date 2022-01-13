import os
from collections.abc import Iterable

import numpy
from hera.toolkit import abstractToolkit, TOOLKIT_SAVEMODE_NOSAVE,TOOLKIT_SAVEMODE_ONLYFILE,TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,TOOLKIT_SAVEMODE_FILEANDDB,TOOLKIT_SAVEMODE_FILEANDDB_REPLACE

#### import the simple module from the paraview
import paraview.simple as pvsimple
from paraview import servermanager
#### disable automatic camera reset on 'Show'
pvsimple._DisableFirstRenderCameraReset()

RECONSRUCTED_CASE = "Reconstructed Case"
DECOMPOSED_CASE = 'Decomposed Case'



class ReaderNode(abstractToolkit):
    """
        A toolkit to handle VTK pipelines.

        The ReaderNode is the rootnode. Other nodes get the data from it.
    """

    _componentsNames = None  # names of components for reading.
    _timeList        = None
    _reader         = None
    _readerName     = None

    def __init__(self,projectName,filesDirectory=None):
        super().__init__(toolkitName="VTKPipelineToolkit",projectName=projectName,filesDirectory=filesDirectory)

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

    def initializeReader(self, readerName, casePath, timeList=None , CaseType=DECOMPOSED_CASE, fieldnames=None,servername=None):
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

        if servername is not None:
            pvsimple.Connect(servername)

        self._readerName  = readerName
        self._reader = pvsimple.OpenFOAMReader(FileName="%s/tmp.foam" % casePath, CaseType=CaseType, guiName=readerName)
        self._reader.MeshRegions.SelectAll()
        self._possibleRegions = list(self._reader.MeshRegions)
        self._reader.MeshRegions = ['internalMesh']
        if fieldnames is not None:
            self._reader.CellArrays = fieldnames
        self._reader.UpdatePipeline()

        return self._reader

    @property
    def readerName(self):
        return self._readerName

    @property
    def reader(self):
        return self._reader

    @property
    def timelist(self):
        return self._timeList

    @timelist.setter
    def timelist(self, value):
        if not isinstance(value,Iterable):
            self._timeList = numpy.atleast_1d(value)
        else:
            self._timeList = value

class VTKNode:

    _filter     = None
    _nodeName   = None
    _reader     = None
    _path       = None

    @property
    def filter(self):
        return self._filter

    @property
    def reader(self):
        return self._reader

    @property
    def path(self):
        return self._path

    @property
    def name(self):
        return self._nodeName

    def __init__(self,nodeName,reader,path,**kwargs):
        """
            Initializes a node.

            The parameters in vtk are given in structured names :

                a.b = 3

                so we changed it to : a_b .

        Parameters
        ----------
        nodeName
        reader
        path: str
            The path to the root (the list of all the nodes from the root to the current node).

        kwargs : dict
            The name->value of the paremters.
            The paremter names are separated with '_'.

        """
        self._nodeName = nodeName
        self._reader = reader
        self._path = path

        self._filter = getattr(pvsimple, filtertype)(Input=father, guiName=filterGuiName)
        for param, pvalue in kwargs.items():
            pvalue = str(pvalue) if isinstance(pvalue, unicode) else pvalue  # python2, will be removed in python3.
            paramnamelist = param.split("_")
            paramobj = self._filter
            for pname in paramnamelist[:-1]:
                paramobj = getattr(paramobj, pname)
            setattr(paramobj, paramnamelist[-1], pvalue)
        self._filter.UpdatePipeline()

    def getData(self,timeList=None,regularize=None,fieldNames=None,saveMode=TOOLKIT_SAVEMODE_FILEANDDB):
        """
            Return the data of the filter.
        Parameters
        ----------
        timeList
        saveMode: str
            Determine whether to cache the result in the database.

            - stored as a file (but not saved to the DB)  - return dask.DataFrame
            - stored and saved in db                      - return dask.DataFrame
            - not store (just read all to the memeory and return pandas.DataFrame).

            The files are stored to the root.filesDirectory/... path ..../ data.parquet/

        regularize : str/bool
            If None - save the data as is in parquet (DataFrame)
            If true - try to convert to xarray. If saved than as an netcdf file.
            if str  - run a regularize function and return an xarray. If saved than as an netcdf file. [Not implemented yet]

        Returns
        -------
            pandas/dask of the results. or xarray.
        """
        timeList = self.reader.TimestepValues if timeList is None else numpy.atleast_1d(timelist)

        # Search if that filter results already exists in the cache.
        doc = None # The document of the file, if the requested time is not in the time list.

        ret = []
        for timeList in timeList:

            self.filter.UpdatePipeline(timeslice)
            rawData = servermanager.Fetch(datasource)
            data = dsa.WrapDataObject(rawData)
            #print(data.PointData)
            #print(type(data.PointData))

            if isinstance(data.Points, dsa.VTKArray):
                points = numpy.array(data.Points).squeeze()
            else:
                points = numpy.concatenate([numpy.array(x) for x in data.Points.GetArrays()]).squeeze()

            curstep = pandas.DataFrame()

            # create index
            curstep['x'] = points[:, 0]
            curstep['y'] = points[:, 1]
            curstep['z'] = points[:, 2]
            curstep['time'] = timeslice

            fieldlist = data.PointData.keys() if fieldNames is None else fieldNames
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
                    L = [slice(None, None, None)] + list(indxiter)
                    try:
                        curstep["%s%s" % (field, self._componentsNames[indxiter])] = arry[L]
                    except ValueError:
                        print("Field %s is problematic... ommiting" % field)

            curstep = curstep.set_index(['time', 'x', 'y', 'z'])

            if saveMode == TOOLKIT_SAVEMODE_NOSAVE:
                ret.append(curstep)
            elif saveMode in [TOOLKIT_SAVEMODE_ONLYFILE, TOOLKIT_SAVEMODE_FILEANDDB,TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,TOOLKIT_SAVEMODE_FILEANDDB_REPLACE]:
                # getting the file name
                if doc is not None:
                    # The file exists but not with that time step.
                    filename = doc.resource
                else:
                    extension = "parquet" if regularize is None else "nc"
                    filepath = os.path.join(self.reader.filesDirectory,self.path)
                    os.makedirs(filepath,exist_ok=True)
                    filename = os.path.join(filepath,f"{self.name}.{extension}")

                if os.path.exists(filename):
                    # append to dask/nc file
                else:
                    if regularize is None:
                        # append to parquet.

                    else:
                        # append to netcdf.
                        data.to_netcdf(curfilename)