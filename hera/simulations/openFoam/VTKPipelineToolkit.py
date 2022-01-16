import os
import sys
import numpy
if sys.version_info > (3, 7):
    from hera.toolkit import abstractToolkit, TOOLKIT_SAVEMODE_NOSAVE,TOOLKIT_SAVEMODE_ONLYFILE,TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,TOOLKIT_SAVEMODE_FILEANDDB,TOOLKIT_SAVEMODE_FILEANDDB_REPLACE
    python27=False
else:
    # Define a dummy class for the 2.7 version.
    python27=True
    class abstractToolkit(object):

        _filesDirectory = None
        _projectName = None

        @property
        def filesDirectory(self):
            return self._filesDirectory

        @property
        def projectName (self):
            return self._projectName

        def __init__(self,projectName,toolkitName,filesDirectory):
            self._filesDirectory = filesDirectory
            self._projectName = projectName



#### import the simple module from the paraview
try:
    import paraview.simple as pvsimple
    from paraview import servermanager

    #### disable automatic camera reset on 'Show'
    pvsimple._DisableFirstRenderCameraReset()
    paraview = True
except ImportError:
    ## Working in 3.8 with hera.
    paraview = False

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

    _children = None # A map of the child filters.

    def __init__(self, projectName, readerName, casePath, filesDirectory=None, timeList=None, caseType=DECOMPOSED_CASE, fieldNames=None, serverName=None):
        if python27:
            super(ReaderNode,self).__init__(toolkitName="VTKPipelineToolkit",projectName=projectName,filesDirectory=filesDirectory)
        else:
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

        self._children = dict()
        self.initializeReader(readerName=readerName,
                              casePath=casePath,
                              timeList=timeList,
                              caseType=caseType,
                              fieldNames=fieldNames,
                              serverName=serverName)

    def initializeReader(self, readerName, casePath, timeList=None, caseType=DECOMPOSED_CASE, fieldNames=None, serverName=None):
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

        if serverName is not None:
            pvsimple.Connect(serverName)

        self._readerName  = readerName
        self._reader = pvsimple.OpenFOAMReader(FileName="%s/tmp.foam" % casePath, CaseType=caseType, guiName=readerName)
        self._reader.MeshRegions.SelectAll()
        self._possibleRegions = list(self._reader.MeshRegions)
        self._reader.MeshRegions = ['internalMesh']
        if fieldNames is not None:
            self._reader.CellArrays = fieldNames
        self._reader.UpdatePipeline()

        return self._reader

    def addNode(self, item):
        """
            Define a new filter under this one.
        Parameters
        ----------
        item

        Returns
        -------

        """
        if paraview:
            if item in self._children.keys():
                ret =  self._children[item]
            else:
                newnode = VTKNode(item, self.reader, self.readerName)
                self._children[item] = newnode
                ret = newnode
            return ret
        else:
            raise ValueError("Cannot define VTK pipeline when paraview.simple is not installed. Install or switch to the right environment")

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
        self._timeList = numpy.atleast_1d(value)

    def __getitem__(self, item):
        return self._children[item]

    def keys(self):
        return self._children.keys()

    def writePipeline(self):
        """
            Writes all the nodes in the pipeline.

            First make a list of all the nodes that are written as parquet and a list of all the node that are written as

        Returns
        -------

        """

class VTKNode:

    _name       = None
    _filter     = None
    _filterType   = None
    _father     = None
    _JSONpath       = None
    _saveNode = None
    _regularize = None

    _children = None

    @property
    def name(self):
        return self._name

    @property
    def filter(self):
        return self._filter

    @property
    def father(self):
        return self._father

    @property
    def JSONpath(self):
        return self._JSONpath

    @property
    def filterType(self):
        return self._filterType

    @property
    def saveNode(self):
        return self._saveNode

    @saveNode.setter
    def saveNode(self, value):
        self._saveNode = value

    @property
    def regularize(self):
        return self._regularize

    @regularize.setter
    def regularize(self, value):
        self._regularize = value

    def __init__(self, name, father, JSONpath,saveNode=False,regularize=None):
        """
            Initializes a node.

            The parameters in vtk are given in structured names :

                a.b = 3

                so we changed it to : a_b .

        Parameters
        ----------
        nodeType : str
            The type of the filter (i.e clip,slice and ect).
        father : filter
            The main fileter.

        JSONpath: str
            The path to the root (the list of all the nodes from the root to the current node).

        saveNode : bool
            If true save the data of the filter.
            The format depends on the regularize data.

        regularize : str/None/True
            If None, save as raw data (in parquet)
            if True, try to convert to xarray and save as netcdf
            if 'str', use the function to regularize the data. [not implemented yet]
        """
        self._name = name
        self._father = father
        self._JSONpath = JSONpath
        self._children = dict
        self._saveNode = saveNode
        self._regularize = regularize

    def __call__(self, filterType, **kwargs):
        self._filterType = filterType
        filtobj = getattr(pvsimple, self.filterType)
        if filtobj is None:
            raise ValueError("%s is not a valid VTK object. The valid objects are: %s " % (item, '\n'.join(dir(pvsimple))))

        self._filter= filtobj(Input=self.father, guiName=self.name)
        for param, pvalue in kwargs.items():
            pvalue = str(pvalue) if isinstance(pvalue, unicode) else pvalue  # python2, will be removed in python3.
            paramnamelist = param.split("_")
            paramobj = self._filter
            for pname in paramnamelist[:-1]:
                paramobj = getattr(paramobj, pname)
            setattr(paramobj, paramnamelist[-1], pvalue)
        self._filter.UpdatePipeline()

    def addNode(self, item):
        """
            Define a new filter under this one.
        Parameters
        ----------
        item

        Returns
        -------

        """
        if paraview:
            if item in self._children:
                ret =  self._children[item]
            else:
                newnode = VTKNode(item, self.filter, self.JSONpath + "_" + self.filterType)
                self._children[item] = newnode
                ret = newnode
            return ret
        else:
                raise ValueError("Cannot define VTK pipeline when paraview.simple is not installed. Install or switch to the right environment")

    def __getitem__(self, item):
        return self._children[item]

    def keys(self):
        return self._children.keys()

    def getData(self):
        """
            Return the data of the filter from the hera datalayer.
            Must write the filter first (use writeData in the root node) and load (use loadData) before using
            this function.

        Parameters
        ----------

        Returns
        -------
            pandas/dask of the results. or xarray.
        """
        if paraview:
            raise NotImplementedError("Get data for a node is not implemented. Use the write at the root node. Then in the conda environment load it, and get it")
        else:
            pass




if __name__=="__main__":
    rootnode = ReaderNode(projectName="test",readerName="root",casePath="channel",filesDirectory="cache")
    rootnode.addNode("CellCenters")(filterType="CellCenters")
    print(rootnode['CellCenters'])
