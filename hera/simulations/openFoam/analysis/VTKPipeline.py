import pandas
import os
import json
import glob

import sys
from .. import DECOMPOSED_CASE,RECONSTRUCTED_CASE
from ....utils import loadJSON,loggedObject
from ....datalayer import datatypes
paraviewExists = False
try:
    import paraview.simple as pvsimple
    from ..datalayer.pvOpenFOAMBase import paraviewOpenFOAM

    paraviewExists = True
except ImportError:
    print("paraview not Found. Cannot use the VTK pipeline.")

class VTKpipeline(loggedObject):
    """

    This is a helper class to handle VTK pipelines.
    It has 2 main functions:

        * Execute the pipelines and save the results to parquet or netcdf file
                Must run in the python 2.7 environment.

        * Loads the results to the database.
                Must run in the python 3+ (with hera) environment.

    Currently works only for the JSON pipeline. The XML (paraview native pipelines) will be built in the future.

    The pipeline is initialized with a reader.

    The VTK pipeline JSON structure.
        {
           "metadata" : {
                  "guiname" : <gui name>,
                  "timeList" : ...
            },
            "pipeline" : {
                    "filterName" : {
                            "type" : The type of the filter. (clip,slice,...).
                            "write"   : None/parquet (pandas)/netcdf (xarray),
                            "params" : [
                                    ("key","value"),
                                          .
                                          .
                                          .
                            ],...
                            "downstream" : [Another pipeline]
                            }
                        }
      }
     The write to writes the requested filters to the disk.
     Each filter is saved to a parquet/netcdf file.  The file name is the filter name
    """

    _VTKpipelineJSON = None  # Holds the json of the VTK pipeline.
    _pvOFBase = None  # Holds the OF base.
    _casePath = None
    _reader = None
    _datalayer = None
    _simulationDocument = None

    @property
    def simulationDocument(self):
        return self._simulationDocument

    @property
    def reader(self):
        return self._reader


    @property
    def datalayer(self):
        return self._datalayer

    @property
    def pvOFBase(self):
        return self._pvOFBase

    @property
    def VTKpipelineJSON(self):
        return self._VTKpipelineJSON

    def __init__(self, datalayer, pipelineJSON, nameOrWorkflowFileOrJSONOrResource, caseType=DECOMPOSED_CASE, serverName=None, fieldNames=None):
        """
            Initializes a VTK pipeline.

            {
                "metadata" : {
                        "fieldNames"   : [optional] The field names to Load.
                },
                "pipelines" : {
                        <filter name> : {
                            "type" : <filter type>,
                            "params" : [
                                {
                                    <parameter name>  : <parameter value>

                                }
                            ],
                            "write" : ['dataFrame','dataArray'],
                            "downstream" : {
                                    .... filters...
                            }
                        }
                }
            }


        Parameters
        ----------
        datalayer : simulation.openFOAM toolkit
                The data layer to use to load the pipeline.

        pipelineJSON : str,dict
                the JSON (filename, or string) or dict.

        nameOrWorkflowFileOrJSONOrResource: str
             - Resource (the directory name)
             - Simulation name
             - Its workflow
             - workfolow dict.

        CaseType:  str
                Either 'Decomposed Case' for parallel cases or 'Reconstructed Case'
                for single processor cases.

        fieldnames: None or list of field names.  default: None.
                The list of fields to load.
                if None, read all fields

        servername: str
                if None, work locally.
                connection string to the paraview server.

                The connection string is printed when the server is initialized.
        """
        super().__init__()
        self._VTKpipelineJSON = loadJSON(pipelineJSON)
        self._datalayer = datalayer
        self._serverName = serverName
        self._caseType = caseType


        self._simulationDocument = self.datalayer.getSimulationDocumentFromDB(nameOrWorkflowFileOrJSONOrResource)

        if os.path.exists(nameOrWorkflowFileOrJSONOrResource) and self._simulationDocument is None:
            # the simulation is not in the DB, but we wish to analyze it.
            casePath = nameOrWorkflowFileOrJSONOrResource
        elif self._simulationDocument is not None:
            casePath = self._simulationDocument.resource
        else:
            raise ValueError(f"Simulation {nameOrWorkflowFileOrJSONOrResource} is not in the DB, and does not represent a valid case directory")

        if paraviewExists:
            self._pvOFBase = paraviewOpenFOAM(casePath=nameOrWorkflowFileOrJSONOrResource,
                                              caseType=caseType,
                                              servername=self._serverName,
                                              fieldNames=fieldNames)

            self._pvOFBase.parquetdir = os.path.abspath(os.path.join(casePath, 'parquet'))
            self._pvOFBase.netcdfdir = os.path.abspath(os.path.join(casePath, "netcdf"))
        else:
            self._parquetdir = os.path.abspath(os.path.join(casePath, 'parquet'))
            self._netcdfdir = os.path.abspath(os.path.join(casePath, "netcdf"))

        self._casePath = casePath
        self._fieldNames = pipelineJSON["metadata"].get("fieldNames", fieldNames)

    def _initializeReader(self, readerName, casePath):
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
        self._readerName  = readerName
        self._reader = pvsimple.OpenFOAMReader(FileName="%s/tmp.foam" % casePath, CaseType=self._caseType, guiName=readerName)
        self._reader.MeshRegions.SelectAll()
        self._possibleRegions = list(self._reader.MeshRegions)
        self._reader.MeshRegions = ['internalMesh']
        if self._fieldNames is not None:
            self._reader.CellArrays = self._fieldNames
        self._reader.UpdatePipeline()

        return self._reader


    def execute(self, sourceOrName=None,timeList=None, tsBlockNum=100, overwrite=False):
        """
            Executes the pipeline from the JSON vtk.
            Saves the output of the requested filters.

        Parameters
        ----------
        sourceOrName: str , None
                A reader, the name of the reader, or none (to create reader).

                If reader, then assume that the server is already connected.

        timeList : list
                The list of timesteps to read.

        tsBlockNum: int
                The block number
        overwrite : bool
                If true, overwrite on the parquet. Else append to it.

        Returns
        -------
            None
        """
        self.logger.info(f"Executing pipeline")
        if isinstance(sourceOrName,str):

            self.logger.debug(f"Getting the reader {sourceOrName}")
            reader = pvsimple.FindSource(sourceOrName)
            if reader is None:
                reader = self._initializeReader(readerName="reader", casePath=self._casePath, timeList=timeList, caseType=self._casePath, fieldNames=self._fieldNames)
        else:
            # assume server is connected.
            reader = sourceOrName

        # build the pipeline.
        filterWrite = {}
        self._buildFilterLayer(father=reader, structureJson=self._VTKpipelineJSON["pipelines"], filterWrite=filterWrite)
        # Now execute the pipeline.

        if timeList is None:
            timelist = self._VTKpipelineJSON["metadata"].get("timelist", None)

        if timelist is not None and isinstance(timelist, str):
            # a bit of parsing.
            readerTL = reader.TimestepValues
            BandA = [readerTL[0], readerTL[-1]]

            for i, val in enumerate(timelist.split(":")):
                BandA[i] = BandA[i] if len(val) == 0 else float(val)

            tl = pandas.Series(readerTL)
            timelist = tl[tl.between(*BandA)].values

        # Get the mesh regions.

        if "MeshRegions" in self._VTKpipelineJSON["metadata"]:
            reader.MeshRegions = self._VTKpipelineJSON["metadata"]["MeshRegions"]

        for frmt, datasourceslist in filterWrite.items():
            writer = getattr(self._pvOFBase, "write_%s" % frmt)
            if writer is None:
                raise ValueError("The write %s is not found" % writer)
            writer(datasourcenamelist=datasourceslist,
                   timelist=timelist,
                   fieldnames=self._fieldNames,
                   tsBlockNum=tsBlockNum,
                   overwrite=overwrite)

    def _buildFilterLayer(self, father, structureJson, filterWrite):
        """
            Recursively builds the structure of the leaf.
            Populates the self._filterWrite map

            Since the order of setting the params might be of importance (for example, setting the
            plane type determine the rest of the parameters), we set it as a list.

        :param father:
                The current filter father of the layer.

        :param structureJson:
                The portion of Json to build.

        :param[output]   filterWrite
                an  dictionary with the names of the filters that are about
                to be printed according to format.

        """
        if structureJson is None:
            return

        for filterGuiName in structureJson:
            paramPairList = structureJson[filterGuiName]['params']  # must be a list to enforce order in setting.
            filtertype = structureJson[filterGuiName]['type']
            filter = getattr(pvsimple, filtertype)(Input=father, guiName=filterGuiName)
            for param, pvalue in paramPairList:
                #pvalue = str(pvalue) if isinstance(pvalue, unicode) else pvalue  # python2, will be removed in python3.
                paramnamelist = param.split(".")
                paramobj = filter
                for pname in paramnamelist[:-1]:
                    paramobj = getattr(paramobj, pname)
                setattr(paramobj, paramnamelist[-1], pvalue)
            filter.UpdatePipeline()
            writeformat = structureJson[filterGuiName].get("write", None)

            if (writeformat is not None) and (str(writeformat) != "None"):
                filterlist = filterWrite.setdefault(writeformat, [])
                filterlist.append(filterGuiName)

            self._buildFilterLayer(filter, structureJson[filterGuiName].get("downstream", None), filterWrite)

    def loadToProject(self,overwrite=False):
        """
            Loads the pipeline results to the database.
            We assume that the pipeline was already executed, and that the

        Parameters
        ----------
        projectName : str
                The project name to load the data to.

        Returns
        -------

        """
        self.logger.info(f"Loading pipeline to project {self.datalayer.projectName}")
        for filterName, filterData in self.VTKpipelineJSON["pipelines"].items():
            self.logger.execution(f"Handling filter {filterName}")
            self._recurseNode(filterName=filterName, filterData=filterData, path="",overwrite=overwrite)

    def _recurseNode(self, filterName, filterData, path,overwrite=False):
        """
        This function is recursively passed on the pipeline tree from json file and loads the filters to the database.

        Parameter
        ----------

        """
        formatName = filterData.get('write', None)
        if (formatName is not None) and (formatName != "None"):
            self.logger.debug("Getting the simulation flow document. Use the path as identifier ")
            # Find the flow document.
            flowDoc = self.datalayer.getSimulationDocumentFromDB(self._casePath)

            if len(flowDoc) == 0:
                self.logger.error(f"P{self._casePath} is not found in the project. Load it using the hera-OF-workflow add ")
                raise ValueError(
                    "The simulation in path %s is not found!. Load it to the database first" % self._casePath)

            import pdb
            pdb.set_trace()

            flowDoc = flowDoc[0]
            # Create the new description
            filterDesc = dict(flowDoc['desc'])

            simqry = {
                "simulationName" : filterDesc['simulationName'],
                "groupName" : filterDesc['groupName'],
                "filterName" : filterName
            }

            # check if it is already loaded.
            docList = self.datalayer.getCacheDocuments(type="vtk_filter",**simqry)

            if len(docList) ==0 or overwrite:
                filterDesc['filterName'] = filterName
                filterDesc['filterParameters'] = filterData.get('parameters',{})
                filterDesc['pipeline'] = self.VTKpipelineJSON
                filterDesc['filterPath'] = path

                if formatName in ["netcdf","dataArray"]:
                    resource = glob.glob(os.path.join(self._netcdfdir, "%s_*.nc" % filterName))
                    currentdatatype = datatypes.NETCDF_XARRAY
                elif formatName in ["parquet","dataFrame"]:
                    resource = os.path.join(self._parquetdir, "%s.parquet" % filterName)
                    currentdatatype = datatypes.PARQUET
                else:
                    raise ValueError(
                        "Format type %s unknown for filter %s. Must be 'netcdf' or 'parquet' (case sensitive!)" % (
                        formatName, filterName))

                self.datalayer.addCacheDocument(desc=filterDesc,type="vtk_filter",resource=resource,dataFormat = currentdatatype)
            else:
                self.logger.warning("Filter %s for simulation %s in group %s already in the database" % (filterName,filterDesc['simulationName'],filterDesc['groupName']))

        # processing the
        ds = filterData.get("downstream", {})
        for filterName, filterData in ds.items():
            path += "." + filterName
            self._recurseNode(filterName=filterName, filterData=filterData, path=path,overwrite=overwrite)


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
            #pvalue = str(pvalue) if isinstance(pvalue, unicode) else pvalue  # python2, will be removed in python3.
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
        if paraviewExists:
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
