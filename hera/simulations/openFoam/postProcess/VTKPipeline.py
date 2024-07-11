import pydoc
from .. import CASETYPE_DECOMPOSED,CASETYPE_RECONSTRUCTED
from .VTKPipelineExecutionContext import VTKpipelineExecutionContext


VTKFILTER_TYPE_PARQUET = "parquet"
VTKFILTER_TYPE_NETCDF = "netcdf"

class VTKPipeLine:

    """
        This class represents a VTK pipeline to handle.

        The execution of the VTK Pipeline is handled with the VTKPipelineExecution.
        Hence, the current class represents only the
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

      Execution context
      ------------------

        The execution context is a VTKpipelineExecutionContext.

      Executing the pipeline requires:
        1. the case directory
        2. server name (if it exists).
        3.

        After the context is set, a openFOAM reader is created.
        This is because it takes time to initialize the openFOAM readers.

        However, ometimes the user would want to get the time lists.
        This is done by defining a reader and returning the time step.

        Hence, we allow the user to set an execution context that is a directory of
        these two parameters.

        Once the context is set, it is possible to get the timesteps, and the fields.

    """

    filters = None
    timeList = None

    def __init__(self):
        self.filters = dict()

    def addFilter(self, vtkFilter):
        self.filters[vtkFilter.name] = vtkFilter

    def setExecutionContext(self,OFDatalayer,nameOrWorkflowFileOrJSONOrResource,serverName,caseType=CASETYPE_DECOMPOSED,fieldNames=None):
        """
            Create the execution context (instance of the VTKpipelineExecutionContext)
        Parameters
        ----------
        OFDatalayer : toolkit
            The handle to the project management.
        nameOrWorkflowFileOrJSONOrResource : str
            The reference to the flow directory.
            Can be:
                - Case directory
                - workflow object
                - workflow paramters. (that result in 1 flow. If more than 1 flow are
        serverName : str
            The name of the paraview server. Use None if no server is used.
        caseType : str
            CASETYPE_DECOMPOSED,CASETYPE_RECONSTRUCTED

        fieldNames : list
            Thenames of the fields to load.

        Returns
        -------

        """

        self._nameOrWorkflowFileOrJSONOrResource = nameOrWorkflowFileOrJSONOrResource

        if os.path.isdir(nameOrWorkflowFileOrJSONOrResource):
            casePath = nameOrWorkflowFileOrJSONOrResource
        else:
                simulationDocument = self.datalayer.getWorkflowDocumentFromDB(nameOrWorkflowFileOrJSONOrResource)
                if simDoc is not None:
                    casePath = simulationDocument.resource
                else:
                    raise ValueError(
                        f"Simulation {nameOrWorkflowFileOrJSONOrResource} is not in the DB, and does not represent a valid case directory")

        self.executionContext = VTKpipelineExecutionContext(self.toJSON(),
                                                            casePath,
                                                            caseType=CASETYPE_DECOMPOSED,
                                                            serverName=serverName,
                                                            fieldNames=fieldNames)

    def loadToProject(self,datalayer,overwrite=False):
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
        logger = get_classMethod_logger(self, "loadToProject")
        logger.info(f"Loading pipeline to project {datalayer.projectName}")
        for filterName, filterData in self.VTKpipelineJSON["pipelines"].items():
            logger.execution(f"Handling filter {filterName}")
            self._recurseNode(filterName=filterName, filterData=filterData, path="",overwrite=overwrite,datalayer=datalayer)

    def _recurseNode(self, filterName, filterData, path,datalayer,overwrite=False):
        """
        This function is recursively passed on the pipeline tree from json file and loads the filters to the database.

        Parameter
        ----------

        """
        logger = get_classMethod_logger(self, "_recurseNode")
        logger.info(f"Adding filter {filterName} in {path}")
        formatName = filterData.get('write', None)
        if (formatName is not None) and (formatName != "None"):
            logger.debug("Getting the simulation flow document. Use the path as identifier ")
            # Find the flow document.
            flowDoc = datalayer.getCaseListDocumentFromDB(self._casePath)

            if len(flowDoc) == 0:
                logger.error(f"P{self._casePath} is not found in the project. Load it using the hera-OF-workflow add ")
                raise ValueError(
                    "The simulation in path %s is not found!. Load it to the database first" % self._casePath)

            # Create the new description

            filterDesc = dict(flowDoc['desc'])

            simqry = {
                "simulationName" : filterDesc['simulationName'],
                "groupName" : filterDesc['groupName'],
                "filterName" : filterName
            }

            # check if it is already loaded.
            logger.debug("Checking if the filter is loaded in the project")
            docList = datalayer.getCacheDocuments(type=TYPE_VTK_FILTER,**simqry)

            if len(docList) ==0 or overwrite:
                logger.debug(f"Updating the filter data. Overwrite flag {overwrite}")

                filterDesc = dict(simqry)
                filterDesc['filterParameters'] = filterData.get('params',{})
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

                if len(docList) > 0:
                    logger.debug("Updating the existing data")
                    doc = docList[0]
                    doc['desc'] = filterDesc
                    doc['resource'] = resource
                    doc['dataFormat'] =currentdatatype
                    doc.update()
                else:
                    logger.debug("Adding a new record")
                    datalayer.addCacheDocument(desc=filterDesc,type=TYPE_VTK_FILTER,resource=resource,dataFormat = currentdatatype)
            else:
                logger.warning("Filter %s for simulation %s in group %s already in the database" % (filterName,filterDesc['simulationName'],filterDesc['groupName']))

        logger.execution("Processing the downstream filters. ")
        ds = filterData.get("downstream", {})
        for filterName, filterData in ds.items():
            path += "." + filterName
            self._recurseNode(filterName=filterName, filterData=filterData, path=path,overwrite=overwrite,datalayer=datalayer)

    def addFilterByName(self,filterName,shouldWrite,**kwargs):
        """
            Adds the filter to the piple line and return it for configuration.
        Parameters
        ----------
        filterName : str
            The filter name
        shouldWrite : bool
            The boolean

        kwargs : dict
            Extra paremters.

        Returns
        -------
            the added iVTKFilter.
        """
        fltr = self.getFilter(filterName=filterName,shouldWrite=shouldWrite,**kwargs)
        self.filters[filterName] = fltr
        return fltr

    def getFilter(self,filterName,shouldWrite,**kwagrs):
        """
            Returns the requested filter.
        Parameters
        ----------
        filterName : str
            The name of the filter
        shouldWrite : bool
            Should we write the output, or leave it.

        Returns
        -------

        """
        logger = get_classMethod_logger(self,"getFilter")
        try:
            fltr = pydoc.locate(f"{Slice}VTKFilter")
        except ImportError:
            err = f"Filter {filterName} not found!"
            logger.error(err)
            raise ValueError(err)

        return fltr(name=filterName,shouldWrite = shouldWrite,**kwargs)

    def execute(self,sourceOrName=None,timeList=None, tsBlockNum=50, overwrite=False,append=False):
        """
            Returns a VTK Execution that will allow the execution of the VTK pipeline.

            We use the execution to let the user access the reader (and thus, the timesteps).

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
               If true, overwrite on the parquet.

        append : bool
               If true, append to existing parquet

        Returns
        -------

        """

        logger = get_classMethod_logger(self,"execute")

        if self.executionContext is None:
            err = f"The execution context is not set. use the method setExecutionContext to set it"
            logger.error(err)
            raise ValueError(err)

        self.executionContext.execute(sourceOrName=sourceOrName,timeList=timeList, tsBlockNum=tsBlockNum, overwrite=overwrite,append=append)

    def toJSON(self):
        """
            Converts the pipeline to a VTK JSON of the executions.

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

        Returns
        -------

        """
        retDict = dict()
        for filterName,filterData in self.filters.items():
            retDict.update(filterData.toJSON())

        return dict(metadata=dict(files=None),pipelines=retDict)


class VTKFilter:
    name = None
    type = None  # The type of the filter, clip/slice and ect.
    write = None  # The writing format : None/parquet/netcdf
    shouldWrite = None  # true or false.
    params = None
    downstream = None

    def __init__(self, name, filterType, writeFormat, shouldWrite, params):
        """
            Init an abstract node.

            "filterName": {
                "type": The type of the filter.(clip, slice, ...).
                "write": None / parquet(pandas) / netcdf(xarray),
                "params": [
                    ("key", "value"),
                    .
                    .
                    .
                ], ...
                    "downstream": [Another pipeline]
            }

        Parameters
        ----------
        filterType  : str
            A VTK filter name : clip/cell centers and ect.
        writeFormat : str or None
            None
            parquet - save as parquet (using pandas dataframe)
            netcdf -  save as netcdf  (using xarray dataframe)
        shouldWrite : bool
            Should we write this filter output?

        params : dict
            A list of key, value for the parameters of the dict.
            The order of the paramters is important, as it might change the behaviour of the object (like in the slice or clip filters).
        """
        self.name = name
        self.filterType = filterType
        self.writeFormat = writeFormat
        self.shouldWrite = shouldWrite
        self.params = params
        self.downstream = []

    def toJSON(self):
        """
            converts the node
        Returns
        -------

        """
        ret = dict()
        ret['type'] = self.filterType
        ret['write'] = self.writeFormat if self.shouldWrite else None
        ret['params'] = self.params

        for dsFilter in self.downstream:
            ret['downstream'] = [x.toJSON() for x in dsFilter]

        return {self.name: ret}

    @property
    def filterType(self):
        return self.write

    @filterType.setter
    def filterType(self, value):
        logger = get_classMethod_logger(self,"filterType")
        if value is not None:
            if value not in [VTKFILTER_TYPE_PARQUET,VTKFILTER_TYPE_NETCDF]:
                err = f"The filter must be either {VTKFILTER_TYPE_PARQUET} or {VTKFILTER_TYPE_NETCDF}"
                logger.error(err)
                raise ValueError(err)
        self.write = value

class SliceVTKFilter(VTKFilter):

    def __init__(self, name, shouldWrite,**kwargs):
        kwargs.setdefault("params",[])
        super().__init__(name=name, filterType="Slice", writeFormat="parquet", shouldWrite=shouldWrite,**kwargs)

    def setPlane(self, origin):
        """
            Setting slice as plane with the origin
        Parameters
        ----------
        origin : list
            a 3-tuple of the origin

        Returns
        -------

        """
        self.params = params = [("SliceType", "Plane"),
                                ("HyperTreeGridSlicer", "Plance"),
                                ("SliceOffsetValues", 0),
                                ("SliceType.Origin", [0, 0, 0]),
                                ("HyperTreeGridSlicer.Origin", [0, 0, 0])
                                ]

class CellCentersVTKFilter(VTKFilter):
    """
        The Cell center filter.
    """
    def __init__(self,name,shouldWrite,**kwargs):
        kwargs.setdefault("params",[])
        super().__init__(name=name, filterType="CellCenters", writeFormat="parquet", shouldWrite=shouldWrite,**kwargs)




















