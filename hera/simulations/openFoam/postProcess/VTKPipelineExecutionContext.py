import pandas
import os
import glob
from .. import CASETYPE_DECOMPOSED, TYPE_VTK_FILTER
from ....utils import loadJSON,get_classMethod_logger
from ....datalayer import datatypes
paraviewExists = False
try:
    import paraview.simple as pvsimple
    paraviewExists = True
except ImportError:
    print("paraview not Found. Cannot execute the VTK pipeline.")

from .pvOpenFOAMBase import paraviewOpenFOAM

class VTKpipelineExecutionContext:
    """

    This is a helper class to execute VTK pipelines.
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

    @property
    def reader(self):
        return self._reader

    @property
    def pvOFBase(self):
        return self._pvOFBase

    @property
    def VTKpipelineJSON(self):
        return self._VTKpipelineJSON

    def __init__(self, pipelineJSON, casePath, caseType=CASETYPE_DECOMPOSED, serverName=None, fieldNames=None):
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
                            "write" : ['parquet','netcdf'],
                            "downstream" : {
                                    .... filters...
                            }
                        }
                }
            }

        use parquet if non-regular data and netcdf for regular data.


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
        logger = get_classMethod_logger(self,"__init__")
        self._VTKpipelineJSON = loadJSON(pipelineJSON)

        self._serverName = serverName
        self._caseType = caseType

        if paraviewExists:
            self._pvOFBase = paraviewOpenFOAM(casePath=casePath,
                                              caseType=caseType,
                                              servername=self._serverName,
                                              fieldNames=fieldNames)

            self._pvOFBase.parquetdir = os.path.abspath(os.path.join(casePath, 'parquet'))
            self._pvOFBase.netcdfdir = os.path.abspath(os.path.join(casePath, "netcdf"))

        self._parquetdir = os.path.abspath(os.path.join(casePath, 'parquet'))
        self._netcdfdir = os.path.abspath(os.path.join(casePath, "netcdf"))

        self._casePath = casePath
        self._fieldNames = self._VTKpipelineJSON["metadata"].get("fieldNames", fieldNames)

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
        self._readerName = readerName
        self._reader = self.getReader(readerName=readerName)


    def getReader(self, readerName="reader"):
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
        #self._readerName  = readerName
        reader = pvsimple.OpenFOAMReader(FileName="%s/tmp.foam" % self._casePath, CaseType=self._caseType, guiName=readerName)
        reader.MeshRegions.SelectAll()
        possibleRegions = list(reader.MeshRegions)
        reader.MeshRegions = ['internalMesh']
        if self._fieldNames is not None:
            reader.CellArrays = self._fieldNames
        reader.UpdatePipeline()

        return reader


    def execute(self, sourceOrName=None,timeList=None, tsBlockNum=50, overwrite=False,append=False):
        """
            Executes the pipeline from the JSON vtk.
            Saves the output of the requested filters.

            IF file exists, and it is without --append or --overwrite then raise exception.

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
            None
        """
        logger = get_classMethod_logger(self, "execute")
        logger.info(f"Executing pipeline")
        if sourceOrName is None:
            reader = self.initializeReader(readerName="reader")
        elif isinstance(sourceOrName,str):
            logger.debug(f"Getting the reader {sourceOrName}")
            reader = pvsimple.FindSource(sourceOrName)
            if reader is None:
                reader = self.initializeReader(readerName="reader")
        else:
            logger.debug(f"Got reader object:  {sourceOrName}")
            # assume server is connected.
            reader = sourceOrName

        # build the pipeline.
        filterWrite = {}
        self._buildFilterLayer(father=reader, structureJson=self._VTKpipelineJSON["pipelines"], filterWrite=filterWrite)
        # Now execute the pipeline.

        if timeList is None:
            timelist = self._VTKpipelineJSON["metadata"].get("timelist", None)

        if timeList is not None and isinstance(timeList, str):
            # a bit of parsing.
            readerTL = reader.TimestepValues
            BandA = [readerTL[0], readerTL[-1]]

            for i, val in enumerate(timeList.split(":")):
                BandA[i] = BandA[i] if len(val) == 0 else float(val)

            tl = pandas.Series(readerTL)
            timeList = tl[tl.between(*BandA)].values

        # Get the mesh regions.
        if "MeshRegions" in self._VTKpipelineJSON["metadata"]:
            reader.MeshRegions = self._VTKpipelineJSON["metadata"]["MeshRegions"]

        for frmt, datasourceslist in filterWrite.items():
            writer = getattr(self._pvOFBase, "write_%s" % frmt)
            if writer is None:
                raise ValueError("The write %s is not found" % writer)
            writer(datasourcenamelist=datasourceslist,
                   timeList=timeList,
                   fieldnames=self._fieldNames,
                   tsBlockNum=tsBlockNum,
                   overwrite=overwrite,
                   append=append)

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
        logger = get_classMethod_logger(self, "_buildFilterLayer")
        logger.info(f"building Filter layer {structureJson}")
        if structureJson is None:
            return
        for filterGuiName in structureJson:

            paramPairList = structureJson[filterGuiName]['params']  # must be a list to enforce order in setting.
            filtertype = structureJson[filterGuiName]['type']
            filter = getattr(pvsimple, filtertype)(Input=father, guiName=filterGuiName)
            logger.debug(f"Adding filter {filterGuiName} of type {filtertype} to {father}")
            for param, pvalue in paramPairList:
                logger.debug(f"...Adding parameters {param} with value {pvalue}")
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





class VTKPipeLineold:
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
                            "write"   : bool True/False
                            "params" : [
                                    ("key","value"),
                                          .
                                          .
                                          .
                            ],...
                            "downstream" : {filterName : ....}
                            }
                        }
      }

      Execution context
      ------------------

        The execution context is a VTKpipelineExecutionContext.

      Executing the pipeline requires:
        1. the case directory
        2. server name (if it exists).

        After the context is set, a openFOAM reader is created.
        This is because it takes time to initialize the openFOAM readers.

        However, ometimes the user would want to get the time lists.
        This is done by defining a reader and returning the time step.

        Hence, we allow the user to set an execution context that is a directory of
        these two parameters.

        Once the context is set, it is possible to get the timesteps, and the fields.

    """

    datalayer = None
    filters = None
    timeList = None
    casePath = None

    def __init__(self, datalayer, nameOrWorkflowFileOrJSONOrResource, serverName=None, caseType=CASETYPE_DECOMPOSED,
                 timeList=None):
        self.datalayer = datalayer
        if os.path.isdir(nameOrWorkflowFileOrJSONOrResource):
            self.casePath = nameOrWorkflowFileOrJSONOrResource
            self.simulationDocument = None
            simName = os.path.basename(self.casePath)
            groupName = simName.split("_")[0] if '_' in simName else ""
            self.simulationParams = {
                "simulationName": os.path.basename(self.casePath),
                "groupName": groupName,
                "simulationProperties": {}
            }

        else:
            simulationDocument = self.datalayer.getWorkflowDocumentFromDB(nameOrWorkflowFileOrJSONOrResource)
            if simDoc is not None:
                self.casePath = simulationDocument.resource
                self.simulationDocument = simulationDocument
                simulationProperties = self.simulationDocument['desc'].copy()
                self.simulationParams = {
                    "simulationName": simulationProperties['simulationName'],
                    "groupName": simulationProperties['groupName'],
                    "simulationProperties": simulationProperties
                }

            else:
                raise ValueError(
                    f"Simulation {nameOrWorkflowFileOrJSONOrResource} is not in the DB, and does not represent a valid case directory")

        self.filters = dict()
        self.timeList = timeList


    def loadToProject(self, datalayer, overwrite=False):
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
            logger.debug(f"Handling filter {filterName}")
            self._buildFiltersToUpdate(filterName=filterName, filterData=filterData, path="", overwrite=overwrite,
                                       datalayer=datalayer)

    def getData(self, filterName=None):
        """
            Returns the data for the requested nodes
            If None, return for all the nodes that were write=True.

        Parameters
        ----------
        vtkpiplines
        casePath
        serverName

        Returns
        -------

        """
        logger = get_classMethod_logger(self, "executePipelines")
        logger.info(f"Getting the resource case path ")


    def getAllFilters(self):
        """
            Flatten all the filter names to be in the list of filters to build.
        Returns
        -------

        """
        pass

    def getFiltersToUpdate(self, filterName, timeList=None):
        """
            This function checks which filters need to be updated and returns a list of their names.

            It checks if the fitlers data is in the db already with the requested timeList.

        Parameters
        ----------
        filterName
        filterData
        overwrite
        timeList

        Returns
        -------

        """
        logger = get_classMethod_logger(self, "_recurseNode")
        filterFullName = filter.fullName
        logger.info(f"Adding filter {filterFullName}")

        filtersToUpdate = []
        if filterData['write']:
            filterdescriptor = dict(filterName=filterFullName,
                                    pipeline=vtkPipelineJSON,
                                    simulationParameters=simulationparams)

            # check if it is already loaded.
            logger.debug("Checking if the filter is loaded in the project")
            docList = datalayer.getCacheDocuments(type=TYPE_VTK_FILTER, **dictToMongoQuery(filterdescriptor))
            logger.debug("Found document" if len(docList) > 0 else "no documents")

            if len(docList) == 0 or overwrite:
                logger.debug(f"Updating the filter data. Overwrite flag {overwrite}")

                if len(docList) > 0:
                    logger.debug("Checking if the time list is similar to the one in the DB")
                    doc = docList[0]
                    if doc.desc['timelist'] != timeList:
                        filtersToUpdate.append(filterFullName)
                else:

            else:
                logger.debug("Adding a new record")
                doc = datalayer.addCacheDocument(type=TYPE_VTK_FILTER,
                                                 dataFormat=currentdatatype,
                                                 resource=resource,
                                                 desc=qryparams)

    logger.debug("Processing the downstream filters. ")
    ds = filterData.get("downstream", {})
    for filterName, filterData in ds.items():
        path += "." + filterName
        self._buildFiltersToUpdate(filterName=filterName, filterData=filterData, path=path, overwrite=overwrite,
                                   datalayer=datalayer)


def execute(self, sourceOrName=None, timeList=None, tsBlockNum=50, overwrite=False, append=False):
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
    logger = get_classMethod_logger(self, "execute")

    if self.executionContext is None:
        err = f"The execution context is not set. use the method setExecutionContext to set it"
        logger.error(err)
        raise ValueError(err)

    self.executionContext.execute(sourceOrName=sourceOrName, timeList=timeList, tsBlockNum=tsBlockNum,
                                  overwrite=overwrite, append=append)


def toJSON(self):
    """
        Converts the pipeline to a VTK JSON of the executions.

        {
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
                    },
         }

    Returns
    -------

    """
    retDict = dict()
    for filterName, filterData in self.filters.items():
        retDict.update(filterData.toJSON())

    return dict(metadata=dict(timeList=self.timeList, casePath=self.casePath),
                pipelines=retDict)


def writeJSON(self, fileName):
    """
        Writes the pipeline to a JSON file.
    Parameters
    ----------
    fileName

    Returns
    -------

    """
    if '.json' not in fileName:
        fileName = f'{fileName}.json'

    with open(fileName, 'w') as outfile:
        json.dump(self.toJSON(), fileName, indent=4)
