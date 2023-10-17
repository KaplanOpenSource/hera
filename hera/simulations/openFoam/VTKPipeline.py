import pdb

import pandas
import os
import json
import glob

import sys
from . import DECOMPOSED_CASE,RECONSTRUCTED_CASE,TYPE_VTK_FILTER
from ...utils import loadJSON
from ...utils.logging import helpers as hera_logging
from ...datalayer import datatypes
paraviewExists = False
try:
    import paraview.simple as pvsimple
    from simulations.openFoam.datalayer.pvOpenFOAMBase import paraviewOpenFOAM

    paraviewExists = True
except ImportError:
    print("paraview not Found. Cannot execute the VTK pipeline.")

class VTKpipeline:
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

        if self._simulationDocument is None:
            if self.datalayer is not None:
                self._simulationDocument = self.datalayer.getWorkflowDocumentFromDB(self._nameOrWorkflowFileOrJSONOrResource)

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
        self.logger = hera_logging.get_logger(self)
        self._VTKpipelineJSON = loadJSON(pipelineJSON)
        self._datalayer = datalayer
        self._serverName = serverName
        self._caseType = caseType
        self._nameOrWorkflowFileOrJSONOrResource = nameOrWorkflowFileOrJSONOrResource

        if os.path.isdir(nameOrWorkflowFileOrJSONOrResource):
            casePath = nameOrWorkflowFileOrJSONOrResource
        else:
            simDoc = self.simulationDocument
            if simDoc is not None:
                casePath = simDoc.resource
            else:
                raise ValueError(f"Simulation {nameOrWorkflowFileOrJSONOrResource} is not in the DB, and does not represent a valid case directory")


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

    def _initializeReader(self, readerName):
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
        self._reader = pvsimple.OpenFOAMReader(FileName="%s/tmp.foam" % self._casePath, CaseType=self._caseType, guiName=readerName)
        self._reader.MeshRegions.SelectAll()
        self._possibleRegions = list(self._reader.MeshRegions)
        self._reader.MeshRegions = ['internalMesh']
        if self._fieldNames is not None:
            self._reader.CellArrays = self._fieldNames
        self._reader.UpdatePipeline()

        return self._reader


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
        self.logger.info(f"Executing pipeline")
        if sourceOrName is None:
            reader = self._initializeReader(readerName="reader")
        elif isinstance(sourceOrName,str):
            self.logger.debug(f"Getting the reader {sourceOrName}")
            reader = pvsimple.FindSource(sourceOrName)
            if reader is None:
                reader = self._initializeReader(readerName="reader")
        else:
            self.logger.debug(f"Gor reader object:  {sourceOrName}")
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
        self.logger.info(f"building Filter layer {structureJson}")
        if structureJson is None:
            return
        for filterGuiName in structureJson:

            paramPairList = structureJson[filterGuiName]['params']  # must be a list to enforce order in setting.
            filtertype = structureJson[filterGuiName]['type']
            filter = getattr(pvsimple, filtertype)(Input=father, guiName=filterGuiName)
            self.logger.execution(f"Adding filter {filterGuiName} of type {filtertype} to {father}")
            for param, pvalue in paramPairList:
                self.logger.debug(f"...Adding parameters {param} with value {pvalue}")
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
        self.logger.info(f"Adding filter {filterName} in {path}")
        formatName = filterData.get('write', None)
        if (formatName is not None) and (formatName != "None"):
            self.logger.debug("Getting the simulation flow document. Use the path as identifier ")
            # Find the flow document.
            flowDoc = self.datalayer.getWorkflowDocumentFromDB(self._casePath)

            if len(flowDoc) == 0:
                self.logger.error(f"P{self._casePath} is not found in the project. Load it using the hera-OF-workflow add ")
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
            self.logger.debug("Checking if the filter is loaded in the project")
            docList = self.datalayer.getCacheDocuments(type=TYPE_VTK_FILTER,**simqry)

            if len(docList) ==0 or overwrite:
                self.logger.debug(f"Updating the filter data. Overwrite flag {overwrite}")

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
                    self.logger.debug("Updating the existing data")
                    doc = docList[0]
                    doc['desc'] = filterDesc
                    doc['resource'] = resource
                    doc['dataFormat'] =currentdatatype
                    doc.update()
                else:
                    self.logger.debug("Adding a new record")
                    self.datalayer.addCacheDocument(desc=filterDesc,type=TYPE_VTK_FILTER,resource=resource,dataFormat = currentdatatype)
            else:
                self.logger.warning("Filter %s for simulation %s in group %s already in the database" % (filterName,filterDesc['simulationName'],filterDesc['groupName']))

        self.logger.execution("Processing the downstream filters. ")
        ds = filterData.get("downstream", {})
        for filterName, filterData in ds.items():
            path += "." + filterName
            self._recurseNode(filterName=filterName, filterData=filterData, path=path,overwrite=overwrite)

