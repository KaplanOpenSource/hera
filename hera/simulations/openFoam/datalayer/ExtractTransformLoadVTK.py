import pandas
import os
import json
import glob

import sys
from .. import DECOMPOSED_CASE, VTK_PIPELINE_FILTER

version = sys.version_info[0]
if version ==2 :
    import paraview.simple as pvsimple
    from .pvOpenFOAMBase import paraviewOpenFOAM
else: # version == 3:
    from hera.datalayer import Project, datatypes
    from ...hermesWorkflowToolkit import workflowToolkit

class VTKpipeline(object):
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
    _mainpath = None
    _casePath = None

    @property
    def pvOFBase(self):
        return self._pvOFBase

    @property
    def VTKpipelineJSON(self):
        return self._VTKpipelineJSON

    def __init__(self, pipelineJSON, casePath, caseType=DECOMPOSED_CASE, servername=None):
        """
            Initializes a VTK pipeline.

        :param pipelineJSON:
            JSON of the pipeline.
        :param name: a name for the files.

        casePath: str
                A full path to the case directory.

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
        self._VTKpipelineJSON = pipelineJSON
        outputdir = pipelineJSON["metadata"].get("CaseDirectory", os.getcwd())

        if version == 2:
            self._pvOFBase = paraviewOpenFOAM(casePath=casePath,
                                              caseType=caseType,
                                              servername=servername)

            self._pvOFBase.parquetdir = os.path.abspath(os.path.join(outputdir, casePath, 'parquet'))
            self._pvOFBase.netcdfdir = os.path.abspath(os.path.join(outputdir, casePath, "netcdf"))
        else:
            self._parquetdir = os.path.abspath(os.path.join(outputdir, casePath, 'parquet'))
            self._netcdfdir = os.path.abspath(os.path.join(outputdir, casePath, "netcdf"))

        self._mainpath = os.path.abspath(outputdir)
        self._casePath = os.path.abspath(casePath)

    def execute(self, source, tsBlockNum=100, overwrite=False):
        """
            Builds the pipeline from the JSON vtk.

        :param source:
            The source filter guiName that the pipeline will be build on.
        :param writeMetadata:
            if True, copies the json file to the results directory.

        """

        if version > 2:
            raise NotImplementedError("The execution of the pipeline must take place in the python 2.7 environment.")

        # build the pipeline.
        reader = pvsimple.FindSource(source)
        filterWrite = {}
        self._buildFilterLayer(father=reader, structureJson=self._VTKpipelineJSON["pipelines"], filterWrite=filterWrite)
        # Now execute the pipeline.
        timelist = self._VTKpipelineJSON["metadata"].get("timelist", None)

        if timelist is not None and (isinstance(timelist, str) or isinstance(timelist, unicode)):
            # a bit of parsing.
            readerTL = reader.TimestepValues
            BandA = [readerTL[0], readerTL[-1]]

            for i, val in enumerate(timelist.split(":")):
                BandA[i] = BandA[i] if len(val) == 0 else float(val)

            tl = pandas.Series(readerTL)
            timelist = tl[tl.between(*BandA)].values

        # else just take it from the json (it should be a list).

        # Get the mesh regions.
        if "MeshRegions" in self._VTKpipelineJSON["metadata"]:
            reader.MeshRegions = self._VTKpipelineJSON["metadata"]["MeshRegions"]

        for frmt, datasourceslist in filterWrite.items():
            writer = getattr(self._pvOFBase, "write_%s" % frmt)
            if writer is None:
                raise ValueError("The write %s is not found" % writer)
            writer(datasourcenamelist=datasourceslist,
                   timelist=timelist,
                   fieldnames=self._VTKpipelineJSON["metadata"].get('fields', None),
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
                pvalue = str(pvalue) if isinstance(pvalue, unicode) else pvalue  # python2, will be removed in python3.
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

    def loadToProject(self,projectName):
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
        if version < 3:
            raise NotImplementedError(
                "Loading the result of the pipeline to the project must take place in the conda-python 3+ environment")

        proj = Project(projectName=projectName)

        for filterName, filterData in self.VTKpipelineJSON["pipelines"].items():
            self._recurseNode(filterName=filterName, filterData=filterData, path="",project=proj)

    def _recurseNode(self, filterName, filterData, path, project):

        """
        This function is recursively passed on the pipeline tree from json file and loads the filters to the database.

        Parameter
        ----------

        Tree : a list of filters applied in the current tree before the current node
        nodeName : the current filter node name from the pipeline
        nodeData : the current filter node properties from the pipeline
        metadata : the metadata from the json file
        pipelines : the pipelines from the json file
        path : path to the main directory
        name : output folder/ hdf files string
        projectName : projectName string
        """
        formatName = filterData.get('write', None)
        if (formatName is not None) and (formatName != "None"):

            # Find the flow document.
            flowDoc = project.getSimulationsDocuments(resource=self._casePath)
            if len(flowDoc) == 0:
                raise ValueError(
                    "The simulation in path %s is not found!. Load it to the database first" % self._casePath)

            flowDoc = flowDoc[0]
            # Create the new description
            filterDesc = dict(flowDoc['desc'])

            simqry = {
                workflowToolkit.DESC_SIMULATIONNAME : filterDesc['simulationName'],
                workflowToolkit.DESC_GROUPNAME : filterDesc['groupName'],
                "filterName" : filterName
            }

            # check if it is already loaded.
            docList = project.getCacheDocuments(type=VTK_PIPELINE_FILTER,**simqry)

            if len(docList) ==0:
                filterDesc['filterName'] = filterName
                filterDesc['filterParameters'] = filterData.get('parameters',{})
                filterDesc['pipeline'] = self.VTKpipelineJSON
                filterDesc['filterPath'] = path

                if formatName == "netcdf":
                    resource = glob.glob(os.path.join(self._netcdfdir, "%s_*.nc" % filterName))
                    currentdatatype = datatypes.NETCDF_XARRAY
                elif formatName == "parquet":
                    resource = os.path.join(self._parquetdir, "%s.parquet" % filterName)
                    currentdatatype = datatypes.PARQUET
                else:
                    raise ValueError(
                        "Format type %s unknown for filter %s. Must be 'netcdf' or 'parquet' (case sensitive!)" % (
                        formatName, filterName))

                project.addCacheDocument(desc=filterDesc,type=VTK_PIPELINE_FILTER,resource=resource,dataFormat = currentdatatype)
            else:
                print("Filter %s for simulation %s in group %s already in the database" % (filterName,filterDesc['simulationName'],filterDesc['groupName']))

        # processing the
        ds = filterData.get("downstream", {})
        for filterName, filterData in ds.items():
            path += "." + filterName
            self._recurseNode(filterName=filterName, filterData=filterData, path=path,project=project)
