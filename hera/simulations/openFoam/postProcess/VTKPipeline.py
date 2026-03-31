import json
from typing import Literal
import numpy
import os.path
import importlib
import tqdm
import pandas
import pydoc
from hera import get_classMethod_logger
import hera.simulations.openFoam.postProcess.VTKPipeline as pipelineModule
from hera.simulations.openFoam import CASETYPE_DECOMPOSED, CASETYPE_RECONSTRUCTED, TYPE_VTK_FILTER
from hera.utils import dictToMongoQuery
from hera.simulations.openFoam.postProcess.pvOpenFOAMBase import paraviewOpenFOAM
import paraview.simple as pvsimple
from deprecated import deprecated
import os
import shutil

class VTKPipeLine:
    """
        Holds a vtk pipline. That is a structure of filters and their parameters
    """
    filters = None

    FILTER_CELLCENTERS = "CellCenters"
    FILTER_SLICE = "Slice"
    FILTER_PLOTOVERLINE = "PlotOverLine"
    FILTER_EXTRACTBLOCK = "ExtractBlock"
    FILTER_INTEGRATEVARIABLES= "IntegrateVariables"

    def __init__(self, datalayer, vtkPipeline=None):
        """Initialize the VTK pipeline with an optional existing filter dict."""
        self.datalayer = datalayer
        self.filters = dict() if vtkPipeline is None else vtkPipeline

    @staticmethod
    def newVTKPipelineFilter(name, filterType, write=True, params=None):
        """
            Initializes a new filter from the list.
        Parameters
        ----------
        name
        filterType
        write
        params : list
            This is because in VTK the order in which the parameters are applied changes the behaviour of the function.

            NOTE: Because we create the modules with getattr, using param=[] and param += ... is a problem.
                  Either set the param in the function, or change the creation mechanism to pydoc.

        father : vtk Node.

        Returns
        -------

        """
        if params is None:
            params = []
        vtkFilterList = [x.split("_")[1] for x in dir(pipelineModule) if x.startswith("vtkFilter")]
        if filterType not in vtkFilterList:
            filterNameList = ",".join(vtkFilterList)
            raise ValueError(f"{filterType} is not Known. Must be one of {filterNameList}")
        newClassPath = f"hera.simulations.openFoam.postProcess.VTKPipeline.vtkFilter_{filterType}"
        newCls = pydoc.locate(newClassPath)
        if newCls is None:
            raise RuntimeError(f"{filterType} does not exist")
        return newCls(name=name, write=write, params=params)

    def addFilter(self, name, filterType, write=True, params=None):
        """Create and add a new filter to the pipeline by type name."""
        if params is None:
            params = []
        newFilter = VTKPipeLine.newVTKPipelineFilter(name=name, filterType=filterType, write=write, params=params)
        self[name] = newFilter
        return newFilter

    def addFilterFromObj(self,newFilter):
        """
        Adds a filter to the pipeline using an already existent instance.

        Parameters

            filter(VTKFilter)
                an instance of a filter

        """
        self[newFilter.name] = newFilter

    @deprecated("Use addFilterFromObj")
    def addExistingFilter(self, newFilter):
        """
        Adds a filter to the pipeline using an already existent instance.

        Parameters
        
            filter(VTKFilter)
                an instance of a filter

        """
        self[newFilter.name] = newFilter

    def __setitem__(self, key, value):
        """Set a filter in the pipeline by key."""
        self.filters[key] = value

    def __getitem__(self, item):
        """
            Return the filter.

            Allows path syntax A.B.C
        Parameters
        ----------
        item : str
            The filter name or a path to the filter.

        Returns
        -------

        """
        try:
            pathList = item.split(".")
            val = self.filters[pathList[0]]
            for it in pathList[1:]:
                val = val[it]
        except KeyError:
            raise KeyError(f"The filter {item} is not found in the current pipeline")
        return val

    def registerPipeline(self, nameOrWorkflowFileOrJSONOrResource, serverName=None, caseType=CASETYPE_DECOMPOSED):
        """Bind this pipeline to a simulation case and return a registered pipeline."""

        return registeredVTKPipeLine(datalayer=self.datalayer,
                                     vtkpipeline=self,
                                     nameOrWorkflowFileOrJSONOrResource=nameOrWorkflowFileOrJSONOrResource,
                                     serverName=serverName,
                                     caseType=caseType)

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
        for _, filterData in self.filters.items():
            retDict.update(filterData.toJSON())

        return dict(filters=retDict)

    def allFilterNames(self,writeOnly=False):
        """
            Return a list of the full name of all the filters.
        Parameters
        ----------
        writeOnly: bool
            If true, return only the filter that are marked to write.

        Returns
        -------

        """

        def recurseAllNames(fatherPath, filtersList):
            """Recursively collect all filter names in the pipeline tree."""
            ret = []
            for filterName, filterObj in filtersList.items():
                # Build the dot-separated path for this filter.  Root-level filters
                # use their bare name; nested filters prepend their parent's path
                # (e.g. "ExtractBlock.CellCenters").
                currentName = filterName if fatherPath is None else f"{fatherPath}.{filterName}"
                # Include this filter in the result if we want all filters, or if
                # writeOnly is requested and this filter is marked for output.
                if (writeOnly and filterObj.write) or not writeOnly:
                        ret.append(currentName)
                # Recurse into downstream children so the full tree is traversed
                # depth-first, producing names in topological order.
                ret += recurseAllNames(currentName, filterObj.downstream)
            return ret

        # Start the recursion from the top-level filters with no parent path.
        return recurseAllNames(None, self.filters)


class registeredVTKPipeLine:
    """
        Represents binding of a vtk pipline to a case.
    """
    vtkpipeline = None
    datalayer = None
    casePath = None
    pvOFBase = None

    def __init__(self, datalayer, vtkpipeline, nameOrWorkflowFileOrJSONOrResource, serverName=None,
                 caseType=CASETYPE_DECOMPOSED):
        """Initialize by binding a VTK pipeline to a simulation case from DB or directory."""
        logger = get_classMethod_logger(self, "__init__")
        self.datalayer = datalayer
        self.vtkpipeline = vtkpipeline
        self.tsBlockNum = 50

        simulationDocumentList = []
        simulationDocumentList += self.datalayer.getWorkflowDocumentFromDB(nameOrWorkflowFileOrJSONOrResource)
        if len(simulationDocumentList) != 0:
            logger.info(f"found {nameOrWorkflowFileOrJSONOrResource} in database")
            simulationDocument = simulationDocumentList[0]
            self.casePath = simulationDocument.resource
            self.simulationDocument = simulationDocument
            simulationProperties = self.simulationDocument['desc'].copy()
            self.simulationParams = {
                "workflowName": simulationProperties['workflowName'],
                "groupName": simulationProperties['groupName'],
                "workflowParameters": simulationProperties['parameters']
            }
        elif os.path.isdir(nameOrWorkflowFileOrJSONOrResource):
            logger.info(f"{nameOrWorkflowFileOrJSONOrResource} is a directory, creating a psuedo document")
            self.casePath = nameOrWorkflowFileOrJSONOrResource
            self.simulationDocument = None
            simName = os.path.basename(self.casePath)
            groupName = simName.split("_")[0] if '_' in simName else ""
            self.simulationParams = {
                "workflowName": os.path.basename(self.casePath),
                "groupName": groupName,
                "workflowParameters": {}
            }

        else:
            raise ValueError(f"Simulation {nameOrWorkflowFileOrJSONOrResource} is not in the DB, and does not represent a valid case directory")

        self.pvOFBase = paraviewOpenFOAM(casePath=self.casePath,
                                         caseType=caseType,
                                         servername=serverName)

    def clearCache(self, regularMesh=None,filterName=None):
        """Remove cached filter output documents and files from disk."""
        # 1. Get the potential filters to process
        logger = get_classMethod_logger(self, "clearCache")
        if filterName is None:
            requestedFiltersToProcess = self.vtkpipeline.allFilterNames()
        else:
            requestedFiltersToProcess = list(numpy.atleast_1d(filterName))
        logger.info(f"Removing the cache for filters {requestedFiltersToProcess}")

        paramDict = dict()
        if regularMesh is not None:
            paramDict['regularMesh'] = regularMesh

        for filterName in requestedFiltersToProcess:
            logger.debug(f"Removing {filterName}")

            qry = self._buildFilterQuery(filterName=filterName,**paramDict)
            docList = self.datalayer.deleteSimulationsDocuments(type=TYPE_VTK_FILTER, **dictToMongoQuery(qry))
            logger.info(f"Found {len(docList)} documents to delete. ")
            for doc in docList:
                logger.debug(f"Deleting resource {doc['desc']['workflowName']} : {doc['desc']['pipeline']['filters']} ")
                outputFile = doc['resource']

                if os.path.exists(outputFile):
                    if os.path.isfile(outputFile):
                        os.remove(outputFile)
                    else:
                        shutil.rmtree(outputFile)

    def getData(self, regularMesh, filterName=None, timeList=None, latestTime=False, fieldNames=None, overwrite=False):
        """
        Return pipeline filter results as a dict keyed by filter name.

        Orchestrates: time resolution, cache lookup, ParaView execution, and
        DB persistence.  Each logical step is delegated to a private helper.

        Parameters
        ----------
        regularMesh : bool
            If True, output as zarr (xarray); otherwise parquet (pandas).
        filterName : str, list of str, or None
            Filter(s) to retrieve.  None means all write-enabled filters.
        timeList : None, str, or list
            Timesteps to process.  None = all; str = "start:end" range; list = explicit.
        latestTime : bool
            If True, restrict to only the last available timestep.
        fieldNames : list or None
            Optional field-name whitelist to limit reader I/O.
        overwrite : bool
            If True, recompute even when cached results exist.
        """
        logger = get_classMethod_logger(self, "getData")
        filext = self.getFilterOutputFileExt(regularMesh)

        # Step 1: Determine which filters to process.
        requestedFilters = self._resolveRequestedFilters(filterName)
        logger.info(f"The requested filters are : {requestedFilters}")

        # Step 2: Resolve the list of timesteps to compute.
        caseTimeList = self.datalayer.getTimeList(self.casePath)
        timeList = self._parseTimeList(timeList, caseTimeList)
        if latestTime:
            timeList = [timeList[-1]]
        logger.debug(f"Getting timeList {timeList}")

        # Step 3: Check cache — remove already-computed timesteps per filter.
        timeList, filtersToProcess, filtersOutputFilename, DBDocumentsDict = \
            self._filterCachedTimesteps(requestedFilters, timeList, regularMesh, filext, overwrite)
        logger.info(f"Computing filters {filtersToProcess}")

        # Step 4: Build ParaView pipeline and execute for uncached timesteps.
        if len(filtersToProcess) > 0:
            filtersToComputeDict = self._buildAndExecuteParaViewPipeline(
                filtersToProcess, filtersOutputFilename, timeList,
                fieldNames, overwrite, regularMesh)

        # Step 5: Persist newly computed timesteps into the DB cache.
        self._updateCacheDB(
            filtersToProcess, timeList, DBDocumentsDict,
            filtersOutputFilename if len(filtersToProcess) > 0 else {},
            regularMesh)

        # Step 6: Load and return cached data for every requested filter.
        ret = {}
        for fName in requestedFilters:
            ret[fName] = DBDocumentsDict[fName].getData()
        return ret

    # ------------------------------------------------------------------
    # Private helpers — each encapsulates one logical step of getData
    # ------------------------------------------------------------------

    def _resolveRequestedFilters(self, filterName):
        """Return the list of filter names to process.

        If *filterName* is None, return every filter in the pipeline that has
        ``write=True``.  Otherwise, coerce *filterName* (str or list) into a
        list.
        """
        if filterName is None:
            return self.vtkpipeline.allFilterNames(writeOnly=True)
        return list(numpy.atleast_1d(filterName))

    def _parseTimeList(self, timeList, caseTimeList):
        """Normalise the caller-supplied *timeList* into a concrete list.

        Handles three forms:
        - ``None``  — use every timestep available in the case.
        - ``str``   — a colon-separated ``"start:end"`` range where either
          bound may be omitted (defaults to first/last case time).
        - ``list``  — returned as-is.
        """
        if timeList is None:
            # No restriction — use the full case time range.
            return caseTimeList

        if isinstance(timeList, str):
            # Parse "start:end" range; missing sides default to case bounds.
            bounds = [caseTimeList[0], caseTimeList[-1]]
            for i, val in enumerate(timeList.split(":")):
                bounds[i] = bounds[i] if len(val) == 0 else float(val)
            tl = pandas.Series(caseTimeList)
            return tl[tl.between(*bounds)].values

        # Explicit list — pass through unchanged.
        return timeList

    def _filterCachedTimesteps(self, requestedFilters, timeList, regularMesh, filext, overwrite):
        """Check the DB cache and strip already-computed timesteps.

        For each requested filter, look up an existing cache document.  If one
        exists, remove its timesteps from *timeList* so only the delta needs
        to be computed (incremental strategy).

        Returns
        -------
        timeList : list
            Timesteps still requiring computation after cache subtraction.
        filtersToProcess : list[str]
            Subset of *requestedFilters* that actually need (re-)computation.
        filtersOutputFilename : dict[str, str]
            Mapping from filter name to its output file path on disk.
        DBDocumentsDict : dict
            Mapping from filter name to its existing cache document (if any).
        """
        logger = get_classMethod_logger(self, "_filterCachedTimesteps")
        filtersToProcess = []
        filtersOutputFilename = dict()
        DBDocumentsDict = dict()

        for fName in requestedFilters:
            qry = self._buildFilterQuery(filterName=fName, regularMesh=regularMesh)
            docList = self.datalayer.getCacheDocuments(type=TYPE_VTK_FILTER, **dictToMongoQuery(qry))

            if len(docList) > 0:
                # Cache hit — reuse output path and subtract known timesteps.
                logger.info("Found existing filter output in cache")
                cached_filter = docList[0]
                filtersOutputFilename[fName] = cached_filter.resource
                dbTimeList = cached_filter['desc']['simulation']['timeList']
                timeList = [ts for ts in timeList if ts not in dbTimeList]
                DBDocumentsDict[fName] = cached_filter
            else:
                # Cache miss — generate a fresh output file path.
                outputFilePath = self.getFilterOutputFilePath(fName, filext, generate_new=True)
                filtersOutputFilename[fName] = outputFilePath

            # A filter needs computation if: forced overwrite, no cache, or
            # there are timesteps not yet in the cache.
            logger.debug(
                "Compute the filter if you need to overwrite the results, it is not in the DB, or there are times not in the DB")
            if overwrite or len(docList) == 0 or len(timeList) > 0:
                logger.debug(f"{fName} added to process because overwrite=True or filter not in DB")
                filtersToProcess.append(fName)

        return timeList, filtersToProcess, filtersOutputFilename, DBDocumentsDict

    def _buildAndExecuteParaViewPipeline(self, filtersToProcess, filtersOutputFilename,
                                          timeList, fieldNames, overwrite, regularMesh):
        """Instantiate the ParaView filter tree and run it over *timeList*.

        Steps:
        1. Create the OpenFOAM reader.
        2. Optionally restrict the reader to *fieldNames*.
        3. Build the full filter tree from the pipeline JSON.
        4. Execute and write results to disk (parquet or zarr).
        5. Clean up ParaView proxies to free server-side memory.

        Returns the dict mapping filter names to their output file paths.
        """
        logger = get_classMethod_logger(self, "_buildAndExecuteParaViewPipeline")
        filtersToComputeDict = dict()

        logger.info(f"Building the vtk objects from the JSON")
        # The reader is the root of the ParaView pipeline graph.
        reader = self.pvOFBase.initializeReader(readerName="reader")

        # Restrict to specific field arrays to reduce memory and I/O.
        if fieldNames is not None:
            reader.CellArrays = fieldNames

        # Recursively create ParaView filter proxies from the pipeline JSON.
        filtersToCompute = self._buildFilterLayer(
            fatherName=None, father=reader,
            structureJson=self.vtkpipeline.toJSON()['filters'])
        logger.info(f"Added all filters to the layer. Computing filters {filtersToCompute}")

        # Map each filter that needs computation to its output file path.
        for fName in filtersToProcess:
            logger.debug(f"\t{fName} will be saved in {filtersOutputFilename[fName]}")
            filtersToComputeDict[fName] = filtersOutputFilename[fName]

        # Execute the pipeline over the remaining timesteps and write output.
        self.pvOFBase.writeCase(filtersDict=filtersToComputeDict,
                                timeList=timeList,
                                fieldnames=fieldNames,
                                tsBlockNum=self.tsBlockNum,
                                overwrite=overwrite, regularMesh=regularMesh)

        # Clean up all ParaView sources/filters to free server-side memory.
        for name, proxy in list(pvsimple.GetSources().items()):
            logger.debug(f"Deleting source {name}")
            pvsimple.Delete(proxy)

        return filtersToComputeDict

    def _updateCacheDB(self, filtersToProcess, timeList, DBDocumentsDict,
                        filtersToComputeDict, regularMesh):
        """Merge newly computed timesteps into the DB cache.

        For filters that already have a cache document, append the new
        timesteps and save.  For first-time filters, create a brand-new
        cache document pointing to the output file on disk.
        """
        logger = get_classMethod_logger(self, "_updateCacheDB")

        for fName in filtersToProcess:
            logger.debug(f"Updating times {timeList} to filter {fName}")

            if fName in DBDocumentsDict:
                # Incremental update: merge new timesteps with existing ones.
                doc = DBDocumentsDict[fName]
                fullTime = sorted(timeList + doc['desc']['simulation']['timeList'])
                doc.desc['simulation']['timeList'] = fullTime
                doc.save()
            else:
                # First computation — create a new cache record in the DB.
                logger.debug("...Adding a new record to the DB")
                recordData = self._buildFilterQuery(filterName=fName, regularMesh=regularMesh)
                recordData['simulation']['timeList'] = timeList
                dataFormat = self.datalayer.datatypes.ZARR_XARRAY if regularMesh else self.datalayer.datatypes.PARQUET
                doc = self.datalayer.addCacheDocument(
                    dataFormat=dataFormat,
                    resource=os.path.abspath(filtersToComputeDict[fName]),
                    type=TYPE_VTK_FILTER,
                    desc=recordData)
                DBDocumentsDict[fName] = doc

            logger.debug(f"Reading filter {fName} data")

    def getFilterOutputFileExt(self, regularMesh):
        """Return the file extension based on mesh regularity."""
        return "zarr" if regularMesh else "parquet"

    def getFilterOutputFilePath(self, filterName, filext, generate_new=True):
        """Generate or retrieve the output file path for a filter."""
        workflow_name = self.simulationParams['workflowName']
        counter_name = f"{workflow_name}_{filterName}_counter" # more consistent to have a counter per filter name
        curr_number = self.datalayer.getCounter(counter_name)
        if curr_number is None and not generate_new: # means there are none
            return None
        counter = self.datalayer.getCounterAndAdd(counter_name) if generate_new else curr_number
        outputFileName = f"{filterName.replace('.','_')}_{counter}.{filext}"
        outputFilePath = os.path.join(os.path.abspath(self.casePath), "vtkpipelinedata", outputFileName)
        return outputFilePath

    def getRegularData(self, filterName=None, timeList=None,  fieldNames=None, overwrite=False):
        """Retrieve pipeline data as regular mesh (xarray/zarr) format."""
        return self.getData(regularMesh=True, filterName=filterName, timeList=timeList,
                            fieldNames=fieldNames,
                            overwrite=overwrite)

    def getNonRegularData(self, filterName=None, timeList=None,  fieldNames=None, overwrite=False):
        """Retrieve pipeline data as non-regular mesh (pandas/parquet) format."""
        return self.getData(regularMesh=False, filterName=filterName, timeList=timeList,
                            fieldNames=fieldNames,
                            overwrite=overwrite)

    def _buildFilterQuery(self, filterName, regularMesh=None):
        """Build a database query dict for a specific filter."""
        qry = dict(simulation=self.simulationParams,
                   pipeline=self.vtkpipeline.toJSON())

        if regularMesh is not None:
            qry['simulation']['regularMesh'] = regularMesh

        qry['filterName'] = filterName
        return qry

    def _buildFilterLayer(self, fatherName, father, structureJson):
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
        logger.debug(f"Initialized logger {logger}")
        logger.info(f"building Filter layer {json.dumps(structureJson, indent=4)}")
        # Accumulates the fully-qualified names (e.g. "ExtractBlock.CellCenters") of
        # every filter created during this recursive traversal, so the caller knows
        # which ParaView sources exist and can be looked up via pvsimple.FindSource().
        ret = []


        if structureJson is not None:
            # Iterate over each sibling filter at this level of the tree.
            for filterGuiName in structureJson:
                # params is deliberately a *list* of (key, value) pairs rather than a
                # dict, because the order in which VTK filter properties are set can
                # change behaviour (e.g. setting SliceType before SliceType.Origin).
                paramPairList = structureJson[filterGuiName]['params']  # must be a list to enforce order in setting.
                filtertype = structureJson[filterGuiName]['filterType']
                # Build the dot-separated full name: root filters use their own name,
                # child filters prepend their parent's name (e.g. "parent.child").
                newFilterName = filterGuiName if fatherName is None else f"{fatherName}.{filterGuiName}"
                # Instantiate the actual ParaView filter object. `father` is the
                # upstream pipeline source (reader or another filter) that feeds data
                # into this filter.
                filter = getattr(pvsimple, filtertype)(Input=father, guiName=newFilterName)
                logger.debug(
                    f"Adding filter {filterGuiName} of type {filtertype} to {'Reader' if fatherName is None else fatherName}")

                # Apply each parameter in order.  Parameters may use dot-notation
                # (e.g. "SliceType.Origin") to reach nested sub-proxy attributes.
                for param, pvalue in paramPairList:
                    logger.debug(f"...Adding parameters {param} with value {pvalue}")
                    # pvalue = str(pvalue) if isinstance(pvalue, unicode) else pvalue  # python2, will be removed in python3.
                    # Split the param name on "." to traverse nested proxy attributes.
                    # For "SliceType.Origin", we first resolve filter.SliceType, then
                    # set the "Origin" attribute on that sub-proxy.
                    paramnamelist = param.split(".")
                    paramobj = filter
                    for pname in paramnamelist[:-1]:
                        paramobj = getattr(paramobj, pname)
                    setattr(paramobj, paramnamelist[-1], pvalue)
                # Force the filter to execute so downstream filters see updated data.
                filter.UpdatePipeline()
                logger.debug(f"Filter {newFilterName} added to the pipeline. Now adding its downstream filters.")
                ret.append(newFilterName)
                # Recurse into the downstream children of this filter, building the
                # tree depth-first.  The returned names are appended so the final list
                # is in creation (topological) order.
                ret += self._buildFilterLayer(newFilterName, filter,
                                              structureJson[filterGuiName].get("downstream", None))

        return ret


class VTKFilter:
    """Base class representing a single VTK filter node in a pipeline tree."""

    name = None
    filterType = None  # The type of the filter, clip/slice and ect.
    write = None  # true or false.
    params = None
    downstream = None
    father = None

    def __init__(self, name, filterType, write, params):
        """Init an abstract node.
        
        "filterName": {
            "type": The type of the filter.(clip, slice, ...).
            "write": True/False
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
            name (str): Name of the filter
            filterType (str): A VTK filter name : clip/cell centers and ect.
            write (bool): Should we write this filter output?
            params (dict): A list of key, value for the parameters of the dict. The order of the paramters is important, as it might change the behaviour of the object (like in the slice or clip filters).
        """
        self.name = name
        self.filterType = filterType
        self.write = write
        self.params = params
        self.downstream = dict()

    @property
    def fullName(self):
        """
        
        Returns
        -------
        full path of the filter from the father.
        """

        def traverse(filter):
            """Traverse up the filter tree to build the full path."""
            if filter.father is not None:
                fatherName = filter.father.fullName()
            return [fatherName, self.name]

        return ".".join(traverse(self))

    def toJSON(self):
        """converts the node
        
        Returns
        -------
        Jsonified filter
        """
        ret = dict()
        ret['filterType'] = self.filterType
        ret['write'] = True if self.write else False
        ret['params'] = self.params
        ret['downstream'] = {}
        for filterName, dsFilter in self.downstream.items():
            ret['downstream'][filterName] = dsFilter.toJSON()[filterName]

        return {self.name: ret}

    def __setitem__(self, key, value):
        """Add a downstream filter by key."""
        self.downstream[key] = value

    def __getitem__(self, item):
        """Return the filter.
        
        Allows path syntax A.B.C

        Parameters
        ----------
        item : str
            The filter name or a path to the filter.

        Returns
        -------

        """
        pathList = item.split(".")
        try:
            val = self.downstream[pathList[0]]
            for it in pathList[1:]:
                val = val[it]
        except KeyError:
            raise KeyError(f"Filter {item} not found")
        return val

    def addFilter(self, name, filterType, write=None, params=None):
        """Create and add a new downstream filter by type name."""
        if params is None:
            params=[] # The dynamic initialization causes the parameter param to get old values.
        newFilter = VTKPipeLine.newVTKPipelineFilter(name=name, filterType=filterType, write=write, params=params)
        self.downstream[name] = newFilter
        return newFilter

    def addFilterFromObj(self,newFilter):
        """
        Adds a filter to the pipeline using an already existent instance.

        Parameters

            filter(VTKFilter)
                an instance of a filter

        """
        self.downstream[newFilter.name] = newFilter


    def set_param(self, param_name:str, new_val):
        """Set or update a parameter value in the params list."""
        # change param value in list, keeping at the same index
        self.params= [(param_name, new_val) if param[0] == param_name else param for param in self.params]
        # in case iteration didn't pass it we append it
        if (param_name, new_val)  not in self.params:
            self.params .append((param_name,new_val))
            
    def enforce_param(self, param_name:str, enforced_param):
        """In case param isn't initialized we set it to the enforced param

        Raises:
            RuntimeError: Raises in case param is already set to something that's not the enforced param 
        """
        if all(param[0] != param_name for param in self.params):
            self.params.append((param_name, enforced_param))
        elif (param_name, enforced_param) not in self.params:
            raise RuntimeError(f"{param_name} must not be different than {enforced_param}")



class vtkFilter_Slice(VTKFilter):
    """VTK Slice filter for cutting data with a plane."""

    def __init__(self, name, write, **kwargs):
        """Initialize a Slice filter."""
        kwargs.setdefault("params",[])
        super().__init__(name=name, filterType="Slice", write=write, **kwargs)
                
    def setPlaneOrigin(self, origin:list[3]):
        """Setting plane slice origin

        Parameters:
            origin (list[3]): a 3-tuple of the origin
        """
        self.enforce_param("SliceType", "Plane")
        self.set_param("SliceType.Origin", origin)
            
    def setPlaneNormal(self, normal:list[3]):
        """Setting plane slice normal

        Parameters:
            normal (list[3]): a 3-tuple of the normal
        """
        self.enforce_param("SliceType", "Plane")
        self.set_param("SliceType.Normal", normal)

class vtkFilter_PlotOverLine(VTKFilter):
    """VTK PlotOverLine filter for sampling data along a line segment."""

    SAMPLE_UNIFORMLY = 'Sample Uniformly'
    SAMPLE_AT_CELL_BOUNDARIES = 'Sample At Cell Boundaries'
    SAMPLE_AT_SEGMENT_CENTERS = 'Sample At Segment Centers'

    def __init__(self, name, write, **kwargs):
        """Initialize a PlotOverLine filter."""
        kwargs.setdefault("params",[])
        super().__init__(name=name, filterType="PlotOverLine", write=write, **kwargs)

    def setSamplePattern(self, pattern: Literal['Sample Uniformly', 'Sample At Cell Boundaries', 'Sample At Segment Centers']):
        """Set the sampling pattern for the line plot."""
        self.set_param("SamplingPattern",  pattern)

    def setUniformSampleResolution(self, res: int):
        """Set the resolution for uniform sampling along the line."""
        self.enforce_param("SamplingPattern", 'Sample Uniformly')
        self.set_param("Resolution", res)

    def setPoints(self, point1:list[3],point2:list[3]):
        """Setting PlotOverLine point1(start point) and point2(end point)
        
        Parameters:
            point1 (list[3]): a 3-tuple of the start point
            point2 (list[3]): a 3-tuple of the end point       
        """
        self.set_param("Point1", point1)
        self.set_param("Point2", point2)

class vtkFilter_CellCenters(VTKFilter):
    """
        The Cell center filter.
    """

    def __init__(self, name, write, **kwargs):
        """Initialize a CellCenters filter."""
        kwargs.setdefault("params", [])
        super().__init__(name=name, filterType="CellCenters", write=write, **kwargs)


class vtkFilter_ExtractBlock(VTKFilter):
    """VTK ExtractBlock filter for extracting specific mesh regions."""

    def __init__(self, name, write, patchList=[],internalMesh=False,params=None):
        """Initialize an ExtractBlock filter with selected patches."""
        if params is None:
            params = []
        selectors = [f'/Root/boundary/{patchName}' for patchName in patchList]
        if internalMesh:
            selectors += ['/Root/internalMesh']
        final_params =  [("Selectors",selectors)] + params
        super().__init__(name=name, filterType="ExtractBlock", write=write,params=final_params)

    def setRegionsToExtract(self,patchList=[], internalMesh=True):
        """Set the patch list and internal mesh flag for block extraction."""
        selectors = [f'/Root/boundary/{patchName}' for patchName in patchList]
        if internalMesh:
            selectors += ['/Root/internalMesh']

        self.set_param("Selectors", selectors)


class vtkFilter_DescriptiveStatistics(VTKFilter):
    """VTK DescriptiveStatistics filter for computing summary statistics."""

    def __init__(self, name, write, **kwargs):
        """Initialize a DescriptiveStatistics filter."""
        kwargs.setdefault("params",[])
        super().__init__(name=name, filterType="DescriptiveStatistics", write=write, **kwargs)

    def setVariablesOfInterest(self, variables=[]):
        """Set the list of variables to compute statistics for."""
        self.set_param("VariablesofInterest", variables)


class vtkFilter_IntegrateVariables(VTKFilter):
    """VTK IntegrateVariables filter for integrating fields over a mesh."""

    def __init__(self, name, write, **kwargs):
        """Initialize an IntegrateVariables filter."""
        super().__init__(name=name, filterType="IntegrateVariables", write=write, **kwargs)
