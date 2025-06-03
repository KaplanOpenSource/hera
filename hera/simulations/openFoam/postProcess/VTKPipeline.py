import json
import numpy
import os.path
import importlib

from hera import get_classMethod_logger
from hera.simulations.openFoam import CASETYPE_DECOMPOSED, CASETYPE_RECONSTRUCTED
from hera.utils import dictToMongoQuery

class VTKPipeLine:
    """
        Holds a vtk pipline. That is a structure of filters and their parameters
    """
    filters = None

    def __init__(self, datalayer, vtkPipeline=None):
        self.datalayer = datalayer
        self.filters = dict() if vtkPipeline is None else vtkPipeline

    @staticmethod
    def newVTKPipelineFilter(name, filterType, write=None, params=[]):
        """
            Initializes a new filter from the list.
        Parameters
        ----------
        name
        filterType
        write
        params : list
            This is because in VTK the order in which the parameters are applied changes the behaviour of the function.
        father : vtk Node.

        Returns
        -------

        """
        pipelineModule = importlib.import_module("hera.simulations.openFoam.postProcess.VTKPipeline")
        vtkFilterList = [x.split("_")[1] for x in dir(pipelineModule) if x.startswith("vtkFilter")]
        if filterType not in vtkFilterList:
            filterNameList = ",".join(vtkFilterList)
            raise ValueError(f"{filterType} is not Known. Must be one of {filterNameList}")
        return getattr(pipelineModule, f"vtkFilter_{filterType}")(name=name, filterType=filterType, write=write,
                                                                  params=params)

    def addFilter(self,name, filterType, write=None, params=[]):
        newFilter = VTKPipeLine.newVTKPipelineFilter(name=name,filterType=filterType,write=write,params=params)
        self.__setitem__(name,newFilter)

    def __setitem__(self, key, value):
        self.filters[key] = item

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
        pathList = item.split(".")
        val = self.filters[pathList[0]]
        for it in pathList[1:]:
            val = val[it]

        return val

    def registerPipeline(self, nameOrWorkflowFileOrJSONOrResource, serverName=None, caseType=CASETYPE_DECOMPOSED,
                         timeList=None):

        registeredVTKPipeLine(datalayer=self.datalayer,
                              vtkpipeline=self,
                              nameOrWorkflowFileOrJSONOrResource=nameOrWorkflowFileOrJSONOrResource,
                              serverName=serverName,
                              caseType=caseType,
                              timeList=timeList)

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

    def getAllFilters(self):
        """
            Return a list of the full name of all the filters.
        Parameters
        ----------
        self

        Returns
        -------

        """

        def recurseAllNames(fatherPath, filtersList):
            ret = []
            nextList = []
            if len(filtersList) > 0:
                currentName = filter if fatherPath is None else f"{fatherPath}.{filterName}"
                for filterName in filtersList:
                    ret.append(currentName)
                    curNode = self.__getitem__(currentName)
                    nextList += curNode.downstream
                ret += recurseAllNames(currentName, nextList)
            return ret

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

        self.datalayer = datalayer
        self.vtkpipeline = vtkpipeline

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

        self.pvOFBase = paraviewOpenFOAM(casePath=casePath,
                                         caseType=caseType,
                                         servername=serverName)

    def clearCache(self):
        pass

    def getData(self, filterName=None, timeList=None, meshRegions=None, nonRegularCase=True):
        """
            Returns the data of the vtkpipeline as a dict.
            The stuctucture is similar to that of a vtkpipline.
        Parameters
        ----------
        filterName
        timeList
        meshRegions
        nonRegularCase

        Returns
        -------

        """
        logger = get_classMethod_logger(self, "execute")
        ret = dict()

        # 1. Get the potential filters to process
        if filterName is None:
            requestedFiltersToProcess = self.vtkpipeline.getAllFilters()
        else:
            requestedFiltersToProcess = list(numpy.atleast_1d(filterName))
        logger.info(f"The potential filters are : {requestedFiltersToProcess}")

        # 2. Initialiaze the reader and the time list.

        if sourceOrName is None:
            reader = self.initializeReader(readerName="reader")
        elif isinstance(sourceOrName, str):
            logger.debug(f"Getting the reader {sourceOrName}")
            reader = pvsimple.FindSource(sourceOrName)
            if reader is None:
                reader = self.getReader(readerName="reader")
        else:
            logger.debug(f"Got reader object:  {sourceOrName}")
            # assume server is connected.
            reader = sourceOrName

        # Now execute the pipeline.
        if timeList is not None and isinstance(timeList, str):
            # a bit of parsing.
            readerTL = reader.TimestepValues
            BandA = [readerTL[0], readerTL[-1]]

            for i, val in enumerate(timeList.split(":")):
                BandA[i] = BandA[i] if len(val) == 0 else float(val)

            tl = pandas.Series(readerTL)
            timeList = tl[tl.between(*BandA)].values
        else:
            timelist = reader.TimestepValues

        # Get the mesh regions.
        reader.MeshRegions = meshRegions

        # 1. Check if the filters are already in the db
        qry = dict(simulations=self.simulationParams,
                   pipeline=self.vtkpipeline.toJSON())
        qry['simulations']['timeList'] = timelist
        qry['simulations']['MeshRegions'] = "" if meshRegions is None else meshRegions

        filtersToProcess = []
        for filterName in requestedFiltersToProcess:
            qry['filterName'] = filterName
            docList = self.datalayer.getSimulationDocuments(TYPE_VTK_FILTER, **dictToMongoQuery(qry))
            if docList > 0:
                ret[filterName] = docList[0].getData()
            else:
                filtersToProcess.append(filterName)

        # build the pipeline.
        filtersToCompute = self._buildFilterLayer(fatherName=None,
                                                  father=reader,
                                                  structureJson=self.vtkpipeline["filters"],
                                                  filtersToProcess=filtersToProcess)

        logger.info(f"Building the filter layer pipeline")

        # 3. Compute the values and save the parquet/zarr.
        methodName = "writeNonRegularCase" if nonRegularCase else "writeRegularCase"
        method = getattr(self.pvOFBase, methodName)

        method(datasourcenamelist=filtersToCompute,
               timeList=timeList,
               fieldnames=self._fieldNames,
               tsBlockNum=tsBlockNum,
               overwrite=overwrite,
               append=append)

        # 4. Update the DB.
        for filterName in filtersToProcess:
            print(f"adding {filterName} to DB")


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
        # self._readerName  = readerName
        reader = pvsimple.OpenFOAMReader(FileName="%s/tmp.foam" % self._casePath, CaseType=self._caseType,
                                         guiName=readerName)
        reader.MeshRegions.SelectAll()
        possibleRegions = list(reader.MeshRegions)
        reader.MeshRegions = ['internalMesh']
        if self._fieldNames is not None:
            reader.CellArrays = self._fieldNames
        reader.UpdatePipeline()

        return reader

    def _buildFilterLayer(self, fatherName, father, structureJson, filtersToProcess):
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
        ret = []

        if structureJson is not None:
            for filterGuiName in structureJson:
                paramPairList = structureJson[filterGuiName]['params']  # must be a list to enforce order in setting.
                filtertype = structureJson[filterGuiName]['type']
                filter = getattr(pvsimple, filtertype)(Input=father, guiName=filterGuiName)
                logger.debug(f"Adding filter {filterGuiName} of type {filtertype} to {fatherName}: {father}")

                for param, pvalue in paramPairList:
                    logger.debug(f"...Adding parameters {param} with value {pvalue}")
                    # pvalue = str(pvalue) if isinstance(pvalue, unicode) else pvalue  # python2, will be removed in python3.
                    paramnamelist = param.split(".")
                    paramobj = filter
                    for pname in paramnamelist[:-1]:
                        paramobj = getattr(paramobj, pname)
                    setattr(paramobj, paramnamelist[-1], pvalue)
                filter.UpdatePipeline()

                newFilterName = filterGuiName if fatherName is None else f"{fatherName}.{filterGuiName}"
                if newFilterName in filtersToProcess:
                    ret.append(filterGuiName)

                ret += self._buildFilterLayer(newFilterName, filter,
                                              structureJson[filterGuiName].get("downstream", None), filtersToProcess,
                                              filtersList)
        return ret


class VTKFilter:
    name = None
    filterType = None  # The type of the filter, clip/slice and ect.
    write = None  # true or false.
    params = None
    downstream = None
    father = None

    def __init__(self, name, filterType, write, params):
        """
            Init an abstract node.

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
        filterType  : str
            A VTK filter name : clip/cell centers and ect.
        write : bool
            Should we write this filter output?

        params : dict
            A list of key, value for the parameters of the dict.
            The order of the paramters is important, as it might change the behaviour of the object (like in the slice or clip filters).
        """
        self.name = name
        self.filterType = filterType
        self.write = shouldWrite
        self.params = params
        self.downstream = dict()

    @property
    def fullName(self):
        """
            Returns the full path of the filter from the father.
        Returns
        -------

        """

        def traverse(filter):
            if filter.father is not None:
                fatherName = filter.father.fullName()
            return [fatherName, self.name]

        return ".".join(traverse(self))

    def toJSON(self):
        """
            converts the node
        Returns
        -------

        """
        ret = dict()
        ret['filterType'] = self.filterType
        ret['write'] = True if self.write else False
        ret['params'] = self.params
        for filterName, dsFilter in self.downstream.items():
            ret['downstream'][filterName] = dsFilter.toJSON()

        return {self.name: ret}

    def __setitem__(self, key, value):
        self.downstream[key] = item

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
        pathList = item.split(".")
        try:
            val = self.downstream[pathList[0]]
            for it in pathList[1:]:
                val = val[it]
        except KeyError:
            raise KeyError(f"Filter {item} not found")

        return val

    def __getitem__(self, item):

        return None

    def append(self, newFilters):
        """
            Adds a new filter to the results.
        Parameters
        ----------
        newFilters : list
            A filter or a filter list.

        Returns
        -------

        """
        if isinstance(newFilter, list):
            self.downstream += newFilter
            for nf in newFilters:
                nf.father = self
        else:
            self.downstream.append(newFilter)
            newFilter.father = self


class vtkFilter_Slice(VTKFilter):

    def __init__(self, name, write, **kwargs):
        kwargs.setdefault("params", [])
        super().__init__(name=name, filterType="Slice", write=write, **kwargs)

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
        self.params = [("SliceType", "Plane"),
                       ("HyperTreeGridSlicer", "Plance"),
                       ("SliceOffsetValues", 0),
                       ("SliceType.Origin", [0, 0, 0]),
                       ("HyperTreeGridSlicer.Origin", [0, 0, 0])
                       ]


class vtkFilter_CellCenters(VTKFilter):
    """
        The Cell center filter.
    """

    def __init__(self, name, write, **kwargs):
        kwargs.setdefault("params", [])
        super().__init__(name=name, filterType="CellCenters", write=write, **kwargs)
