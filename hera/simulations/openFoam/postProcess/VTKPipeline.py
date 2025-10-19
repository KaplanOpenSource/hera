import json
import numpy
import os.path
import importlib
import tqdm

from hera import get_classMethod_logger
from hera.simulations.openFoam import CASETYPE_DECOMPOSED, CASETYPE_RECONSTRUCTED, TYPE_VTK_FILTER
from hera.utils import dictToMongoQuery
from hera.simulations.openFoam.postProcess.pvOpenFOAMBase import paraviewOpenFOAM
import paraview.simple as pvsimple
import os
import shutil

class VTKPipeLine:
    """
        Holds a vtk pipline. That is a structure of filters and their parameters
    """
    filters = None

    FILTER_CELLCENTERS = "CellCenters"
    FILTER_SLICE = "Slice"

    def __init__(self, datalayer, vtkPipeline=None):
        self.datalayer = datalayer
        self.filters = dict() if vtkPipeline is None else vtkPipeline

    @staticmethod
    def newVTKPipelineFilter(name, filterType, write=True, params=[]):
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
        return getattr(pipelineModule, f"vtkFilter_{filterType}")(name=name, write=write, params=params)

    def addFilter(self, name, filterType, write=True, params=[]):
        newFilter = VTKPipeLine.newVTKPipelineFilter(name=name, filterType=filterType, write=write, params=params)
        self.__setitem__(name, newFilter)

    def __setitem__(self, key, value):
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
        for filterName, filterData in self.filters.items():
            retDict.update(filterData.toJSON())

        return dict(filters=retDict)

    def allFilterNames(self):
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
            for filterName, filterObj in filtersList.items():
                currentName = filterName if fatherPath is None else f"{fatherPath}.{filterName}"
                ret.append(currentName)
                ret += recurseAllNames(currentName, filterObj.downstream)
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
        self.tsBlockNum = 50

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
            if simulationDocument is not None:
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

        self.pvOFBase = paraviewOpenFOAM(casePath=self.casePath,
                                         caseType=caseType,
                                         servername=serverName)

    def clearCache(self, regularMesh,filterName=None):
        # 1. Get the potential filters to process
        logger = get_classMethod_logger(self, "clearCache")
        if filterName is None:
            requestedFiltersToProcess = self.vtkpipeline.allFilterNames()
        else:
            requestedFiltersToProcess = list(numpy.atleast_1d(filterName))
        logger.info(f"Removing the cache for filters {requestedFiltersToProcess}")

        for filterName in requestedFiltersToProcess:
            logger.debug(f"Removing {filterName}")
            qry = self._buildFilterQuery(filterName=filterName, meshRegions=None,regularMesh=regularMesh)
            docList = self.datalayer.deleteSimulationsDocuments(type=TYPE_VTK_FILTER, **dictToMongoQuery(qry))
            logger.info(f"Found {len(docList)} documents to delete. ")
            for doc in docList:
                logger.debug(f"Deleting resource {doc['desc']['simulation']['simulationName']} : {doc['desc']['pipeline']['filters']} ")
                outputFile = doc['resource']

                if os.path.exists(outputFile):
                    if os.path.isfile(outputFile):
                        os.remove(outputFile)
                    else:
                        shutil.rmtree(outputFile)

    def getData(self, regularMesh, filterName=None, timeList=None, meshRegions=None , fieldNames=None,overwrite=False):
        """
            Returns the data of the vtkpipeline as a dict.
            The stuctucture is similar to that of a vtkpipline.
        Parameters
        ----------
        filterName : str, list of str, None
            The name of the filter to get, a list of filters or get all the filters that write=True in the pipeline.
        timeList : None, list
            The list of timestep
        meshRegions
        nonRegularCase
        sourceOrName  : str
            A name of a filter that is already in the pipeline that will be used as a base for the pipeline.
            If NOne, then initialize a reader.

        Returns
        -------

        """
        logger = get_classMethod_logger(self, "execute")
        ret = dict()
        filext = "zarr" if regularMesh else "parquet"

        # 1. Get the potential filters to process
        if filterName is None:
            requestedFiltersToProcess = self.vtkpipeline.allFilterNames()
        else:
            requestedFiltersToProcess = list(numpy.atleast_1d(filterName))
        logger.info(f"The requested filters are : {requestedFiltersToProcess}")

        # 2. Initialiaze the reader and the time list.
        reader = self.pvOFBase.initializeReader(readerName="reader")
        if fieldNames is not None:
            reader.CellArrays = fieldNames

        # Now execute the pipeline.
        if timeList is not None:
            if isinstance(timeList, str):
                # a bit of parsing.
                readerTL = reader.TimestepValues
                BandA = [readerTL[0], readerTL[-1]]

                for i, val in enumerate(timeList.split(":")):
                    BandA[i] = BandA[i] if len(val) == 0 else float(val)

                tl = pandas.Series(readerTL)
                timeList = tl[tl.between(*BandA)].values
        else:
            timeList = reader.TimestepValues

        logger.debug(f"Getting timeList {timeList}")

        # Get the mesh regions.
        if meshRegions is not None:
            reader.MeshRegions = meshRegions

        filtersToProcess = []
        filtersOutputFilename = dict()
        DBDocumentsDict = dict()

        for filterName in requestedFiltersToProcess:
            qry = self._buildFilterQuery(filterName=filterName, meshRegions=meshRegions,regularMesh=regularMesh)
            docList = self.datalayer.getSimulationsDocuments(type=TYPE_VTK_FILTER, **dictToMongoQuery(qry))
            if len(docList) > 0:
                filtersOutputFilename[filterName] = docList[0].resource
                dbTimeList = docList[0]['desc']['simulation']['timeList']

                # Filtering the timesteps.
                timeList = [ts for ts in timeList if ts not in dbTimeList]
                DBDocumentsDict[filterName] = docList[0]

            else:
                # adding counter to prevent overwriting similar filter from a different pipeline.
                countr = self.datalayer.getCounterAndAdd("OpenFOAMData")
                outputFileName = f"{filterName}_{countr}.{filext}"
                filtersOutputFilename[filterName] = os.path.join(os.path.abspath(self.casePath), "vtkpipelinedata", outputFileName)

            logger.debug(
                "Compute the filter if you need to overwrite the results, it is not in the DB, or there are times not in the DB")
            if overwrite or len(docList) == 0 or len(timeList) > 0:
                logger.debug(f"{filterName} added to process because overwrite=True or filter not in DB")
                filtersToProcess.append(filterName)

        logger.info(f"Computing filters {filtersToProcess}")
        # Compute the filters.
        if len(filtersToProcess) > 0:
            filtersToComputeDict = dict()  #
            logger.info(f"Building the vtk objects from the JSON")
            # Remember that the buildFilterlayer adds the real objects to the pipeline,
            # and just return the guinames that can be used to retrieve the object with findSource.
            # Hence it is a list of strings.
            filtersToCompute = self._buildFilterLayer(fatherName=None,
                                                      father=reader,
                                                      structureJson=self.vtkpipeline.toJSON()['filters'])
            logger.info(f"Added all filters to the layer. Computing filters {filtersToCompute}")

            # Build the ourput file name.
            for filterName in filtersToCompute:
                logger.debug(f"\t{filterName} will saved in {filtersOutputFilename[filterName]}")
                filtersToComputeDict[filterName] = filtersOutputFilename[filterName]

            # Compute the values and save the parquet/zarr.
            self.pvOFBase.writeCase(filtersDict=filtersToComputeDict,
                                    timeList=timeList,
                                    fieldnames=fieldNames,
                                    tsBlockNum=self.tsBlockNum,
                                    overwrite=overwrite, regularMesh=regularMesh)

        # 4. Update the DB.
        for filterName in filtersToProcess:
            logger.debug(f"Updating times {timeList} to filter {filterName}")
            if filterName in DBDocumentsDict:
                doc = DBDocumentsDict[filterName]
                fullTime = sorted(timeList + doc['desc']['simulation']['timeList'])
                doc.desc['simulation']['timeList'] = fullTime
                doc.save()
            else:
                logger.debug("...Adding a new record to the DB")

                recordData = self._buildFilterQuery(filterName=filterName, meshRegions=meshRegions,regularMesh=regularMesh)
                recordData['simulation']['timeList'] = timeList
                dataFormat = self.datalayer.datatypes.ZARR_XARRAY if regularMesh else self.datalayer.datatypes.PARQUET
                doc = self.datalayer.addSimulationsDocument(dataFormat=dataFormat,
                                                            resource=os.path.abspath(filtersToComputeDict[filterName]),
                                                            type=TYPE_VTK_FILTER,
                                                            desc=recordData)

            logger.debug(f"Reading filter {filterName} data")

        for filterName in requestedFiltersToProcess:
            qry = self._buildFilterQuery(filterName=filterName, meshRegions=meshRegions,regularMesh=regularMesh)
            docList = self.datalayer.getSimulationsDocuments(type=TYPE_VTK_FILTER, **dictToMongoQuery(qry))
            ret[filterName] = docList[0].getData()

        return ret

    def getRegularData(self, filterName=None, timeList=None, meshRegions=None, fieldNames=None, overwrite=False):
        return self.getData(regularMesh=True, filterName=filterName, timeList=timeList, meshRegions=meshRegions,
                            fieldNames=fieldNames,
                            overwrite=overwrite)

    def getNonRegularData(self, filterName=None, timeList=None, meshRegions=None, fieldNames=None, overwrite=False):
        return self.getData(regularMesh=False, filterName=filterName, timeList=timeList, meshRegions=meshRegions,
                            fieldNames=fieldNames,
                            overwrite=overwrite)

    def _buildFilterQuery(self, filterName, meshRegions,regularMesh=None):
        qry = dict(simulation=self.simulationParams,
                   pipeline=self.vtkpipeline.toJSON())
        if meshRegions is not None:
            qry['simulation']['MeshRegions'] = meshRegions

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
        ret = []

        if structureJson is not None:
            for filterGuiName in structureJson:
                paramPairList = structureJson[filterGuiName]['params']  # must be a list to enforce order in setting.
                filtertype = structureJson[filterGuiName]['filterType']
                filter = getattr(pvsimple, filtertype)(Input=father, guiName=filterGuiName)

                logger.debug(
                    f"Adding filter {filterGuiName} of type {filtertype} to {'Reader' if fatherName is None else fatherName}")

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
                logger.debug(f"Filter {newFilterName} added to the pipeline. Now adding its downstream filters.")
                ret.append(filterGuiName)
                ret += self._buildFilterLayer(newFilterName, filter,
                                              structureJson[filterGuiName].get("downstream", None))

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
        self.write = write
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
        ret['downstream'] = {}
        for filterName, dsFilter in self.downstream.items():
            ret['downstream'][filterName] = dsFilter.toJSON()

        return {self.name: ret}

    def __setitem__(self, key, value):
        self.downstream[key] = value

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

    def addFilter(self, name, filterType, write=None, params=[]):
        newFilter = VTKPipeLine.newVTKPipelineFilter(name=name, filterType=filterType, write=write, params=params)
        self.__setitem__(name, newFilter)
        return newFilter


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
                       ("SliceType.Origin", origin),
                       ("HyperTreeGridSlicer.Origin", [0, 0, 0])
                       ]


class vtkFilter_CellCenters(VTKFilter):
    """
        The Cell center filter.
    """

    def __init__(self, name, write, **kwargs):
        kwargs.setdefault("params", [])
        super().__init__(name=name, filterType="CellCenters", write=write, **kwargs)


class vtkFilter_ExtractBlock(VTKFilter):
    def __init__(self, name, write, **kwargs):
        selectors = [f'/Root/boundary/{patchName}' for patchName in patchList]
        if internalMesh:
            selectors += ['/Root/internalMesh']
        params=  [
            ("Selectors",selectors)
        ]
        super().__init__(name=name, filterType="ExtractBlock", write=write, **kwargs)


class vtkFilter_IntegrateVariables(VTKFilter):
    def __init__(self, name, write, patchList=[],internalMesh=False):
        super().__init__(name=name, filterType="IntegrateVariables", write=write, params=params)
