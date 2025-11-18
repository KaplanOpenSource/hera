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
        if params is None:
            params = []
        newFilter = VTKPipeLine.newVTKPipelineFilter(name=name, filterType=filterType, write=write, params=params)
        self.filters[name] = newFilter
        return newFilter

    def addFilterFromObj(self,newFilter):
        """
        Adds a filter to the pipeline using an already existent instance.

        Parameters

            filter(VTKFilter)
                an instance of a filter

        """
        self.filters[newFilter.name] = newFilter

    @deprecated("Use addFilterFromObj")
    def addExistingFilter(self, filter):
        """
        Adds a filter to the pipeline using an already existent instance.

        Parameters
        
            filter(VTKFilter)
                an instance of a filter

        """
        self.filters[newFilter.name] = newFilter

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
            ret = []
            nextList = []
            for filterName, filterObj in filtersList.items():
                currentName = filterName if fatherPath is None else f"{fatherPath}.{filterName}"
                if (writeOnly and filterObj.write) or not writeOnly:
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
                "workflowName": os.path.basename(self.casePath),
                "groupName": groupName,
                "workflowProperties": {}
            }

        else:
            simulationDocumentList = self.datalayer.getWorkflowDocumentFromDB(nameOrWorkflowFileOrJSONOrResource)
            if simulationDocumentList is not None:
                simulationDocument = simulationDocumentList[0]
                self.casePath = simulationDocument.resource
                self.simulationDocument = simulationDocument
                simulationProperties = self.simulationDocument['desc'].copy()
                self.simulationParams = {
                    "workflowName": simulationProperties['workflowName'],
                    "groupName": simulationProperties['groupName'],
                    "workflowParameters": simulationProperties['parameters']
                }

            else:
                raise ValueError(
                    f"Simulation {nameOrWorkflowFileOrJSONOrResource} is not in the DB, and does not represent a valid case directory")

        self.pvOFBase = paraviewOpenFOAM(casePath=self.casePath,
                                         caseType=caseType,
                                         servername=serverName)

    def clearCache(self, regularMesh=None,filterName=None):
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

    def getData(self, regularMesh, filterName=None, timeList=None , fieldNames=None,overwrite=False):
        """
            Returns the data of the vtkpipeline as a dict.
            The stuctucture is similar to that of a vtkpipline.
        Parameters
        ----------
        filterName : str, list of str, None
            The name of the filter to get, a list of filters or get all the filters that write=True in the pipeline.
        timeList : None, list
            The list of timestep
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
            requestedFiltersToProcess = self.vtkpipeline.allFilterNames(writeOnly=True)
        else:
            requestedFiltersToProcess = list(numpy.atleast_1d(filterName))
        logger.info(f"The requested filters are : {requestedFiltersToProcess}")

        CaseTimeList = self.datalayer.getTimeList(self.casePath)

        # Now execute the pipeline.
        if timeList is not None:
            if isinstance(timeList, str):
                # a bit of parsing.
                readerTL = CaseTimeList
                BandA = [readerTL[0], readerTL[-1]]

                for i, val in enumerate(timeList.split(":")):
                    BandA[i] = BandA[i] if len(val) == 0 else float(val)

                tl = pandas.Series(readerTL)
                timeList = tl[tl.between(*BandA)].values
            else:
                timeList = timeList
        else:
            timeList = CaseTimeList

        logger.debug(f"Getting timeList {timeList}")

        filtersToProcess = []
        filtersOutputFilename = dict()
        DBDocumentsDict = dict()
        for filterName in requestedFiltersToProcess:
            qry = self._buildFilterQuery(filterName=filterName,regularMesh=regularMesh)
            docList = self.datalayer.getSimulationsDocuments(type=TYPE_VTK_FILTER, **dictToMongoQuery(qry))
            # attempt to extract cached filter results
            if len(docList) > 0:
                # There should only be one filter result cached
                cached_filter = docList[0]
                filtersOutputFilename[filterName] = cached_filter.resource
                dbTimeList = cached_filter['desc']['simulation']['timeList']

                # Filtering the timesteps.
                timeList = [ts for ts in timeList if ts not in dbTimeList]
                DBDocumentsDict[filterName] = cached_filter

            else:
                counter = self.datalayer.getCounterAndAdd("OpenFOAMData")
                outputFileName = f"{filterName.replace('.','_')}_{counter}.{filext}"
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
            reader = self.pvOFBase.initializeReader(readerName="reader")

            # 2. Initialiaze the reader and the time list.
            if fieldNames is not None:
                reader.CellArrays = fieldNames

            filtersToCompute = self._buildFilterLayer(fatherName=None,
                                                      father=reader,
                                                      structureJson=self.vtkpipeline.toJSON()['filters'])
            logger.info(f"Added all filters to the layer. Computing filters {filtersToCompute}")

            # Build the ourput file name.
            for filterName in filtersToProcess:
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

                recordData = self._buildFilterQuery(filterName=filterName, regularMesh=regularMesh)
                recordData['simulation']['timeList'] = timeList
                dataFormat = self.datalayer.datatypes.ZARR_XARRAY if regularMesh else self.datalayer.datatypes.PARQUET
                doc = self.datalayer.addSimulationsDocument(dataFormat=dataFormat,
                                                            resource=os.path.abspath(filtersToComputeDict[filterName]),
                                                            type=TYPE_VTK_FILTER,
                                                            desc=recordData)
                DBDocumentsDict[filterName] = doc

            logger.debug(f"Reading filter {filterName} data")

        for filterName in requestedFiltersToProcess:
            ret[filterName] = DBDocumentsDict[filterName].getData()

        return ret

    def getRegularData(self, filterName=None, timeList=None,  fieldNames=None, overwrite=False):
        return self.getData(regularMesh=True, filterName=filterName, timeList=timeList,
                            fieldNames=fieldNames,
                            overwrite=overwrite)

    def getNonRegularData(self, filterName=None, timeList=None,  fieldNames=None, overwrite=False):
        return self.getData(regularMesh=False, filterName=filterName, timeList=timeList,
                            fieldNames=fieldNames,
                            overwrite=overwrite)

    def _buildFilterQuery(self, filterName, regularMesh=None):
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
        ret = []


        if structureJson is not None:
            for filterGuiName in structureJson:
                paramPairList = structureJson[filterGuiName]['params']  # must be a list to enforce order in setting.
                filtertype = structureJson[filterGuiName]['filterType']
                newFilterName = filterGuiName if fatherName is None else f"{fatherName}.{filterGuiName}"
                filter = getattr(pvsimple, filtertype)(Input=father, guiName=newFilterName)
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
                logger.debug(f"Filter {newFilterName} added to the pipeline. Now adding its downstream filters.")
                ret.append(newFilterName)
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

    def __init__(self, name, write, **kwargs):
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
    SAMPLE_UNIFORMLY = 'Sample Uniformly'
    SAMPLE_AT_CELL_BOUNDARIES = 'Sample At Cell Boundaries'
    SAMPLE_AT_SEGMENT_CENTERS = 'Sample At Segment Centers'

    def __init__(self, name, write, **kwargs):
        kwargs.setdefault("params",[])
        super().__init__(name=name, filterType="PlotOverLine", write=write, **kwargs)

    def setSamplePattern(self, pattern: Literal['Sample Uniformly', 'Sample At Cell Boundaries', 'Sample At Segment Centers']):
        self.set_param("SamplingPattern",  pattern)

    def setUniformSampleResolution(self, res: int):
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
        kwargs.setdefault("params", [])
        super().__init__(name=name, filterType="CellCenters", write=write, **kwargs)


class vtkFilter_ExtractBlock(VTKFilter):
    def __init__(self, name, write, patchList=[],internalMesh=False,params=None):
        if params is None:
            params = []
        selectors = [f'/Root/boundary/{patchName}' for patchName in patchList]
        if internalMesh:
            selectors += ['/Root/internalMesh']
        final_params =  [("Selectors",selectors)] + params
        super().__init__(name=name, filterType="ExtractBlock", write=write,params=final_params)

    def setRegionsToExtract(self,patchList, internalMesh=True):
        selectors = [f'/Root/boundary/{patchName}' for patchName in patchList]
        if internalMesh:
            selectors += ['/Root/internalMesh']

        self.set_param("Selectors", selectors)



class vtkFilter_IntegrateVariables(VTKFilter):
    def __init__(self, name, write, **kwargs):
        super().__init__(name=name, filterType="IntegrateVariables", write=write, **kwargs)
