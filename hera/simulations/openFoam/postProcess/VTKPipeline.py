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
