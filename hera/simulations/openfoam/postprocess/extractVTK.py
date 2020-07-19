import pandas
import os
import json
from pvOpenFOAMBase import paraviewOpenFOAM
import paraview.simple as pvsimple

class VTKpipeline(object):
    """This class executes a pipeline (runs and saves the outputs).
    It also holds the metadata.

    Currently works only for the JSON pipeline. The XML (paraview native pipelines) will be built in the future.

    The pipeline is initialized with a reader.
    The metadata holds the names of all the filters that has to be executed.

    The VTK pipeline JSON structure.
        {
           "metadata" : {
                  "guiname" : <gui name>,
                   ... all other meta data ..........

            },
            "pipeline" : {
                            "type" : The type of the filter. (clip,slice,...).
                            "write"   : None/hdf (pandas)/netcdf (xarray),
                            "params" : [
                                    ("key","value"),
                                          .
                                          .
                                          .
                            ],...
                            "downstream" : [Another pipeline]
                        }
    }
     The write to writes the requested filters to the disk.
     The files are saved to a single .nc/hdf file with keys/fields in it. The file name is the
     pipelinename.

     As before, the hdf/nc is for one timestep.
    """

    _VTKpipelineJSON = None  # Holds the json of the VTK pipeline.
    _pvOFBase = None  # Holds the OF base.
    _mainpath = None
    _name = None

    @property
    def name(self):
        return self._name

    @property
    def pvOFBase(self):
        return self._pvOFBase

    def __init__(self, name, pipelineJSON, casePath, caseType='Decomposed Case', servername=None):
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
        self._pvOFBase = paraviewOpenFOAM(casePath=casePath, caseType=caseType, servername=servername)
        self._VTKpipelineJSON = pipelineJSON
        self._name = name

        outputdir = pipelineJSON["metadata"].get("datadir", "None")
        if outputdir != "None":
            self._pvOFBase.hdfdir = "%s/%s/hdf" % (outputdir, name)
            self._pvOFBase.netcdfdir = os.path.join(outputdir, name, "netcdf")
            self._mainpath = "%s/%s" % (outputdir, name)
        else:
            self._pvOFBase.hdfdir = os.path.join(name, "hdf")
            self._mainpath = os.path.join(name)

    def execute(self, source, writeMetadata=True):
        """
            Builds the pipeline from the JSON vtk.

        :param source:
            The source filter guiName that the pipeline will be build on.
        :param writeMetadata:
            if True, copies the json file to the results directory.

        """

        # build the pipeline.
        reader = pvsimple.FindSource(source)

        filterWrite = {}
        self._buildFilterLayer(father=reader, structureJson=self._VTKpipelineJSON["pipelines"], filterWrite=filterWrite)

        # Now execute the pipeline.
        timelist = self._VTKpipelineJSON["metadata"].get("timelist", "None")
        if (timelist == "None"):
            timelist = None

        elif isinstance(timelist, str) or isinstance(timelist, unicode):
            # a bit of parsing.
            BandA = [0, 1e6]

            for i, val in enumerate(timelist.split(":")):
                BandA[i] = BandA[i] if len(val) == 0 else float(val)

            tl = pandas.Series(reader.TimestepValues)
            timelist = tl[tl.between(*BandA)].values

        # else just take it from the json (it should be a list).

        # Get the mesh regions.
        if "MeshRegions" in self._VTKpipelineJSON["metadata"]:
            reader.MeshRegions = self._VTKpipelineJSON["metadata"]["MeshRegions"]

        for frmt, datasourceslist in filterWrite.items():
            writer = getattr(self._pvOFBase, "write_%s" % frmt)
            if writer is None:
                raise ValueError("The write %s is not found" % writer)
            writer(readername=source, datasourcenamelist=datasourceslist, timelist=timelist,
                   fieldnames=self._VTKpipelineJSON["metadata"].get('fields', None), outfile=self.name)

        if writeMetadata:
            with open('%s/meta.json' % self._mainpath, 'w') as outfile:
                json.dump(self._VTKpipelineJSON, outfile)

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

    def get_slice_height(filename, itteration, axis_index, axis_value, offset=0.0,
                         clipxmin=0.0, clipxmax=0.0, clipymin=0.0, clipymax=0.0, clipzmin=0.0, clipzmax=0.0,
                         reader=None):
        if reader is None:
            reader = bse.ReadCase('case name', filename, CaseType='Reconstructed Case')  # 'Reconstructed Case')
            print('Initial timelist print:', reader.TimestepValues)

        contour1 = pvsimple.Contour(reader)
        contour1.ContourBy = ['POINTS', 'HeightFromTopo']
        try:
            realoffset = offset + axis_value
            print('realoffset easy')
        except:
            print('realoffset except', len(offset))
            realoffset = []
            for i in range(len(offset)):
                realoffset.append(offset[i] + axis_value)

        contour1.Isosurfaces = realoffset
        contour1.PointMergeMethod = 'Uniform Binning'
        # // "downstream": {
        #                  // "pipeline": {
        #                                 // "type": "Contour",
        # // "guiname": "Contour1",
        # // "write": "hdf",
        # // "params": [["ContourBy", ["POINTS", "HeightFromTopo"]], ["Isosurfaces", 30],
        #               ["PointMergeMethod", "Uniform Binning"]]
        #              //}

# if __name__ == "__main__":
    # bse = pvOFBase()
    # reader = bse.ReadCase("Test", "AC4_3Da.foam", CaseType='Decomposed Case')  # 'Reconstructed Case')
    #
    # R = pipelineFactory_JSON().getPipeline("VTKPipe.json")
    # with open('test.json') as json_file:
    #     data = json.load(json_file)
    # vtkpipe = VTKpipeline(name="test", pipelineJSON=data, casePath="/home/ofir/Projects/openFoamUsage/askervein", caseType="Reconstructed Case")
    # vtkpipe.execute("mainReader")
