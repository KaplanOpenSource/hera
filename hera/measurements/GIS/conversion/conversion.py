import geopandas
from shapely.geometry import LineString
from scipy.interpolate import griddata
from numpy import array, cross, sqrt
import numpy
import pandas
from ..datalayer.datalayer import GIS_datalayer
from hera.measurements.GIS.utils import dataManipulations
from .... import datalayer

# import FreeCAD

class converter():

    _projectName = None
    _FilesDirectory = None
    _Measurments = None
    _GISdatalayer = None
    _manipulator = None

    def __init__(self, projectName, FilesDirectory):

        self._FilesDirectory = FilesDirectory
        self._projectName = projectName
        self._GISdatalayer = GIS_datalayer(projectName=projectName, FilesDirectory=FilesDirectory)
        self._Measurments = datalayer.Measurements
        self._manipulator = dataManipulations()

    def addSTLpath(self, path, NewFileName, **kwargs):
        """
        Adds a path to the dataframe under the type 'stlType'.

        Parameters:

            path: The path (string)
            NewFileName: A name for the file (string)
            kwargs: Additional parameters for the document.

        Returns: A list that contains the stl string and a dict that holds information about it.
        -------

        """

        self._Measurments.addDocument(projectName=self._projectName,
                                      desc=dict(name = NewFileName, **kwargs),
                                      type="stlFile",
                                      resource=path,
                                      dataFormat="string")

    def convert_to_stl(self, data, NewFileName, dxdy=None, save=True, flat=None, path=None, **kwargs):

        """
        Converts a geopandas dataframe data to an stl file.

        Parameters:

            data: The data that should be converted to stl. May be a dataframe or a name of a saved polygon in the database.
            NewFileName: A name for the new stl file, also used in the stl string. (string)
            dxdy: the dimention of each cell in the mesh in meters, the default is 50.
            save: Default is True. If True, the new stl string is saved as a file and the path to the file is added to the database.
            flat: Default is None. Else, it assumes that the area is flat and the value of flat is the height of the mesh cells.
            path: Default is None. Then, the path in which the data is saved is the given self.FilesDirectory. Else, the path is path. (string)
            kwargs: Any additional metadata to be added to the new document in the database.

        Returns
        -------

        """

        if type(data) == str:
            polygon = self._GISdatalayer.getGeometry(data)
            dataframe = self._GISdatalayer.getGISData(Geometry=data, GeometryMode="contains")[0]
            geodata = self._manipulator.PolygonDataFrameIntersection(polygon=polygon, dataframe=dataframe)
        elif type(data) == geopandas.geodataframe.GeoDataFrame:
            geodata = data
        else:
            raise KeyError("data should be geopandas dataframe or a polygon.")

        stl = STL(NewFileName,dxdy)
        stlstr, newdata = stl.Convert_geopandas_to_stl(gpandas=geodata, flat=flat)

        if save:
            p = self._FilesDirectory if path is None else path
            new_file_path = p + "/" + NewFileName + ".stl"
            new_file = open(new_file_path, "w")
            new_file.write(stlstr)
            self.addSTLpath(p, NewFileName, **kwargs)

        return stlstr, newdata

class STL():

    _dxdy = None
    _NewFileName = None

    @property
    def dxdy(self):
        return self._dxdy

    @dxdy.setter
    def dxdy(self, value):
        if value is not None:
            self._dxdy = value

    def __init__(self, NewFileName, dxdy=None):

        self._dxdy = 50 if dxdy is None else dxdy  # m
        self._NewFileName = NewFileName

    def _make_facet_str(self, n, v1, v2, v3):
        facet_str = 'facet normal ' + ' '.join(map(str, n)) + '\n'
        facet_str += '  outer loop\n'
        facet_str += '      vertex ' + ' '.join(map(str, v1)) + '\n'
        facet_str += '      vertex ' + ' '.join(map(str, v2)) + '\n'
        facet_str += '      vertex ' + ' '.join(map(str, v3)) + '\n'
        facet_str += '  endloop\n'
        facet_str += 'endfacet\n'
        return facet_str

    def _makestl(self, X, Y, elev):
        """
            Takes a mesh of x,y and elev and convert it to stl file.

            X - matrix of x coordinate. [[ like meshgrid ]]
            Y - matrix of y coordinate. [[ like meshgrid ]]
            elev - matrix of elevation.

        """
        base_elev = 0
        stl_str = 'solid ' + self._NewFileName + '\n'
        for i in range(elev.shape[0] - 1):
            for j in range(elev.shape[1] - 1):

                x = X[i, j];
                y = Y[i, j]
                v1 = [x, y, elev[i, j]]

                x = X[i + 1, j];
                y = Y[i, j]
                v2 = [x, y, elev[i + 1, j]]

                x = X[i, j];
                y = Y[i, j + 1]
                v3 = [x, y, elev[i, j + 1]]

                x = X[i + 1, j + 1];
                y = Y[i + 1, j + 1]
                v4 = [x, y, elev[i + 1, j + 1]]

                # dem facet 1
                n = cross(array(v1) - array(v2), array(v1) - array(v3))
                n = n / sqrt(sum(n ** 2))
                stl_str += self._make_facet_str(n, v1, v2, v3)

                # dem facet 2
                n = cross(array(v2) - array(v3), array(v2) - array(v4))
                n = n / sqrt(sum(n ** 2))
                # stl_str += self._make_facet_str( n, v2, v3, v4 )
                stl_str += self._make_facet_str(n, v2, v4, v3)

                # base facets
                v1b = list(v1)
                v2b = list(v2)
                v3b = list(v3)
                v4b = list(v4)

                v1b[-1] = base_elev
                v2b[-1] = base_elev
                v3b[-1] = base_elev
                v4b[-1] = base_elev

                n = [0.0, 0.0, -1.0]

                stl_str += self._make_facet_str(n, v1b, v2b, v3b)
                stl_str += self._make_facet_str(n, v2b, v3b, v4b)

                vlist = [v1, v2, v3, v4]
                vblist = [v1b, v2b, v3b, v4b]

                # Now the walls.
                for k, l in [(0, 1), (0, 2), (1, 3), (2, 3)]:
                    # check if v[i],v[j] are on boundaries.
                    kboundary = False
                    if vlist[k][0] == X.min() or vlist[k][0] == X.max():
                        kboundary = True

                    lboundary = False
                    if vlist[l][1] == Y.min() or vlist[l][1] == Y.max():
                        lboundary = True

                    if (kboundary or lboundary):
                        # Add i,j,j-base.
                        n = cross(array(vlist[k]) - array(vlist[l]), array(vblist[l]) - array(vlist[l]))
                        n = n / sqrt(sum(n ** 2))
                        stl_str += self._make_facet_str(n, vlist[k], vblist[l], vlist[l])

                        # add j-base,i-base,i
                        n = cross(array(vlist[k]) - array(vblist[k]), array(vlist[k]) - array(vblist[l]))
                        n = n / sqrt(sum(n ** 2))
                        stl_str += self._make_facet_str(n, vlist[k], vblist[k], vblist[l])

        stl_str += 'endsolid ' + self._NewFileName + '\n'
        return stl_str

    def Convert_geopandas_to_stl(self, gpandas, flat=None):
        """
            Gets a shape file of topography.
            each contour line has property 'height'.
            Converts it to equigrid xy mesh and then build the STL.
        """

        # 1. Convert contour map to regular height map.
        # 1.1 get boundaries
        xmin = gpandas['geometry'].bounds['minx'].min()
        xmax = gpandas['geometry'].bounds['maxx'].max()

        ymin = gpandas['geometry'].bounds['miny'].min()
        ymax = gpandas['geometry'].bounds['maxy'].max()

        print("Mesh boundaries x=(%s,%s) ; y=(%s,%s)" % (xmin, xmax, ymin, ymax))
        # 1.2 build the mesh.
        grid_x, grid_y = numpy.mgrid[(xmin):(xmax):self.dxdy, (ymin):(ymax):self.dxdy]

        # 2. Get the points from the geom
        Height = []
        XY = []

        for i, line in enumerate(gpandas.iterrows()):
            if isinstance(line[1]['geometry'], LineString):
                linecoords = [x for x in line[1]['geometry'].coords]
                lineheight = [line[1]['HEIGHT']] * len(linecoords)
                XY += linecoords
                Height += lineheight
            else:
                for ll in line[1]['geometry']:
                    linecoords = [x for x in ll.coords]
                    lineheight = [line[1]['HEIGHT']] * len(linecoords)
                    XY += linecoords
                    Height += lineheight
        if flat is not None:
            for i in range(len(Height)):
                Height[i] = flat

        grid_z2 = griddata(XY, Height, (grid_x, grid_y), method='cubic')

        if numpy.isnan(grid_z2).any():
            print("Found some NaN in cubic iterpolation. consider increasing the boundaries of the interior")

        stlstr = self._makestl(grid_x, grid_y, grid_z2)

        data = pandas.DataFrame({"XY": XY, "Height": Height})

        return stlstr, data