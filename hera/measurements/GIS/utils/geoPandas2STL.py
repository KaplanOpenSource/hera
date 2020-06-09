import geopandas
from shapely.geometry import MultiLineString, LineString
from scipy.interpolate import griddata
from numpy import array, cross, sqrt
import numpy
import matplotlib.pyplot as plt
from itertools import product
import argparse

#import FreeCAD

class Topography(object):
    """
        Shape file utils for openfoam.

    """

    _dxdy = None

    @property
    def dxdy(self):
        return self._dxdy

    @dxdy.setter
    def dxdy(self, value):
        if value is not None:
            self._dxdy = value

    def __init__(self, dxdy=None):
        self._dxdy = 10 if dxdy is None else dxdy  # m
        self._skipinterior = 100 #100 # m, the meters to remove from the interpolation to make sure all exists.

    def _make_facet_str(self, n, v1, v2, v3):
        facet_str = 'facet normal ' + ' '.join(map(str, n)) + '\n'
        facet_str += '  outer loop\n'
        facet_str += '      vertex ' + ' '.join(map(str, v1)) + '\n'
        facet_str += '      vertex ' + ' '.join(map(str, v2)) + '\n'
        facet_str += '      vertex ' + ' '.join(map(str, v3)) + '\n'
        facet_str += '  endloop\n'
        facet_str += 'endfacet\n'
        return facet_str

    def _makestl(self, X, Y, elev, solidname):
        """
            Takes a mesh of x,y and elev and convert it to stl file.

            X - matrix of x coordinate. [[ like meshgrid ]]
            Y - matrix of y coordinate. [[ like meshgrid ]]
            elev - matrix of elevation.

        """
        base_elev = 0
        stl_str = 'solid ' + solidname + '\n'
        for i in range(elev.shape[0] - 1):
            print(i)
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

        stl_str += 'endsolid ' + solidname + '\n'
        return stl_str

    def Convert_shp_to_stl(self, shpfile, solidname, flat=None):
        """
            Gets a shape file of topography.
            each contour line has property 'height'.
            Converts it to equigrid xy mesh and then build the STL.
        """

        # 1. read the shp file.
        gpandas = geopandas.read_file(shpfile)

        # 2. Convert contour map to regular height map.
        # 2.1 get boundaries
        xmin = gpandas['geometry'].bounds['minx'].min()
        xmax = gpandas['geometry'].bounds['maxx'].max()

        ymin = gpandas['geometry'].bounds['miny'].min()
        ymax = gpandas['geometry'].bounds['maxy'].max()

        print("Mesh boundaries x=(%s,%s) ; y=(%s,%s)" % (xmin, xmax, ymin, ymax))
        # 2.2 build the mesh.
        inter = self._skipinterior
        grid_x, grid_y = numpy.mgrid[(xmin + inter):(xmax - inter):self.dxdy, (ymin + inter):(ymax - inter):self.dxdy]

        # 3. Get the points from the geom
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
                Height[i]=flat
            
        grid_z2 = griddata(XY, Height, (grid_x, grid_y), method='cubic')

        if numpy.isnan(grid_z2).any():
            print("Found some NaN in cubic iterpolation. consider increasing the boundaries of the interior")

        stlstr = self._makestl(grid_x, grid_y, grid_z2, solidname)

        data = {"grid_x": grid_x, "grid_y": grid_y, "grid_z": grid_z2, "XY": XY, "Height": Height, "geo": gpandas}
        print("X min: %s , X max: %s " % (numpy.min(grid_x), numpy.min(grid_x)))
        print("Y min: %s , Y max: %s " % (numpy.min(grid_y), numpy.min(grid_y)))
        return stlstr, data


class Buildings(object):
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Write stlfile')
    parser.add_argument('--shpfile', help='The shp file name')
    parser.add_argument('--output', help='The stl file name')
    parser.add_argument('--dxdy', help='The dxdy resolution (meters)', default=50)

    args = parser.parse_args()

    fname = args.shpfile
    outfilename = args.output
    dxdy = args.dxdy
    # "/home/yehudaa/Dropbox/Projects/data/Galim/GALIM-CONTOUR.shp"
    # fname = "/home/yehudaa/Dropbox/EyalYehuda/Jerusalem/BLDG1.shp"
    fname = "/home/nirb/test/testme1/testme1-CONTOUR.shp"
    outfilename = "/home/nirb/test/testme1/testme1-CONTOUR.STL"
    dxdy = 50

    utls = Topography(dxdy=float(dxdy))
    print("Converting...")
    stlstr, data = utls.Convert_shp_to_stl(shpfile=fname, solidname="Topography")
    f = open(outfilename, 'w')
    f.write(stlstr)
    f.close()

