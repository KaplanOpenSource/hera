import os
import logging
import numpy
from .abstractLocation import datalayer as locationDatalayer
from ....datalayer import datatypes
import matplotlib.pyplot as plt

import requests
import random
from osgeo import gdal
import numpy as np


from shapely.geometry import Point,box,MultiLineString, LineString

try:
    from freecad import app as FreeCAD
    import Part
    import Mesh
except ImportError as e:
    logging.warning("FreeCAD not Found, cannot convert to STL")


class datalayer(locationDatalayer):
    """
        Holds a polygon with description. Allows querying on the location of the shapefile.

        The projectName in the public DB is 'locationGIS'
    """
    _dxdy = None
    _publicProjectName = None

    @property
    def dxdy(self):
        return self._dxdy

    @dxdy.setter
    def dxdy(self, value):
        if value is not None:
            self._dxdy = value

    @property
    def skipinterior(self):
        return self._skipinterior

    def __init__(self, projectName, FilesDirectory="", databaseNameList=None, useAll=False,publicProjectName="Topography",Source="BNTL", dxdy=None):

        self.publicProjectName = publicProjectName
        super().__init__(projectName=projectName,publicProjectName=self.publicProjectName,FilesDirectory=FilesDirectory,databaseNameList=databaseNameList,useAll=useAll,Source=Source)

        self._dxdy = 10 if dxdy is None else dxdy  # m
        self._skipinterior = 100  # 100 # m, the meters to remove from the interpolation to make sure all exists.


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
                n = numpy.cross(numpy.array(v1) - numpy.array(v2), numpy.array(v1) - numpy.array(v3))
                n = n / numpy.sqrt(sum(n ** 2))
                stl_str += self._make_facet_str(n, v1, v2, v3)

                # dem facet 2
                n = numpy.cross(numpy.array(v2) - numpy.array(v3), numpy.array(v2) - numpy.array(v4))
                n = n / numpy.sqrt(sum(n ** 2))
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
                        n = numpy.cross(numpy.array(vlist[k]) - numpy.array(vlist[l]), numpy.array(vblist[l]) - numpy.array(vlist[l]))
                        n = n / numpy.sqrt(sum(n ** 2))
                        stl_str += self._make_facet_str(n, vlist[k], vblist[l], vlist[l])

                        # add j-base,i-base,i
                        n = numpy.cross(numpy.array(vlist[k]) - numpy.array(vblist[k]), numpy.array(vlist[k]) - numpy.array(vblist[l]))
                        n = n / numpy.sqrt(sum(n ** 2))
                        stl_str += self._make_facet_str(n, vlist[k], vblist[k], vblist[l])

        stl_str += 'endsolid ' + solidname + '\n'
        return stl_str

    def toSTL(self, doc, solidname, flat=None):
        """
            Gets a shape file of topography.
            each contour line has property 'height'.
            Converts it to equigrid xy mesh and then build the STL.
        """

        # 1. read the shp file.
        gpandas = doc.getData()

        # 2. Convert contour map to regular height map.
        # 2.1 get boundaries
        xmin = gpandas['geometry'].bounds['minx'].min()
        xmax = gpandas['geometry'].bounds['maxx'].max()

        ymin = gpandas['geometry'].bounds['miny'].min()
        ymax = gpandas['geometry'].bounds['maxy'].max()

        print("Mesh boundaries x=(%s,%s) ; y=(%s,%s)" % (xmin, xmax, ymin, ymax))
        # 2.2 build the mesh.
        inter = self.skipinterior
        grid_x, grid_y = numpy.mgrid[(xmin + inter):(xmax - inter):self.dxdy,
                         (ymin + inter):(ymax - inter):self.dxdy]

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
                Height[i] = flat

        # adding fill values for places outside the map, e.g. inside the sea.
        grid_z2 = numpy.griddata(XY, Height, (grid_x, grid_y), method='cubic', fill_value=0.)
        # replace zero height with small random values so the stl file won't contain NaNs
        for i in range(grid_z2.shape[0]):
            for j in range(grid_z2.shape[1]):
                if (grid_z2[i, j] == 0.):
                    grid_z2[i, j] = numpy.random.random()

        if numpy.isnan(grid_z2).any():
            print("Found some NaN in cubic iterpolation. consider increasing the boundaries of the interior")

        stlstr = self._makestl(grid_x, grid_y, grid_z2, solidname)

        data = {"grid_x": grid_x, "grid_y": grid_y, "grid_z": grid_z2, "XY": XY, "Height": Height, "geo": gpandas}
        print("X min: %s , X max: %s " % (numpy.min(grid_x), numpy.min(grid_x)))
        print("Y min: %s , Y max: %s " % (numpy.min(grid_y), numpy.min(grid_y)))
        return stlstr, data


def get_altitdue_ip(lat, lon):
    """
    returning the altitude of the point, uses free mapquest data that is limited in the amount of calls per month, it uses Nir BAMBA Benmoshe key

    param lat - the path where we save the stl
    param lon - the width of the domain in the x direction

    return:
    altitude - meters above sea level
    """

    resp = requests.get(
        'http://open.mapquestapi.com/elevation/v1/profile?key=D5z9RSebQJLbUs4bohANIB4TzJdbvyvm&shapeFormat=raw&latLngCollection=' + str(
            lat) + ',' + str(lon))

    height = resp.json()['elevationProfile'][0]['height']

    return height


def get_altitdue_gdal(lat, lon):
    #        if lat<=30:
    # USGS EROS Archive - Digital Elevation - Global Multi-resolution Terrain Elevation Data 2010 (GMTED2010)
    fheight = r'/data3/nirb/10N030E_20101117_gmted_med075.tif'
    fheight = r'/data3/nirb/30N030E_20101117_gmted_med075.tif'
    # https://dds.cr.usgs.gov/srtm/version2_1/SRTM3/Africa/   # 90m resolution
    if lat > 29 and lat < 30 and lon > 34 and lon < 35:
        fheight = r'/data3/nirb/N29E034.hgt'
    elif lat > 29 and lat < 30 and lon > 35 and lon < 36:
        fheight = r'/data3/nirb/N29E035.hgt'
    elif lat > 30 and lat < 31 and lon > 34 and lon < 35:
        fheight = r'/data3/nirb/N30E034.hgt'
    elif lat > 30 and lat < 31 and lon > 35 and lon < 36:
        fheight = r'/data3/nirb/N30E035.hgt'
    elif lat > 31 and lat < 32 and lon > 34 and lon < 35:
        fheight = r'/data3/nirb/N31E034.hgt'
    elif lat > 31 and lat < 32 and lon > 35 and lon < 36:
        fheight = r'/data3/nirb/N31E035.hgt'
    elif lat > 32 and lat < 33 and lon > 34 and lon < 35:
        fheight = r'/data3/nirb/N32E034.hgt'
    elif lat > 32 and lat < 33 and lon > 35 and lon < 36:
        fheight = r'/data3/nirb/N32E035.hgt'
    elif lat > 33 and lat < 33 and lon > 35 and lon < 36:
        fheight = r'/data3/nirb/N33E035.hgt'
    else:
        print('!!!!NOT in Israel !!!!!!!!')
        # taken from https://earthexplorer.usgs.gov/
        fheight = r'/ibdata2/nirb/gt30e020n40.tif'

    ds = gdal.Open(fheight)
    myarray = np.array(ds.GetRasterBand(1).ReadAsArray())
    myarray[myarray < -1000] = 0
    gt = ds.GetGeoTransform()
    rastery = (lon - gt[0]) / gt[1]
    rasterx = (lat - gt[3]) / gt[5]
    height11 = myarray[int(rasterx), int(rastery)]
    height12 = myarray[int(rasterx) + 1, int(rastery)]
    height21 = myarray[int(rasterx), int(rastery) + 1]
    height22 = myarray[int(rasterx) + 1, int(rastery) + 1]
    height1 = (1. - (rasterx - int(rasterx))) * height11 + (rasterx - int(rasterx)) * height12
    height2 = (1. - (rasterx - int(rasterx))) * height21 + (rasterx - int(rasterx)) * height22
    height = (1. - (rastery - int(rastery))) * height1 + (rastery - int(rastery)) * height2

    return height


if __name__ == "__main__":
    lon = random.randint(35750, 35800) / 1000.0  # Hermon
    lat = random.randint(33250, 33800) / 1000.0
    #    lon = 34.986008  // elevation should be 273m according to amud anan
    #    lat = 32.808486  // elevation should be 273m according to amud anan
    #    lon = 35.755  // elevation should be ~820 according to amud anan
    #    lat = 33.459  // elevation should be ~820 according to amud anan
    lon = 35.234987  # 744m
    lat = 31.777978  # 744m
    alt1 = get_altitdue_ip(lat, lon)
    alt2 = get_altitdue_gdal(lat, lon)
    print("the altitude at position: ", lat, lon, " is ", alt1, alt2)

