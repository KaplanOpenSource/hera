import numpy
import geopandas
import shapely
import pandas
import geopandas
import io
import scipy
import numpy
import math
import os
from scipy.interpolate import griddata
import numpy as np
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
import xarray as xr


# ESPG codes
WGS84 = 4326
WSG84 = 4326 # A mistake, keep it for previuos versions that use the code.
ITM   = 2039 # Israeli
ED50_ZONE36N = 23036
BETA = 0.2
KARMAN = 0.41

def convertCRS(points,inputCRS,outputCRS,**kwargs):
    """
        Converts the points in data, assuming to be in input CRS to points in output CRS.

    Parameters
    ----------
    data : numpy.array (2D), pandas.DataFrame, list of 2-tuples
            An array of (x,y) points.

    inputCRS : integer,
            An EPSG code of the original points

    outputCRS : integer, an EPSG
            An EPSG code of the output points.

    kwargs :
            Additional information.
            if data is DataFrame:
                - 'x' : The name of the x column, default "x"
                - 'y' : The name of the y column, default "y"

    Returns
    -------
            pandas with columns x,y in the correct CRS.

    """
    if isinstance(points,numpy.ndarray):
        if len(points.shape) == 1:
            origpoints = geopandas.points_from_xy([points[0]],[points[1]])
        else:
            origpoints = geopandas.points_from_xy(points[:,0],points[:,1])

    elif isinstance(points,pandas.DataFrame):
        origpoints = geopandas.points_from_xy(points[kwargs.get("x","x")],
                                               points[kwargs.get("y","y")])
    elif isinstance(points,list):
        if len(points) == 1:
            origpoints = geopandas.points_from_xy([points[0][0]],[points[0][1]])
        else:
            origpoints = geopandas.points_from_xy([x[0] for x in points],[x[1] for x in points])

    else:
        raise ValueError(f"points must be numpy.array,pandas.DataFrame, or list. Not {type(points)}")


    origpoints.crs = inputCRS
    return origpoints.to_crs(outputCRS)


def create_xarray(minx,miny,maxx,maxy,dxdy=30, inputCRS=WSG84):
    if inputCRS == WSG84:
        min_pp = convertCRS(points=[[miny, minx]], inputCRS=inputCRS, outputCRS=ITM)[0]
        max_pp = convertCRS(points=[[maxy, maxx]], inputCRS=inputCRS, outputCRS=ITM)[0]
    else:
        min_pp = Point(minx, miny)
        max_pp = Point(maxx, maxy)

    x = numpy.arange(min_pp.x, max_pp.x, dxdy)
    y = numpy.arange(min_pp.y, max_pp.y, dxdy)
    xx, yy = np.meshgrid(x, y[::-1])
    grid_points = pd.DataFrame({
        'x': xx.ravel(),
        'y': yy.ravel()
    })

    gdf = gpd.GeoDataFrame(
        grid_points,
        geometry=gpd.points_from_xy(grid_points['x'], grid_points['y']),
        crs=ITM
    )

    gdf_transformed = gdf.to_crs(WSG84)

    gdf_transformed['lat'] = gdf_transformed.geometry.y
    gdf_transformed['lon'] = gdf_transformed.geometry.x

    lat_grid = gdf_transformed['lat'].values.reshape(xx.shape)
    lon_grid = gdf_transformed['lon'].values.reshape(xx.shape)

    i = np.arange(xx.shape[0])
    j = np.arange(xx.shape[1])

    return xr.DataArray(
        coords={
            'i': i,
            'j': j,
            'lat': (['i', 'j'], lat_grid),
            'lon': (['i', 'j'], lon_grid),
            'dxdy': dxdy
        },
        dims=['i', 'j']
    )

class stlFactory:
    """
        Helper class to convert topography geopandas to STL
    """

    _heightColumnsNames=None

    def __init__(self):
        self._heightColumnsNames= "HEIGHT"

    @property
    def heightColumnsNames(self):
        return self._heightColumnsNames

    @heightColumnsNames.setter
    def heightColumnsNames(self, value):
        self._heightColumnsNames = value

    def _make_facet_str(self, n, v1, v2, v3):
        facet_str = 'facet normal ' + ' '.join(map(str, n)) + '\n'
        facet_str += '  outer loop\n'
        facet_str += '      vertex ' + ' '.join(map(str, v1)) + '\n'
        facet_str += '      vertex ' + ' '.join(map(str, v2)) + '\n'
        facet_str += '      vertex ' + ' '.join(map(str, v3)) + '\n'
        facet_str += '  endloop\n'
        facet_str += 'endfacet\n'
        return facet_str

    def rasterizeGeopandas(self, gpandas, dxdy=50):
        """
            Convert the height-contour in geopandas to raster (regular mesh)


        Parameters
        -----------
        gpandas: geopandas.geodataframe
                The height contours as geopandas.

                The structure of the geodataframe is geometry
                which are height lines and a height column (its name given in self.heightColumnsNames, default is HEIGHT).

        dxdy: float
            The resolution of the regularization


        Returns
        -------
            A dict after interpolating to regular mesh.

            Return dict with the keys:
              x:    A 2D map of the x coordinate
              y :   A 2D map of the y coordinate
              height: A 2D map of the height


        """
        # 1. Convert contour map to regular height map.
        # 1.1 get boundaries
        xmin = gpandas.bounds.min()["minx"]
        ymin = gpandas.bounds.min()["miny"]
        xmax = gpandas.bounds.max()["maxx"]
        ymax = gpandas.bounds.max()["maxy"]
        #print("Mesh boundaries x=(%s,%s) ; y=(%s,%s)" % (xmin, xmax, ymin, ymax))
        # 1.2 build the mesh.

        grid_x, grid_y = numpy.mgrid[(xmin):(xmax):dxdy, (ymin):(ymax):dxdy]

        # 2. Get the points from the geom
        Height = []
        XY = []
        for i, line in enumerate(gpandas.iterrows()):
            if isinstance(line[1]['geometry'], LineString):
                linecoords = [x for x in line[1]['geometry'].coords]
                lineheight = [line[1][self.heightColumnsNames]] * len(linecoords)
                XY += linecoords
                Height += lineheight
            else:
                for ll in line[1]['geometry']:
                    linecoords = [x for x in ll.coords]
                    lineheight = [line[1][self.heightColumnsNames]] * len(linecoords)
                    XY += linecoords
                    Height += lineheight

        grid_z2 = griddata(XY, Height, (grid_x, grid_y), method='linear')
        grid_z2 = self._organizeGrid(grid_z2)

        return dict(x=grid_x,y=grid_y,height=grid_z2)

    def rasterToSTL(self, grid_x, grid_y, grid_z, solidName):
        """

        Parameters
        ----------
        grid_x: numpy.array

            2D array with the x grid
        grid_y: numpy.array

            2D array with the y grid
        grid_z: numpy.array

            2D array with the z (height) grid

        solidName: str
            The name of the topography solid

        Returns
        -------
            The STL as string.
        """
        base_elev = grid_z.min() - 10
        stl_str = f"solid {solidName}\n"
        for i in range(grid_z.shape[0] - 1):
            for j in range(grid_z.shape[1] - 1):

                x = grid_x[i, j]
                y = grid_y[i, j]
                v1 = [x, y, grid_z[i, j]]

                x = grid_x[i + 1, j];
                y = grid_y[i, j]
                v2 = [x, y, grid_z[i + 1, j]]

                x = grid_x[i, j];
                y = grid_y[i, j + 1]
                v3 = [x, y, grid_z[i, j + 1]]

                x = grid_x[i + 1, j + 1];
                y = grid_y[i + 1, j + 1]
                v4 = [x, y, grid_z[i + 1, j + 1]]

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
                    if vlist[k][0] == grid_x.min() or vlist[k][0] == grid_x.max():
                        kboundary = True

                    lboundary = False
                    if vlist[l][1] == grid_y.min() or vlist[l][1] == grid_y.max():
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

        stl_str += f"endsolid {solidName}\n"
        return stl_str

    def vectorToSTL(self,gpandas,dxdy=50,solidName="Topography"):
        """
            Convert the vector to topography

        Parameters
        ----------
        gpandas: geopandas.geodataframe
            The height contours as geopandas.

        dxdy: float
            The resolution of the conversion

        solidName: str
            The name of the solid. The default is 'Topogrphay'.

        Returns
        -------
            string of the STL
        """
        if isinstance(gpandas, geopandas.GeoDataFrame):
            rasterMap = self.rasterizeGeopandas(gpandas, dxdy=dxdy)

        elif isinstance(gpandas, pandas.DataFrame) or isinstance(gpandas, dask.dataframe.DataFrame):
            rasterMap = self.rasterizePandas(gpandas,dxdy=dxdy)


        return self.rasterToSTL(grid_x=rasterMap['x'],
                                    grid_y=rasterMap['y'],
                                    grid_z=rasterMap['height'],
                                    solidName=solidName)

    def rasterizePandas(self, gpandas, dxdy=50.,xColumn="x",yColumn="y", heightColumn="height"):
        """
            Gets a shape file of topography.
            each contour line has property 'height'.
            Converts it to equigrid xy mesh and then build the STL.
        """

        # 1. Convert contour map to regular height map.
        # 1.1 get boundaries
        xmin = gpandas[xColumn].min()
        xmax = gpandas[xColumn].max()

        ymin = gpandas[yColumn].min()
        ymax = gpandas[yColumn].max()

        # 1.2 build the mesh.
        grid_x, grid_y = numpy.mgrid[(xmin):(xmax):dxdy, (ymin):(ymax):dxdy]
        # 2. Get the points from the geom
        Nx = int(((xmax - xmin) / dxdy))
        Ny = int(((ymax - ymin) / dxdy))
        grid_z2 = coordinateHandler.regularizeTimeSteps(data=gpandas, fieldList=[heightColumn],
                                              coord1=xColumn,
                                              coord2=yColumn,
                                              n=(Nx, Ny), addSurface=False, toPandas=False)[0][heightColumn]
        grid_z2 = self._organizeGrid(grid_z2)
        return {"x":grid_x,"y":grid_y,"height":grid_z2.values}

    def _organizeGrid(self, grid):

        for row in grid:
            for i in range(len(row)):
                if math.isnan(row[i]):
                    pass
                else:
                    break
            for n in range(i):
                row[n] = row[i]
            for i in reversed(range(len(row))):
                if math.isnan(row[i]):
                    pass
                else:
                    break
            for n in range(len(row)-i):
                row[-n-1] = row[i]
        return grid

    # def ITMtolatlong(point1):
    #     """
    #         change the coordinate system, from x y (ITM) to lat long
    #
    #     Parameters
    #     ----------
    #     x
    #     y
    #
    #     Returns
    #     -------
    #     lat, long
    #
    #     How to use
    #     ----------
    #     from hera.measurements.GIS.raster.topography import TopographyToolkit
    #     TopographyToolkit.ITMtolatlong([200000, 740000])
    #
    #     What to fix
    #     -----------
    #     """
    #     ITM = "epsg:2039"
    #     WGS84 = "epsg:4327"
    #     x = point1[0]
    #     y = point1[1]
    #     coordTransformer2 = Transformer.from_crs(ITM, WGS84, always_xy=True)
    #     lat, long = coordTransformer2.transform(x, y)
    #     return lat, long


# def PolygonDataFrameIntersection(dataframe, polygon):
#     """
#     Creates a new dataframe based on the intersection of a dataframe and a polygon.
#
#     Parameters:
#     ----------
#     dataframe: A geopandas dataframe.
#     polygon: A shapely polygon
#
#     Returns:
#     --------
#
#     geopandas.dataframe.
#
#     A new geopandas dataframe
#
#     """
#
#     newlines = []
#     for line in dataframe["geometry"]:
#         newline = polygon.intersection(line)
#         newlines.append(newline)
#     dataframe["geometry"] = newlines
#     dataframe = dataframe[~dataframe["geometry"].is_empty]
#
#     return dataframe
#
# def getBoundaries(doc):
#     """
#     Returns a dict with  the boundaries of the document.
#
#     Parameters:
#     -----------
#
#     doc:
#
#     Returns:
#     ---------
#
#     dict with the boundaries.
#     """
#
#     dataframe = doc.getDocFromDB()
#     points = doc.desc["points"]
#
#     boundaries = dict(xmin=points[0], xmax=points[2], ymin=points[1], ymax=points[3], zmin=dataframe["HEIGHT"].min(), zmax=dataframe["HEIGHT"].max())
#
#     return boundaries
#
# def makePolygonFromEndPoints(points):
#
#     polygon = shapely.geometry.Polygon([[points[0],points[1]],
#                                         [points[0],points[3]],
#                                         [points[2],points[3]],
#                                         [points[2],points[1]]])
#     return polygon
#
# def ConvexPolygons(data, buffer=100):
#     """
#     Returns polygons of groups of buildings.
#     """
#     data = data.reset_index()
#     d = data.buffer(buffer)
#     indicelist=[[0]]
#     for i in range(1,len(data)):
#         found = False
#         for g in range(len(indicelist)):
#             for n in indicelist[g]:
#                 if d[i].intersection(d[n]).is_empty:
#                     continue
#                 else:
#                     indicelist[g].append(i)
#                     found = True
#                     break
#             if found:
#                 break
#             if g==len(indicelist)-1:
#                 indicelist.append([i])
#
#     geo = data.loc[indicelist[0]].unary_union.convex_hull
#     gpd = geopandas.GeoDataFrame.from_dict([{"geometry":geo,"area":geo.area}])
#     for indice in indicelist[1:]:
#         geo = data.loc[indice].unary_union.convex_hull
#         gpd = pandas.concat([gpd,geopandas.GeoDataFrame.from_dict([{"geometry":geo,"area":geo.area}])])
#
#     gpd = gpd.sort_values(by="area", ascending=False).reset_index()
#     found=False
#     for i in range(len(gpd)):
#         for j in range(i+1,len(gpd)):
#             if gpd.loc[i].geometry.intersection(gpd.loc[j].geometry).is_empty:
#                 continue
#             else:
#                 found = True
#                 break
#         if found:
#             break
#     if found:
#         gpd = ConvexPolygons(gpd,buffer=1)
#
#     return gpd
#
#
