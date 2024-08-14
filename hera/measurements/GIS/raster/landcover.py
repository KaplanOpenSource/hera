from .... import toolkit
from ..utils import stlFactory,convertCRS,ITM,WSG84,ED50_ZONE36N
from ....utils.logging import get_classMethod_logger
import numpy
import math
from osgeo import gdal
from pyproj import Transformer
from itertools import product
import pandas
import os
import xarray
import geopandas
import xarray as xr
import numpy as np

class LandCoverToolkit(toolkit.abstractToolkit):
    # """
    # We should use different roughness to buildings and to topo
    #
    #  MCD12Q1 - MODIS/Terra+Aqua Land Cover Type Yearly L3 Global 500m SIN Grid
    #
    #   data is from https://zenodo.org/records/8367523
    #   type are in https://ladsweb.modaps.eosdis.nasa.gov/filespec/MODIS/6/MCD12Q1
    #
    #   Land Cover Type 1: IGBP global vegetation classification scheme
    #   Land Cover Type 2: University of Maryland (UMD) scheme
    #
    #   DataField    Name                    Data      Dimensions
    #                                        Type
    #
    #   DataField_1  Land_Cover_Type_1       UINT8     Dimension_1
    #                                                  Dimension_2
    #
    #           Description:    land cover types (IGBP)
    #
    #           DataField_1 HDF Attributes:
    #
    #              long_name      STRING  1   PGE    "Land_Cover_Type_1"
    #              units          STRING  1   PGE     "class number"
    #              valid_range    uint8    2   PGE     0 16
    #              _FillValue     uint8    1   PGE     255
    #              Water          uint8    1   PGE     0
    #              Evergreen needleleaf forest
    #                             uint8    1   PGE     1
    #              Evergreen broadleaf forest
    #                             uint8    1   PGE     2
    #              Deciduous needleleaf forest
    #                             uint8    1   PGE     3
    #              Deciduous broadleaf forest
    #                             uint8    1   PGE     4
    #              Mixed forests  unit8    1   PGE     5
    #              Closed shrubland
    #                             unit8    1   PGE     6
    #              Open shrublands
    #                             unit8    1   PGE     7
    #              Woody savannas unit8    1   PGE     8
    #              Savannas       unit8    1   PGE     9
    #              Grasslands     unit8    1   PGE    10
    #              Permanent wetlands
    #                             unit8    1   PGE    11
    #              Croplands      unit8    1   PGE    12
    #              Urban and built-up
    #                             unit8    1   PGE    13
    #              Cropland/natural vegetation mosaic
    #                             unit8    1   PGE    14
    #              Snow and ice   unit8    1   PGE    15
    #              Barren or sparsely vegetated
    #                             unit8    1   PGE    16
    #
    #   DataField_2  Land_Cover_Type_2       UINT8     Dimension_1
    #                                                  Dimension_2
    #
    #           Description:    land cover types (UMD)
    #
    #           DataField_2 HDF Attributes:
    #              long_name      STRING   1   PGE    "Land_Cover_Type_2"
    #              units          STRING   1   PGE     "class number"
    #              valid_range    uint8    2   PGE     0 16
    #              _FillValue     uint8    1   PGE     255
    #              Water          uint8    1   PGE     0
    #              Evergreen needleleaf forest
    #                             uint8    1   PGE     1
    #              Evergreen broadleaf forest
    #                             uint8    1   PGE     2
    #              Deciduous needleleaf forest
    #                             uint8    1   PGE     3
    #              Decidous broadleaf forest
    #                             uint8    1   PGE     4
    #              Mixed forests  uint8    1   PGE     5
    #              Closed shrublands
    #                              unit8   1   PGE     6
    #              Open shrubland  unit8   1   PGE     7
    #              Woody savannas  unit8   1   PGE     8
    #              Savannas       unit8    1   PGE     9
    #              Grasslands     uint8    1   PGE    10
    #              Croplands      unit8    1   PGE    12
    #              Urban and built-up
    #                             unit8    1   PGE    13
    #              Barren or sparsely vegetated
    #                             unit8    1   PGE    16
    #
    #   Land cover may change between winter and summer
    #   urban area may change over the years
    # """

    def __init__(self, projectName, filesDirectory=None):
        """
            Initializes land cover data toolkit.

        Parameters
        ----------
        projectName: str
            The project Name that the toolkit is initialized on
        toolkitName: str
            the specific toolkit, getting from the child.

        FilesDirectory: str or None
                The path to save a regions files when they are created.

                if str then represents a path (relative or absolute) to save the files in. The directory is created automatically.

                if None, then tries to get the default path of the project from the config. if it does not
                exist, then use the current directory.

        """
        super().__init__(projectName=projectName, toolkitName = 'LandCoverToolkit', filesDirectory=filesDirectory)



    def getLandCoverAtPoint(self,lon,lat,inputCRS=WSG84, dataSourceName=None):
        """
        Get the landcover type integer value in a specific point.

        Parameters
        ----------
        lon : float
            The longitude coodinates.

        lat : float
            The latitute coordinates.

        inputCRS : int , default=WSG84
            The ESPRG of the coordinates.

        dataSourceName : string , default=None
            The name of the data source to use.

        Returns
        -------
            int

        """
        dataSourceName = self.getConfig()['defaultLandCover'] if dataSourceName is None else dataSourceName
        datasourceDocument = self.getDataSourceDocument(dataSourceName)
        ds = self.getDataSourceData(dataSourceName)
        img = ds.GetRasterBand(1).ReadAsArray()
        gt = ds.GetGeoTransform()

        width = ds.RasterXSize
        height = ds.RasterYSize
        minx = gt[0]
        maxx = gt[0] + width * gt[1] + height * gt[2]
        miny = gt[3] + width * gt[4] + height * gt[5]
        maxy = gt[3]
        x = math.floor((lat - gt[3]) / gt[5])
        y = math.floor((lon - gt[0]) / gt[1])
        return img[x, y]

    def getLandCover(self,minlon,minlat,maxlon,maxlat,dxdy = 30, inputCRS=WSG84, dataSourceName=None):
        """
        Get Xarray LandCover map.

        Parameters
        ----------
        minlon : float
            Minimum value of longitude coodinates.

        minlat : float
            Minimum value of latitute coodinates.

        maxlon: float
            Maximum value of longitude coodinates.

        maxlat: float
            Maximum value of latitute coodinates.

        dxdy: int, default=30
            Spatial resolution of the output land cover map. Defines the step size for the grid points
            in both the x (longitude) and y (latitude) directions within the specified bounding box.
            Smaller values result in finer resolution.

        inputCRS : int , default=WSG84
            The ESPRG of the coordinates.

        dataSourceName : string , default=None
            The name of the data source to use.

        Returns
        -------
            xarray.DataArray

        """

        def vectorizeLandCoverCalc(lat, lon, img, lonUpperLeft, lonResolution, latUpperLeft, latResolution):
            ilat = math.floor((lat - latUpperLeft) / latResolution)
            ilon = math.floor((lon - lonUpperLeft) / lonResolution)
            return img[ilat, ilon]

        min_pp = convertCRS(points=[[minlon, minlat]], inputCRS=WSG84, outputCRS=ITM)[0]
        max_pp = convertCRS(points=[[maxlon, maxlat]], inputCRS=WSG84, outputCRS=ITM)[0]
        x = numpy.arange(min_pp.x, max_pp.x, dxdy)
        y = numpy.arange(min_pp.y, max_pp.y, dxdy)
        xx = numpy.zeros((len(x), len(y)))
        yy = numpy.zeros((len(x), len(y)))
        for ((i, vx), (j, vy)) in product([(i, vx) for (i, vx) in enumerate(x)], [(j, vy) for (j, vy) in enumerate(y)]):
            print((i, j), end="\r")
            newpp = convertCRS(points=[[vx, vy]], inputCRS=ITM, outputCRS=inputCRS)[0]
            xx[i, j] = newpp.x
            yy[i, j] = newpp.y

        ds = self.getDataSourceData(dataSourceName)
        img = ds.GetRasterBand(1).ReadAsArray()
        lonUpperLeft, lonResolution, lonRotation, latUpperLeft, latRotation, latResolution = ds.GetGeoTransform()
        vectorizedLandCover = numpy.vectorize(lambda tlat, tlon: vectorizeLandCoverCalc(lat=tlat,
                                                                                        lon=tlon,
                                                                                        img=img,
                                                                                        lonUpperLeft=lonUpperLeft,
                                                                                        lonResolution=lonResolution,
                                                                                        latUpperLeft=latUpperLeft,
                                                                                        latResolution=latResolution))
        landcover = vectorizedLandCover(xx, yy)
        ### Transform to XArray
        i = np.arange(landcover.shape[0])
        j = np.arange(landcover.shape[1])
        xarray = xr.DataArray(
            landcover,
            coords={
                'i': i,
                'j': j,
                'lat': (['i', 'j'], xx),
                'lon': (['i', 'j'], yy),
                'landcover': (['i', 'j'], landcover)
                },
            dims=['i', 'j']
            )
        xarray.attrs['landcover_description'] = self.getCodingMap(dataSourceName)
        return xarray

    def getRoughnessAtPoint(self,lon,lat,inputCRS=WSG84, dataSourceName=None):
        """
        Get the roughness value of a specific point in the map.

        Parameters
        ----------
        lon : float
            The longitude coodinates.

        lat : float
            The latitute coordinates.

        inputCRS : int , default=WSG84
            The ESPRG of the coordinates.

        dataSourceName : string , default=None
            The name of the data source to use.

        Returns
        -------
            float

        """
        dataSourceName = self.getConfig()['defaultLandCover'] if dataSourceName is None else dataSourceName
        datasourceDocument = self.getDataSourceDocument(dataSourceName)
        landcover = self.getLandCoverAtPoint(lon=lon,lat=lat,inputCRS=inputCRS,dataSourceName=dataSourceName)

        handlerFunction = getattr(self, f"handleType{datasourceDocument['desc']['type']}")
        return handlerFunction(landcover)


    def getRoughnessFromLandcover(self,landcover,isBuilding=False,dataSourceName=None):
        """
        Adds Roughness field (z0) to landcover Xarray.

        Parameters
        ----------
        landcover : Xarray
            Landcover Xarray map in which the Roughness field will be added to.

        isBuilding : bool, default=False
            Is the landcover containts building area.

        dataSourceName : string , default=None
            The name of the data source to use.

        Returns
        -------
            xarray.DataArray

        """
        dataSourceName = self.getConfig()['defaultLandCover'] if dataSourceName is None else dataSourceName
        datasourceDocument = self.getDataSourceDocument(dataSourceName)
        if not isBuilding:
            handlerFunction = getattr(self, f"handleType{datasourceDocument['desc']['type']}")
            roughness_values = np.vectorize(handlerFunction)(landcover['landcover'])
            landcover = landcover.assign_coords(z0=(['i', 'j'], roughness_values))
        return landcover

    def getRoughness(self,minlon,minlat,maxlon,maxlat,dxdy = 30, inputCRS=WSG84, dataSourceName=None,isBuilding=False):
        """
        Returns Xarray LandCover map with Roughness (zo) field. Just as applying getLandCover and getRoughnessFromLandcover at the same time.

        Parameters
        ----------
        minlon : float
            Minimum value of longitude coodinates.

        minlat : float
            Minimum value of latitute coodinates.

        maxlon: float
            Maximum value of longitude coodinates.

        maxlat: float
            Maximum value of latitute coodinates.

        dxdy: int, default=30
            Spatial resolution of the output land cover map. Defines the step size for the grid points
            in both the x (longitude) and y (latitude) directions within the specified bounding box.
            Smaller values result in finer resolution.

        inputCRS : int , default=WSG84
            The ESPRG of the coordinates.

        dataSourceName : string , default=None
            The name of the data source to use.

        isBuilding : bool, default=False
            Is the landcover containts building area.

        Returns
        -------
            xarray.DataArray
        """
        landcover = self.getLandCover(minlon,minlat,maxlon,maxlat,dxdy = dxdy, inputCRS=inputCRS, dataSourceName=dataSourceName)
        landcover = self.getRoughnessFromLandcover(landcover,isBuilding=isBuilding,dataSourceName=dataSourceName)
        return landcover

    def handleType1(self,landcover):
        """
        Converting land type of Type-1 to roughness.
        Based on the paper:
            * https://wes.copernicus.org/articles/6/1379/2021/ table a2
            * https://doi.org/10.5194/wes-6-1379-2021 Satellite-based estimation of roughness lengths and displacement heights for wind resource modelling, Rogier Floors, Merete Badger, Ib Troen, Kenneth Grogan, and Finn-Hendrik Permien

        Parameters
        ----------
        landcover : int
            Landcover type value.

        Returns
        -------
            float
        """
        roughnessDict = {
            0 : 0.0001, # Water,  depends on waves that depends on wind
            1 : 1,       # Evergreen needleleaf forest
            2:  1,         # Evergreen broadleaf forest
            3:  1,  # Deciduous needleleaf forest
            4:  1, # Deciduous broadleaf forest
            5:  1, # Mixed forests
            6:  0.05, # Closed shrubland
            7:  0.06, # Open shrublands
            8: 0.05, # Woody savannas
            9: 0.15,  # Savannas
            10: 0.12,  # Grasslands
            11: 0.3, # Permanent wetlands
            12: 0.15, # Croplands
            13: 0.8, # Urban and built-up
            14: 0.14, # Cropland/natural vegetation mosaic
            15: 0.001, # Snow and ice
            16: 0.01, # Barren or sparsely vegetated
        }
        return roughnessDict.get(landcover,0.05)


    def getCodingMap(self,datasourceName):
        """
        Returns dictionary that maps landcover int value to string of landcover.

        Parameters
        ----------
        datasourceName: str
            Datasource type (Type-1, Type-2 etc.)

        Returns
        -------
            dict
        """
        dict = {}
        if datasourceName=='Type-1':
            dict = {
                    0: "Water",
                    1: "Evergreen needleleaf forest",
                    2: "Evergreen broadleaf forest",
                    3: "Deciduous needleleaf forest",
                    4: "Deciduous broadleaf forest",
                    5: "Mixed forests",
                    6: "Closed shrubland",
                    7: "Open shrublands",
                    8: "Woody savannas",
                    9: "Savannas",
                    10: "Grasslands",
                    11: "Permanent wetlands",
                    12: "Croplands",
                    13: "Urban and built-up",
                    14: "Cropland/natural vegetation mosaic",
                    15: "Snow and ice",
                    16: "Barren or sparsely vegetated"
                }

        return dict

# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Tue May 21 18:01:43 2024
#
# @author: nirb
# """
#
#
# from osgeo import gdal
# import math
# import matplotlib.pyplot as plt
#
# def getlc(filepath, lat, long):
#
#     ds = gdal.Open(filepath)
#     img = ds.GetRasterBand(1).ReadAsArray()
#     gt = ds.GetGeoTransform()
#
#     width = ds.RasterXSize
#     height = ds.RasterYSize
#     minx = gt[0]
#     maxx = gt[0] + width*gt[1] + height*gt[2]
#     miny = gt[3] + width*gt[4] + height*gt[5]
#     maxy = gt[3]
#     # plt.figure()
#     # plt.imshow(img)
#     # plt.show()
#     x = math.floor((lat - gt[3])/gt[5])
#     y = math.floor((long - gt[0])/gt[1])
#     return (img[x,y])
#
# def lc2roughnesslength(lc, lctype=1):
#     # https://wes.copernicus.org/articles/6/1379/2021/ table a2
#     # https://doi.org/10.5194/wes-6-1379-2021 Satellite-based estimation of roughness lengths and displacement heights for wind resource modelling, Rogier Floors, Merete Badger, Ib Troen, Kenneth Grogan, and Finn-Hendrik Permien
#
#     if lctype == 1:
#         if lc == 0: # Water
#            rl = 0.0001 # depends on waves that depends on wind
#         elif lc == 1: # Evergreen needleleaf forest
#            rl = 1.0
#         elif lc == 2: # Evergreen broadleaf forest
#            rl = 1.0
#         elif lc == 3: # Deciduous needleleaf forest
#            rl = 1.0
#         elif lc == 4: # Deciduous broadleaf forest
#            rl = 1.0
#         elif lc == 5: # Mixed forests
#            rl = 1.0
#         elif lc == 6: # Closed shrubland
#            rl = 0.05
#         elif lc == 7: # Open shrublands
#            rl = 0.06
#         elif lc == 8: # Woody savannas
#            rl = 0.05
#         elif lc == 9: # Savannas
#            rl = 0.15
#         elif lc == 10: # Grasslands
#            rl = 0.12
#         elif lc == 11: # Permanent wetlands
#            rl = 0.3
#         elif lc == 12: # Croplands
#            rl = 0.15
#         elif lc == 13: # Urban and built-up
#            rl = 0.8
#         elif lc == 14: # Cropland/natural vegetation mosaic
#            rl = 0.14
#         elif lc == 15: # Snow and ice
#            rl = 0.001
#         elif lc == 16: # Barren or sparsely vegetated
#            rl = 0.01
#         else:
#            rl = 0.05 # Bamba choice
#
#     return rl
#
# def roughnesslength2sandgrainroughness(rl):
# #Desmond, C. J., Watson, S. J., & Hancock, P. E. (2017). Modelling the wind energy resource in complex terrain and atmospheres. Numerical simulation and wind tunnel investigation of non-neutral forest canopy flow. Journal of wind engineering and industrial aerodynamics, 166, 48-60.‏
# # https://www.sciencedirect.com/science/article/pii/S0167610516300083#bib12
# # eq. 5: Equivalent sand grain roughness (m) is z0*30
#
# # we can you it for "nutkRoughWallFunction" boundary condition for Ks (sand grain roughness) parameter
# # Cs value can be set as 0.5
#
#     return rl*30.0 # return Ks value
#
# if __name__ == "__main__":
#     filename = r'lc_mcd12q1v061.t1_c_500m_s_20210101_20211231_go_epsg.4326_v20230818.tif' # 500m resolution
#     # filename = r'lc_mcd12q1v061.t2_c_500m_s_20210101_20211231_go_epsg.4326_v20230818.tif'
#     # filename = r'lc_mcd12q1v061.t5_c_500m_s_20210101_20211231_go_epsg.4326_v20230818.tif'
#     # filename = r'lc_mcd12q1v061.t1_c_500m_s_20200101_20201231_go_epsg.4326_v20230818.tif'
#
#     filepath = r'/data3/GIS_Data/LC/'+filename
#
#     lat = 31.88
#     long = 34.743
#
#     mylc = getlc(filepath, lat, long)
#     print (mylc, lc2roughnesslength(mylc))
#
