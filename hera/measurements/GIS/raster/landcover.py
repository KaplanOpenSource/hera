#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 18:01:43 2024

@author: nirb
"""
# MCD12Q1 - MODIS/Terra+Aqua Land Cover Type Yearly L3 Global 500m SIN Grid

# data is from https://zenodo.org/records/8367523
# type are in https://ladsweb.modaps.eosdis.nasa.gov/filespec/MODIS/6/MCD12Q1

# Land Cover Type 1: IGBP global vegetation classification scheme
# Land Cover Type 2: University of Maryland (UMD) scheme


        # DataField    Name                    Data      Dimensions
        #                                      Type

        # DataField_1  Land_Cover_Type_1       UINT8     Dimension_1
        #                                                Dimension_2

        #         Description:    land cover types (IGBP)


        #         DataField_1 HDF Attributes:

        #            long_name      STRING  1   PGE    "Land_Cover_Type_1"
        #            units          STRING  1   PGE     "class number"
        #            valid_range    uint8    2   PGE     0 16
        #            _FillValue     uint8    1   PGE     255
        #            Water          uint8    1   PGE     0
        #            Evergreen needleleaf forest 
        #                           uint8    1   PGE     1
        #            Evergreen broadleaf forest 
        #                           uint8    1   PGE     2
        #            Deciduous needleleaf forest 
        #                           uint8    1   PGE     3
        #            Deciduous broadleaf forest 
        #                           uint8    1   PGE     4
        #            Mixed forests  unit8    1   PGE     5
        #            Closed shrubland
        #                           unit8    1   PGE     6
        #            Open shrublands
        #                           unit8    1   PGE     7
        #            Woody savannas unit8    1   PGE     8
        #            Savannas       unit8    1   PGE     9
        #            Grasslands     unit8    1   PGE    10
        #            Permanent wetlands
        #                           unit8    1   PGE    11
        #            Croplands      unit8    1   PGE    12
        #            Urban and built-up
        #                           unit8    1   PGE    13
        #            Cropland/natural vegetation mosaic
        #                           unit8    1   PGE    14
        #            Snow and ice   unit8    1   PGE    15
        #            Barren or sparsely vegetated
        #                           unit8    1   PGE    16

        # DataField_2  Land_Cover_Type_2       UINT8     Dimension_1
        #                                                Dimension_2

        #         Description:    land cover types (UMD)


        #         DataField_2 HDF Attributes:
        #            long_name      STRING   1   PGE    "Land_Cover_Type_2"
        #            units          STRING   1   PGE     "class number"
        #            valid_range    uint8    2   PGE     0 16
        #            _FillValue     uint8    1   PGE     255
        #            Water          uint8    1   PGE     0
        #            Evergreen needleleaf forest 
        #                           uint8    1   PGE     1
        #            Evergreen broadleaf forest 
        #                           uint8    1   PGE     2
        #            Deciduous needleleaf forest 
        #                           uint8    1   PGE     3
        #            Decidous broadleaf forest 
        #                           uint8    1   PGE     4
        #            Mixed forests  uint8    1   PGE     5
        #            Closed shrublands
        #                            unit8   1   PGE     6
        #            Open shrubland  unit8   1   PGE     7
        #            Woody savannas  unit8   1   PGE     8	
        #            Savannas       unit8    1   PGE     9   
        #            Grasslands     uint8    1   PGE    10
        #            Croplands      unit8    1   PGE    12
        #            Urban and built-up
        #                           unit8    1   PGE    13
        #            Barren or sparsely vegetated
        #                           unit8    1   PGE    16
                                  

# Land cover may change between winter and summer
# urban area may change over the years

from osgeo import gdal
import math
import matplotlib.pyplot as plt
    
def getlc(filepath, lat, long):

    ds = gdal.Open(filepath)
    width = ds.RasterXSize
    height = ds.RasterYSize
    gt = ds.GetGeoTransform()
    minx = gt[0]
    maxx = gt[0] + width*gt[1] + height*gt[2]
    miny = gt[3] + width*gt[4] + height*gt[5] 
    maxy = gt[3] 
    img = ds.GetRasterBand(1).ReadAsArray()
    # plt.figure()
    # plt.imshow(img)
    # plt.show()
    x = math.floor((lat - gt[3])/gt[5])
    y = math.floor((long - gt[0])/gt[1])
    return (img[x,y])

def lc2roughnesslength(lc, lctype=1):
    # https://wes.copernicus.org/articles/6/1379/2021/ table a2
    if lctype == 1:
        if lc == 0: # Water
           rl = 0.0001 # depends on waves that depends on wind
        elif lc == 1: # Evergreen needleleaf forest
           rl = 1.0
        elif lc == 2: # Evergreen broadleaf forest 
           rl = 1.0 
        elif lc == 3: # Deciduous needleleaf forest 
           rl = 1.0
        elif lc == 4: # Deciduous broadleaf forest 
           rl = 1.0 
        elif lc == 5: # Mixed forests
           rl = 1.0 
        elif lc == 6: # Closed shrubland
           rl = 0.05
        elif lc == 7: # Open shrublands
           rl = 0.06
        elif lc == 8: # Woody savannas
           rl = 0.05
        elif lc == 9: # Savannas
           rl = 0.15
        elif lc == 10: # Grasslands
           rl = 0.12
        elif lc == 11: # Permanent wetlands
           rl = 0.3
        elif lc == 12: # Croplands
           rl = 0.15
        elif lc == 13: # Urban and built-up
           rl = 0.8
        elif lc == 14: # Cropland/natural vegetation mosaic
           rl = 0.14
        elif lc == 15: # Snow and ice
           rl = 0.001
        elif lc == 16: # Barren or sparsely vegetated
           rl = 0.01
        else:
           rl = 0.05 # Bamba choice
    
    return rl
    

if __name__ == "__main__":
    filename = r'lc_mcd12q1v061.t1_c_500m_s_20210101_20211231_go_epsg.4326_v20230818.tif' # 500m resolution
    # filename = r'lc_mcd12q1v061.t2_c_500m_s_20210101_20211231_go_epsg.4326_v20230818.tif'
    # filename = r'lc_mcd12q1v061.t5_c_500m_s_20210101_20211231_go_epsg.4326_v20230818.tif'
    # filename = r'lc_mcd12q1v061.t1_c_500m_s_20200101_20201231_go_epsg.4326_v20230818.tif'
    
    filepath = r'/data3/GIS_Data/LC/'+filename
    
    lat = 31.88
    long = 34.743
    
    mylc = getlc(filepath, lat, long)
    print (mylc, lc2roughnesslength(mylc))
    print (getlc(filepath, lat, long))
        
