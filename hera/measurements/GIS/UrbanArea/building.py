
import geopandas as gpd
from shapely.geometry import Point
import shapely.affinity
import numpy as np
import numpy
import math
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from shapely.geometry import LineString



class Building(gpd.GeoDataFrame):

    polygon = None
    bldHeight = 0



    def __init__(self, _polygon,_bldHeight ):
        self.polygon = _polygon
        self.bldHeight = _bldHeight




    def PlotBuilding(self):
        g = gpd.GeoSeries([self.polygon])
        gdf = gpd.GeoDataFrame(geometry=g)
        gdf.plot()



    def A_f2(self,MeteoAngle):
        A_f = 0
        p1 = self.polygon
        g = gpd.GeoSeries([p1])

        gdf = gpd.GeoDataFrame(geometry=g)
        gdf['angle'] = [MeteoAngle]

        for index, row in gdf.iterrows():
            rotated = Polygon(shapely.affinity.rotate(row['geometry'], row['angle']))
            gdf.loc[index, 'geometry'] = rotated
            bounds = rotated.bounds
            xMin = bounds[0]
            xMax = bounds[2]
            A_f = (xMax - xMin) * self.bldHeight
        return A_f
        gdf.plot()











