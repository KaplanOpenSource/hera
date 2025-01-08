"""
    Just pieces of code from the vector.topography.

    Should be organized when needed.

"""
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
import numpy as np
from hera import toolkitHome
import pandas as pd
import xarray as xr


class TopographyToolkit(toolkit.abstractToolkit):

    def __init__(self, projectName, filesDirectory=None):
        """
            Initializes vector data toolkit.

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
        super().__init__(projectName=projectName, toolkitName = 'TopographyToolkit', filesDirectory=filesDirectory)
        self._analysis = topographyAnalysis(self)

    def getPointElevation(self,lat, long,dataSourceName=None):
        """
            get the elevation above sea level in meters from DEM for a point
            SRTM30 means 30 meters resolution

        Parameters
        ----------
        lat : float
            Latitiute in the WSG projection
        long : float
            Longitude in the WSG projection

        dataSourceName: string
            The name of the datasources loaded. Use getDataSourceList to see which datasource were loaded.

        Returns
        -------
        eleveation above sea level

        How to use
        ----------
        from hera import toolkitHome
        tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_RASTER)
        tk.getPointElevation(lat = 33.459,long = 35.755)

        This should return ~  ~820 according to amud anan.
        """
        logger = get_classMethod_logger(self,"getPointElevation")

        filename = 'N'+str(int(lat))+'E'+str(int(long)).zfill(3)+'.hgt'
        logger.info(f"Getting elevation from file {filename}")

        if dataSourceName is None:
            dataSourceName = self.getConfig()['defaultSRTM']

        if dataSourceName is None:
            raise ValueError(f"The datasource is not found!")

        logger.debug(f"Getting the data source {dataSourceName}")
        fheight = os.path.join(self.getDataSourceData(dataSourceName),filename)

        ds =  gdal.Open(fheight)
        myarray = numpy.array(ds.GetRasterBand(1).ReadAsArray())
        myarray[myarray < -1000] = 0  # deal with problematic points
        gt = ds.GetGeoTransform()
        rastery = (long - gt[0]) / gt[1]
        rasterx = (lat - gt[3]) / gt[5]
        height11 = myarray[int(rasterx), int(rastery)]
        if (int(rasterx)+1>=myarray.shape[0]) or (int(rastery)+1>=myarray.shape[1]):
            height12 = height11
            height21 = height11
            height22 = height11
        else:
            height12 = myarray[int(rasterx) + 1, int(rastery)]
            height21 = myarray[int(rasterx), int(rastery) + 1]
            height22 = myarray[int(rasterx) + 1, int(rastery) + 1]
        height1 = (1. - (rasterx - int(rasterx))) * height11 + (rasterx - int(rasterx)) * height12
        height2 = (1. - (rasterx - int(rasterx))) * height21 + (rasterx - int(rasterx)) * height22
        elevation = (1. - (rastery - int(rastery))) * height1 + (rastery - int(rastery)) * height2

        return elevation

    def getPointListElevation(self, pointList, dataSourceName=None, inputCRS=WSG84):
        """
            Return the elevation of the point list.

            For now we assume that pointList is in WSG84 coordinates.

        Parameters
        ----------
        latList
        longList
        dataSourceName

        Returns
        -------

        """
        totalNumberOfPoints = pointList.shape[0]
        logger = get_classMethod_logger(self, "getPointListElevation")
        logger.info(f"Computing the elevation for a list of {totalNumberOfPoints} points")
        logger.debug("Computing the file name for each point")


        if isinstance(pointList,geopandas.geoseries.GeoSeries):
            pointList = pointList.to_frame("point").assign(filename=pointList.apply(lambda row: 'N' + str(int(row.y)) + 'E' + str(int(row.x)).zfill(3) + '.hgt' ),
                                                           lon=pointList.x,
                                                           lat=pointList.y,
                                                           elevation=0)

        elif isinstance(pointList,geopandas.geodataframe.GeoDataFrame):
            if 'point' not in pointList:
                error = "GeoDataFrame must contain a field 'point' that contain the points of interest"
                logger.error(error)
                raise ValueError(error)
            pointList = pointList.assign(filename=pointList.apply(lambda row: 'N' + str(int(row.point.y)) + 'E' + str(int(row.point.x)).zfill(3) + '.hgt', axis=1),
                                         lon=pointList.point.x,
                                         lat=pointList.point.y,
                                         elevation=0)
        else:
            pointList = pointList.assign(filename=pointList.apply(lambda x: 'N' + str(int(x.lat)) + 'E' + str(int(x.lon)).zfill(3) + '.hgt' ,axis=1),elevation=0)


        if dataSourceName is None:
            dataSourceName = self.getConfig()['defaultSRTM']
        logger.info(f"Getting the elevation for the points. Using datasource {dataSourceName}")

        if dataSourceName is None:
            err = "The datasource is not found!"
            logger.error(err)
            raise ValueError(err)

        processed = 0

        for grpid,grp in pointList.groupby("filename"):
            logger.info(f"\tProcessing the group {grpid} with {grp.shape[0]}. Processed so far {processed}/{totalNumberOfPoints}")
            fheight = os.path.join(self.getDataSourceData(dataSourceName), grpid)
            logger.debug(f"Getting data from the file {fheight}")
            ds = gdal.Open(fheight)
            myarray = numpy.array(ds.GetRasterBand(1).ReadAsArray())
            myarray[myarray < -1000] = 0  # deal with problematic points
            gt = ds.GetGeoTransform()
            for itemindx,item in grp.iterrows():
                rastery = (item.lon - gt[0]) / gt[1]
                rasterx = (item.lat - gt[3]) / gt[5]
                height11 = myarray[int(rasterx), int(rastery)]
                if (int(rasterx)+1>=myarray.shape[0]) or (int(rastery)+1>=myarray.shape[1]):
                    height12 = height11
                    height21 = height11
                    height22 = height11
                else:
                    height12 = myarray[int(rasterx) + 1, int(rastery)]
                    height21 = myarray[int(rasterx), int(rastery) + 1]
                    height22 = myarray[int(rasterx) + 1, int(rastery) + 1]
                height1 = (1. - (rasterx - int(rasterx))) * height11 + (rasterx - int(rasterx)) * height12
                height2 = (1. - (rasterx - int(rasterx))) * height21 + (rasterx - int(rasterx)) * height22
                elevation = (1. - (rastery - int(rastery))) * height1 + (rastery - int(rastery)) * height2
                pointList.loc[itemindx,'elevation'] = elevation
                processed += grp.shape[0]
        return pointList

    def getElevation(self,minx,miny,maxx,maxy,dxdy=30, inputCRS=WSG84, dataSourceName=None):
        if inputCRS == WSG84:
            min_pp = convertCRS(points=[[miny, minx]], inputCRS=WSG84, outputCRS=ITM)[0]
            max_pp = convertCRS(points=[[maxy, maxx]], inputCRS=WSG84, outputCRS=ITM)[0]
        else:
            min_pp = convertCRS(points=[[minx, miny]], inputCRS=ITM, outputCRS=ITM)[0]
            max_pp = convertCRS(points=[[maxx, maxy]], inputCRS=ITM, outputCRS=ITM)[0]

        x = numpy.arange(min_pp.x, max_pp.x, dxdy)
        y = numpy.arange(min_pp.y, max_pp.y, dxdy)
        xx = numpy.zeros((len(x), len(y)))
        yy = numpy.zeros((len(x), len(y)))
        for ((i, vx), (j, vy)) in product([(i, vx) for (i, vx) in enumerate(x)],
                                          [(j, vy) for (j, vy) in enumerate(y[::-1])]):
            print((i, j), end="\r")
            newpp = convertCRS(points=[[vx, vy]], inputCRS=ITM, outputCRS=WSG84)[0]
            lat = newpp.y
            lon = newpp.x
            xx[i, j] = lat
            yy[i, j] = lon

        pointList = pd.DataFrame({'lat': xx.flatten(), 'lon': yy.flatten()})
        elevation_df = self.getPointListElevation(pointList,dataSourceName)
        i = np.arange(xx.shape[0])
        j = np.arange(xx.shape[1])

        return xr.DataArray(
            coords={
                'i': i,
                'j': j,
                'lat': (['i', 'j'], xx),
                'lon': (['i', 'j'], yy),
                'elevation': (['i', 'j'], elevation_df['elevation'].values.reshape(xx.shape[0], xx.shape[1])),
                'dxdy': dxdy
            },
            dims=['i', 'j']
        )

    def getElevationOfXarray(self,xarray_dataset,dataSourceName=None):
        pointList = pd.DataFrame({'lat': xarray_dataset['lat'].values.flatten(), 'lon': xarray_dataset['lon'].values.flatten()})
        elevation_df = self.getPointListElevation(pointList, dataSourceName)
        xarray_dataset = xarray_dataset.assign_coords(elevation=(['i', 'j'], elevation_df['elevation'].values.reshape(xarray_dataset.shape[0], xarray_dataset.shape[1])))
        return xarray_dataset

    def createElevationSTL(self, minx, miny, maxx, maxy, dxdy = 30,shiftx=0,shifty=0,inputCRS=WSG84, dataSourceName=None, solidName="Topography"):
        """
            Return the STL string from xarray dataset with the following fields:
        Parameters
        ----------
        lowerleft_point : float
                The lower left corner
        upperright_point: float
                The upper right corner

        dxdy : float
                The resolution in m (default m).
        inputCRS : The ESPG code of the input projection.

        outputCRS : The ESPG code of the output projection.
                    [Default ITM]

        shiftx : Used when one wants to set another point as origin center

        shifty : Used when one wants to set another point as origin center

        Returns
        -------

        """
        elevation = self.getElevation(minx=minx,miny=miny, maxx=maxx, maxy=maxy, dxdy=dxdy, inputCRS=inputCRS, dataSourceName=dataSourceName)

        return self.getElevationSTL(elevation,shiftx,shifty,solidName)

    def getElevationSTL(self,elevation,shiftx=0,shifty=0,solidName="Topography"):
        # Convert to ITM coordinates for STL file
        xx = numpy.zeros((elevation.shape[0],elevation.shape[1]))
        yy = numpy.zeros((elevation.shape[0],elevation.shape[1]))
        for i in range(elevation.shape[0]):
            for j in range(elevation.shape[1]):
                ITM_pp = convertCRS([[elevation[i, j].lon.values, elevation[i, j].lat.values]], inputCRS=WSG84, outputCRS=ITM)[0]
                xx[i, j] = ITM_pp.x
                yy[i, j] = ITM_pp.y

        stlstr = stlFactory().rasterToSTL(xx - shiftx, yy - shifty, elevation['elevation'].values, solidName=solidName)
        return stlstr

class topographyAnalysis:

    datalayer = None

    def __init__(self,datalayer):
        self.datalayer = datalayer

    def calculateStastics(self,elevation):
        """
            calculate the domain elevation statistics
        Parameters
        ----------
        elevation - elevation data including coordinates

        Returns
        -------
        domain size (in meters for ITM, assuming small domain  (rectangle)
        max elevation + location
        min elevation + location
        mean elevation
        std of elevation

        example
        -------
        after creating a project (hera-project project create project-name)
        from hera import toolkitHome
        tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_RASTER_TOPOGRAPHY,projectName="topotest")
        elevation = tk.getDomainElevation(xmin=34.529,xmax=34.531,ymin=31.160,ymax=31.162)
        print(tk.analysis.calculateStastics(elevation))
        """
        xmin = elevation['X'].values.min()
        xmax = elevation['X'].values.max()
        ymin = elevation['Y'].values.min()
        ymax = elevation['Y'].values.max()
        domainsize = (xmax - xmin) * (ymax - ymin)
        domainmean = elevation['Elevation'].values.mean()
        domainstd = elevation['Elevation'].values.std()
        domainmax = elevation['Elevation'].values.max()
        maxpos = numpy.argmax(elevation['Elevation'].values.ravel())
        domainmaxlocation = (elevation['X'].values.ravel()[maxpos],elevation['Y'].values.ravel()[maxpos])
        domainmin = elevation['Elevation'].values.min()
        minpos = numpy.argmin(elevation['Elevation'].values.ravel())
        domainminlocation = (elevation['X'].values.ravel()[minpos],elevation['Y'].values.ravel()[minpos])
        elevation_statistics = {
                "xmin": xmin,
                "xmax": xmax,
                "ymin": ymin,
                "ymax": ymax,
                "size": domainsize,
                "mean": domainmean,
                "std": domainstd,
                "domainmax": domainmax,
                "domainmaxlocation": domainmaxlocation,
                "domainmin": domainmin,
                "domainminlocation": domainminlocation,
                }
        return elevation_statistics
