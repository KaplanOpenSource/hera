"""
    Just pieces of code from the vector.topography.

    Should be organized when needed.

"""
from hera import toolkit
from hera.utils import stlFactory, convertCRS, ITM, WSG84, ED50_ZONE36N, create_xarray
from hera.utils.logging import get_classMethod_logger

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
import geopandas as gpd
from shapely.geometry import Point


WSG84 = 4326



class TopographyToolkit(toolkit.abstractToolkit):

    def __init__(self, projectName, filesDirectory=None):
        """
        Initializes the TopographyToolkit.

        Parameters
        ----------
        projectName : str
            The project Name that the toolkit is initialized on.
        filesDirectory : str or None
            The path to save region files when they are created.

            - If str, it represents a path (relative or absolute) to save the files. The directory is created automatically.
            - If None, tries to get the default path of the project from the config. If it does not
              exist, then uses the current working directory.
        """

        # Important change:
        # Instead of passing a full Project object to the parent class (Toolkit),
        # we now pass only the project name (as a string).
        # This is necessary because MongoDB expects simple types like strings,
        # and cannot serialize full complex Python objects (like Project instances).
        super().__init__(projectName=projectName, toolkitName='TopographyToolkit', filesDirectory=filesDirectory)

        # Initialize the analysis module for topography calculations
        self._analysis = topographyAnalysis(self)
    def findElevationFile(self, filename, dataSourceName):
        """
        Attempts to find the .hgt file in one of the registered resource folders.
        Supports both single path or list of paths.
        """
        resources = self.getDataSourceData(dataSourceName)
        if isinstance(resources, str):
            resources = [resources]

        for folder in resources:
            candidate = os.path.join(folder, filename)
            if os.path.exists(candidate):
                return candidate

        raise FileNotFoundError(f"{filename} not found in any of: {resources}")


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
        fheight = self.findElevationFile(filename, dataSourceName)

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
            fheight = self.findElevationFile(grpid, dataSourceName)
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



    def convertPointsCRS(self, points, inputCRS, outputCRS, **kwargs):
        """
        Convert list/array/DataFrame of points from one CRS to another, using GeoPandas.

        Parameters
        ----------
        points : list of tuples, numpy array, or pandas.DataFrame
            Points to convert.

        inputCRS : int
            EPSG code of input coordinate system.

        outputCRS : int
            EPSG code of output coordinate system.

        Returns
        -------
        geopandas.GeoDataFrame
            Converted points with 'geometry' column.
        """

        if isinstance(points, np.ndarray):
            if points.ndim == 1:
                gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy([points[0]], [points[1]]))
            else:
                gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(points[:, 0], points[:, 1]))

        elif isinstance(points, pd.DataFrame):
            x_col = kwargs.get("x", "x")
            y_col = kwargs.get("y", "y")
            gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(points[x_col], points[y_col]))

        elif isinstance(points, list):
            if len(points) == 1:
                gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy([points[0][0]], [points[0][1]]))
            else:
                gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy([x[0] for x in points], [x[1] for x in points]))

        else:
            raise ValueError(f"Unsupported type: {type(points)}")

        gdf.set_crs(inputCRS, inplace=True)
        return gdf.to_crs(outputCRS)

    def create_xarray(self, minx, miny, maxx, maxy, dxdy=30, inputCRS=WSG84):
        """
        Create an xarray grid of lat/lon points within the given bounding box.

        Parameters
        ----------
        minx, miny, maxx, maxy : float
            Bounding box coordinates.
        dxdy : float
            Grid spacing (in meters for ITM or degrees for WGS84).
        inputCRS : int
            EPSG code of the coordinate system (default is WGS84).

        Returns
        -------
        xarray.Dataset
            Dataset containing lat/lon coordinates over the i, j grid.

        Notes
        -----
        ✅ Added memory guard: throws an error if grid is larger than 1,000,000 points.
        """

        if inputCRS == WSG84:
            # Convert lat/lon to ITM using internal conversion method
            min_pp = self.convertPointsCRS(points=[[miny, minx]], inputCRS=inputCRS, outputCRS=ITM).geometry.iloc[0]
            max_pp = self.convertPointsCRS(points=[[maxy, maxx]], inputCRS=inputCRS, outputCRS=ITM).geometry.iloc[0]
        else:
            min_pp = Point(minx, miny)
            max_pp = Point(maxx, maxy)

        x = np.arange(min_pp.x, max_pp.x, dxdy)
        y = np.arange(min_pp.y, max_pp.y, dxdy)

        # 🔒 Memory guard
        if len(x) * len(y) > 1_000_000:
            raise MemoryError(f"Too many grid points: {len(x)} x {len(y)} = {len(x) * len(y)}. "
                              f"Increase dxdy or reduce area.")

        xx, yy = np.meshgrid(x, y[::-1])
        grid_points = pd.DataFrame({'x': xx.ravel(), 'y': yy.ravel()})

        gdf = gpd.GeoDataFrame(
            grid_points,
            geometry=gpd.points_from_xy(grid_points['x'], grid_points['y']),
            crs=ITM
        )

        # Convert back to WGS84
        gdf_transformed = gdf.to_crs(WSG84)
        gdf_transformed['lat'] = gdf_transformed.geometry.y
        gdf_transformed['lon'] = gdf_transformed.geometry.x

        lat_grid = gdf_transformed['lat'].values.reshape(xx.shape)
        lon_grid = gdf_transformed['lon'].values.reshape(xx.shape)

        i = np.arange(xx.shape[0])
        j = np.arange(xx.shape[1])

        return xr.Dataset(
            coords={
                'i': i,
                'j': j,
            },
            data_vars={
                'lat': (['i', 'j'], lat_grid),
                'lon': (['i', 'j'], lon_grid)
            }
        )

    def getElevation(self, minx, miny, maxx, maxy, dxdy=30, inputCRS=WSG84, dataSourceName=None):
        """
        Generates elevation data over a rectangular area using given resolution.

        Parameters
        ----------
        minx, miny, maxx, maxy : float
            Bounding box coordinates of the area to analyze.
        dxdy : float
            Resolution (spacing) in coordinate units (default is 30).
        inputCRS : CRS
            Coordinate reference system of the input coordinates.
        dataSourceName : str, optional
            Name of the data source to fetch elevations from.

        Returns
        -------
        xarray.Dataset
            Dataset with 'lat', 'lon' and calculated 'elevation' layer.

        Notes
        -----
        🔧 FIXED: Replaced xarray_dataset.shape access with .shape from 'lat' variable,
        since xarray.Dataset has no shape attribute.
        """

        # Create initial lat/lon grid
        xarray_dataset = self.create_xarray(minx, miny, maxx, maxy, dxdy, inputCRS)

        # Flatten to point list
        pointList = pd.DataFrame({
            'lat': xarray_dataset['lat'].values.flatten(),
            'lon': xarray_dataset['lon'].values.flatten()
        })

        # Get elevations
        elevation_df = self.getPointListElevation(pointList, dataSourceName)

        # 🔧 FIX: Use shape of lat array instead of dataset shape
        i_dim, j_dim = xarray_dataset['lat'].shape

        # Add elevation coordinate to dataset
        xarray_dataset = xarray_dataset.assign_coords(
            elevation=(('i', 'j'), elevation_df['elevation'].values.reshape(i_dim, j_dim))
        )

        return xarray_dataset

    def getElevationOfXarray(self, xarray_dataset, dataSourceName=None):
        """
        Computes elevation values for each (lat, lon) point in an xarray dataset
        and returns the same dataset with an added 'elevation' coordinate.

        Parameters
        ----------
        xarray_dataset : xarray.Dataset
            Dataset with 'lat' and 'lon' variables defined over dimensions ['i', 'j']
        dataSourceName : str, optional
            Name of the data source to use. If not given, will be extracted from config.

        Returns
        -------
        xarray.Dataset
            Same dataset with added 'elevation' coordinate over ['i', 'j']

        Notes
        -----
        🔧 FIXED: previously tried to access `xarray_dataset.shape` which doesn't exist on Dataset.
        Now correctly gets shape from 'lat' variable.
        """

        # Convert the xarray dataset into a flat DataFrame of points
        pointList = pd.DataFrame({
            'lat': xarray_dataset['lat'].values.flatten(),
            'lon': xarray_dataset['lon'].values.flatten()
        })

        # Get elevation values for the points
        elevation_df = self.getPointListElevation(pointList, dataSourceName)

        # 🔧 FIX: Replace xarray_dataset.shape with actual shape from 'lat' variable
        i_dim, j_dim = xarray_dataset['lat'].shape

        # Assign elevation as a new coordinate in the dataset
        xarray_dataset = xarray_dataset.assign_coords(
            elevation=(['i', 'j'], elevation_df['elevation'].values.reshape(i_dim, j_dim))
        )

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
        """
        Generates STL string from elevation dataset.

        Parameters
        ----------
        elevation : xarray.Dataset
            Dataset with 'lat', 'lon' and 'elevation' coordinates.
        shiftx, shifty : float
            Optional shifts in X and Y.
        solidName : str
            Name of the STL solid.

        Returns
        -------
        str
            STL content as string.

        Notes
        -----
        🔧 FIXED: Accessed shape using elevation['elevation'].shape instead of elevation.shape
        because xarray.Dataset has no attribute 'shape'.
        """
        grid_points = pd.DataFrame({
            'x': elevation['lat'].values.flatten(),
            'y': elevation['lon'].values.flatten()
        })
        gdf = gpd.GeoDataFrame(
            grid_points,
            geometry=gpd.points_from_xy(grid_points['y'], grid_points['x']),
            crs=WSG84
        )
        gdf_transformed = gdf.to_crs(ITM)
        gdf_transformed['y'] = gdf_transformed.geometry.y
        gdf_transformed['x'] = gdf_transformed.geometry.x

        # 🔧 FIX: Use actual shape from elevation['elevation'] instead of elevation.shape
        i_dim, j_dim = elevation['elevation'].shape
        xx = gdf_transformed['x'].values.reshape(i_dim, j_dim)
        yy = gdf_transformed['y'].values.reshape(i_dim, j_dim)
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
