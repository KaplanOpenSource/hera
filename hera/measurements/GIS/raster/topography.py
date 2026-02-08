"""
Vector Topography Toolkit (GIS).

This module provides utilities for reading DEM (HGT/SRTM) elevation data and building
derived products such as:
- Point and point-list elevation queries
- Regular lat/lon grids (xarray) for a bounding box
- Elevation layers over grids
- STL generation from elevation rasters (for 3D printing / CAD workflows)

Notes:
- DEM files are resolved from toolkit datasources (ToolkitDataSource documents).
- Several functions depend on optional libraries (GDAL, GeoPandas, etc.).
- Coordinate systems are handled via EPSG codes (default: WGS84 / EPSG:4326).
"""

from hera import toolkit
from hera.measurements.GIS.utils import stlFactory, convertCRS, ITM, WSG84, ED50_ZONE36N, create_xarray
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


WSG84 = 4326



class TopographyToolkit(toolkit.abstractToolkit):
    """
    Toolkit for managing raster topography and elevation data.

    The TopographyToolkit provides comprehensive support for Digital Elevation Model (DEM)
    data, primarily in HGT/SRTM format. It enables elevation queries, grid generation,
    and derived products for terrain analysis.

    Key Features
    ------------
    - Point and point-list elevation queries with bilinear interpolation
    - Regular lat/lon grid generation (xarray) for bounding boxes
    - Elevation layer creation over custom grids
    - STL generation from elevation rasters for 3D printing/CAD workflows
    - Support for multiple DEM file formats (HGT, SRTM)
    - Coordinate system transformations

    Data Sources
    ------------
    DEM data sources are stored as ToolkitDataSource documents pointing to directories
    containing HGT files. Each datasource can reference one or more resource directories
    where DEM tiles are stored.

    File Naming Convention
    -----------------------
    HGT files follow the naming pattern: N{lat}E{lon}.hgt
    - lat: Integer degrees of latitude (e.g., 33 for 33°N)
    - lon: Integer degrees of longitude, zero-padded to 3 digits (e.g., 035 for 35°E)
    Example: N33E035.hgt for tile covering 33°N-34°N, 35°E-36°E

    Examples
    --------
    >>> from hera import toolkitHome
    >>> 
    >>> # Get the Topography toolkit
    >>> topo_tk = toolkitHome.getToolkit("GIS_Raster_Topography", projectName="my_project")
    >>> 
    >>> # Get elevation at a single point
    >>> elevation = topo_tk.getPointElevation(
    ...     lat=31.7683,
    ...     long=35.2137,
    ...     dataSourceName="srtm_data"
    ... )
    >>> print(f"Elevation: {elevation} meters")
    >>> 
    >>> # Get elevations for multiple points
    >>> import geopandas as gpd
    >>> from shapely.geometry import Point
    >>> points = gpd.GeoDataFrame({
    ...     'geometry': [Point(35.2, 31.7), Point(35.3, 31.8)]
    ... })
    >>> elevations = topo_tk.getPointListElevation(points)
    >>> 
    >>> # Create elevation grid
    >>> grid = topo_tk.analysis.getElevationGrid(
    ...     minlat=31.0, maxlat=32.0,
    ...     minlon=34.0, maxlon=35.0,
    ...     dxdy=0.001
    ... )
    """

    def __init__(self, projectName, filesDirectory=None):
        """
            Initialize the TopographyToolkit in a given project context.

            The toolkit is a Project-backed component (inherits from abstractToolkit) and uses
            the project's data layer for configuration and datasource discovery. The toolkit also
            attaches an analysis helper object (topographyAnalysis) under `self._analysis`.

            Parameters
            ----------
            projectName : str
                Project name to bind the toolkit to (used for config and DB-backed documents).
            filesDirectory : str or None, optional
                Output directory for files produced by toolkit operations.
                - If a string path is provided (relative or absolute), it will be created if needed.
                - If None, the toolkit will try to use the project's configured default directory.
                  If missing, it falls back to the current working directory.

            Notes
            -----
            This toolkit passes only `projectName` (not a full Project object) to the parent class.
            This avoids storing non-serializable Python objects in DB-backed documents.
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
        Resolve a DEM/HGT filename within registered datasource resource folders.

        This helper looks up the datasource payload (via `getDataSourceData`) and searches
        for the requested file under each configured resource path.

        Parameters
        ----------
        filename : str
            DEM filename to locate (e.g., "N33E035.hgt").
        dataSourceName : str
            Name of the toolkit datasource that contains one or more resource folders.

        Returns
        -------
        str
            Absolute path to the located DEM file.

        Raises
        ------
        FileNotFoundError
            If the file cannot be found under any of the configured resource paths.
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
        Return DEM elevation (meters above sea level) for a single WGS84 point.

        The function reads the relevant HGT tile, loads it via GDAL, and performs a bilinear
        interpolation at the requested (lat, lon) position.

        Parameters
        ----------
        lat : float
            Latitude in WGS84 (EPSG:4326).
        long : float
            Longitude in WGS84 (EPSG:4326).
        dataSourceName : str, optional
            Toolkit datasource name that provides the DEM resource path(s).
            If None, the datasource is taken from project config under 'defaultSRTM'.

        Returns
        -------
        float
            Elevation in meters above sea level.

        Raises
        ------
        ValueError
            If no datasource can be resolved (explicitly or from config).
        FileNotFoundError
            If the required HGT file is not found in the datasource resources.

        Notes
        -----
        - HGT naming convention is inferred from integer degrees of (lat, lon).
        - Values below -1000 are treated as invalid and replaced with 0.
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
        Compute elevation values for a list of points and return an enriched table.

        The function groups points by their required HGT tile file, loads each tile once,
        and assigns elevations to all points in that group using bilinear interpolation.

        Parameters
        ----------
        pointList : pandas.DataFrame or geopandas.GeoDataFrame or geopandas.GeoSeries
            Supported inputs:
            - GeoSeries of Points: uses point.x / point.y
            - GeoDataFrame with a 'point' column containing shapely Points
            - DataFrame with columns 'lat' and 'lon' (or attributes accessible as x.lat/x.lon)
        dataSourceName : str, optional
            Toolkit datasource name that provides the DEM resource path(s).
            If None, uses project config key 'defaultSRTM'.
        inputCRS : int, optional
            EPSG code of the input points. Currently the implementation assumes WGS84 for
            filename resolution (EPSG:4326). (Kept for API compatibility.)

        Returns
        -------
        pandas.DataFrame or geopandas.GeoDataFrame
            A table containing the original points plus:
            - filename: DEM tile filename per point
            - elevation: computed elevation value (meters)

        Raises
        ------
        ValueError
            If the datasource cannot be resolved or GeoDataFrame structure is invalid.
        FileNotFoundError
            If a required HGT file is missing.

        Notes
        -----
        Values below -1000 are treated as invalid and replaced with 0.
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
        Convert point coordinates from one CRS to another using GeoPandas.

        This helper accepts multiple input formats (array, list, DataFrame) and returns
        a GeoDataFrame re-projected into the requested output CRS.

        Parameters
        ----------
        points : list[tuple] or numpy.ndarray or pandas.DataFrame
            Points to convert:
            - list of (x, y) tuples
            - numpy array shaped (2,) or (N, 2)
            - DataFrame with x/y columns (defaults to 'x' and 'y', configurable via kwargs)
        inputCRS : int
            EPSG code of the input CRS.
        outputCRS : int
            EPSG code of the output CRS.
        **kwargs : dict
            Optional column name overrides for DataFrame inputs:
            - x: name of x column (default 'x')
            - y: name of y column (default 'y')

        Returns
        -------
        geopandas.GeoDataFrame
            GeoDataFrame in output CRS with a 'geometry' column.

        Raises
        ------
        ValueError
            If `points` type is not supported.
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
        Create an xarray grid of (lat, lon) coordinates for a bounding box.

        The bounding box is interpreted in the provided CRS. Internally, when input is WGS84,
        the corners are converted to ITM to generate a meter-based grid spacing (dxdy), then
        converted back to WGS84 lat/lon arrays.

        Parameters
        ----------
        minx, miny, maxx, maxy : float
            Bounding box coordinates. Interpretation depends on inputCRS:
            - If WGS84: (lon/lat) style values are converted to ITM for gridding
            - If ITM: values are treated as meter coordinates directly
        dxdy : float, optional
            Grid spacing:
            - meters if gridding in ITM
            - degrees if gridding directly in WGS84
        inputCRS : int, optional
            EPSG code of input CRS (default: EPSG:4326 / WGS84).

        Returns
        -------
        xarray.Dataset
            Dataset with coordinates:
            - i, j
            and variables:
            - lat[i, j], lon[i, j]

        Raises
        ------
        MemoryError
            If the generated grid exceeds 1,000,000 points (safety guard).

        Notes
        -----
        The memory guard prevents accidental creation of massive grids.
        Increase dxdy or reduce the bounding box to lower memory usage.
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
        Generate an elevation layer over a rectangular area.

        This method:
        1) Builds a (lat, lon) grid using `create_xarray(...)`
        2) Flattens the grid to a point list
        3) Computes elevation for all points via `getPointListElevation(...)`
        4) Reshapes and attaches elevation back onto the returned xarray.Dataset

        Parameters
        ----------
        minx, miny, maxx, maxy : float
            Bounding box coordinates of the area.
        dxdy : float, optional
            Grid spacing (default: 30).
        inputCRS : int, optional
            EPSG code of the input coordinates (default: EPSG:4326).
        dataSourceName : str, optional
            DEM datasource name. If None, falls back to project config 'defaultSRTM'.

        Returns
        -------
        xarray.Dataset
            Dataset containing:
            - lat[i, j]
            - lon[i, j]
            - elevation[i, j]

        Notes
        -----
        This method uses the shape of `xarray_dataset['lat']` (not Dataset.shape) because
        xarray.Dataset does not define a `.shape` attribute.
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
        Compute and attach elevation values for an existing (lat, lon) xarray grid.

        Parameters
        ----------
        xarray_dataset : xarray.Dataset
            Dataset with 'lat' and 'lon' variables defined over dimensions ['i', 'j'].
        dataSourceName : str, optional
            DEM datasource name. If None, falls back to project config 'defaultSRTM'.

        Returns
        -------
        xarray.Dataset
            The same dataset with an added 'elevation' coordinate over ['i', 'j'].

        Notes
        -----
        Uses `xarray_dataset['lat'].shape` for reshaping because xarray.Dataset does not
        have a `.shape` attribute.
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
        Create an STL representation of the elevation surface for a bounding box.

        This is a convenience wrapper that generates an elevation grid via `getElevation(...)`
        and then converts it to STL via `getElevationSTL(...)`.

        Parameters
        ----------
        minx, miny, maxx, maxy : float
            Bounding box coordinates of the area.
        dxdy : float, optional
            Grid spacing/resolution (default: 30).
        shiftx, shifty : float, optional
            Optional XY shifts in the output ITM coordinate space (useful for re-centering).
        inputCRS : int, optional
            EPSG code of input CRS (default: WGS84 / EPSG:4326).
        dataSourceName : str, optional
            DEM datasource name.
        solidName : str, optional
            Name to embed as the STL "solid" identifier.

        Returns
        -------
        str
            STL content as a string.
        """
        elevation = self.getElevation(minx=minx,miny=miny, maxx=maxx, maxy=maxy, dxdy=dxdy, inputCRS=inputCRS, dataSourceName=dataSourceName)

        return self.getElevationSTL(elevation,shiftx,shifty,solidName)


    def getElevationSTL(self,elevation,shiftx=0,shifty=0,solidName="Topography"):
        """
        Convert an elevation xarray dataset into an STL string.

        The method converts the (lat, lon) grid into ITM coordinates, reshapes X/Y arrays to
        match the elevation raster shape, and uses stlFactory().rasterToSTL(...) to generate
        the final STL content.

        Parameters
        ----------
        elevation : xarray.Dataset
            Dataset containing 'lat', 'lon', and 'elevation' variables.
        shiftx, shifty : float, optional
            Optional XY shifts (applied in ITM space).
        solidName : str, optional
            STL solid name.

        Returns
        -------
        str
            STL content as a string.

        Notes
        -----
        Uses `elevation['elevation'].shape` for reshaping because xarray.Dataset does not
        define `.shape`.
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
        Compute basic elevation statistics over a domain.

        Parameters
        ----------
        elevation : xarray.Dataset or similar
            Elevation data including X/Y coordinates and an Elevation layer.
            Expected keys/variables (as used by this implementation):
            - 'X', 'Y', 'Elevation'

        Returns
        -------
        dict
            Dictionary with:
            - xmin, xmax, ymin, ymax
            - size (domain area, assuming rectangular domain)
            - mean, std
            - domainmax, domainmaxlocation
            - domainmin, domainminlocation

        Notes
        -----
        This function assumes the elevation object exposes the expected keys and that the
        domain is small enough for simple rectangular area estimation.
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
