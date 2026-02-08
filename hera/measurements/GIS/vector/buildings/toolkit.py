import numpy
import os
import geopandas
from ..toolkit import VectorToolkit
from .analysis import analysis

try:
    import FreeCAD
    import Part
    import Mesh
except ImportError as e:
    print("FreeCAD module is not installed in this environment. Cannot convert to STL")

import matplotlib.pyplot as plt
from .....utils.logging import get_classMethod_logger
from ...utils import WGS84, ITM, ED50_ZONE36N
from ..... import toolkitHome
import geopandas as gpd
from shapely.geometry import *


class BuildingsToolkit(VectorToolkit):
    """
    Toolkit for managing building geometry and height data.

    The BuildingsToolkit handles vector building data in GeoPandas format, providing
    capabilities for spatial queries, height extraction, STL conversion, and integration
    with topography data.

    Key Features
    ------------
    - Building geometry management (polygons, footprints)
    - Height data extraction and processing
    - Integration with raster topography for elevation data
    - STL file generation for 3D visualization (requires FreeCAD)
    - Spatial filtering and region extraction
    - Building height estimation from number of levels

    Data Structure
    --------------
    Building data is stored as GeoPandas GeoDataFrame with the following expected columns:
    - **geometry**: Polygon geometry representing building footprint
    - **Building height column**: Column name specified in BuildingHeightColumn.
      Default column name is 'BLDG_HT' (height in meters)
    - **Land height column**: Column name specified in LandHeightColumn.
      Default column name is 'HT_LAND' (ground elevation in meters)

    Data Sources
    ------------
    Building data sources are stored as GeoPandas GeoDataFrames. Each datasource
    can contain multiple buildings with their geometric and attribute information.

    Examples
    --------
    >>> from hera import toolkitHome
    >>> from shapely.geometry import Polygon
    >>> 
    >>> # Get the Buildings toolkit
    >>> buildings_tk = toolkitHome.getToolkit("GIS_Buildings", projectName="my_project")
    >>> 
    >>> # Get buildings from a rectangular region
    >>> buildings = buildings_tk.getBuildingsFromRectangle(
    ...     minx=34.0, miny=31.0,
    ...     maxx=35.0, maxy=32.0,
    ...     dataSourceName="city_buildings",
    ...     withElevation=True
    ... )
    >>> 
    >>> # Get building heights from topography
    >>> buildings_with_elevation = buildings_tk.getBuildingHeightFromRasterTopographyToolkit(
    ...     buildingData=buildings,
    ...     topographyDataSource="srtm_data"
    ... )
    >>> 
    >>> # Convert buildings to STL format (requires FreeCAD)
    >>> buildings_tk.buildingsGeopandasToSTLRasterTopography(
    ...     buildingData=buildings,
    ...     buildingHeightColumn="BLDG_HT",
    ...     buildingElevationColumn="HT_LAND",
    ...     outputFileName="/path/to/output.stl",
    ...     flatTerrain=False
    ... )
    """
    def __init__(self, projectName, filesDirectory=None):

        super().__init__(projectName=projectName, toolkitName="Buildings", filesDirectory=filesDirectory)
        self._analysis = analysis(dataLayer=self)

    def getBuildingHeightFromRasterTopographyToolkit(self, buildingData, topographyDataSource=None):
        """
        Extract topography elevation for each building's centroid from raster topography data.

        This method uses the raster topography toolkit to query elevation values at the
        centroid of each building polygon. The elevation is added as a new column to the
        building GeoDataFrame.

        Parameters
        ----------
        buildingData : geopandas.GeoDataFrame
            GeoDataFrame containing building polygons. Must have a 'geometry' column
            with Polygon geometries. The geometries are converted to WGS84 for
            elevation queries.

        topographyDataSource : str, optional
            Name of the datasource in the raster topography toolkit to use for
            elevation queries. If None, uses the default datasource configured in
            the topography toolkit. Default is None.

        Returns
        -------
        geopandas.GeoDataFrame
            Input GeoDataFrame with an additional 'evaluation' column containing
            elevation values (meters above sea level) for each building's centroid.

        Examples
        --------
        >>> # Get buildings
        >>> buildings = buildings_tk.getBuildingsFromRectangle(
        ...     minx=34.0, miny=31.0, maxx=35.0, maxy=32.0
        ... )
        >>> 
        >>> # Add elevation data
        >>> buildings_with_elevation = buildings_tk.getBuildingHeightFromRasterTopographyToolkit(
        ...     buildingData=buildings,
        ...     topographyDataSource="srtm_elevation"
        ... )
        >>> 
        >>> # Check elevation values
        >>> print(buildings_with_elevation['evaluation'].head())
        """
        topotk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_RASTER_TOPOGRAPHY, projectName=self.projectName)
        elevations = topotk.getPointListElevation(buildingData.centroid.to_crs(WGS84))
        return buildingData.join(elevations)

    def buildingsGeopandasToSTLRasterTopography(self,
                                                buildingData,
                                                buildingHeightColumn,
                                                buildingElevationColumn,
                                                outputFileName,
                                                flatTerrain=False,
                                                referenceTopography=0,
                                                nonFlatTopographyShift=10):
        """
        Convert building polygons to STL 3D mesh format using FreeCAD.

        This method creates 3D solid models of buildings by extruding their footprints
        to their specified heights. Buildings are positioned according to topography
        elevation data or a flat reference plane.

        Parameters
        ----------
        buildingData : geopandas.GeoDataFrame
            GeoDataFrame containing building polygons with height and elevation columns.
            Each building must have a valid Polygon geometry.

        buildingHeightColumn : str
            Name of the column containing building heights in meters. This value
            determines the vertical extent of each building.

        buildingElevationColumn : str
            Name of the column containing ground elevation values in meters. Used
            to position buildings vertically when flatTerrain=False.

        outputFileName : str
            Absolute path to the output STL file. The directory must exist or be
            writable. File will be overwritten if it exists.

        flatTerrain : bool, optional
            If True, all buildings are placed on a flat reference plane at
            referenceTopography elevation. If False, uses buildingElevationColumn
            for individual building placement. Default is False.

        referenceTopography : float, optional
            Reference elevation in meters when flatTerrain=True. All buildings
            are placed with their base at this elevation. Default is 0.

        nonFlatTopographyShift : float, optional
            Vertical offset in meters applied to building base elevation when
            flatTerrain=False. Positive values shift buildings upward. Default is 10.

        Returns
        -------
        None
            The STL file is written to disk. No return value.

        Raises
        ------
        ValueError
            If FreeCAD is not installed or not found in PYTHONPATH.
        KeyError
            If required columns are missing from buildingData.

        Notes
        -----
        - Requires FreeCAD to be installed and accessible via Python import
        - This is a low-level procedure. For easier usage, consider using
          regionToSTL methods from the analysis layer
        - Large datasets may take significant time to process
        - Progress is logged every 100 buildings

        Examples
        --------
        >>> # Convert buildings with topography
        >>> buildings_tk.buildingsGeopandasToSTLRasterTopography(
        ...     buildingData=buildings,
        ...     buildingHeightColumn="BLDG_HT",
        ...     buildingElevationColumn="HT_LAND",
        ...     outputFileName="/path/to/buildings.stl",
        ...     flatTerrain=False,
        ...     nonFlatTopographyShift=10
        ... )
        >>> 
        >>> # Convert buildings on flat terrain
        >>> buildings_tk.buildingsGeopandasToSTLRasterTopography(
        ...     buildingData=buildings,
        ...     buildingHeightColumn="height",
        ...     buildingElevationColumn="elevation",
        ...     outputFileName="/path/to/flat_buildings.stl",
        ...     flatTerrain=True,
        ...     referenceTopography=100.0
        ... )
        """
        logger = get_classMethod_logger(self, "geoPandasToSTL")
        logger.info(f"Converting {len(buildingData)} to STL. Using {'flat' if flatTerrain else 'topography'} settings")

        try:
            FreeCADDOC = FreeCAD.newDocument("Unnamed")
        except:
            err = "FreeCAD not found. Install before using this function and add to the PYTHONPATH"
            raise ValueError(err)

        for indx, building in buildingData.iterrows():  # converting al the buildings

            try:
                walls = building['geometry'].exterior.xy
                walls = numpy.stack(walls).T
            except:
                continue

            if indx % 100 == 0:
                logger.debug(f"{indx}/{len(buildingData)} shape file is executed")

            wallsheight = building[buildingHeightColumn]
            altitude = referenceTopography if flatTerrain else building[buildingElevationColumn] - nonFlatTopographyShift
            logger.debug(f" Building -- {indx} --")
            logger.debug("newSketch = FreeCADDOC.addObject('Sketcher::SketchObject', 'Sketch" + str(indx) + "')")
            newSketch = FreeCADDOC.addObject('Sketcher::SketchObject', 'Sketch' + str(indx))

            logger.debug(
                f"newSketch.Placement = FreeCAD.Placement(FreeCAD.Vector(0.000000, 0.000000, {altitude}), FreeCAD.Rotation(0.000000, 0.000000, 0.000000, 1.000000))")
            newSketch.Placement = FreeCAD.Placement(FreeCAD.Vector(0.000000, 0.000000, altitude),  # 2*k-1
                                                    FreeCAD.Rotation(0.000000, 0.000000, 0.000000, 1.000000))

            for xy0, xy1 in zip(walls[:-1], walls[1:]):
                logger.debug(
                    f"newSketch.addGeometry(Part.LineSegment(FreeCAD.Vector({xy0[0]}, {xy0[1]}, {altitude}),FreeCAD.Vector({xy1[0]}, {xy1[1]}, {altitude})))")
                newSketch.addGeometry(
                    Part.LineSegment(FreeCAD.Vector(xy0[0], xy0[1], altitude),
                                     FreeCAD.Vector(xy1[0], xy1[1], altitude)))

            logger.debug("FreeCADDOC.addObject('Part::Extrusion', 'building" + str(indx) + "')")
            newPad = FreeCADDOC.addObject("Part::Extrusion", "building" + str(indx))
            logger.debug("newPad.Base = newSketch")
            newPad.Base = newSketch
            buildingTopAltitude = wallsheight + nonFlatTopographyShift # the paddign is from the bottom of the buildings, which is nonFlatTopographyShift lower.
            logger.debug(f"newPad.LengthFwd = {buildingTopAltitude}")
            newPad.LengthFwd = buildingTopAltitude
            logger.debug("newPad.Solid = True")
            newPad.Solid = True
            logger.debug("newPad.Symmetric = False")
            newPad.Symmetric = False
            FreeCADDOC.recompute()

        logger.info(f"Writing the STL {outputFileName}")
        Mesh.export(FreeCADDOC.Objects, outputFileName)


    def getBuildingsFromRectangle(self, minx, miny, maxx, maxy, dataSourceName=None, inputCRS=WGS84, withElevation=False):
        """
        Extract buildings within a rectangular bounding box.

        This method queries the building datasource and returns all buildings whose
        geometries intersect with or are contained within the specified rectangular
        region. Optionally adds elevation data from raster topography.

        Parameters
        ----------
        minx : float
            Minimum x-coordinate (longitude) of the bounding box.

        miny : float
            Minimum y-coordinate (latitude) of the bounding box.

        maxx : float
            Maximum x-coordinate (longitude) of the bounding box.

        maxy : float
            Maximum y-coordinate (latitude) of the bounding box.

        dataSourceName : str, optional
            Name of the building datasource to query. If None, uses the default
            datasource configured in the project (from 'defaultBuildingDataSource'
            config key). Default is None.

        inputCRS : int, optional
            EPSG code of the input coordinate system for minx, miny, maxx, maxy.
            Default is WGS84 (4326).

        withElevation : bool, optional
            If True, automatically adds elevation data to each building's centroid
            using the raster topography toolkit. The elevation is added as an
            'evaluation' column. Default is False.

        Returns
        -------
        geopandas.GeoDataFrame
            GeoDataFrame containing buildings within the specified region. Includes
            all original columns from the datasource plus optionally an 'evaluation'
            column if withElevation=True.

        Examples
        --------
        >>> # Get buildings in a region (WGS84 coordinates)
        >>> buildings = buildings_tk.getBuildingsFromRectangle(
        ...     minx=34.75, miny=31.75,
        ...     maxx=34.85, maxy=31.85,
        ...     dataSourceName="tel_aviv_buildings"
        ... )
        >>> 
        >>> # Get buildings with elevation data
        >>> buildings_with_elevation = buildings_tk.getBuildingsFromRectangle(
        ...     minx=34.75, miny=31.75,
        ...     maxx=34.85, maxy=31.85,
        ...     withElevation=True
        ... )
        >>> 
        >>> # Get buildings in ITM coordinates
        >>> buildings = buildings_tk.getBuildingsFromRectangle(
        ...     minx=180000, miny=660000,
        ...     maxx=190000, maxy=670000,
        ...     inputCRS=2039  # ITM
        ... )
        """
        if dataSourceName is None:
            dataSourceName = self.getConfig()["defaultBuildingDataSource"]

        doc = self.getDataSourceDocument(dataSourceName)
        buildings = self.cutRegionFromSource(doc, shape=[minx, miny, maxx, maxy], inputCRS=inputCRS)

        if withElevation:
            buildings = self.getBuildingHeightFromRasterTopographyToolkit(buildings)

        return buildings
    @staticmethod
    def get_buildings_height(gdf):
        """
        Extract building names, geometries, and height information from a GeoDataFrame.

        This utility method processes a GeoDataFrame to extract building information,
        including name, geometry, and height. If height is not directly available,
        it estimates height from the number of building levels (levels * 3 meters).

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            Input GeoDataFrame containing building geometries and properties. Expected
            columns include:
            - 'name': Building name or identifier
            - 'height': Building height in meters (optional)
            - 'building:levels': Number of building levels/floors (optional)
            - 'geometry': Polygon geometry of building footprint

        Returns
        -------
        geopandas.GeoDataFrame
            GeoDataFrame with columns:
            - 'name': Building name
            - 'geometry': Building polygon geometry
            - 'height': Building height in meters (from 'height' column or
              estimated from 'building:levels' * 3)

        Notes
        -----
        - Height estimation assumes 3 meters per floor/level
        - If neither 'height' nor 'building:levels' is available, height will be None
        - Original geometry is preserved

        Examples
        --------
        >>> import geopandas as gpd
        >>> from shapely.geometry import Polygon
        >>> 
        >>> # Create sample building data
        >>> buildings = gpd.GeoDataFrame({
        ...     'name': ['Building A', 'Building B'],
        ...     'building:levels': [5, None],
        ...     'height': [None, 20.0],
        ...     'geometry': [Polygon([(0,0), (10,0), (10,10), (0,10)]),
        ...                  Polygon([(20,20), (30,20), (30,30), (20,30)])]
        ... })
        >>> 
        >>> # Extract height information
        >>> result = BuildingsToolkit.get_buildings_height(buildings)
        >>> print(result[['name', 'height']])
        """
        building_info = []

        for index, row in gdf.iterrows():
            name = row.get('name', 'Unnamed')  # Extract the building name

            # Retain the original geometry
            geometry = row.geometry

            # Try to get the height or number of levels
            height = row.get('height')
            levels = row.get('building:levels')

            # Determine the height value
            if height is None and levels is not None:
                height = levels * 3  # Estimate height based on number of levels

            # Append the building information
            building_info.append({
                'name': name,
                'geometry': geometry,  # Keep the original geometry
                'height': height
            })

        return gpd.GeoDataFrame(building_info)
    @staticmethod
    def filter_buildings_in_area(buildings_data, min_longitude, min_latitude, max_longitude, max_latitude):
        """
        Filter building features within a geographic bounding box.

        This static method filters building features from a GeoJSON-like dictionary
        based on spatial intersection with a rectangular bounding box defined by
        latitude/longitude coordinates.

        Parameters
        ----------
        buildings_data : dict
            GeoJSON-like dictionary containing building features. Must have a
            'features' key containing a list of feature dictionaries, each with:
            - 'geometry': GeoJSON geometry object (Point, Polygon, MultiPolygon, etc.)
            - 'properties': Dictionary of feature attributes

        min_longitude : float
            Minimum longitude (x-coordinate) of the bounding box in WGS84.

        min_latitude : float
            Minimum latitude (y-coordinate) of the bounding box in WGS84.

        max_longitude : float
            Maximum longitude (x-coordinate) of the bounding box in WGS84.

        max_latitude : float
            Maximum latitude (y-coordinate) of the bounding box in WGS84.

        Returns
        -------
        geopandas.GeoDataFrame
            GeoDataFrame containing only buildings whose geometries intersect with
            or are contained within the specified bounding box. Includes all
            properties from the original features plus a 'geometry' column with
            Shapely geometry objects.

        Examples
        --------
        >>> import json
        >>> 
        >>> # Load GeoJSON data
        >>> with open('buildings.geojson') as f:
        ...     buildings_geojson = json.load(f)
        >>> 
        >>> # Filter buildings in a region (Tel Aviv area)
        >>> filtered = BuildingsToolkit.filter_buildings_in_area(
        ...     buildings_data=buildings_geojson,
        ...     min_longitude=34.75,
        ...     min_latitude=31.75,
        ...     max_longitude=34.85,
        ...     max_latitude=31.85
        ... )
        >>> 
        >>> print(f"Found {len(filtered)} buildings in the area")
        """
        # Create a polygon from the bounding box
        bbox_polygon = Polygon([(min_longitude, min_latitude),
                                (max_longitude, min_latitude),
                                (max_longitude, max_latitude),
                                (min_longitude, max_latitude),
                                (min_longitude, min_latitude)])  # Close the polygon

        # List to hold filtered building data
        filtered_buildings = []

        for feature in buildings_data['features']:
            geometry = feature['geometry']
            properties = feature['properties']

            # Create a geometry shape from the geometry data
            geom = shape(geometry)  # Handles Point, Polygon, MultiPolygon, etc.

            # Check if the geometry intersects with the bounding box polygon
            if geom.intersects(bbox_polygon):
                # Add all properties and the geometry
                properties['geometry'] = geom
                filtered_buildings.append(properties)

        # Create a GeoDataFrame
        gdf = gpd.GeoDataFrame(filtered_buildings)

        return gdf