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
    Toolkit to manage the buildings. Reading the shapefile with geoDataFrame will result in dataframe
    with the following columns:
        - Geometry: the polygon of the building.
        - Building height column : the column name is in BuildingHeightColumn. Default value=BLDG_HT.
        - Land height  : the columns name is in LandHeightColumn. Default value=HT_LAND.
    """
    def __init__(self, projectName, filesDirectory=None,connectionName=None):

        super().__init__(projectName=projectName, toolkitName="Buildings", filesDirectory=filesDirectory,connectionName=connectionName)
        self._analysis = analysis(dataLayer=self)

    def getBuildingHeightFromRasterTopographyToolkit(self, buildingData, topographyDataSource=None):
        """
        Get the topography height of each building (at its center) in the building data using the topography toolkit. Return data frame wtih 'evaluation' as a column.

        Parameters
        ----------
        buildingData : geopandas.geoDataFrame.
            The building.

        topographyDataSource : string,default=None.
            The name of the datasource in the topography toolkit. If None, use the default datasource there.

        Returns
        -------
            geopandas.DataFrame
        """
        topotk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_RASTER_TOPOGRAPHY, projectName=self.projectName)
        elevations = topotk.getPointListElevation(buildingData.centroid.to_crs(WGS84))
        return buildingData.join(elevations)

    def buildingsGeopandasToSTLRasterTopography(self,
                                                buildingData,
                                                buildingHeightColumn,
                                                buildingElevationColumn,
                                                outputFileName,
                                                flatTerrain = False,
                                                referenceTopography = 0,
                                                nonFlatTopographyShift=10):
        """
        Converts a building data (in geopandas format) to STL using the FreeCAD module.
        Using the raster topography to estimate the height of each building.
        This is a low level procedure. It can be used, but the easier way to use the toolkit is to generate the buildings from an area using the regionToSTL procedure.
        We must save the file to the disk, as it is the current implementation of FreeCAD.

        Parameters
        ----------
        buildingData : geopandas.DataFrame
            The buildings data.

        buildingHeightColumn : string
            The name of the column that holds the height of the buildings in [m].

        buildingElevationColumn: string
            The name of the column that holds the elevation of the building.

        outputFileName : string
            The absolute path of the output STL.

        flatTerrain : bool
            If true, use a flat terrain.

        nonFlatTopographyShift : float
            Shift the house with respect to its height in the topography.

        referenceTopography : float [default 0]
            If flatTerrain, use this as the reference height for the buildings.

        Returns
        -------

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
        Return the buildings geopandas for the rectangle region.

        Parameters
        ----------
        minx: float
            Minimum value of x-axis.

        miny: float
            Minimum value of y-axis.

        maxx: float
            Maximum value of x-axis.

        may: float
            Maximum value of y-axis.

        dataSourceName: str,default=None.
            The datasource name. If none, will use the default datasource.

        withElevation : bool,default=False.
            If True, use the topography (raster) toolkit to get the heights.

        Returns
        -------
            geopandas.DataFrame
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
        Extract building names, geometries(coordination), and height information from a GeoDataFrame.
        if there is no height available - it will calculate the number of levels in the building * 3/
        else: none

        Parameters:
        gdf (GeoDataFrame): A GeoPandas DataFrame containing building geometries
                            and associated properties.

        Returns:
        geopandas.DataFrame

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
        Filter building features by a specified geographic area using a bounding box.

        Parameters:
        ----------
        buildings_data : dict
            A GeoJSON-like dictionary containing building features to filter.

        min_longitude : float
            Minimum longitude defining the bounding box.

        min_latitude : float
            Minimum latitude defining the bounding box.

        max_longitude : float
            Maximum longitude defining the bounding box.

        max_latitude : float
            Maximum latitude defining the bounding box.

        Returns:
        -------
        gpd.GeoDataFrame
            A GeoDataFrame containing buildings that are located within the specified area.
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