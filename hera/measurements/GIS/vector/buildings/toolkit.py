import numpy
import os
import geopandas
from ..toolkit import VectorToolkit
from .analysis import analysis

try:
#    logger.execution("Trying to Load the FreeCAD module")
    import FreeCAD
    import Part
    import Mesh
except ImportError as e:
#    logger.error(f"Loading the Building Toolkit. FreeCAD not Found, cannot convert to STL: {e}")
    raise ImportError("FreeCAD module is not installed in this environment. Cannot convert to STL")

import matplotlib.pyplot as plt
from .....utils.logging import get_classMethod_logger
from ...utils import WSG84, ITM, ED50_ZONE36N
from ..... import toolkitHome


class BuildingsToolkit(VectorToolkit):
    """
        Toolkit to manage the buildings.

        Reading the shapefile with geoDataFrame will result in dataframe
        with the following columns:

        - geometry: the polygon of the building
        - Building height column : the column name is in BuildingHeightColumn
                                    default value: BLDG_HT
        - Land height  : the columns name is in LandHeightColumn
                                default value: HT_LAND
    """

    def __init__(self, projectName, filesDirectory=None):

        super().__init__(projectName=projectName, toolkitName="Buildings", filesDirectory=filesDirectory)
        self._analysis = analysis(dataLayer=self)

    def getBuildingHeightFromRasterTopographyToolkit(self, buildingData, topographyDataSource=None):
        """
            Get the topography height of each building (at its center) in the building data using the topography toolkit.

        Parameters
        ----------
        buildingData : geopandas.geoDataFrame.
                A building

        topographyDataSource : string  default None
            The name of the datasource in the topography toolkit. If None, use the default datasource there.

        Returns
        -------
            geopandas.DataFrame with 'elevation' as a column.
        """
        topotk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_RASTER_TOPOGRAPHY, projectName=self.projectName)
        elevations = topotk.getPointListElevation(buildingData.centroid.to_crs(WSG84))
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

                This is a low level procedure. It can be used, but the easier way to use the toolkit is
                to generate the buildings from an area using the regionToSTL procedure.

                We must save the file to the disk, as it is the current implementation of FreeCAD.

        Parameters
        ----------
        buildingData : geopandas.DataFrame
        buildingHeightColumn : string
            The name of the column that holds the height of the buildings in [m].
        buildingElevationColumn: string
            The name of the column that holds the elevation of the building.

        outputFileName : string
            The absolute path of the output STL

        flatTerrain : bool
            If true, use a flat terrain.
        nonFlatTopographyShift : float
            shift the house with respect to its height in the topography.

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
                logger.execution(f"{indx}/{len(buildingData)} shape file is executed")

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

        logger.execution(f"Writing the STL {outputFileName}")
        Mesh.export(FreeCADDOC.Objects, outputFileName)


    def getBuildingsFromRectangle(self, minx, miny, maxx, maxy, dataSourceName=None, inputCRS=WSG84,withElevation=False):
        """
            Return the buildings geopandas for the region.

        Parameters
        ----------
        minx : float
        miny
        maxx
        may
        dataSourceName
        withElevation : bool
            If True, use the topography (raster) toolkit to get the heghts.
            for now, uses the default datasource in the raster toolkit.

        Returns
        -------

        """
        if dataSourceName is None:
            dataSourceName = self.getConfig()["defaultBuildingDataSource"]

        doc = self.getDataSourceDocument(dataSourceName)
        buildings = self.cutRegionFromSource(doc, shape=[minx, miny, maxx, maxy], inputCRS=inputCRS)

        if withElevation:
            buildings = self.getBuildingHeightFromRasterTopographyToolkit(buildings)

        return buildings

    def regionToSTL(self, regionNameOrData, outputFileName, dataSourceName, flat=None, saveMode=None, crs=None):
        """
            Converts the document to the stl and saves it to the disk.
            Adds the stl file to the DB.


            Parameters
            ----------

            regionNameOrData: str or geopandas .
                The name of the datasource or the geopandas file to convert.

                If geopandas has the following columns:

            dataSourceName: string
                The name of the datasource that contains the database of the buildings

            flat: None or float.
                The base of the building.
                If None, use the buildings height.

            outputfile: str
                a path to the output file.

            saveMode: str
                - None, does not add to the DB.
                - TOOLKIT_SAVEMODE_FILEANDDB_REPLACE replace the document if exists
                - TOOLKIT_SAVEMODE_FILEANDDB         throws exception.

            crs : int [optional]
                Used if the CRS of the shapeData is different than the CRS of the datasource.


            Returns
            -------
                The maximal height

        """
        self.logger.info("---- Start ----")

        if not isinstance(regionNameOrData, geopandas.GeoDataFrame):
            shp = self.cutRegionFromSource(regionNameOrData, datasourceName=dataSourceName, isBounds=True, crs=crs)
        else:
            shp = regionNameOrData

        self.logger.info(f"Region {regionNameOrData} consists of {len(shp)} buildings")

        maxheight = self.buildingsGeopandasToSTL(buildingData=shp, outputFileName=outputFileName, flat=flat)

        if saveMode in [TOOLKIT_SAVEMODE_FILEANDDB_REPLACE,
                        TOOLKIT_SAVEMODE_FILEANDDB]:

            regionNameSTL = outputFileName.split(".")[0] if "." in outputFileName else outputFileName

            doc = self.getSTL(regionNameSTL)

            if doc is not None and saveMode == TOOLKIT_SAVEMODE_FILEANDDB:
                raise ValueError(f"STL {regionNameSTL} exists in project {self.projectName}")

            desc = {
                toolkit.TOOLKIT_VECTOR_REGIONNAME: regionNameSTL,
                toolkit.toolkitExtension.TOOLKIT_TOOLKITNAME_FIELD: self.toolkitName
            }

            if doc is None:
                self.addCacheDocument(type=self.doctype,
                                      resource=os.path.abspath(outputFileName),
                                      dataFormat=self.datatypes.STRING,
                                      desc=desc)
            else:
                doc.resource = os.path.abspath(outputFileName)
                doc.desc = desc
                doc.save()

        return maxheight

    def buildingsToSTL(self, point1, point2, domainname, outputfile):
        """
            write stl file of the buildings in a domain

        Parameters
        ----------
        point1 lower left corner (ITM)
        point2 upper right corner
        domainname - domain name
        outputfile - where to write the file

        Returns
        -------

        How to use
        ----------
        from hera.measurements.GIS.vector.buildings.toolkit import BuildingsToolkit
        BuildingsToolkit.buildingstorToSTL([200000, 740000], [201000, 741000], 'test', 'test2.stl')

        What to fix
        -----------

        """
        # print("poinat1")
        bounding = [point1[0], point1[1], point2[0], point2[1]]
        # bounding = point1

        self.addRegion(bounding, domainname, crs=2039)
        reg = self.cutRegionFromSource(domainname, datasourceName='BNTL', isBounds=True, crs=2039)

        bt.regionToSTL(domainname, outputfile, 'BNTL')
        return

    if __name__ == "__main__":
        fig, ax = plt.subplots()
        bt = BuildingsToolkit('kaplan', filesDirectory='/home/hadas/Hadas_S/Kaplan/')
        reg = [177933, 663923, 178933, 664423]
        # reg = bt.cutRegionFromSource(reg, 'BNTL', True)
        # reg.to_file('/home/hadas/Hadas_S/Kaplan/BNTL.shp')
        reg.plot(ax=ax)
        bt.addRegion(reg, 'buildings', crs=2039)
        # reg =[177933,663923,178933,664923]
        # lis = vt.getRegionNameList()
        reg = bt.cutRegionFromSource('new9', dataSourceName='BNTL', isBounds=True, crs=2039)
        data = gps.GeoDataFrame.from_file(
            '/mnt/public/omri_hadas/Production_Mode/Dispersion_Model/Haifa09_aerosols/LSM_for_SOURCE_ESTIMATION_epsilon_version/Lambda_Inputs/Haifa_Krayot_202323_741796/290_rez_200_afterBLD_correction/BLD_krayot_after_correction.shp')
        lm = bt.analysis.LambdaOfDomain(270, 200, buildingsDataSourceNameOrData=data, crs=2039)
        lm = bt.LambdaFromBuildingData(270, 200, data)
        # lm = bt._analysis.LambdaFromDatasource(270, 200, 'test', exteriorBlockNameOrData=reg, crs=2039)
        # p=1

    # newSketch = FreeCADDOC.addObject('Sketcher::SketchObject', 'Sketch3180')
    # newSketch.Placement = FreeCAD.Placement(FreeCAD.Vector(0.000000, 0.000000, 19.86), FreeCAD.Rotation(0.000000, 0.000000, 0.000000, 1.000000))
    # newSketch.addGeometry(Part.LineSegment(FreeCAD.Vector(179062.00002653955, 662887.8099938995, 19.86),FreeCAD.Vector(179028.70002653956, 662904.1199938995, 19.86)))
    # newSketch.addGeometry(Part.LineSegment(FreeCAD.Vector(179028.70002653956, 662904.1199938995, 19.86),FreeCAD.Vector(179038.50002653955, 662923.9999938995, 19.86)))
    # newSketch.addGeometry(Part.LineSegment(FreeCAD.Vector(179038.50002653955, 662923.9999938995, 19.86),FreeCAD.Vector(179071.80002653957, 662907.6199938995, 19.86)))
    # newSketch.addGeometry(Part.LineSegment(FreeCAD.Vector(179071.80002653957, 662907.6199938995, 19.86),FreeCAD.Vector(179062.00002653955, 662887.8099938995, 19.86)))
    # FreeCADDOC.addObject('Part::Extrusion', 'building3180')
    # newPad.Base = newSketch
    # newPad.LengthFwd = 21.62
