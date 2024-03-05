import numpy
import os
import logging
from hera import toolkitHome
from hera.measurements.GIS.vector import toolkit
from hera.measurements.GIS.vector.buildings.analysis import analysis
from hera.toolkit import TOOLKIT_SAVEMODE_NOSAVE,TOOLKIT_SAVEMODE_ONLYFILE,TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,TOOLKIT_SAVEMODE_FILEANDDB,TOOLKIT_SAVEMODE_FILEANDDB_REPLACE
from hera.datalayer import datatypes
import geopandas
import matplotlib.pyplot as plt

#FREECADPATH = '/usr/lib/freecad-python3/lib/'  # Or add to PythonPath
#import sys
#sys.path.append(FREECADPATH)

import FreeCAD
import Part
import Mesh


class BuildingsToolkit(toolkit.VectorToolkit):
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

    _BuildingHeightColumn = None

    @property
    def BuildingHeightColumn(self):
        return self._BuildingHeightColumn

    @BuildingHeightColumn.setter
    def BuildingHeightColumn(self, value):
        self._BuildingHeightColumn = value

    @property
    def LandHeightColumn(self):
        return self._LandHeightColumns

    @LandHeightColumn.setter
    def LandHeightColumn(self, value):
        self._LandHeightColumns = value


    def __init__(self, projectName, filesDirectory=None):

        super().__init__(projectName=projectName, toolkitName="Buildings", filesDirectory=filesDirectory)
        self._analysis = analysis(dataLayer=self)

        self._BuildingHeightColumn = "BLDG_HT"
        self._LandHeightColumns   = 'HT_LAND'
        self.name   = projectName 

    @property
    def doctype(self):
        return f"{self.name}_STL"


    def geoPandasToSTL(self,buildingData,outputFileName,flat=None):
        """
                Converts a building data (in geopandas format) to STL using the Freecad module.

        Parameters
        ----------
        buildingData : geopandas.GeoDataFrame
                The building data

        outputFileName : str
                Full path to the output file.


        flat: None or float.
                The base of the building.
                If None, use the buildings height.

        Returns
        -------
            float
                The maximal height that was found
        """
        self.logger.info("---- Start ---")
        try:
            self.logger.execution("Loading Freecad")
         #   from freecad import app as FreeCAD
         #   import Part
         #   import Mesh
        except ImportError as e:
            logging.error("Loading the Building Toolkit. FreeCAD not Found, cannot convert to STL")
            raise ImportError("FreeCAD module is not installed in this environment. Cannot convert to STL")

        maxheight = -500
        try:
                FreeCADDOC = FreeCAD.newDocument("Unnamed")
        except:
            err = "FreeCAD not found. Install before using this function and add to the PYTHONPATH"
            raise ValueError(err)

        shp = buildingData # olv version compatability

        for indx,building in shp.iterrows():  # converting al the buildings

            try:
                walls = building['geometry'].exterior.xy
                walls = numpy.stack(walls).T
            except:
                continue

            if indx % 100 == 0:
                self.logger.execution(f"{indx}/{len(shp)} shape file is executed")


            wallsheight = building[self.BuildingHeightColumn]
            if flat is None:
                altitude = building[self.LandHeightColumn]
            else:
                altitude = flat

            self.logger.debug(f" Building -- {indx} --")
            self.logger.debug("newSketch = FreeCADDOC.addObject('Sketcher::SketchObject', 'Sketch"+ str(indx) +"')")
            newSketch = FreeCADDOC.addObject('Sketcher::SketchObject', 'Sketch' + str(indx))

            self.logger.debug(f"newSketch.Placement = FreeCAD.Placement(FreeCAD.Vector(0.000000, 0.000000, {altitude}), FreeCAD.Rotation(0.000000, 0.000000, 0.000000, 1.000000))")
            newSketch.Placement = FreeCAD.Placement(FreeCAD.Vector(0.000000, 0.000000, altitude),  # 2*k-1
                                                             FreeCAD.Rotation(0.000000, 0.000000, 0.000000, 1.000000))

            for xy0,xy1 in zip(walls[:-1],walls[1:]):
                self.logger.debug(f"newSketch.addGeometry(Part.LineSegment(FreeCAD.Vector({xy0[0]}, {xy0[1]}, {altitude}),FreeCAD.Vector({xy1[0]}, {xy1[1]}, {altitude})))")
                newSketch.addGeometry(Part.LineSegment(FreeCAD.Vector(xy0[0], xy0[1], altitude),FreeCAD.Vector(xy1[0], xy1[1], altitude)))

            self.logger.debug("FreeCADDOC.addObject('Part::Extrusion', 'building" + str(indx) + "')")
            newPad = FreeCADDOC.addObject("Part::Extrusion", "building" + str(indx))
            self.logger.debug("newPad.Base = newSketch")
            newPad.Base = newSketch
            buildingTopAltitude = wallsheight # + altitude  # wallsheight + altitude
            maxheight = max(maxheight, buildingTopAltitude)
            self.logger.debug(f"newPad.LengthFwd = {buildingTopAltitude}")
            newPad.LengthFwd = wallsheight

            self.logger.debug("newPad.Solid = True")
            newPad.Solid = True

            self.logger.debug("newPad.Symmetric = False")
            newPad.Symmetric = False



        FreeCADDOC.recompute() # maybe it should be in the loop.

        outputfileFull = os.path.abspath(os.path.join(self.FilesDirectory,outputFileName))
        self.logger.execution(f"Writing the STL {outputfileFull}")
        Mesh.export(FreeCADDOC.Objects, outputfileFull)

        self.logger.info(f"geoPandasToSTL - End. Max height found {maxheight}")
        return maxheight


    def regionToSTL(self, regionNameOrData, outputFileName,dataSourceName,flat=None,saveMode=TOOLKIT_SAVEMODE_FILEANDDB_REPLACE, crs=None):
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

        if not isinstance(regionNameOrData,geopandas.GeoDataFrame):
            shp = self.cutRegionFromSource(regionNameOrData, datasourceName=dataSourceName, isBounds=True, crs=crs)
        else:
            shp = regionNameOrData

        self.logger.info(f"Region {regionNameOrData} consists of {len(shp)} buildings")

        maxheight = self.geoPandasToSTL(buildingData=shp, outputFileName=outputFileName, flat=flat)

        if saveMode in [TOOLKIT_SAVEMODE_FILEANDDB_REPLACE,
                        TOOLKIT_SAVEMODE_FILEANDDB]:

            regionNameSTL = outputFileName.split(".")[0] if "." in outputFileName else outputFileName

            doc = self.getSTL(regionNameSTL)

            if doc is not None and saveMode == TOOLKIT_SAVEMODE_FILEANDDB:
                raise ValueError(f"STL {regionNameSTL} exists in project {self.projectName}")

            desc = {
                    toolkit.TOOLKIT_VECTOR_REGIONNAME: regionNameSTL,
                    toolkit.toolkitExtension.TOOLKIT_TOOLKITNAME_FIELD : self.toolkitName
            }

            if doc is None:
                self.addCacheDocument(type=self.doctype,
                                      resource=os.path.abspath(outputFileName),
                                      dataFormat=datatypes.STRING,
                                      desc=desc)
            else:
                doc.resource = os.path.abspath(outputFileName)
                doc.desc = desc
                doc.save()

        return maxheight


    def getSTL(self,regionNameSTL):
        """
            Retrive the STL from the DB

        Parameters
        ----------
        regionNameSTL: str
            The name of the regiona name STL.

        Returns
        -------
            The document of the STL.
        """
        desc = {
                "vector":regionNameSTL,
                "type" : self.doctype
                }
        docList = self.getCacheDocuments(**desc)
        return None if len(docList)==0 else docList[0]

    def buildingsToSTL(self,point1, point2, domainname, outputfile):
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
        #print("poinat1")
        bounding = [point1[0], point1[1], point2[0], point2[1]]
        # bounding = point1

        self.addRegion(bounding, domainname, crs=2039)
        reg = self.cutRegionFromSource(domainname, datasourceName='BNTL', isBounds=True, crs=2039)

        bt.regionToSTL(domainname, outputfile, 'BNTL')
        return


if __name__ == "__main__":
    fig, ax = plt.subplots()
    bt = BuildingsToolkit('kaplan',filesDirectory ='/home/hadas/Hadas_S/Kaplan/')
    reg = [177933, 663923, 178933, 664423]
    # reg = bt.cutRegionFromSource(reg, 'BNTL', True)
    # reg.to_file('/home/hadas/Hadas_S/Kaplan/BNTL.shp')
    reg.plot(ax=ax)
    bt.addRegion(reg, 'buildings',crs = 2039)
    # reg =[177933,663923,178933,664923]
    # lis = vt.getRegionNameList()
    reg = bt.cutRegionFromSource('new9',dataSourceName='BNTL',isBounds = True, crs = 2039)
    data = gps.GeoDataFrame.from_file('/mnt/public/omri_hadas/Production_Mode/Dispersion_Model/Haifa09_aerosols/LSM_for_SOURCE_ESTIMATION_epsilon_version/Lambda_Inputs/Haifa_Krayot_202323_741796/290_rez_200_afterBLD_correction/BLD_krayot_after_correction.shp')
    lm = bt.analysis.LambdaOfDomain(270,200,buildingsDataSourceNameOrData=data,crs = 2039)
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
