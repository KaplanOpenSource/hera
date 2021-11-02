import os
import logging
from hera.measurements.GIS.vector import toolkit
from hera.measurements.GIS.vector.buildings.analysis import analysis
from hera.toolkit import TOOLKIT_SAVEMODE_NOSAVE,TOOLKIT_SAVEMODE_ONLYFILE,TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,TOOLKIT_SAVEMODE_FILEANDDB,TOOLKIT_SAVEMODE_FILEANDDB_REPLACE

from hera.datalayer import datatypes



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
        self.logger.info("geoPandasToSTL - Start")
        try:
            from freecad import app as FreeCAD
            import Part
            import Mesh
        except ImportError as e:
            logging.error("Loading the Building Toolkit. FreeCAD not Found, cannot convert to STL")
            raise ImportError("FreeCAD module is not installed in this environment. Cannot convert to STL")

        maxheight = -500
        FreeCADDOC = FreeCAD.newDocument("Unnamed")

        shp = buildingData # olv version compatability

        k = -1
        for j in range(len(shp)):  # converting al the buildings
            try:
                walls = shp['geometry'][j].exterior.xy
            except:
                continue

            if j % 100 == 0:
                self.logger.execution(f"{j} shape file is executed. Length Shape: {len(shp)}")

            k = k + 1
            wallsheight = shp[self.BuildingHeightColumn][j]
            if flat is None:
                altitude = shp[self.LandHeightColumn][j]
            else:
                altitude = flat

            FreeCADDOC.addObject('Sketcher::SketchObject', 'Sketch' + str(j))
            FreeCADDOC.Objects[2 * k].Placement = FreeCAD.Placement(FreeCAD.Vector(0.000000, 0.000000, 0.000000),  # 2*k-1
                                                             FreeCAD.Rotation(0.000000, 0.000000, 0.000000, 1.000000))

            for i in range(len(walls[0]) - 1):
                FreeCADDOC.Objects[2 * k].addGeometry(Part.Line(FreeCAD.Vector(walls[0][i], walls[1][i], altitude),
                                                         FreeCAD.Vector(walls[0][i + 1], walls[1][i + 1], altitude)))

            FreeCADDOC.addObject("PartDesign::Pad", "Pad" + str(j))
            FreeCADDOC.Objects[2 * k + 1].Sketch = FreeCADDOC.Objects[2 * k]
            buildingTopAltitude = wallsheight + altitude  # wallsheight + altitude
            maxheight = max(maxheight, buildingTopAltitude)
            FreeCADDOC.Objects[2 * k + 1].Length = buildingTopAltitude  # 35.000000
            FreeCADDOC.Objects[2 * k + 1].Reversed = 0
            FreeCADDOC.Objects[2 * k + 1].Midplane = 0
            FreeCADDOC.Objects[2 * k + 1].Length2 = wallsheight  # 100.000000
            FreeCADDOC.Objects[2 * k + 1].Type = 0
            FreeCADDOC.Objects[2 * k + 1].UpToFace = None


        FreeCADDOC.recompute() # maybe it should be in the loop.

        outputfileFull = os.path.abspath(os.path.join(self.FilesDirectory,outputFileName))
        Mesh.export(FreeCADDOC.Objects, outputfileFull)

        self.logger.info(f"geoPandasToSTL - End. Max height found {maxheight}")
        return maxheight


    def regionToSTL(self, regionNameOrData, outputFileName,dataSourceName,flat=None,saveMode=TOOLKIT_SAVEMODE_FILEANDDB_REPLACE):
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

            Returns
            -------
                The maximal height

        """

        if isinstance(regionNameOrData, str):
            shp = self.cutRegionFromSource(regionNameOrData,dataSourceName=dataSourceName)
        else:
            shp = regionNameOrData

        maxheight = self.geoPandasToSTL(buildingData=shp, outputFileName=outputFileName, flat=flat)

        if saveMode in [TOOLKIT_SAVEMODE_FILEANDDB_REPLACE,
                        TOOLKIT_SAVEMODE_FILEANDDB]:

            regionNameSTL = outputFileName.split(".")[0] if "." in outputFileName else outputFileName

            doc = self.getSTL(regionNameSTL)

            if doc is not None and saveMode == TOOLKIT_SAVEMODE_FILEANDDB:
                raise ValueError(f"STL {regionNameSTL} exists in project {self.projectName}")

            desc = {
                    toolkit.TOOLKIT_LOCATION_REGIONNAME: regionNameSTL,
                    toolkit.TOOLKIT_TOOLKITNAME_FIELD : self.toolkitName
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

if __name__ == "__main__":

    bt = BuildingsToolkit('test4')
    bt.addRegion([177933,663923,187933,673923], 'new9',crs = 2039)
    # lis = vt.getRegionNameList()
    #reg = bt.cutRegionFromSource('new9',dataSourceName='BNTL',isBounds = True, crs = 2039)
    #data = gps.GeoDataFrame.from_file('/mnt/public/omri_hadas/Production_Mode/Dispersion_Model/Haifa09_aerosols/LSM_for_SOURCE_ESTIMATION_epsilon_version/Lambda_Inputs/Haifa_Krayot_202323_741796/290_rez_200_afterBLD_correction/BLD_krayot_after_correction.shp')
    #lm = bt._analysis.LambdaOfDomain(270,200,buildingsDataSourceNameOrData=data,crs = 2039)
    lm = bt._analysis.LambdaOfDomain(270, 200, exteriorBlockNameOrData='new9', crs=2039)
    p=1

