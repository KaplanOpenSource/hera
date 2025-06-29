import geopandas
import io
import os
from hera import toolkit
from hera.measurements.GIS.vector import toolkit
from hera.datalayer import datatypes, nonDBMetadataFrame
from hera import toolkitHome
from hera.toolkit import TOOLKIT_SAVEMODE_NOSAVE,TOOLKIT_SAVEMODE_ONLYFILE,TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,TOOLKIT_SAVEMODE_FILEANDDB,TOOLKIT_SAVEMODE_FILEANDDB_REPLACE, get_classMethod_logger  # ✅ FIXED: added missing import

class DemographyToolkit(toolkit.VectorToolkit):
    """
        A toolkit to manage demography data

    """

    _populationTypes = None # A dictionary with the names of the population types.

    _buildings = None
    _analysis = None

    @property
    def buildings(self):
        return self._buildings

    @property
    def shapes(self):
        return self._shapes

    @property
    def analysis(self):
        return self._analysis

    @property
    def populationTypes(self):
        return self._populationTypes

    def __init__(self, projectName, filesDirectory=None):
        """
            Initializes the demography tool.

        Parameters
        ----------
        projectName: str
            The project name

        SourceOrData: str or geopandas or None.
            None:   Try to load the default data source. If does not exist, set data to None.
            str:    The name of the source in the DB that will be used as a data.
            geoDataFrame: Use this as the demography data. See documentation of the structure of the demographic dataframe.

        dataSourceVersion : tuple of integers
            If specified load this version of the data.

        """
        self._projectName = projectName
        super().__init__(projectName=projectName, toolkitName="Demography", filesDirectory=filesDirectory)
        self._analysis = analysis(self)

        self._populationTypes = {"All":"total_pop","Children":"age_0_14","Youth":"age_15_19",
                           "YoungAdults":"age_20_29","Adults":"age_30_64","Elderly":"age_65_up"}

        self._buildings = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_BUILDINGS,projectName=projectName)

    def setDefaultDirectory(self, fileDirectory, create=True):
        """
            Set the default directory for the project.

        Parameters
        ----------
        fileDirectory: str
                The path to save the regions in.
                The directory is created if create flag is true (and directory does not exist).

        create: bool
            If false and directory does not exist, raise a NotADirectoryError exception.

        Returns
        -------
            str, the path.
        """
        fllpath = os.path.abspath(fileDirectory)
        logger = get_classMethod_logger(self, "setDefaultDirectory")

        if not os.path.exists(fllpath):
            if create:
                logger.debug(f"Directory {fllpath} does not exist. create")
                # ✅ Changed: safer than using os.system("mkdir -p ...")
                os.makedirs(fllpath, exist_ok=True)
            else:
                raise NotADirectoryError(f"{fllpath} is not a directory, and create is False.")

        # ✅ Changed: removed call to self.setConfig(...) since it's not defined in base class
        # ✅ Replaced with direct assignment to self._FilesDirectory as done elsewhere in abstractToolkit
        self._FilesDirectory = fllpath

    def projectPolygonOnPopulation(self, shapelyPolygon, dataSourceOrData, populationTypes="All", dataSourceVersion=None):
        import warnings
        warnings.warn("Depracted in Version 2.0.0+. Use datalayer.calculatePopulationInPolygon",
                      category=DeprecationWarning,
                      stacklevel=2)
        self.analysis.calculatePopulationInPolygon(shapelyPolygon=shapelyPolygon, dataSourceOrData=dataSourceOrData, dataSourceVersion=dataSourceVersion, populationTypes=populationTypes)  # 🔧 FIXED: Shape → shapelyPolygon

    @property
    def FilesDirectory(self):
        return self._FilesDirectory


    def loadData(self, regionaName, dataOrFileName, version=(0,0,1), metadata=dict(), overwrite=False):
        """
            Loading an demographic data to the DB

            Currently only an existing shp files/geojson files and replaces

            TODO:
                - Should be extended in the future to get a string and save it.
                - Extend to different save modes (i.e just save a file or warn if already exists).

        Parameters
        ----------
        regionaName: str
                The name of the region
        dataOrFileName: str
                Currently only a file name
        version: list
                A list of number to state the version of the data.
        metadata:
                any additional
        :return:
        """
        # ✅ We no longer delete the document manually – handled in addDataSource if needed
        self.addDataSource(
            dataSourceName=regionaName,
            resource=dataOrFileName,
            dataFormat=datatypes.GEOPANDAS,
            overwrite=overwrite,  # ✅ מוסיפים את זה
            **metadata
        )


class analysis:
    """
        Analysis of the demography toolkit.
    """

    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, dataLayer):
        self._datalayer = dataLayer

    def createNewArea(self,
                      shapeNameOrData,
                      dataSourceOrData,
                      dataSourceVersion=None,
                      populationTypes=None,
                      convex=True,
                      saveMode=TOOLKIT_SAVEMODE_NOSAVE,
                      regionName=None, metadata=dict()):
        """
            Make a geoDataFrame with a polygon as the geometry,
            and the sum of the population in the polygons that intersect it as its population.

            If saveMode is set to save to file (with or without DB) the regionName
            is used as the file name.

        Parameters
        -----------
        shapeNameOrData: str, geopandas
            A shape name, geopandas dataframe, geoJSON str

        dataSourceOrData: str, geopandas
            A demography data source name, a geoJSON (shapes with the population) or geopandas.dataframe

        dataSourceVersion: 3-tuple of int
            A version of the demography data source

        convex: bool

        saveMode: str
                Can be either:

                    - TOOLKIT_SAVEMODE_NOSAVE   : Just load the data from file and return the datafile

                    - TOOLKIT_SAVEMODE_ONLYFILE : Loads the data from file and save to a file.
                                                  raise exception if file exists.

                    - TOOLKIT_SAVEMODE_ONLYFILE_REPLACE: Loads the data from file and save to a file.
                                                  Replace the file if it exists.

                    - TOOLKIT_SAVEMODE_FILEANDDB : Loads the data from file and save to a file and store to the DB as a source.
                                                    Raise exception if the entry exists.

                    - TOOLKIT_SAVEMODE_FILEANDDB_REPLACE: Loads the data from file and save to a file and store to the DB as a source.
                                                    Replace the entry in the DB if it exists.

        regionName: str
            optional. If saved to the DB, use this as a region.
        metadata: dict
            Metadata to be saved to the DB if needed.

        Returns
        -------
            Document with the new data
        """
        # 🛡️ Ensure regionName is provided if saving is requested
        if saveMode in [TOOLKIT_SAVEMODE_ONLYFILE,
                        TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,
                        TOOLKIT_SAVEMODE_FILEANDDB,
                        TOOLKIT_SAVEMODE_FILEANDDB_REPLACE] and regionName is None:
            raise ValueError("Must specify regionName if saveMode is set to save data")

        # 📦 Load population layer
        if isinstance(dataSourceOrData, str):
            Data = self.datalayer.getDataSourceData(dataSourceOrData, dataSourceVersion)
            if Data is None:
                Data = geopandas.read_file(io.StringIO(dataSourceOrData))
        else:
            Data = dataSourceOrData

        # 📦 Load the polygon
        if isinstance(shapeNameOrData, str):
            polydoc = self.datalayer.shapes.getShape(shapeNameOrData)
            if polydoc is None:
                polydoc = geopandas.read_file(io.StringIO(shapeNameOrData))
        elif isinstance(shapeNameOrData, geopandas.geodataframe.GeoDataFrame):
            polydoc = shapeNameOrData
        else:
            poly = shapeNameOrData  # already a shapely geometry

        # 🧠 Construct final geometry if needed
        if isinstance(shapeNameOrData, str) or isinstance(shapeNameOrData, geopandas.geodataframe.GeoDataFrame):
            if convex:
                polys = self.datalayer.buildings.analysis.ConvexPolygons(polydoc)
                poly = polys.loc[polys.area == polys.area.max()].geometry[0]
            else:
                poly = polydoc.unary_union

        # 📐 Intersect population data with polygon (including topology fix)
        try:
            res_intersect_poly = Data.loc[Data["geometry"].intersection(poly).is_empty == False]
        except Exception as e:
            from shapely.errors import TopologicalError
            if isinstance(e, TopologicalError):
                if not poly.is_valid:
                    poly = poly.buffer(0)
                Data["geometry"] = Data["geometry"].apply(lambda g: g if g.is_valid else g.buffer(0))
                res_intersect_poly = Data.loc[Data["geometry"].intersection(poly).is_empty == False]
            else:
                raise e

        # 🧾 Construct resulting GeoDataFrame
        newData = geopandas.GeoDataFrame.from_dict([{"geometry": poly}])
        newData.crs = Data.crs

        # 🧮 Sum population attributes
        populationTypes = list(
            self.datalayer.populationTypes.values()) if populationTypes is None else populationTypes
        for populationType in populationTypes:
            if populationType in res_intersect_poly:
                newData[populationType] = res_intersect_poly.sum()[populationType]

        doc = None  # will hold DB document if stored

        # 💾 Save to file or DB if required
        if saveMode != TOOLKIT_SAVEMODE_NOSAVE:
            filename = regionName if "." in regionName else f"{regionName}.shp"
            fullname = os.path.join(self.datalayer.filesDirectory, filename)
            newData.to_file(fullname)

            if saveMode == TOOLKIT_SAVEMODE_FILEANDDB:
                desc = {
                    toolkit.TOOLKIT_DATASOURCE_NAME: regionName,
                    toolkit.TOOLKIT_TOOLKITNAME_FIELD: self.datalayer.toolkitName
                }
                desc.update(**metadata)
                doc = self.datalayer.addCacheDocument(
                    desc=desc,
                    resource=fullname,
                    type=toolkit.TOOLKIT_DATASOURCE_TYPE,
                    dataFormat=datatypes.GEOPANDAS
                )

        # 📤 Return result: either wrapped metadata frame or DB document
        return nonDBMetadataFrame(newData) if doc is None else doc

    def calculatePopulationInPolygon(self,
                                     shapelyPolygon,
                                     dataSourceOrData,
                                     dataSourceVersion=None,
                                    populationTypes=None):
        """
            Finds the population in a polygon.

        Parameters:
        -----------

            shapeNameOrData: str, shapely.Polygon, geopandas
                    The polygon to calculate the poulation in.

                    Can be either:
                        - the polygon itself (shapely.Polygon)
                        - shape name in the DB (str)
                        - geoJSON (str)
                        - geopandas.

            dataSourceOrData: str or geopandas
                    The demographic data.

                    Can be either:
                        - demography data source name
                        - geoJSON (str)
                        - geopandas


            dataSourceVersion: 3-tuple of int
                    If dataSourceOrData is demography data source name
                    then dataSourceVersion is a possible version.

            populationTypes: str or list of str
                    Additional population columns that will be calculated.

        Returns
        -------
            geopandas
            The intersection of the demography and the polygon.

        """


        # if isinstance(shapeNameOrData,str):
        #     poly = self.datalayer.shapes.getShape(shapeNameOrData)
        #     if poly is None:
        #         poly = geopandas.read_file(io.StringIO(shapeNameOrData))
        # else:

        poly = shapelyPolygon

        if isinstance(dataSourceOrData,str):
            demography = self.datalayer.getDataSourceData(dataSourceOrData, dataSourceVersion)
            if demography is None:
                demography = geopandas.read_file(io.StringIO(dataSourceOrData))
        else:
            demography = dataSourceOrData

        populationTypes = self.datalayer.populationTypes if populationTypes is None else populationTypes

        if isinstance(populationTypes,str):
            populationTypes = [populationTypes]

        res_intersect_poly = demography.loc[demography["geometry"].intersection(poly).is_empty == False]
        intersection_poly = res_intersect_poly["geometry"].intersection(poly)

        res_intersection = geopandas.GeoDataFrame.from_dict(
            {"geometry": intersection_poly.geometry,
             "areaFraction": intersection_poly.area / res_intersect_poly.area})

        for populationType in populationTypes:
            populationType = self.datalayer.populationTypes.get(populationType,populationType)
            if populationType in res_intersect_poly:
                res_intersection[populationType] = intersection_poly.area / res_intersect_poly.area * res_intersect_poly[populationType]

        return res_intersection
