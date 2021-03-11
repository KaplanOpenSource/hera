import geopandas
from .locations.shapes import datalayer as shapeDatalayer
from ...toolkit import toolkit
from .locations.buildings import datalayer as buildingsDatalayer
import shapely

SAVEMODE_NOSAVE = None
SAVEMODE_ONLYFILE = "File"
SAVEMODE_FILEANDDV = "DB"


class DemographyToolkit(toolkit):
    """
        A toolkit to manage demography data


    """
    _Data = None

    _populationTypes = None # A dictionary with the names of the population types.

    _shapes = None
    _buildings = None

    @property
    def buildings(self):
        return self._buildings

    @property
    def shapes(self):
        return self._shapes

    def __init__(self, projectName,SourceOrData=None):
        """
            Initializes the demography tool.

        :param projectName: str
            The project name

        :param SourceOrData: str or geopandas or None. 
            None:   Does not have an inn
            str:    The name of the source in the DB that will be used as a data.
            pandas: Use this as the demography data. See documentation of the structure of the demographic dataframe.
        """
        self._projectName = projectName
        super().__init__(projectName=projectName,toolkitName="Demography")


        if isinstance(SourceOrData,str):
            if SourceOrData is None:
                raise ValueError(f"Source must be one of: {','.join(self.listSources())}")

            sourceDoc = self.getSource(SourceOrData)
            if sourceDoc is None:
                raise ValueError(f"Please load the f{SourceOrData} to the database for the project {projectName}. Use the hera-data package")

            self._Data = sourceDoc.getData()
        elif isinstance(SourceOrData,geopandas.GeoDataFrame):
            self._Data = SourceOrData
        else:
            raise ValueError("SourceOrData must be the name of the source (str) or the geopandas dataframe of the demography")


        self._populationTypes = {"All":"total_pop","Children":"age_0_14","Youth":"age_15_19",
                           "YoungAdults":"age_20_29","Adults":"age_30_64","Elderly":"age_65_up"}

        self._shapes    = shapeDatalayer(projectName=projectName)
        self._buildings = buildingsDatalayer(projectName=projectName)

    def projectPolygonOnPopulation(self, Shape, projectName=None, populationTypes="All", Data=None):
        import warnings
        warnings.warn("projectPolygonOnPopulation was changed to analysis.calculatePopulationInPolygon",DeprecationWarning)
        self.analysis.calculatePopulationInPolygon(Shape=Shape, projectName=projectName, populationTypes=populationTypes, Data=Data)

    @property
    def data(self):
        return self._Data

    @data.setter
    def data(self, value):
        if not isinstance(value,geopandas.GeoDataFrame):
            raise ValueError("Must be a geodataframe.")
        self._Data = value
        return self

class analysis:

    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, dataLayer):
        self._datalayer = dataLayer



    def populateNewArea(self, Shape, populationTypes=None, convex=True, saveMode=SAVEMODE_NOSAVE, path=None, name=None, **kwargs):
        """
        Make a geoDataFrame with a selected polygon as the geometry,
        and the sum of the population in the polygons that intersect it as its population.

        :param Shape:
        :param populationTypes:
        :param convex:
        :param saveToFile:
        :param addToDB:
        :param path:
        :param name:
        :param Data:
        :param kwargs:
        :return:
        """
        Data = self.datalayer._Data if Data is None else Data
        if isinstance(Shape,str):
            poly = self.datalayer.shapes.getShape(Shape)
            if convex:
                polys = self.datalayer.buildings.analysis.ConvexPolygons(poly)
                poly = polys.loc[polys.area==polys.area.max()].geometry[0]
            else:
                poly = documents[0].getData().unary_union
        else:
            poly = Shape



        res_intersect_poly = Data.loc[Data["geometry"].intersection(poly).is_empty == False]
        populationTypes = self.datalayer.getConfig()["populationTypes"].values() if populationTypes is None else populationTypes

        newData = geopandas.GeoDataFrame.from_dict([{"geometry": poly}])
        for populationType in populationTypes:
            newData[populationType] = res_intersect_poly.sum()[populationType]

        if saveMode!=SAVEMODE_NOSAVE:
            if path is None:
                raise KeyError("Select a path for the new file")
            newData.to_file(path)
            if addToDB:
                if name is None:
                    if type(Shape) == str:
                        name = Shape
                    else:
                        raise KeyError("Select a name for the new area")
                self.datalayer.addMeasurementsDocument(desc=(dict(name=name, **kwargs)),
                                                   resource=path, type="Demography", dataFormat="geopandas")
        return newData

    def calculatePopulationInPolygon(self, Shape, projectName=None, populationTypes="All", Data=None):
        """
            Finds the population in a polygon.

        Params:
            Shape: shapely.Polygon or str.
                    The polygon to calculate the poulation in.
                    Can be either the polygon itself (shapely.Polygon) or a name of a saved geometry in the database.
                    The saved geometries in the DB are obtained using the shapeToolkit.

            populationTypes: str or list of str
                    Additional population columns that will be calcualted.
        """

        Data = self._Data if Data is None else Data
        if type(Shape) == str:
            sDatalayer = shapeDatalayer(projectName=projectName)
            poly = sDatalayer.getShape(Shape)
            if poly is None:
                documents = sDatalayer.getMeasurementsDocuments(CutName=Shape)
                if len(documents) == 0:
                    raise KeyError("Shape %s was not found" % Shape)
                else:
                    points = documents[0].asDict()["desc"]["points"]
                    poly = shapely.geometry.Polygon([[points[0], points[1]],
                                                     [points[0], points[3]],
                                                     [points[2], points[3]],
                                                     [points[2], points[1]]])
        else:
            poly = Shape
        if type(populationTypes) == str:
            populationTypes = [populationTypes]
        res_intersect_poly = Data.loc[Data["geometry"].intersection(poly).is_empty == False]
        intersection_poly = res_intersect_poly["geometry"].intersection(poly)
        res_intersection = geopandas.GeoDataFrame.from_dict(
            {"geometry": intersection_poly.geometry,
             "areaFraction": intersection_poly.area / res_intersect_poly.area})
        for populationType in populationTypes:
            if "populationTypes" in self.getConfig():
                if populationType in self.getConfig()["populationTypes"]:
                    populationType = self.getConfig()["populationTypes"][populationType]
            res_intersection[populationType] = intersection_poly.area / res_intersect_poly.area * res_intersect_poly[
                populationType]

        return res_intersection

    # def setConfig(self, Source="Lamas", units="WGS84", populationTypes = None, **kwargs):
    #     """
    #     Create a config documnet or updates an existing config document.
    #     """
    #     populationTypes = {"All":"total_pop","Children":"age_0_14","Youth":"age_15_19",
    #                        "YoungAdults":"age_20_29","Adults":"age_30_64","Elderly":"age_65_up"} if populationTypes is None else populationTypes
    #     config = dict(source=Source,units=units,populationTypes=populationTypes, **kwargs)
    #     super().setConfig(config=config)
    #
    #     datalist = self.getMeasurementsDocuments(source=config["source"])
    #
    #     if len(datalist) > 0:
    #         self._Data = datalist[0].getData()
    #     else:
    #         self._Data = None
