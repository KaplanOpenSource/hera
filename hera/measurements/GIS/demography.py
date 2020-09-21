import geopandas
from .locations.shapes import datalayer as shapeDatalayer
from ...datalayer import project
from .locations.topography import datalayer as topoDatalayer
import shapely

class datalayer(project.ProjectMultiDBPublic):

    _publicProjectName = None
    _projectName = None

    @property
    def publicProjectName(self):
        return self._publicProjectName

    def __init__(self, projectName,Source="Lamas",publicProjectName="Demography",databaseNameList=None, useAll=False,):
        self._projectName = projectName
        self._publicProjectName = publicProjectName
        super().__init__(projectName=projectName, publicProjectName=publicProjectName,
                         databaseNameList=databaseNameList, useAll=useAll)
        self.setConfig({"source":Source, "units":"WGS84",
                        "populationTypes":{"All":"total_pop","Children":"age_0_14","Youth":"age_15_19","YoungAdults":"age_20_29","Adults":"age_30_64","Elderly":"age_65_up"}})

    def projectPolygonOnPopulation(self, Shape, populationTypes="All"):
        """
        Finds the population in a polygon.
        Params:
            Shape: The polygon, either a shapely polygon or a name of a saved geometry in the database.
            populationTypes: A string or a list of strings with options from config/population or column names
                             from the data.
        """

        Data = self.getMeasurementsDocuments(source=self.getConfig()["source"])[0].getData()
        if type(Shape) == str:
            poly = shapeDatalayer(projectName=self._projectName, databaseNameList=self._databaseNameList,
                                   useAll=self._useAll).getShape(Shape)
            if poly is None:
                documents = topoDatalayer(projectName=self._projectName, databaseNameList=self._databaseNameList,
                                   useAll=self._useAll).getDocuments(CutName=Shape)
                if len(documents) == 0:
                    raise KeyError("Shape %s was not found" % Shape)
                else:
                    points = documents[0].asDict()["desc"]["points"]
                    poly = shapely.geometry.Polygon([[points[0],points[1]],
                                            [points[0],points[3]],
                                            [points[2],points[3]],
                                            [points[2],points[1]]])
        else:
            poly = Shape
        if type(populationTypes) == str:
            populationTypes = [populationTypes]

        res_intersect_poly = Data.loc[Data["geometry"].intersection(poly).is_empty == False]
        intersection_poly = res_intersect_poly["geometry"].intersection(poly)
        res_intersection = geopandas.GeoDataFrame.from_dict(
            {"geometry": intersection_poly.geometry,
             "areaFraction": intersection_poly.area/res_intersect_poly.area})
        for populationType in populationTypes:
            print(self._publicProjectName)
            if "populationTypes" in self.getConfig():
                if populationType in self.getConfig()["populationTypes"]:
                    populationType = self.getConfig()["populationTypes"][populationType]
            import pdb
            pdb.set_trace()
            res_intersection[populationType] = intersection_poly.area / res_intersect_poly.area * res_intersect_poly[populationType]

        return res_intersection

    # def populateNewArea(self, Geometry, data=None, populationTypes=None, convex=True,save=False, addToDB=False, path=None, name=None, **kwargs):
    #     """
    #     make a geodataframe with a selected polygon as the geometry, and the sum of the population in the polygons that intersect it as its population.
    #     """
    #
    #     Data = self._publicMeasure.getDocuments(projectName="PublicData", type="Population")[0].getDocFromDB() if data is None else data
    #
    #     if type(Geometry) == str:
    #         poly = GIS_datalayer(projectName=self._projectName, FilesDirectory="").getGeometry(name=Geometry)
    #         if poly is None:
    #             documents = GIS_datalayer(projectName=self._projectName, FilesDirectory="").getExistingDocuments(CutName=Geometry)
    #             if len(documents) == 0:
    #                 raise KeyError("Geometry %s was not found" % Geometry)
    #             else:
    #                 if convex:
    #                     polys = ConvexPolygons(documents[0].getDocFromDB())
    #                     poly = polys.loc[polys.area==polys.area.max()].geometry[0]
    #                 else:
    #                     poly = documents[0].getDocFromDB().unary_union
    #     else:
    #         poly = Geometry
    #     res_intersect_poly = Data.loc[Data["geometry"].intersection(poly).is_empty == False]
    #     populationTypes = ["total_pop","age_0_14","age_15_19","age_20_29","age_30_64","age_65_up"] if populationTypes is None else populationTypes
    #
    #     newData = geopandas.GeoDataFrame.from_dict([{"geometry": poly}])
    #     for populationType in populationTypes:
    #         newData[populationType] = res_intersect_poly.sum()[populationType]
    #
    #     if save:
    #         if path is None:
    #             raise KeyError("Select a path for the new file")
    #         newData.to_file(path)
    #         if addToDB:
    #             if name is None:
    #                 if type(Geometry) == str:
    #                     name = Geometry
    #                 else:
    #                     raise KeyError("Select a name for the new area")
    #             datalayer.Measurements.addDocument(projectName=self._projectName, desc=(dict(name=name, **kwargs)),
    #                                                resource=path, type="Population", dataFormat="geopandas")
    #     return newData