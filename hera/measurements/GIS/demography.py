import geopandas
from .locations.shapes import datalayer as shapeDatalayer
from ...datalayer import project
from .locations.buildings import datalayer as buildingsDatalayer
import shapely

class datalayer(project.ProjectMultiDBPublic):

    _publicProjectName = None
    _projectName = None
    _analysis = None
    _Data = None


    @property
    def analysis(self):
        return self._analysis

    @property
    def publicProjectName(self):
        return self._publicProjectName

    def __init__(self, projectName,publicProjectName="Demography",databaseNameList=None, useAll=False,):
        self._projectName = projectName
        self._publicProjectName = publicProjectName
        super().__init__(projectName=projectName, publicProjectName=publicProjectName,
                         databaseNameList=databaseNameList, useAll=useAll)
        self.setConfig()
        self._analysis = analysis(projectName=projectName, dataLayer=self)

        datalist = self.getMeasurementsDocuments(source=self.getConfig()["source"])
        if len(datalist) > 0:
            self._Data = datalist[0].getData()
        else:
            self._Data = None

    def setConfig(self, Source="Lamas", units="WGS84", populationTypes = None, dbName=None, **kwargs):
        """
        Create a config documnet or updates an existing config document.
        """
        populationTypes = {"All":"total_pop","Children":"age_0_14","Youth":"age_15_19",
                           "YoungAdults":"age_20_29","Adults":"age_30_64","Elderly":"age_65_up"} if populationTypes is None else populationTypes
        config = dict(source=Source,units=units,populationTypes=populationTypes, **kwargs)
        super().setConfig(config=config, dbName=dbName)

        datalist = self.getMeasurementsDocuments(source=config["source"])

        if len(datalist) > 0:
            self._Data = datalist[0].getData()
        else:
            self._Data = None

    def projectPolygonOnPopulation(self, Shape, projectName=None, populationTypes="All", Data=None):
        """
        Finds the population in a polygon.
        Params:
            Shape: The polygon, either a shapely polygon or a name of a saved geometry in the database.
            populationTypes: A string or a list of strings with options from config/population or column names
                             from the data.
        """

        Data = self._Data if Data is None else Data
        if type(Shape) == str:
            sDatalayer = shapeDatalayer(projectName=projectName, databaseNameList=self._databaseNameList, useAll=self._useAll)
            poly = sDatalayer.getShape(Shape)
            if poly is None:
                documents = sDatalayer.getMeasurementsDocuments(CutName=Shape)
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
            if "populationTypes" in self.getConfig():
                if populationType in self.getConfig()["populationTypes"]:
                    populationType = self.getConfig()["populationTypes"][populationType]
            res_intersection[populationType] = intersection_poly.area / res_intersect_poly.area * res_intersect_poly[populationType]

        return res_intersection

class analysis():

    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, projectName, dataLayer=None, databaseNameList=None, useAll=False,
                 publicProjectName="Demography", Source="Lamas"):

        self._datalayer = datalayer(projectName=projectName, publicProjectName=publicProjectName,
                         databaseNameList=databaseNameList, useAll=useAll, Source=Source) if datalayer is None else dataLayer

    def populateNewArea(self, Shape, projectName, populationTypes=None, convex=True,save=False, addToDB=False, path=None, name=None, Data=None, **kwargs):
        """
        make a geodataframe with a selected polygon as the geometry, and the sum of the population in the polygons that intersect it as its population.
        """

        Data = self.datalayer._Data if Data is None else Data
        if type(Shape) == str:
            sDatalayer = shapeDatalayer(projectName=projectName, databaseNameList=self.datalayer._databaseNameList, useAll=self.datalayer._useAll)
            poly = sDatalayer.getShape(Shape)
            if poly is None:
                documents = sDatalayer.getMeasurementsDocuments(CutName=Shape)
                if len(documents) == 0:
                    raise KeyError("Shape %s was not found" % Shape)
                else:
                    if convex:
                        polys = buildingsDatalayer(projectName=projectName).analysis.ConvexPolygons(documents[0].getData())
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

        if save:
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
