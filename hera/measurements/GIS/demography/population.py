import geopandas
from .... import datalayer
from ..datalayer.datalayer import GIS_datalayer
from ..analytics.dataManipulations import dataManipulations

class population():

    _publicMeasure = None
    _projectName = None
    _populationDict = None

    @property
    def agesDefinition(self):
        return self._populationDict

    def __init__(self, projectName):

        self._publicMeasure = datalayer.Measurements_Collection(user="public")
        self._projectName = projectName
        self._populationDict = {"All":"total_pop","Children":"AGE_0_14","Youth":"AGE_15_19","YoungAdults":"age_20_29","Adults":"age_30_64","Elderly":"AGE_65_up"}

    def projectPolygonOnPopulation(self, Geometry, data=None, populationTypes="All"):
        """
        Finds the population in a polygon.
        Params:
            Geometry: The polygon, either a shapely polygon or a name of a saved geometry in the database.
            data: Default is None. In that case, the function uses the population data saved in PublicData.
                  Otherwise, a geodataframe.
            populationTypes: A string or a list of strings with options from self._populationDict or column names
                             from the data.
        """

        if data is None:
            Data = self._publicMeasure.getDocuments(projectName="PublicData", type="Population")[0].getData()
            Data.crs= {'init' :'epsg:4326'}
            Data = Data.to_crs(epsg=2039)
        else:
            Data = data
        if type(populationTypes) == str:
            populationTypes = [populationTypes]

        if type(Geometry) == str:
            poly = GIS_datalayer(projectName=self._projectName, FilesDirectory="").getGeometry(name=Geometry)
            if poly is None:
                documents=GIS_datalayer(projectName=self._projectName, FilesDirectory="").getExistingDocuments(CutName=Geometry)
                if len(documents)==0:
                    raise KeyError("Geometry %s was not found" % Geometry)
                else:
                    poly = dataManipulations().makePolygonFromEndPoints(documents[0].asDict()["desc"]["points"])
        else:
            poly=Geometry
        poly = poly.buffer(0)
        res_intersect_poly = Data.loc[Data["geometry"].intersection(poly).is_empty == False]
        intersection_poly = res_intersect_poly["geometry"].intersection(poly)
        res_intersection = geopandas.GeoDataFrame.from_dict(
            {"geometry": intersection_poly.geometry,
             "areaFraction": intersection_poly.area/res_intersect_poly.area})
        for populationType in populationTypes:
            if data is None:
                populationType = self._populationDict[populationType]
                print(populationType)
            res_intersection[populationType] = intersection_poly.area / res_intersect_poly.area * res_intersect_poly[populationType]

        return res_intersection