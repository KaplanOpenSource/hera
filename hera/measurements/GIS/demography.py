import geopandas
import io
from ... import toolkit
from ...datalayer import datatypes
from .shapes import ShapesToolKit
from .locations.buildings import BuildingsToolkit
import shapely

from ...toolkit import TOOLKIT_SAVEMODE_NOSAVE,TOOLKIT_SAVEMODE_ONLYFILE,TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,TOOLKIT_SAVEMODE_FILEANDDB,TOOLKIT_SAVEMODE_FILEANDDB_REPLACE



class DemographyToolkit(toolkit.abstractToolkit):
    """
        A toolkit to manage demography data

    """

    _populationTypes = None # A dictionary with the names of the population types.

    _shapes = None
    _buildings = None

    @property
    def buildings(self):
        return self._buildings

    @property
    def shapes(self):
        return self._shapes


    def __init__(self, projectName):
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
            If sepefied load this version of the data.

        """
        self._projectName = projectName
        super().__init__(projectName=projectName,toolkitName="Demography")

        self._populationTypes = {"All":"total_pop","Children":"age_0_14","Youth":"age_15_19",
                           "YoungAdults":"age_20_29","Adults":"age_30_64","Elderly":"age_65_up"}

        self._shapes    = ShapesToolKit(projectName=projectName)
        self._buildings = BuildingsToolkit(projectName=projectName)


    def projectPolygonOnPopulation(self, Shape, projectName=None, populationTypes="All", Data=None):
        import warnings
        warnings.warn("Depracted in Version 2.0.0+. Use analysis.calculatePopulationInPolygon",
                      category=DeprecationWarning,
                      stacklevel=2)
        self.analysis.calculatePopulationInPolygon(Shape=Shape, projectName=projectName, populationTypes=populationTypes, Data=Data)

    @property
    def populationTypes(self):
        return self._populationTypes


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


    def populateNewArea(self,
                        Shape,
                        dataSourceOrData,
                        dataSourceVersion=None,
                        populationTypes=None,
                        convex=True,
                        saveMode=TOOLKIT_SAVEMODE_NOSAVE,
                        path=None,
                        regionName=None, **kwargs):
        """
            Make a geoDataFrame with a selected polygon as the geometry,
            and the sum of the population in the polygons that intersect it as its population.

        Parameters
        -----------
        Shape:
        dataSourceOrData:
        dataSourceVersion: 3-tuple of int

        populationTypes:
        convex:
        saveMode:
        path:
        regionName:
        kwargs:

        Returns
        -------

        """
        if isinstance(dataSourceOrData,str):
            Data = self.datalayer.getDataSourceData(dataSourceOrData,dataSourceVersion)
        else:
            Data = dataSourceOrData

        if isinstance(Shape,str):
            polydoc = self.datalayer.shapes.getShape(Shape)
            if convex:
                polys = self.datalayer.buildings.analysis.ConvexPolygons(polydoc)
                poly = polys.loc[polys.area==polys.area.max()].geometry[0]
            else:
                poly = polydoc.unary_union
        else:
            poly = Shape

        res_intersect_poly = Data.loc[Data["geometry"].intersection(poly).is_empty == False]

        newData = geopandas.GeoDataFrame.from_dict([{"geometry": poly}])

        populationTypes = self.datalayer.populationTypes if populationTypes is None else populationTypes
        for populationType in populationTypes:
            if populationType in res_intersect_poly:
                newData[populationType] = res_intersect_poly.sum()[populationType]

        if saveMode !=toolkit.TOOLKIT_SAVEMODE_NOSAVE:
            if path is None:
                raise KeyError("Select a path for the new file")
            newData.to_file(path)
            if saveMode == toolkit.TOOLKIT_SAVEMODE_FILEANDDB
                if regionName is None:
                    if type(Shape) == str:
                        regionName = Shape
                    else:
                        raise KeyError("Select a regionName for the new area")


                desc = {toolkit.TOOLKIT_DATASOURCE_NAME : regionName, toolkit.TOOLKIT_TOOLKITNAME_FIELD : self.name}
                desc.update(**kwargs)
                self.datalayer.addMeasurementsDocument(desc=desc,
                                                       resource=path,
                                                       type=toolkit.TOOLKIT_DATASOURCE_TYPE,
                                                       dataFormat=datatypes.GEOPANDAS)
        return newData

    def calculatePopulationInPolygon(self,
                                     shapeNameOrData,
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


        if isinstance(shapeNameOrData,str):
            poly = self.datalayer.shapes.getShape(shapeNameOrData)
            if poly is None:
                poly = geopandas.read_file(io.StringIO(shapeNameOrData))
        else:
            poly = shapeNameOrData

        if isinstance(dataSourceOrData,str):
            demography = self.datalayer.getDataSourceData(dataSourceOrData,dataSourceVersion)
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

