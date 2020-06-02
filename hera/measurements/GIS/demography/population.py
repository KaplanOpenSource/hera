import geopandas

class population():

    def projectPolygonOnPopulation(self, data,poly,populationColumn="population",cityName="cityName"):

        poly = poly.buffer(0)
        res_intersect_poly = data[data["geometry"].intersect(poly)]