import geopandas
import shapely

class dataManipulations():

    def PolygonDataFrameIntersection(self, dataframe, polygon):
        """
        Creates a new dataframe based on the intersection of a dataframe and a polygon.
        Parameters:
        ----------
        dataframe: A geopandas dataframe.
        polygon: A shapely polygon

        Returns: A new geopandas dataframe

        """

        newlines = []
        for line in dataframe["geometry"]:
            newline = polygon.intersection(line)
            newlines.append(newline)
        dataframe["geometry"] = newlines
        dataframe = dataframe[~dataframe["geometry"].is_empty]

        return dataframe

    def getBoundaries(self, doc):

        dataframe = doc.getData()
        points = doc.desc["points"]

        boundaries = dict(xMin=points[0], xMax=points[2], yMin=points[1], yMax=points[3], zMin=dataframe["HEIGHT"].min(), zMax=dataframe["HEIGHT"].max())

        return boundaries

    def makePolygonFromEndPoints(self,points):

        polygon = shapely.geometry.Polygon([[points[0],points[1]],
                                            [points[0],points[3]],
                                            [points[2],points[3]],
                                            [points[2],points[1]]])
        return polygon