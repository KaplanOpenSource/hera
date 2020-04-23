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



        # xMin
        # 205000; // L = 350
        # xMax
        # 209000; // 4
        # yMin
        # 736000; // L = 280
        # yMax
        # 740000; // 7
        # zMin - 5;
        # zMax
        # 650;