from .... import datalayer
import pandas
from shapely import geometry
import os

class GIS_datalayer():

    _Measurments = None
    _projectName = None
    _FilesDirectory = None

    def __init__(self, projectName, FilesDirectory):

        self._FilesDirectory = FilesDirectory
        self._projectName = projectName
        self._Measurments = datalayer.Measurements

        os.system("mkdir -p %s" % self._FilesDirectory)

    def getExistingDocuments(self, **kwargs):
        """
        Loads data using an existing document in the database.

        kwargs: any value of any parameter in the database, by which to select the data.

        Returns: The data
        """

        data = self._Measurments.getDocuments(projectName=self._projectName, **kwargs)

        return data

    def makeData(self, points, CutName, mode="Contour", additional_data=None):
        """
        Generates a new document that holds the path of a GIS shapefile.

        Parameters:
            points: Holds the ITM coordinates of a rectangle. It is a list, from the structure [minimum x, minimum y, maximum x, maximum y]\n
            CutName: Used as part of a new file's name. (string)\n
            mode: The data type of the desired data. Recieves "Contour", "Buildings" or "Roads".\n
            additional_data: A dictionary with any additional parameters and their values.

        """

        fullfilesdirect = {"Contour": "RELIEF/CONTOUR.shp", "Buildings": "BUILDINGS/BLDG.shp",
                           "Roads": "TRANSPORTATION/MAIN_ROAD.shp"}

        if additional_data is not None:
            additional_data["CutName"] = CutName
            additional_data["points"] = points
            additional_data["mode"] = mode
        else:
            additional_data = {"CutName": CutName, "points": points, "mode": mode}

        documents = self._Measurments.getDocuments(projectName=self._projectName, points=points, mode=mode)
        if len(documents) == 0:

            FileName = "%s//%s%s-%s.shp" % (self._FilesDirectory, self._projectName, CutName, mode)

            os.system("ogr2ogr -clipsrc %s %s %s %s %s /mnt/public/New-MAPI-data/BNTL_MALE_ARZI/BNTL_MALE_ARZI/%s" % (points[0],points[1],points[2],points[3], FileName, fullfilesdirect[mode]))
            datalayer.Measurements.addDocument(projectName=self._projectName, desc=dict(**additional_data), type="GIS",
                                               resource = FileName, dataFormat = "geopandas")
        else:
            resource = documents[0].asDict()["resource"]
            datalayer.Measurements.addDocument(projectName=self._projectName, desc=dict(**additional_data), type="GIS",
                                               resource = resource, dataFormat = "geopandas")

    def check_data(self, **kwargs):
        """
        Checks whether there is a document that fulfills desired requirements.
        Parameters:
            kwargs: Any desired requirements.

        Returns: True if there is a data that fulfills the requirement. False if there isn't.

        """

        check = self._Measurments.getDocuments(projectName=self._projectName, **kwargs)

        if 0 == len(check):
            result = False
        else:
            result = True

        return result

    def getFilesPointList(self, **kwargs):
        """
        Returns a list of all the coordinates used for defining areas in existing documents.

        Parameters:
            kwargs: kwargs: Any desired requirements for the documents.

        Returns: List of lists, each contains 4 coordinates - [minimum x, minimum y, maximum x, maximum y].

        """

        documents = self._Measurments.getDocuments(projectName=self._projectName, **kwargs)
        points = []
        for document in documents:
            data = document.asDict()
            if "points" in data["desc"].keys():
                if data["desc"]["points"] not in points:
                    points.append(data["desc"]["points"])
        return points

    def getFilesPolygonList(self, **kwargs):
        """
        Returns a list of all the polygons used for defining areas in existing documents.

        Parameters:
            kwargs: kwargs: Any desired requirements for the documents.

        Returns: List of polygons.
        """

        points = self.getFilesPointList(**kwargs)
        polygons = []
        for coordinates in points:
            polygons.append(geometry.Polygon([[coordinates[0], coordinates[1]], [coordinates[2], coordinates[1]],
                                              [coordinates[2], coordinates[3]], [coordinates[0], coordinates[3]]]))
        return polygons

    def getGeometry(self, name):
        """
        Returns the geometry shape of a given name from the database.
        Parameters:
            name: THe shape's name (string)
        Returns: The geometry (shapely Point or Polygon)

        """

        geo, geometry_type = self.getGeometryPoints(name)
        if geometry_type == "Polygon":
            geo = geometry.Polygon(geo)
        elif geometry_type == "Point":
            geo = geometry.Point(geo[0])
        return geo

    def getGeometryPoints(self, name):
        """
        Returns the coordinates (list) and shape type ("Point" or "Polygon") of a geometry shape for a given name from the database.
        Parameters:
            name: THe shape's name (string)
        Returns: The geometry ([[ccoordinates], geometry_type])

        """

        document = self._Measurments.getDocuments(projectName=self._projectName, name=name)
        geo = document[0].asDict()["desc"]["geometry"]
        geometry_type = document[0].asDict()["desc"]["geometry_type"]

        return geo, geometry_type

    def addGeometry(self, Geometry, name):
        """
        This function is used to add a new geometry shape to the database.

        Parameters:
            Geometry: The geometry shape to add to the database. Geometry must be given as one of the following structurs.\n
                      Shapely polygon or point, point coordinates ([x,y]), list of point coordinates ([[x1,y1],[x2,y2],...]),\n
                      list of x coordinates and y coordinates ([[x1,x2,...],[y1,y2,...]]) \n
            name: The name of the shape. (string)
        """

        check = self.check_data(name=name)
        KeyErrorText = "Geometry must be given as one of the following structurs.\n" \
                        "Shapely polygon or point, point coordinates ([x,y]), list of point coordinates ([[x1,y1],[x2,y2],...]),\n" \
                        "list of x coordinates and y coordinates ([[x1,x2,...],[y1,y2,...]])"
        if check:
            raise KeyError("Name is already used.")
        else:
            if type(Geometry)==geometry.polygon.Polygon:
                geopoints = list(zip(*Geometry.exterior.coords.xy))
                geometry_type = "Polygon"
            elif type(Geometry)==geometry.point.Point:
                geopoints = list(Geometry.coords)
                geometry_type = "Point"
            elif type(Geometry)==list:
                if type(Geometry[0])==list:
                    if len(Geometry)>=3:
                        for geo in Geometry:
                            if len(geo)!=2:
                                raise KeyError(KeyErrorText)
                        geometry_type = "Polygon"
                        geopoints = Geometry
                    elif len(Geometry)==2:
                        if len(Geometry[0])==len(Geometry[1])>=3:
                            geometry_type = "Polygon"
                            geopoints=[]
                            for i in range(len(Geometry[0])):
                                geopoints.append([Geometry[0][i], Geometry[1][i]])
                        else:
                            raise KeyError(KeyErrorText)
                    else:
                        raise KeyError(KeyErrorText)
                else:
                    if len(Geometry)!=2:
                        raise KeyError(KeyErrorText)
                    geometry_type = "Point"
                    geopoints = [Geometry]
            else:
                raise KeyError(KeyErrorText)
            datalayer.Measurements.addDocument(projectName=self._projectName,
                                               desc=dict(geometry=geopoints, geometry_type=geometry_type, name=name),
                                               type="GeometryShape",
                                               resource="/mnt/public/New-MAPI-data/BNTL_MALE_ARZI/BNTL_MALE_ARZI/RELIEF/CONTOUR.shp",
                                               dataFormat="geopandas")

    def getGISDocuments(self, points=None, CutName=None, mode="Contour", GeometryMode="contains", Geometry=None, **kwargs):
        """
        This function is used to load GIS data.
        One may use it to get all data that corresponds to any parameters listed in a document,
        or to add a new document that relates to a file that holds GIS data in an area defined by a rectangle.
        Can also be used to perform geometrical queries.

        parameters:
            points: optional, for adding new data. Holds the ITM coordinates of a rectangle. It is a list, from the structure [minimum x, minimum y, maximum x, maximum y]\n
            CutName: optional, for adding new data. Used as part of a new file's name. (string)\n
            mode: The data type of the desired data. Recieves "Contour", "Buildings" or "Roads".\n
            GeometryMode: The mode of a geomtrical queries. Recieves "contains" or "intersects".\n
            Geometry: A shapely geometry or a string with the name of a saved shapely geometry. Used to perform geometrical queries.\n
            **kwargs: any additional parameters that describe the data.
            return: The data.
        """

        if Geometry is not None:
            if type(Geometry)==str:
                try:
                    Geometry = self.getGeometry(Geometry)
                except IndexError:
                    raise IndexError("Geometry isn't defined.")
            containPoints = []
            points = self.getFilesPointList(**kwargs)
            polygons = self.getFilesPolygonList(**kwargs)
            for i in range(len(points)):
                if GeometryMode == "contains":
                    if polygons[i].contains(Geometry):
                        containPoints.append(points[i])
                elif GeometryMode == "intersects":
                    if polygons[i].intersects(Geometry):
                        containPoints.append(points[i])
                else:
                    raise KeyError("GeometryMode incorrectly called. Choose 'contains' or 'intersects'.")
            if 1 == len(containPoints):
                data = self.getExistingDocuments(points=containPoints[0], **kwargs)
            else:
                data = []
                for p in containPoints:
                    data.append(self.getExistingDocuments(points=p, **kwargs))

        else:
            if points==None and CutName==None:
                check = self.check_data(**kwargs)
            elif points==None and CutName!=None:
                check = self.check_data(CutName=CutName, mode=mode, **kwargs)
            elif CutName==None and points!=None:
                check = self.check_data(points=points, mode=mode, **kwargs)
            else:
                check = self.check_data(points=points, CutName=CutName, mode=mode, **kwargs)

            if check:
                data = self.getExistingDocuments(mode=mode, **kwargs)
            else:
                if points == None or CutName == None:
                    raise KeyError("Could not find data. Please insert points and CutName for making new data.")
                else:
                    self.makeData(points=points, CutName=CutName, mode=mode, additional_data=kwargs)
                    data = self.getExistingDocuments(points=points, CutName=CutName, mode=mode)

        return data

    def loadImage(self, path, locationName, extents):
        """
        Loads an image to the database.

        :param projectName: The project name
        :param path: The image path
        :param locationName: The location name
        :param extents: The extents of the image [left, right, bottom, top]
        :return:
        """
        doc = dict(projectName=self._projectName,
                   resource=path,
                   dataFormat='image',
                   type='GIS',
                   desc=dict(locationName=locationName,
                             left=extents[0],
                             right=extents[1],
                             bottom=extents[2],
                             top=extents[3]
                             )
                   )
        self._Measurments.addDocument(**doc)

