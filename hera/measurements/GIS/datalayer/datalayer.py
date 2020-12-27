from .... import datalayer
import pandas
from shapely import geometry
import os

class GIS_datalayer:

    _projectMultiDB = None
    _projectName = None
    _FilesDirectory = None

    @property
    def project(self):
        return self._projectMultiDB

    def __init__(self, projectName, FilesDirectory, users=[None], useAll=False):

        self._FilesDirectory = FilesDirectory
        self._projectName = projectName
        self._projectMultiDB = datalayer.ProjectMultiDB(projectNameDict=projectName, users=users, useAll=useAll)

        os.system("mkdir -p %s" % self._FilesDirectory)

    def getExistingDocuments(self, **kwargs):
        """
        Loads data using an existing document in the database.

        kwargs: any value of any parameter in the database, by which to select the data.

        Returns: The data
        """
        data = self._projectMultiDB.getMeasurementsDocuments(**kwargs)

        return data

    def makeData(self, points, CutName, mode="Contour", additional_data=None, useOwn=False):
        """
        Generates a new document that holds the path of a GIS shapefile.
        Parameters:
            points: Holds the ITM coordinates of a rectangle. It is a list, from the structure [minimum x, minimum y, maximum x, maximum y]\n
            CutName: Used as part of a new file's name. (string)\n
            mode: The data type of the desired data. Recieves any mode specified in the GISOrigin document.\n
            additional_data: A dictionary with any additional parameters and their values.

        """
        if useOwn:
            fullfilesdirect = self._projectMultiDB.getMeasurementsDocumentsAsDict(type="GISOrigin")["documents"][0]["desc"]["modes"]
            path = self._projectMultiDB.getMeasurementsDocumentsAsDict(type="GISOrigin")["documents"][0]["resource"]
            fullPath = "%s/%s" % (path, fullfilesdirect[mode])
        else:
            publicproject = datalayer.ProjectMultiDB(projectNameDict="PublicData", users=["public"])
            fullPath = publicproject.getMeasurementsDocumentsAsDict(type="GIS",mode=mode)["documents"][0]["resource"]

        if additional_data is not None:
            additional_data["CutName"] = CutName
            additional_data["points"] = points
            additional_data["mode"] = mode
        else:
            additional_data = {"CutName": CutName, "points": points, "mode": mode}

        documents = self._projectMultiDB.getMeasurementsDocumentsAsDict(points=points, mode=mode)
        if len(documents) == 0:

            FileName = "%s//%s%s-%s.shp" % (self._FilesDirectory, self._projectName, CutName, mode)

            os.system("ogr2ogr -clipsrc %s %s %s %s %s %s" % (points[0],points[1],points[2],points[3], FileName,fullPath))
            self._projectMultiDB.addMeasurementsDocument(desc=additional_data, type="GIS",
                                               resource = FileName, dataFormat = "geopandas")
        else:
            resource = documents["documents"][0]["resource"]
            self._projectMultiDB.addMeasurementsDocument(desc=dict(**additional_data), type="GIS",
                                               resource = resource, dataFormat = "geopandas")


    def check_data(self, **kwargs):
        """
        Checks whether there is a document that fulfills desired requirements.
        Parameters:
            kwargs: Any desired requirements.

        Returns: True if there is a data that fulfills the requirement. False if there isn't.

        """

        check = self._projectMultiDB.getMeasurementsDocuments(**kwargs)

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

        documents = self._projectMultiDB.getMeasurementsDocumentsAsDict(**kwargs)["documents"]
        points = []
        for document in documents:
            if "points" in document["desc"].keys():
                if document["desc"]["points"] not in points:
                    points.append(document["desc"]["points"])
        return points

    def getFilesPolygonList(self, **kwargs):
        """
        Returns a list of all the polygons used for defining areas in existing documents.

        Parameters:
        -------
            kwargs: MongoDB Query language.

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
        -----------
            name: str
                The shape's name.

        Returns:
        --------
            shapely Point or Polygon

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
        -----------
            name: str
                The shape's name (string)
        Returns:
        -----------
            A tuple: ([ccoordinates], geometry_type)

        """
        document = self._projectMultiDB.getMeasurementsDocumentsAsDict(name=name, type="GeometryShape")
        if len(document) ==0:
            geo=None
            geometry_type=None
        else:
            geo = document["documents"][0]["desc"]["geometry"]
            geometry_type = document["documents"][0]["desc"]["geometry_type"]

        return geo, geometry_type

    def addGeometry(self, geometry, name):
        """
        This function is used to add a new geometry shape to the database.

        Parameters:
        -----------
            geometry:
                      The geometry shape to add to the database.
                      geometry must be given as one of the following structurs.\n

                      * Shapely polygon or point
                      * point coordinates ([x,y]), list of point coordinates ([[x1,y1],[x2,y2],...])
                      * list of x coordinates and y coordinates ([[x1,x2,...],[y1,y2,...]])
            name: str
                The name of the shape.

        Returns
        -------
            None
        """

        check = self.check_data(name=name)
        KeyErrorText = "geometry must be given as one of the following structurs.\n" \
                        "Shapely polygon or point, point coordinates ([x,y]), list of point coordinates ([[x1,y1],[x2,y2],...]),\n" \
                        "list of x coordinates and y coordinates ([[x1,x2,...],[y1,y2,...]])"
        if check:
            raise KeyError("Name is already used.")
        else:
            if type(geometry)==geometry.polygon.Polygon:
                geopoints = list(zip(*geometry.exterior.coords.xy))
                geometry_type = "Polygon"
            elif type(geometry)==geometry.point.Point:
                geopoints = list(geometry.coords)
                geometry_type = "Point"
            elif type(geometry)==list:
                if type(geometry[0])==list:
                    if len(geometry)>=3:
                        for geo in geometry:
                            if len(geo)!=2:
                                raise KeyError(KeyErrorText)
                        geometry_type = "Polygon"
                        geopoints = geometry
                    elif len(geometry)==2:
                        if len(geometry[0])==len(geometry[1])>=3:
                            geometry_type = "Polygon"
                            geopoints=[]
                            for i in range(len(geometry[0])):
                                geopoints.append([geometry[0][i], geometry[1][i]])
                        else:
                            raise KeyError(KeyErrorText)
                    else:
                        raise KeyError(KeyErrorText)
                else:
                    if len(geometry)!=2:
                        raise KeyError(KeyErrorText)
                    geometry_type = "Point"
                    geopoints = [geometry]
            else:
                raise KeyError(KeyErrorText)
            self._projectMultiDB.addMeasurementsDocument(desc=dict(geometry=geopoints, geometry_type=geometry_type, name=name),
                                               type="GeometryShape",
                                               resource="/mnt/public/New-MAPI-data/BNTL_MALE_ARZI/BNTL_MALE_ARZI/RELIEF/CONTOUR.shp",
                                               dataFormat="geopandas")

    def getGISDocuments(self, points=None, cutName=None, mode="Contour", geometryMode="contains", geometry=None, useOwn=False, **kwargs):
        """
        This function is used to load GIS data.
        One may use it to get all data that corresponds to any parameters listed in a document,
        or to add a new document that relates to a file that holds GIS data in an area defined by a rectangle.
        Can also be used to perform geometrical queries.

        Parameters:
        -----------
            points:
                Holds the ITM coordinates of a rectangle.
                If is a list, from the structure [minimum x, minimum y, maximum x, maximum y]

            cutName:
                optional, for adding new data. Used as part of a new file's name. (string)\n

            mode:
                The data type of the desired data. Recieves "Contour", "Buildings" or "Roads".\n

            geometryMode:
                The mode of a geomtrical queries. Recieves "contains" or "intersects".\n

            geometry:
                A shapely geometry or a string with the name of a saved shapely geometry. Used to perform geometrical queries.\n

            kwargs: any additional parameters that describe the data.

            Returns:
            --------
                The data.
        """

        if geometry is not None:
            if type(geometry)==str:
                try:
                    geometry = self.getGeometry(geometry)
                except IndexError:
                    raise IndexError("geometry isn't defined.")
            containPoints = []
            points = self.getFilesPointList(**kwargs)
            polygons = self.getFilesPolygonList(**kwargs)
            for i in range(len(points)):
                if geometryMode == "contains":
                    if polygons[i].contains(geometry):
                        containPoints.append(points[i])
                elif geometryMode == "intersects":
                    if polygons[i].intersects(geometry):
                        containPoints.append(points[i])
                else:
                    raise KeyError("geometryMode incorrectly called. Choose 'contains' or 'intersects'.")
            if 1 == len(containPoints):
                data = self.getExistingDocuments(points=containPoints[0], **kwargs)
            else:
                data = []
                for p in containPoints:
                    data.append(self.getExistingDocuments(points=p, **kwargs))

        else:
            if points==None and cutName==None:
                check = self.check_data(**kwargs)
            elif points==None and cutName!=None:
                check = self.check_data(CutName=cutName, mode=mode, **kwargs)
            elif cutName==None and points!=None:
                check = self.check_data(points=points, mode=mode, **kwargs)
            else:
                check = self.check_data(points=points, CutName=cutName, mode=mode, **kwargs)

            if check:
                data = self.getExistingDocuments(mode=mode, points=points, CutName=cutName, **kwargs)
            else:
                if points == None or cutName == None:
                    raise KeyError("Could not find data. Please insert points and cutName for making new data.")
                else:
                    self.makeData(points=points, CutName=cutName, mode=mode, useOwn=useOwn, additional_data=kwargs)
                    data = self.getExistingDocuments(points=points, CutName=cutName, mode=mode)

        return data


    def loadImage(self, path, locationName, extents):
        """
        Loads an image to the local database.

        Parameters:
        -----------

        projectName: str
                    The project name
        path:  str
                    The image path
        locationName: str
                    The location name
        extents: list or dict
                list: The extents of the image [xmin, xmax, ymin, ymax]
                dict: A dict with the keys xmin,xmax,ymin,ymax


        Returns
        -------
        """

        if isinstance(extents,dict):
            extentList = [extents['xmin'],extents['xmax'],extents['ymin'],extents['ymax']]
        elif isinstance(extents,list):
            extentList = extents
        else:
            raise ValueError("extents is either a list(xmin, xmax, ymin, ymax) or dict(xmin=, xmax=, ymin=, ymax=) ")


        doc = dict(resource=path,
                   dataFormat='image',
                   type='GIS',
                   desc=dict(locationName=locationName,
                             xmin=extentList[0],
                             xmax=extentList[1],
                             ymin=extentList[2],
                             ymax=extentList[3]
                             )
                   )
        self._projectMultiDB.addMeasurementsDocument(**doc)
