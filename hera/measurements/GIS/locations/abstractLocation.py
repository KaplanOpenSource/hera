from shapely import geometry
import os
from ....datalayer import project
from .shapes import datalayer as shapeDatalayer
import numpy
import xarray
import dask
import geopandas

class datalayer(project.ProjectMultiDBPublic):

    _projectName = None
    _publicProjectName = None
    _FilesDirectory = None

    @property
    def FilesDirectory(self):
        return self._FilesDirectory

    @property
    def projectName(self):
        return self._projectName

    @property
    def publicProjectName(self):
        return self._publicProjectName

    def __init__(self, projectName, FilesDirectory="", databaseNameList=None, useAll=False,publicProjectName="Topography",Source="BNTL"):

        self._projectName = projectName
        self._publicProjectName = publicProjectName
        super().__init__(projectName=projectName, publicProjectName=publicProjectName,useAll=useAll)
        self.setConfig({"source":Source})
        if FilesDirectory == "":
            self._FilesDirectory = os.getcwd()
        else:
            os.system("mkdir -p %s" % FilesDirectory)
            self._FilesDirectory = FilesDirectory


    def makeData(self, points, CutName, additional_data=None, Source=None):
        """
        Generates a new document that holds the path of a GIS shapefile.

        Parameters:
            points: Holds the ITM coordinates of a rectangle. It is a list, from the structure [minimum x, minimum y, maximum x, maximum y]\n
            CutName: Used as part of a new file's name. (string)\n
            additional_data: A dictionary with any additional parameters and their values.

        """
        fullPath = self.getMeasurementsDocumentsAsDict(source=Source)["documents"][0]["resource"]

        if additional_data is not None:
            additional_data["CutName"] = CutName
            additional_data["points"] = points
        else:
            additional_data = {"CutName": CutName, "points": points}

        documents = self.getMeasurementsDocumentsAsDict(points=points,type=self._publicProjectName)
        if len(documents) == 0:
            getattr(self,f"makeData_{Source}")(CutName=CutName,points=points,fullPath=fullPath,additional_data=additional_data)
        else:
            resource = documents["documents"][0]["resource"]
            dataFormat = documents["documents"][0]["dataFormat"]
            if self._databaseNameList[0] == "public" or self._databaseNameList[0] == "Public" and len(
                    self._databaseNameList) > 1:
                userName = self._databaseNameList[1]
            else:
                userName = self._databaseNameList[0]
            self.addMeasurementsDocument(desc=dict(**additional_data), type=self._publicProjectName,
                                               resource = resource, dataFormat = dataFormat,users=[userName])

    def makeData_BNTL(self, CutName,points,fullPath,additional_data,**kwargs):
        """
        Generates a new shapefile, used for BNTL source
        """

        FileName = "%s/%s.shp" % (self._FilesDirectory, CutName)
        os.system(
            "ogr2ogr -clipsrc %s %s %s %s %s %s" % (points[0], points[1], points[2], points[3], FileName, fullPath))
        self.addMeasurementsDocument(desc=additional_data, type=self._publicProjectName,
                                     resource=FileName, dataFormat="geopandas")

    def makeData_SRTM(self,CutName,points,additional_data,**kwargs):
        """
        Generates a new dask dataframe, used for SRTM source
        """

        FileName = "%s/%s.parquet" % (self._FilesDirectory, CutName)
        allData = self.getMeasurementsDocuments(source="SRTM",type=self._publicProjectName)[0].getData()
        dataArray = allData.read(1)
        xs = numpy.linspace(allData.bounds.left,allData.bounds.right,dataArray.shape[1])
        ys = numpy.linspace(allData.bounds.top,allData.bounds.bottom,dataArray.shape[0])
        xmin = numpy.where(xs>points[0])[0].min()
        ymin = numpy.where(ys<points[1])[0].min()
        xmax = numpy.where(xs<points[2])[0].max()
        ymax = numpy.where(ys>points[3])[0].max()
        xColumn = self.getConfig()["xColumn"] if "xColumn" in self.getConfig().keys() else "x"
        yColumn = self.getConfig()["yColumn"] if "yColumn" in self.getConfig().keys() else "y"
        heightColumn = self.getConfig()["heightColumn"] if "heightColumn" in self.getConfig().keys() else "height"
        data = xarray.DataArray(data=dataArray[ymax:ymin,xmin:xmax],dims=["yWGS84","xWGS84"],coords=[ys[ymax:ymin],xs[xmin:xmax]]).to_dataframe(heightColumn).reset_index()
        gdf = geopandas.GeoDataFrame(data, geometry=geopandas.points_from_xy(data["xWGS84"], data["yWGS84"]))
        gdf.crs = {"init": "epsg:4326"}
        gdf = gdf.to_crs({"init": "epsg:2039"})
        gdf[xColumn] = gdf.geometry.x
        gdf[yColumn] = gdf.geometry.y
        data = dask.dataframe.from_pandas(gdf.drop(columns="geometry"),npartitions=1)
        data.to_parquet(FileName,compression="GZIP")
        self.addMeasurementsDocument(desc=additional_data, type=self._publicProjectName,
                                     resource=FileName, dataFormat="parquet")

    def check_data(self, **kwargs):
        """
        Checks whether there is a document that fulfills desired requirements.
        Parameters:
            kwargs: Any desired requirements.

        Returns: True if there is a data that fulfills the requirement. False if there isn't.

        """

        check = self.getMeasurementsDocuments(**kwargs)

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
        if self._databaseNameList[0] == "public" or self._databaseNameList[0] == "Public" and len(
                self._databaseNameList) > 1:
            userName = self._databaseNameList[1]
        else:
            userName = self._databaseNameList[0]
        documents = self.getMeasurementsDocumentsAsDict(**kwargs,users=[userName])["documents"]
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
            kwargs: kwargs: Any desired requirements for the documents.

        Returns: List of polygons.
        """

        points = self.getFilesPointList(**kwargs)
        polygons = []
        for coordinates in points:
            polygons.append(geometry.Polygon([[coordinates[0], coordinates[1]], [coordinates[2], coordinates[1]],
                                              [coordinates[2], coordinates[3]], [coordinates[0], coordinates[3]]]))
        return polygons

    def getDocuments(self, points=None, CutName=None, ShapeMode="contains", Shape=None, Source=None, **kwargs):
        """
        This function is used to load GIS data.
        One may use it to get all data that corresponds to any parameters listed in a document,
        or to add a new document that relates to a file that holds GIS data in an area defined by a rectangle.
        Can also be used to perform geometrical queries.

        parameters:
            points: optional, for adding new data. Holds the ITM coordinates of a rectangle. It is a list, from the structure [minimum x, minimum y, maximum x, maximum y]\n
            CutName: optional, for adding new data. Used as part of a new file's name. (string)\n
            mode: The data type of the desired data. Recieves "Contour", "Buildings" or "Roads".\n
            ShapeMode: The mode of a geomtrical queries. Recieves "contains" or "intersects".\n
            Shape: A shapely geometry or a string with the name of a saved shapely geometry. Used to perform geometrical queries.\n
            **kwargs: any additional parameters that describe the data.
            return: The data.
        """
        Source = self.getConfig()["source"] if Source is None else Source
        if Shape is not None:
            if type(Shape)==str:
                try:
                    Shape = shapeDatalayer(projectName=self._projectName, databaseNameList=self._databaseNameList, useAll=self._useAll).getShape(Shape)
                except IndexError:
                    raise IndexError("Shape isn't defined.")
            containPoints = []
            points = self.getFilesPointList(**kwargs)
            polygons = self.getFilesPolygonList(**kwargs)
            for i in range(len(points)):
                if ShapeMode == "contains":
                    if polygons[i].contains(Shape):
                        containPoints.append(points[i])
                elif ShapeMode == "intersects":
                    if polygons[i].intersects(Shape):
                        containPoints.append(points[i])
                else:
                    raise KeyError("ShapeMode incorrectly called. Choose 'contains' or 'intersects'.")
            if 1 == len(containPoints):
                data = self.getMeasurementsDocuments(points=containPoints[0], type=self._publicProjectName, **kwargs)
            else:
                data = []
                for p in containPoints:
                    data += self.getMeasurementsDocuments(points=p, type=self._publicProjectName, **kwargs)

        else:
            if points==None and CutName==None:
                check = self.check_data(type=self._publicProjectName,**kwargs)
                if check:
                    data = self.getMeasurementsDocuments(type=self._publicProjectName,**kwargs)
            elif points==None and CutName!=None:
                check = self.check_data(CutName=CutName,type=self._publicProjectName, **kwargs)
                if check:
                    data = self.getMeasurementsDocuments(CutName=CutName,type=self._publicProjectName, **kwargs)
            elif CutName==None and points!=None:
                check = self.check_data(points=points,type=self._publicProjectName, **kwargs)
                if check:
                    data = self.getMeasurementsDocuments(CutName=CutName, type=self._publicProjectName,**kwargs)
            else:
                check = self.check_data(points=points, CutName=CutName,type=self._publicProjectName, **kwargs)
                if check:
                    data = self.getMeasurementsDocuments(points=points, CutName=CutName,type=self._publicProjectName, **kwargs)
            if not check:
                if points == None or CutName == None:
                    raise KeyError("Could not find data. Please insert points and CutName for making new data.")
                else:
                    self.makeData(points=points, CutName=CutName, Source=Source, additional_data=kwargs)
                    data = self.getMeasurementsDocuments(points=points, CutName=CutName)

        return data