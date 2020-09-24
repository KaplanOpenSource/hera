from ....datalayer import project
from shapely import geometry
import matplotlib.pyplot as plt

class datalayer(project.ProjectMultiDBPublic):

    _projectName = None
    _presentation = None

    @property
    def presentation(self):
        return self._presentation

    def __init__(self, projectName, databaseNameList=None, useAll=False,publicProjectName="Shapes"):

        self._projectName = projectName
        super().__init__(projectName=projectName, publicProjectName=publicProjectName,databaseNameList=databaseNameList,useAll=useAll)
        self._presentation = presentation(projectName=projectName,dataLayer=self)

    def getShape(self, name):
        """
        Returns the geometry shape of a given name from the database.
        Parameters:
            name: THe shape's name (string)
        Returns: The geometry (shapely Point or Polygon)

        """
        geo, shapeType = self.getShapePoints(name)
        if shapeType == "Polygon":
            geo = geometry.Polygon(geo)
        elif shapeType == "Point":
            geo = geometry.Point(geo[0])
        return geo

    def getShapePoints(self, name):
        """
        Returns the coordinates (list) and shape type ("Point" or "Polygon") of a geometry shape for a given name from the database.
        Parameters:
            name: THe shape's name (string)
        Returns: The geometry ([[ccoordinates], shapeType])
        """
        document = self.getMeasurementsDocumentsAsDict(name=name, type="Shape")
        if len(document) ==0:
            geo=None
            shapeType=None
        else:
            geo = document["documents"][0]["desc"]["geometry"]
            shapeType = document["documents"][0]["desc"]["shapeType"]

        return geo, shapeType

    def addShape(self, Shape, name):
        """
        This function is used to add a new geometry shape to the database.

        Parameters:
            Shape: The geometry shape to add to the database. Geometry must be given as one of the following structurs.\n
                      Shapely polygon or point, point coordinates ([x,y]), list of point coordinates ([[x1,y1],[x2,y2],...]),\n
                      list of x coordinates and y coordinates ([[x1,x2,...],[y1,y2,...]]) \n
            name: The name of the shape. (string)
        """

        check = self.getMeasurementsDocuments(name=name)
        KeyErrorText = "Shape must be given as one of the following structurs.\n" \
                        "Shapely polygon or point, point coordinates ([x,y]), list of point coordinates ([[x1,y1],[x2,y2],...]),\n" \
                        "list of x coordinates and y coordinates ([[x1,x2,...],[y1,y2,...]])"
        if len(check)>0:
            raise KeyError("Name is already used.")
        else:
            if type(Shape)==geometry.polygon.Polygon:
                geopoints = list(zip(*Shape.exterior.coords.xy))
                shapeType = "Polygon"
            elif type(Shape)==geometry.point.Point:
                geopoints = list(Shape.coords)
                shapeType = "Point"
            elif type(Shape)==list:
                if type(Shape[0])==list:
                    if len(Shape)>=3:
                        for geo in Shape:
                            if len(geo)!=2:
                                raise KeyError(KeyErrorText)
                        shapeType = "Polygon"
                        geopoints = Shape
                    elif len(Shape)==2:
                        if len(Shape[0])==len(Shape[1])>=3:
                            shapeType = "Polygon"
                            geopoints=[]
                            for i in range(len(Shape[0])):
                                geopoints.append([Shape[0][i], Shape[1][i]])
                        else:
                            raise KeyError(KeyErrorText)
                    else:
                        raise KeyError(KeyErrorText)
                else:
                    if len(Shape)!=2:
                        raise KeyError(KeyErrorText)
                    shapeType = "Point"
                    geopoints = [Shape]
            else:
                raise KeyError(KeyErrorText)
            if self._databaseNameList[0] == "public" or self._databaseNameList[0] == "Public" and len(
                    self._databaseNameList) > 1:
                userName = self._databaseNameList[1]
            else:
                userName = self._databaseNameList[0]
            self.addMeasurementsDocument(desc=dict(geometry=geopoints, shapeType=shapeType, name=name),
                                               type="Shape",
                                               resource="",
                                               dataFormat="string",users=userName)

class presentation():

    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, projectName, dataLayer=None, databaseNameList=None, useAll=False,
                 publicProjectName="Shapes"):

        self._datalayer = datalayer(projectName=projectName, publicProjectName=publicProjectName,
                         databaseNameList=databaseNameList, useAll=useAll) if datalayer is None else dataLayer


    def plot(self, names, color="black", marker="*", ax=None):
        """
        Plots saved geometry shapes.

        Parameters:
            names: The name/s of the shape/s (string or list of strings) \n
            color: The color of the shape (string) n\
            marker: The marker type for points. (string) \n
            ax: The ax of the plot. \n
            return: ax
        """
        if ax is None:
            fig, ax = plt.subplots(1,1)
        else:
            plt.sca(ax)
        if type(names) == list:
            for name in names:
                self._plotSingleShape(name, color, marker, ax)
        else:
            self._plotSingleShape(names, color, marker, ax)

        return ax

    def _plotSingleShape(self, name, color, marker, ax=None):

        if ax is None:
            fig, ax = plt.subplots(1,1)
        else:
            plt.sca(ax)

        geo, geometry_type = self.datalayer.getShapePoints(name)
        if geometry_type == "Point":
            plt.scatter(*geo[0], color=color, marker=marker)
        elif geometry_type == "Polygon":
            geo = self.datalayer.getShape(name)
            x, y = geo.exterior.xy
            ax = plt.plot(x, y, color=color)
        return ax