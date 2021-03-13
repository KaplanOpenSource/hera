from . import abstractLocation
from shapely import geometry
import matplotlib.pyplot as plt
from ....toolkit import TOOLKIT_SAVEMODE_NOSAVE,TOOLKIT_SAVEMODE_ONLYFILE,TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,TOOLKIT_SAVEMODE_FILEANDDB,TOOLKIT_SAVEMODE_FILEANDDB_REPLACE
from ....datalayer import datatypes
import os

class ShapesToolKit(abstractLocation.AbstractLocationToolkit):
    """
        Holds geoJSON shape files.

    """

    _presentation = None

    @property
    def presentation(self):
        return self._presentation

    @property
    def doctype(self):
        return f"{self.name}_GeoJSON"



    def __init__(self, projectName,FilesDirectory=None ):


        super().__init__(projectName=projectName,
                         toolkitName="Shapes",
                         FilesDirectory=FilesDirectory)

        self._presentation = presentation(dataLayer=self)

    def getShape(self, regionName):
        """
            Returns the geometry shape of a given name from the database.

        Parameters
        -----------
            regionName: The shape's name (string)

        Returns
        --------
            The geometry (shapely Point or Polygon)

        """
        geo, shapeType = self.getShapePoints(regionName)
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

    def loadData(self, fileNameOrData, extents, saveMode=TOOLKIT_SAVEMODE_NOSAVE, regionName=None, additionalData=dict()):
        """
            Loading a data from file, and possibly store the region in the database.

        Parameters
        ----------
        fileNameOrData: str
                If str , the datafile to load
                If other objects - convert the
        parser: str
                The name of the parser to use

        saveMode: str
                Can be either:

                    - TOOLKIT_SAVEMODE_NOSAVE   : Just load the data from file and return the datafile

                    - TOOLKIT_SAVEMODE_ONLYFILE : Loads the data from file and save to a file.
                                                  raise exception if file exists.

                    - TOOLKIT_SAVEMODE_ONLYFILE_REPLACE: Loads the data from file and save to a file.
                                                  Replace the file if it exists.

                    - TOOLKIT_SAVEMODE_FILEANDDB : Loads the data from file and save to a file and store to the DB as a source.
                                                    Raise exception if the entry exists.

                    - TOOLKIT_SAVEMODE_FILEANDDB_REPLACE: Loads the data from file and save to a file and store to the DB as a source.
                                                    Replace the entry in the DB if it exists.


            kwargs: Contains:
                regionName: If fileNameOrData is an object, required if saveMode is not NOSAVE.
                additionalData: additional metadata to add if adding to the DB.

        Returns
        -------
            The data or the doc.

            Return the data if the saveMode is either [ TOOLKIT_SAVEMODE_NOSAVE, TOOLKIT_SAVEMODE_ONLYFILE, TOOLKIT_SAVEMODE_ONLYFILE_REPLACE].
            Return the DB document is the saveMode is either  [TOOLKIT_SAVEMODE_FILEANDDB, TOOLKIT_SAVEMODE_FILEANDDB_REPLACE].
        """

        if isinstance(fileNameOrData,str):


            if os.path.exists(os.path.abspath(fileNameOrData)):

                regionName = os.path.basename(fileNameOrData).split(".")[0] if regionName is None else regionName

                data = mpimg.imread(os.path.abspath(fileNameOrData))
            else:
                raise FileNotFoundError(f"The {fileNameOrData} does not exist.")

        else:
            regionName = None
            data = fileNameOrData

        doc = None

        if saveMode in [TOOLKIT_SAVEMODE_ONLYFILE,
                        TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,
                        TOOLKIT_SAVEMODE_FILEANDDB,
                        TOOLKIT_SAVEMODE_FILEANDDB_REPLACE]:

            outputFileName = os.path.join(self.FilesDirectory, f"{regionName}.png")

            if saveMode in [TOOLKIT_SAVEMODE_FILEANDDB,TOOLKIT_SAVEMODE_ONLYFILE]:
                if os.path.exists(outputFileName):
                    raise FileExistsError(f"{outputFileName} exists in project {self.projectName}")

            mpimg.imsave(outputFileName,data)

            if saveMode in [TOOLKIT_SAVEMODE_FILEANDDB,TOOLKIT_SAVEMODE_FILEANDDB_REPLACE]:

                doc = self.getDatasourceData(regionName)
                if doc is not None and saveMode==TOOLKIT_SAVEMODE_FILEANDDB:
                    raise ValueError(f"{regionName} exists in DB for project {self.projectName}")

                if isinstance(extents, dict):
                    extentList = [extents['xmin'], extents['xmax'], extents['ymin'], extents['ymax']]
                elif isinstance(extents, list):
                    extentList = extents
                else:
                    raise ValueError(
                        "extents is either a list(xmin, xmax, ymin, ymax) or dict(xmin=, xmax=, ymin=, ymax=) ")

                additionalData.update({abstractLocation.TOOLKIT_LOCATION_REGIONNAME: regionName,
                                       abstractLocation.toolkit.TOOLKIT_TOOLKITNAME_FIELD: self.toolkitName,
                                       "xmin": extentList[0],
                                       "xmax": extentList[1],
                                       "ymin": extentList[2],
                                       "ymax": extentList[3]
                                       })

                if doc is None:
                    self.addCacheDocument(
                        type = self.doctype,
                        resource=outputFileName,
                        dataFormat=datatypes.IMAGE,
                        desc = additionalData
                    )

                else:
                    doc['resource'] = outputFileName
                    doc.desc = additionalData
                    doc.save()

        return data if doc is None else doc


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

    def __init__(self,dataLayer):

        self._datalayer = datalayer


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