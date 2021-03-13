from shapely import geometry
import os
from .... import toolkit


TOOLKIT_LOCATION_REGIONNAME = "regionName"
TOOLKIT_LOCATION_POINTS     = "points"

class AbstractLocationToolkit(toolkit.abstractToolkit):
    """
        Absract locationToolkit

        This is the father of all the location toolkits. A location toolkit is a
        geographical part (usually a rectangle) of topography, buildings or other shapes.

        When these shapes are given in national shapefile, the abstractToolkit provides tools
        to crop a region, give it a name and load it to the DB.

    """

    _FilesDirectory = None

    @property
    def FilesDirectory(self):
        return self._FilesDirectory


    def __init__(self, projectName, toolkitName, FilesDirectory=None):
        """
            Initializes an abstract location toolkit.

        Parameters
        ----------
        projectName: str
            The project Name that the toolkit is initialized on
        toolkitName: str
            the specific toolkit, getting from the child.

        FilesDirectory: str or None
                The path to save a regions files when they are created.

                if str then represents a path (relative or absolute) to save the files in. The directory is created automatically.

                if None, then tries to get the default path of the project from the config. if it does not
                exist, then use the current directory.

        dataSourceOrData: str or None
                The name of the datasource to use.
                This datasource will be the main region that other subregions will be created from.

        dataSourceVersion: tuple (of int)
                The version of the source.
                if None, and dataSourceOrData is a datasource name (i.e a str) then use the latest version.
        """
        self.logger.info(f"Toolkit {toolkitName} in project {projectName} - Initializing")
        super().__init__(projectName=projectName,toolkitName=toolkitName)

        if FilesDirectory is None:
            self.logger.execution("Directory is not given, tries to load from default or using the current directory")
            self._FilesDirectory = self.getConfig().get("filesDirectory",os.getcwd())
            self.logger.execution(f"Using {self._FilesDirectory}")
        else:
            self.logger.execution(f"Using {os.path.abspath(FilesDirectory)}. Creating if does not exist")
            os.system("mkdir -p %s" % os.path.abspath(FilesDirectory))
            self._FilesDirectory = FilesDirectory

        self.logger.info(f"Toolkit {toolkitName} in project {projectName} - Done")


    def setDefaultDirectory(self,fileDirectory,create=True):
        """
            Set the default directory for the project.

        Parameters
        ----------
        fileDirectory: str
                The path to save the regions in.
                The directory is created if create flag is true (and directory does not exist).

        create: bool
            If false and directory does not exist, raise a NotADirectoryError exception.

        Returns
        -------
            str, the path.
        """
        fllpath = os.path.abspath(fileDirectory)

        if not os.path.exists(fllpath):
            if create:
                self.logger.execution(f"Directory {fllpath} does not exist. create")
                os.system(f"mkdir -p {fllpath}")
            else:
                raise NotADirectoryError(f"{fllpath} is not a directory, and create is False.")

        self.setConfig("filesDirectory",fllpath)


    def makeRegion(self, points, regionName, saveMode=toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE, dataSourceOrFile=None, dataSourceVersion=None, additional_data=dict()):
        """
        Create a new region from the larger data source.

        The region can be created from a file on the disk or from a region that was stored to the DB

        Parameters:
            points: list or dict
                Holds the ITM coordinates of a rectangle.
                 list  - [minimum x, minimum y, maximum x, maximum y]
                 dict  - {minX : , minY : maxX : , maxY : }

            saveMode: str (constants
                The mode of the creation of the region.
                can be one of the following:

                - toolkit.TOOLKIT_SAVEMODE_FILE                 : Write the region to the file. Raise exception if file exists.
                - toolkit.TOOLKIT_SAVEMODE_ONLYFILE_REPLACE     : Write the region to the file. Overwrite file  exists.
                - toolkit.TOOLKIT_SAVEMODE_FILEANDDB            : Write the region to the file, add to the DB as dataresource of that
                                                                  toolkit. If the db record exists, raise exception.
                - toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE    : Write the region to the file, add to the DB as dataresource of that
                                                                  toolkit. If the db record exists, raise exception.

                If true, then update the DB the region exists

            dataSourceOrFile: str
                    The name of the resource, or the path to the  dataset file that is being cropped.

            regionName: str
                The new region name

            additional_data: dict
                A dictionary with any additional metadata parameters for the region.

        """
        doc = None
        self.logger.info(f"Project {self.projectName}, Toolkit {self.toolkitName}: Start")
        if os.path.exists(dataSourceOrFile):
            self.logger.execution(f"Using existing file {dataSourceOrFile}")
            inputData = dataSourceOrFile
        else:
            self.logger.execution(f"Using DB datasource {dataSourceOrFile}")
            resourceDoc = self.getDatasourceDocument(dataSourceOrFile,version=dataSourceVersion)
            inputData    = resourceDoc.resource

        outputFileName = os.path.join(self.FilesDirectory, regionName)

        if saveMode in [toolkit.TOOLKIT_SAVEMODE_FILE,toolkit.TOOLKIT_SAVEMODE_FILEANDDB]:
            if os.path.exists(outputFileName):
                raise ValueError(f"The outputfile {outputFileName} exists. Either remove it or run a saveMode that ends with _REPLACE")

            # Check if the points in the DB are the same as the input. if not and the saveMode is toolkit.TOOLKIT_SAVEMODE_FILEANDDB,
            # then raise exception.

        if saveMode in [toolkit.TOOLKIT_SAVEMODE_FILEANDDB,toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE]:
            doc  = self.getDatasourceDocument(datasourceName=regionName,**additional_data)
            if (doc is not None) and (saveMode==toolkit.TOOLKIT_SAVEMODE_FILEANDDB):
                raise ValueError(f"{regionName} with parameters {additional_data} already exists for project {self.projectName} in toolkit {self.toolkitName}")

        if isinstance(points,list):
            points = dict(minX= points[0], minY=points[1], maxX=points[2] , maxY=points[3])
        elif isinstance(points,dict):
            for field in ["minX","minY","maxX","maxY"]:
                if field not in points:
                    raise ValueError(f"points dict must have the following fields: minX,minY,maxX,maxY")


        try:
            cmd = f"ogr2ogr -clipsrc {points['minX']} {points['minY']} {points['maxX']} {points['maxY']} {outputFileName} {inputData}"
            self.logger.execution(f"Clipping the file: {cmd}")
            os.system(cmd)
        except as e:
            self.logger.error(f"General error cutting the file {e}.")
            raise RuntimeError(f"Error clipping the file {inputData}. Exception {e}")


        if saveMode in [toolkit.TOOLKIT_SAVEMODE_FILEANDDB,toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE]:
            if doc is None:
                # Adding to the DB.
                self.logger.execution(f"Adding {regionName} to the DB in a new record")
                additional_data[toolkit.TOOLKIT_DATASOURCE_NAME] = regionName
                additional_data[toolkit.TOOLKIT_LOCATION_REGIONNAME] = regionName
                additional_data[toolkit.TOOLKIT_LOCATION_POINTS] = points
                additional_data[toolkit.TOOLKIT_TOOLKITNAME_FIELD] = self.toolkitName

                self.addDataSource(dataSourceName=regionName,
                                   resource=outputFileName,
                                   dataFormat=toolkit.datatypes.GEOPANDAS,
                                   version=None,
                                   **additional_data)

            else:
                # Updating DB. the case of existing record and no update was resolved above.
                self.logger.execution(f"Updating the record of {regionName} in the DB")
                docParams = {TOOLKIT_LOCATION_POINTS : points}
                docParams.update(additional_data)

                doc.desc=docParams
                doc.save()

        return outputFileName if doc is None else doc



    def getRegions(self,**kwargs):
        """
            Return all the regions associated with this toolkit.

        Returns
        -------
            List of docs.
        """
        queryDict = {"type":toolkit.TOOLKIT_DATASOURCE_TYPE,
                     toolkit.TOOLKIT_TOOLKITNAME_FIELD : self.toolkitName}

        queryDict.update(**kwargs)
        return self.getMeasurementsDocuments(**queryDict)


    def getLocationByPoints(self,point):
        """

        Parameters
        -----------
        point: tuple
            The point to search

        Returns
        -------
            list with documents that the contain the point
        """
        docList = self.getRegions()

        ret = []
        for doc in docList:
            pol = geometry.Polygon(*doc.desc[TOOLKIT_LOCATION_POINTS])
            if pol.contains(point):
                ret.append(doc)

        return ret


    def getLocationByRegion(self,regionName):
        """

        Parameters
        ---------

            regionName: str
                    The name of the region
        Returns
        -------
            Return the locations with the region anme
        """
        regionQry = { toolkit.TOOLKIT_LOCATION_REGIONNAME : regionName}
        return self.getRegions(**regionQry)



    def getRegionNameList(self):
        return [doc.desc[toolkit.TOOLKIT_LOCATION_REGIONNAME] for doc in self.getRegions()]







    #
    # def getDocuments(self, points=None, CutName=None, ShapeMode="contains", Shape=None, Source=None, **kwargs):
    #     """
    #     This function is used to load GIS data.
    #     One may use it to get all data that corresponds to any parameters listed in a document,
    #     or to add a new document that relates to a file that holds GIS data in an area defined by a rectangle.
    #     Can also be used to perform geometrical queries.
    #
    #     parameters:
    #         points: optional, for adding new data. Holds the ITM coordinates of a rectangle. It is a list, from the structure [minimum x, minimum y, maximum x, maximum y]\n
    #         CutName: optional, for adding new data. Used as part of a new file's name. (string)\n
    #         mode: The data type of the desired data. Recieves "Contour", "Buildings" or "Roads".\n
    #         ShapeMode: The mode of a geomtrical queries. Recieves "contains" or "intersects".\n
    #         Shape: A shapely geometry or a string with the name of a saved shapely geometry. Used to perform geometrical queries.\n
    #         **kwargs: any additional parameters that describe the data.
    #         return: The data.
    #     """
    #
    #     if Shape is not None:
    #         if type(Shape)==str:
    #             try:
    #                 Shape = shapeDatalayer(projectName=self._projectName, databaseNameList=self._databaseNameList, useAll=self._useAll).getShape(Shape)
    #             except IndexError:
    #                 raise IndexError("Shape isn't defined.")
    #         containPoints = []
    #         points = self.getFilesPointList(**kwargs)
    #         polygons = self.getFilesPolygonList(**kwargs)
    #         for i in range(len(points)):
    #             if ShapeMode == "contains":
    #                 if polygons[i].contains(Shape):
    #                     containPoints.append(points[i])
    #             elif ShapeMode == "intersects":
    #                 if polygons[i].intersects(Shape):
    #                     containPoints.append(points[i])
    #             else:
    #                 raise KeyError("ShapeMode incorrectly called. Choose 'contains' or 'intersects'.")
    #         if 1 == len(containPoints):
    #             data = self.getMeasurementsDocuments(points=containPoints[0], type=self._publicProjectName, **kwargs)
    #         else:
    #             data = []
    #             for p in containPoints:
    #                 data += self.getMeasurementsDocuments(points=p, type=self._publicProjectName, **kwargs)
    #
    #     else:
    #         if points==None and CutName==None:
    #             check = self.check_data(type=self._publicProjectName,**kwargs)
    #             if check:
    #                 data = self.getMeasurementsDocuments(type=self._publicProjectName,**kwargs)
    #         elif points==None and CutName!=None:
    #             check = self.check_data(CutName=CutName,type=self._publicProjectName, **kwargs)
    #             if check:
    #                 data = self.getMeasurementsDocuments(CutName=CutName,type=self._publicProjectName, **kwargs)
    #         elif CutName==None and points!=None:
    #             check = self.check_data(points=points,type=self._publicProjectName, **kwargs)
    #             if check:
    #                 data = self.getMeasurementsDocuments(CutName=CutName, type=self._publicProjectName,**kwargs)
    #         else:
    #             check = self.check_data(points=points, CutName=CutName,type=self._publicProjectName, **kwargs)
    #             if check:
    #                 data = self.getMeasurementsDocuments(points=points, CutName=CutName,type=self._publicProjectName, **kwargs)
    #         if not check:
    #             if points == None or CutName == None:
    #                 raise KeyError("Could not find data. Please insert points and CutName for making new data.")
    #             else:
    #                 self.makeRegion(points=points, regionName=CutName, Source=Source, additional_data=kwargs)
    #                 data = self.getMeasurementsDocuments(points=points, CutName=CutName)
    #
    #     return data

    # def check_data(self, **kwargs):
    #     """
    #     Checks whether there is a document that fulfills desired requirements.
    #     Parameters:
    #         kwargs: Any desired requirements.
    #
    #     Returns: True if there is a data that fulfills the requirement. False if there isn't.
    #
    #     """
    #
    #     check = self.getMeasurementsDocuments(**kwargs)
    #
    #     if 0 == len(check):
    #         result = False
    #     else:
    #         result = True
    #
    #     return result
    #
    # def getFilesPointList(self, **kwargs):
    #     """
    #     Returns a list of all the coordinates used for defining areas in existing documents.
    #
    #     Parameters:
    #         kwargs: kwargs: Any desired requirements for the documents.
    #
    #     Returns: List of lists, each contains 4 coordinates - [minimum x, minimum y, maximum x, maximum y].
    #
    #     """
    #     if self._databaseNameList[0] == "public" or self._databaseNameList[0] == "Public" and len(
    #             self._databaseNameList) > 1:
    #         userName = self._databaseNameList[1]
    #     else:
    #         userName = self._databaseNameList[0]
    #     documents = self.getMeasurementsDocumentsAsDict(**kwargs, users=[userName])["documents"]
    #     points = []
    #     for document in documents:
    #         if "points" in document["desc"].keys():
    #             if document["desc"]["points"] not in points:
    #                 points.append(document["desc"]["points"])
    #     return points
    #
    # def getFilesPolygonList(self, **kwargs):
    #     """
    #     Returns a list of all the polygons used for defining areas in existing documents.
    #
    #     Parameters:
    #     -----------
    #         kwargs: kwargs: Any desired requirements for the documents.
    #
    #     Returns:
    #     --------
    #         List of polygons.
    #     """
    #
    #     points = self.getFilesPointList(**kwargs)
    #     polygons = []
    #     for coordinates in points:
    #         polygons.append(geometry.Polygon([[coordinates[0], coordinates[1]], [coordinates[2], coordinates[1]],
    #                                           [coordinates[2], coordinates[3]], [coordinates[0], coordinates[3]]]))
    #     return polygons