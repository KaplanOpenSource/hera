from shapely import geometry
import os
import io
from hera import toolkit
import geopandas
from hera.measurements.GIS.shapes import ShapesToolKit
import pandas
import geojson



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
        #self.logger.info(f"Toolkit {toolkitName} in project {projectName} - Initializing")
        super().__init__(projectName=projectName,toolkitName=toolkitName,FilesDirectory=FilesDirectory)

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

        The region can be created from a file on the disk or from a region that was stored to the DB.

        If stored to DB, they are stored as datasources and therefore can be accessed through
        the datasource retrieval functions of the toolkit class.

        Parameters:
            points: list or dict, or geopandas.
                Holds the ITM coordinates of a rectangle.
                 list  - [minimum x, minimum y, maximum x, maximum y]
                 dict  - {minX : , minY : maxX : , maxY : }

                 if a geopandas, use the bounds of the polygon.

            saveMode: str (constants)
                The mode of the creation of the region.
                can be one of the following:

                - toolkit.TOOLKIT_SAVEMODE_NOSAVE               : dont save.
                - toolkit.TOOLKIT_SAVEMODE_FILE                 : Write the region to the file. Raise exception if file exists.
                - toolkit.TOOLKIT_SAVEMODE_ONLYFILE_REPLACE     : Write the region to the file. Overwrite file  exists.
                - toolkit.TOOLKIT_SAVEMODE_FILEANDDB            : Write the region to the file, add to the DB as dataresource of that
                                                                  toolkit. If the db record exists, raise exception.
                - toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE    : Write the region to the file, add to the DB as dataresource of that
                                                                  toolkit. If the db record exists, raise exception.

                If true, then update the DB the region exists

            dataSourceOrFile: str
                    The name of the resource, or the path to the  dataset file that is being cropped.

            dataSourceVersion: tuple of ints
                    The version of the dataset. If None, take the latest version.

            regionName: str
                The new region name

            additional_data: dict
                A dictionary with any additional metadata parameters for the region.

        """
        doc = None
        #self.logger.info(f"Project {self.projectName}, Toolkit {self.toolkitName}: Start")
        if dataSourceOrFile is not None and os.path.exists(dataSourceOrFile):
            self.logger.execution(f"Using existing file {dataSourceOrFile}")
            inputData = dataSourceOrFile
        else:
            self.logger.execution(f"Using DB datasource {dataSourceOrFile}")
            resourceDoc = self.getDatasourceDocument(dataSourceOrFile,version=dataSourceVersion)
            inputData    = resourceDoc.resource

        outputFileName = os.path.join(self.FilesDirectory, regionName)
        if saveMode in [toolkit.TOOLKIT_SAVEMODE_ONLYFILE,toolkit.TOOLKIT_SAVEMODE_FILEANDDB]:
            if os.path.exists(outputFileName):
                raise ValueError(f"The outputfile {outputFileName} exists. Either remove it or run a saveMode that ends with _REPLACE")

            # Check if the points in the DB are the same as the input. if not and the saveMode is toolkit.TOOLKIT_SAVEMODE_FILEANDDB,
            # then raise exception.

        if saveMode in [toolkit.TOOLKIT_SAVEMODE_FILEANDDB,toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE]:
            doc  = self.getDatasourceDocument(datasourceName=regionName,**additional_data)
            if (doc is not None) and (saveMode==toolkit.TOOLKIT_SAVEMODE_FILEANDDB):
                raise ValueError(f"{regionName} with parameters {additional_data} already exists for project {self.projectName} in toolkit {self.toolkitName}")
        getattr(self,f"makeRegion_{inputData.split('.')[-1]}")(points=points,inputData=inputData,outputFileName=outputFileName)
        if saveMode in [toolkit.TOOLKIT_SAVEMODE_FILEANDDB,toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE]:
            if doc is None:
                # Adding to the DB.
                self.logger.execution(f"Adding {regionName} to the DB in a new record")
                additional_data[toolkit.TOOLKIT_DATASOURCE_NAME] = regionName
                additional_data[TOOLKIT_LOCATION_REGIONNAME] = regionName
                additional_data[TOOLKIT_LOCATION_POINTS] = points
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

    def makeRegion_shp(self,points,inputData,outputFileName):

        if isinstance(points,list):
            points = dict(minX= points[0], minY=points[1], maxX=points[2] , maxY=points[3])
        elif isinstance(points,dict):
            for field in ["minX","minY","maxX","maxY"]:
                if field not in points:
                    raise ValueError(f"points dict must have the following fields: minX,minY,maxX,maxY")
        elif isinstance(points,geopandas.geodataframe.GeoDataFrame):
            bounds = points.unary_union.bounds
            points = dict(minX=bounds[0], minY=bounds[1], maxX=bounds[2], maxY=bounds[3])
        try:
            cmd = f"ogr2ogr -clipsrc {points['minX']} {points['minY']} {points['maxX']} {points['maxY']} {outputFileName} {inputData}"
            self.logger.execution(f"Clipping the file: {cmd}")
            os.system(cmd)
        except Exception as e:
            self.logger.error(f"General error cutting the file {e}.")
            raise RuntimeError(f"Error clipping the file {inputData}. Exception {e}")

    def makeRegionByShapeName(self, shapeNameOrArea, saveMode=toolkit.TOOLKIT_SAVEMODE_FILEANDDB,regionName=None, dataSourceOrFile=None, dataSourceVersion=None, additional_data=dict()):
        """

        Creates the area from the shape given as input.
        The shape can be either the shape name in the DB, a geoJSON
        or a geopandas.

        If stored to DB, they are stored as datasources and therefore can be accessed through
        the datasource retrieval functions of the toolkit class.


        Parameters
        ----------
        shapeNameOrArea: str or geopandas.geoDataFrame
            Either shape name (in the DB), a geoJSON str or the geoDataframe.

        regionName: str
            The new region name

        saveMode: str (constants)
                The mode of the creation of the region.
                TOOLKIT_SAVEMODE_FILEANDDB is default.

                can be one of the following:

                - toolkit.TOOLKIT_SAVEMODE_NOSAVE               : dont save.
                - toolkit.TOOLKIT_SAVEMODE_FILE                 : Write the region to the file. Raise exception if file exists.
                - toolkit.TOOLKIT_SAVEMODE_ONLYFILE_REPLACE     : Write the region to the file. Overwrite file  exists.
                - toolkit.TOOLKIT_SAVEMODE_FILEANDDB            : Write the region to the file, add to the DB as dataresource of that
                                                                  toolkit. If the db record exists, raise exception.
                - toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE    : Write the region to the file, add to the DB as dataresource of that
                                                                  toolkit. If the db record exists, raise exception.


        dataSourceOrFile: str
                The name of the resource, or the path to the  dataset file that is being cropped.

        dataSourceVersion: tuple of ints
                The version of the dataset. If None, take the latest version.

        additional_data: dict
                If saved to DB, use this as additional fields to the meta-data


        :return:
        """
        if isinstance(shapeNameOrArea, str):
            points = ShapesToolKit(projectName=self.projectName).getShape(shapeNameOrArea)
            if points is None:
                points = geopandas.GeoDataFrame.from_features(pandas.read_json(shapeNameOrArea)["features"])
            else:
                points = points.getData()

        return self.makeRegion(points=points,
                                regionName=regionName,
                                saveMode=saveMode,
                                dataSourceOrFile=dataSourceOrFile,
                                dataSourceVersion=dataSourceVersion,
                                additional_data=additional_data)

    def getRegionDocumentByPoints(self,point):
        """

        Parameters
        -----------
        point: tuple
            The point to search

        Returns
        -------
            list with documents that the contain the point
        """
        if isinstance(point, list) or isinstance(point, tuple):
            point = geometry.Point(point)
        elif isinstance(point, geometry.Point):
            pass
        else:
            raise TypeError("point should be list, tuple or shapely.geometry.Point")
        docList = self.getDatasourceDocumentsList()

        ret = []
        for doc in docList:
            points = doc.desc[TOOLKIT_LOCATION_POINTS]
            pol = geometry.Polygon([[points["minX"],points["minY"]],[points["minX"],points["maxY"]],[points["maxX"],points["maxY"]],
                                    [points["maxX"],points["minY"]],[points["minX"],points["minY"]]])
            if pol.contains(point):
                ret.append(doc)

        return ret

    def getRegionDocumentByName(self,regionName):
        """

        Parameters
        ---------

            regionName: str
                    The name of the region
        Returns
        -------
            Return the vectorData with the region anme
        """

        return self.getDatasourceDocument(regionName)

    def getRegionByName(self,regionName):
        """

        Parameters
        ---------

            regionName: str
                    The name of the region
        Returns
        -------
            Return the vectorData with the region anme
        """

        return self.getDatasourceData(regionName)


    def getRegionNameList(self):
        """
            Return the list of region names.

        :return:
        """
        return [doc.desc[TOOLKIT_LOCATION_REGIONNAME] for doc in self.getDatasourceDocumentsList()]



