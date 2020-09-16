import os
import logging
import numpy
from ...datalayer import project
from ...datalayer import datatypes

import matplotlib.pyplot as plt

from shapely.geometry import Point,box,MultiLineString, LineString

try:
    from freecad import app as FreeCAD
    import Part
    import Mesh
except ImportError as e:
    logging.warning("FreeCAD not Found, cannot convert to STL")



class locationImage(project.ProjectMultiDBPublic):
    """
        A class to handle an image that represents a location.

        Looks up the location in the public database in project 'imageLocation'.

    """

    def __init__(self, projectName, useAll=False):
        """
                Initialize

        Parameters
        -----------
        projectName: str
            The name of the project name in the local DB.

        useAll: bool
            whether to return union of all the image locations or just for one DB.
        """

        super.__init__(projectName=projectName,publicProjectName='imageLocation')


    @staticmethod
    def plot(self,data,ax=None,**query):
        """
            Plots the image from the document or from a document.

        Parameters
        -----------

        data: imageLocation doc, str
            Plot the image from the document or query
            the DB for the image with data name and plot the first.
            if more than 1 image exists, raise error.

        ax: matplotlib.Axes
            The axes to plot on, if None create
            a new axes.

        **query: mongoDB
            query map.

        Returns
        --------
            matplotlib.Axes

            The axes of the image
        """
        if ax is None:
            fig, ax = plt.subplots()


        if isinstance(data,str):
            doc = self.getMeasurementsDocuments(dataFormat='image',
                                                type='GIS',
                                                locationName=data,
                                                **query)

            if len(doc) > 1:
                raise ValueError('More than 1 documents fills those requirements')

            doc = doc[0]

        image = doc.getDocFromDB()
        extents = [doc.desc['xmin'], doc.desc['xmax'], doc.desc['ymin'], doc.desc['ymax']]
        ax = plt.imshow(image, extent=extents)


        return ax

    def loadImage(self, path, locationName, extents,sourceName):
        """
        Make region from the

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

        sourceName: str
            The source

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
        self.addMeasurementsDocument(**doc)


    def query(self,imageName=None,point=None,**query):
        """
                get the images.

        Parameters
        ----------

        imageName:  str
                image name.
        point: tuple
            a point inside the domain.
        query:
            querying the data.
        :return:

        """
        docsList = self.getMeasurementsDocuments(imeageName=imageName,**query)
        if point is not None:
            point = point if isinstance(point,Point) else Point(point[0],point[1])
            ret = []
            for doc in docsList:
                bb = box(doc.desc['xmin'],
                         doc.desc['xmax'],
                         doc.desc['ymin'],
                         doc.desc['ymax'])
                if point in bb:
                    ret.append(doc)
        else:
            ret =  docsList
        return ret

class topography(project.ProjectMultiDBPublic):
    """
        Holds a polygon with description. Allows querying on the location of the shapefile.

        The projectName in the public DB is 'locationGIS'



    """


    _dxdy = None

    @property
    def dxdy(self):
        return self._dxdy

    @dxdy.setter
    def dxdy(self, value):
        if value is not None:
            self._dxdy = value

    @def skipinterior(self):
        return self._skipinterior


    def __init__(self, projectName, useAll=False, dxdy=None):
        """
                Initialize

        Parameters
        -----------
        projectName: str
            The name of the project name in the local DB.

        useAll: bool
            whether to return union of all the image locations or just for one DB.
        """
        super().__init__(projectName=projectName,publicProjectName='locationGIS')

        self._dxdy = 10 if dxdy is None else dxdy  # m
        self._skipinterior = 100  # 100 # m, the meters to remove from the interpolation to make sure all exists.


    def _make_facet_str(self, n, v1, v2, v3):
        facet_str = 'facet normal ' + ' '.join(map(str, n)) + '\n'
        facet_str += '  outer loop\n'
        facet_str += '      vertex ' + ' '.join(map(str, v1)) + '\n'
        facet_str += '      vertex ' + ' '.join(map(str, v2)) + '\n'
        facet_str += '      vertex ' + ' '.join(map(str, v3)) + '\n'
        facet_str += '  endloop\n'
        facet_str += 'endfacet\n'
        return facet_str

    def _makestl(self, X, Y, elev, solidname):
        """
            Takes a mesh of x,y and elev and convert it to stl file.

            X - matrix of x coordinate. [[ like meshgrid ]]
            Y - matrix of y coordinate. [[ like meshgrid ]]
            elev - matrix of elevation.

        """
        base_elev = 0
        stl_str = 'solid ' + solidname + '\n'
        for i in range(elev.shape[0] - 1):
            print(i)
            for j in range(elev.shape[1] - 1):

                x = X[i, j];
                y = Y[i, j]
                v1 = [x, y, elev[i, j]]

                x = X[i + 1, j];
                y = Y[i, j]
                v2 = [x, y, elev[i + 1, j]]

                x = X[i, j];
                y = Y[i, j + 1]
                v3 = [x, y, elev[i, j + 1]]

                x = X[i + 1, j + 1];
                y = Y[i + 1, j + 1]
                v4 = [x, y, elev[i + 1, j + 1]]

                # dem facet 1
                n = numpy.cross(numpy.array(v1) - numpy.array(v2), numpy.array(v1) - numpy.array(v3))
                n = n / numpy.sqrt(sum(n ** 2))
                stl_str += self._make_facet_str(n, v1, v2, v3)

                # dem facet 2
                n = numpy.cross(numpy.array(v2) - numpy.array(v3), numpy.array(v2) - numpy.array(v4))
                n = n / numpy.sqrt(sum(n ** 2))
                # stl_str += self._make_facet_str( n, v2, v3, v4 )
                stl_str += self._make_facet_str(n, v2, v4, v3)

                # base facets
                v1b = list(v1)
                v2b = list(v2)
                v3b = list(v3)
                v4b = list(v4)

                v1b[-1] = base_elev
                v2b[-1] = base_elev
                v3b[-1] = base_elev
                v4b[-1] = base_elev

                n = [0.0, 0.0, -1.0]

                stl_str += self._make_facet_str(n, v1b, v2b, v3b)
                stl_str += self._make_facet_str(n, v2b, v3b, v4b)

                vlist = [v1, v2, v3, v4]
                vblist = [v1b, v2b, v3b, v4b]

                # Now the walls.
                for k, l in [(0, 1), (0, 2), (1, 3), (2, 3)]:
                    # check if v[i],v[j] are on boundaries.
                    kboundary = False
                    if vlist[k][0] == X.min() or vlist[k][0] == X.max():
                        kboundary = True

                    lboundary = False
                    if vlist[l][1] == Y.min() or vlist[l][1] == Y.max():
                        lboundary = True

                    if (kboundary or lboundary):
                        # Add i,j,j-base.
                        n = numpy.cross(numpy.array(vlist[k]) - numpy.array(vlist[l]), numpy.array(vblist[l]) - numpy.array(vlist[l]))
                        n = n / numpy.sqrt(sum(n ** 2))
                        stl_str += self._make_facet_str(n, vlist[k], vblist[l], vlist[l])

                        # add j-base,i-base,i
                        n = numpy.cross(numpy.array(vlist[k]) - numpy.array(vblist[k]), numpy.array(vlist[k]) - numpy.array(vblist[l]))
                        n = n / numpy.sqrt(sum(n ** 2))
                        stl_str += self._make_facet_str(n, vlist[k], vblist[k], vblist[l])

        stl_str += 'endsolid ' + solidname + '\n'
        return stl_str

    def toSTL(self, doc, solidname, flat=None):
        """
            Gets a shape file of topography.
            each contour line has property 'height'.
            Converts it to equigrid xy mesh and then build the STL.
        """

        # 1. read the shp file.
        gpandas = doc.getData()

        # 2. Convert contour map to regular height map.
        # 2.1 get boundaries
        xmin = gpandas['geometry'].bounds['minx'].min()
        xmax = gpandas['geometry'].bounds['maxx'].max()

        ymin = gpandas['geometry'].bounds['miny'].min()
        ymax = gpandas['geometry'].bounds['maxy'].max()

        print("Mesh boundaries x=(%s,%s) ; y=(%s,%s)" % (xmin, xmax, ymin, ymax))
        # 2.2 build the mesh.
        inter = self.skipinterior
        grid_x, grid_y = numpy.mgrid[(xmin + inter):(xmax - inter):self.dxdy,
                         (ymin + inter):(ymax - inter):self.dxdy]

        # 3. Get the points from the geom
        Height = []
        XY = []

        for i, line in enumerate(gpandas.iterrows()):
            if isinstance(line[1]['geometry'], LineString):
                linecoords = [x for x in line[1]['geometry'].coords]
                lineheight = [line[1]['HEIGHT']] * len(linecoords)
                XY += linecoords
                Height += lineheight
            else:
                for ll in line[1]['geometry']:
                    linecoords = [x for x in ll.coords]
                    lineheight = [line[1]['HEIGHT']] * len(linecoords)
                    XY += linecoords
                    Height += lineheight
        if flat is not None:
            for i in range(len(Height)):
                Height[i] = flat

        # adding fill values for places outside the map, e.g. inside the sea.
        grid_z2 = numpy.griddata(XY, Height, (grid_x, grid_y), method='cubic', fill_value=0.)
        # replace zero height with small random values so the stl file won't contain NaNs
        for i in range(grid_z2.shape[0]):
            for j in range(grid_z2.shape[1]):
                if (grid_z2[i, j] == 0.):
                    grid_z2[i, j] = numpy.random.random()

        if numpy.isnan(grid_z2).any():
            print("Found some NaN in cubic iterpolation. consider increasing the boundaries of the interior")

        stlstr = self._makestl(grid_x, grid_y, grid_z2, solidname)

        data = {"grid_x": grid_x, "grid_y": grid_y, "grid_z": grid_z2, "XY": XY, "Height": Height, "geo": gpandas}
        print("X min: %s , X max: %s " % (numpy.min(grid_x), numpy.min(grid_x)))
        print("Y min: %s , Y max: %s " % (numpy.min(grid_y), numpy.min(grid_y)))
        return stlstr, data

    def makeRegion(self, points, regionName, mode="Contour", desc=None):
        """
        Load a polygon document that holds the path of a GIS shapefile.

        Parameters:
        -----------
            points: list, dict
                Holds the ITM coordinates of a rectangle.
                If it is a list: [minimum x, maximum x, minimum y, maximum y]
                If it is a dic : dict(xmin=,xmax=,ymin=,ymax=).

            regionName: str
                Used as part of a new file's name.
            mode: str
                The data type of the desired data.
                Recieves any mode specified in the GISOrigin document.\n

            desc: dict
                A dictionary with any additional data for the location.

        Returns:
        --------
            The new document.

        """
        if useOwn:
            fullfilesdirect = self._projectMultiDB.getMeasurementsDocumentsAsDict(type="GISOrigin")["documents"][0]["desc"]["modes"]
            path = self._projectMultiDB.getMeasurementsDocumentsAsDict(type="GISOrigin")["documents"][0]["resource"]
            fullPath = "%s/%s" % (path, fullfilesdirect[mode])
        else:
            publicproject = hera.datalayer.project.ProjectMultiDB(projectName="PublicData", databaseNameList=["public"])
            fullPath = publicproject.getMeasurementsDocumentsAsDict(type="GIS",mode=mode)["documents"][0]["resource"]

        descDict = dict(CutName=regionName, points=points, mode=mode)

        if desc is not None:
            descDict.update(desc)

        documents = self.getMeasurementsDocumentsAsDict(points=points, mode=mode)
        if len(documents) == 0:

            FileName = "%s//%s%s-%s.shp" % (self._FilesDirectory, self._projectName, regionName, mode)

            os.system("ogr2ogr -clipsrc %s %s %s %s %s %s" % (points[0],points[1],points[2],points[3], FileName,fullPath))
            self.addMeasurementsDocument(desc=desc,
                                         type="GIS",
                                         resource = FileName,
                                         dataFormat = "geopandas")
        else:
            resource = documents["documents"][0]["resource"]
            self.addMeasurementsDocument(desc=dict(**desc),
                                         type="GIS",
                                         resource = resource,
                                         dataFormat = "geopandas")


    def query(self, imageName=None, point=None, **query):
        """
            query the existing topography.

        :param imageName:
        :param point:
        :param query:
        :return:
        """
        pass

class buildings(project.ProjectMultiDBPublic):
    """
        Holds the list of buildings.
    """

    @property
    def doctype(self):
        return 'BuildingSTL'


    def toSTL(self, doc, outputfile,flat=None):
        """
            Converts the document to the stl and saves it to the disk.
            Adds the stl file to the DB.


            Parameters
            ----------

            doc: hera.datalayer.document.MetadataFrame, hera.datalayer.
                The document with the data to convert.

            flat: None or float.
                The base of the building.
                If None, use the buildings height.

            outputfile: str
                a path to the output file.

            Returns
            -------
            float
            The string with the STL format.

        """
        self.logger.info(f"begin with doc={doc},outputfile={outputfile},flat={flat}")

        maxheight = -500

        FreeCADDOC = FreeCAD.newDocument("Unnamed")


        shp = doc.getData() if hasattr(doc,"getData") else doc

        k = -1
        for j in range(len(shp)):  # converting al the buildings
            try:
                walls = shp['geometry'][j].exterior.xy
            except:
                continue

            if j % 100 == 0:
                self.logger.execution(f"{j} shape file is executed. Length Shape: {len(shp)}")

            k = k + 1
            wallsheight = shp['BLDG_HT'][j]
            if flat is None:
                altitude = shp['HT_LAND'][j]
            else:
                altitude = flat

            FreeCADDOC.addObject('Sketcher::SketchObject', 'Sketch' + str(j))
            FreeCADDOC.Objects[2 * k].Placement = FreeCAD.Placement(FreeCAD.Vector(0.000000, 0.000000, 0.000000),  # 2*k-1
                                                             FreeCAD.Rotation(0.000000, 0.000000, 0.000000, 1.000000))

            for i in range(len(walls[0]) - 1):
                FreeCADDOC.Objects[2 * k].addGeometry(Part.Line(FreeCAD.Vector(walls[0][i], walls[1][i], altitude),
                                                         FreeCAD.Vector(walls[0][i + 1], walls[1][i + 1], altitude)))

            FreeCADDOC.addObject("PartDesign::Pad", "Pad" + str(j))
            FreeCADDOC.Objects[2 * k + 1].Sketch = FreeCADDOC.Objects[2 * k]
            buildingTopAltitude = wallsheight + altitude  # wallsheight + altitude
            maxheight = max(maxheight, buildingTopAltitude)
            FreeCADDOC.Objects[2 * k + 1].Length = buildingTopAltitude  # 35.000000
            FreeCADDOC.Objects[2 * k + 1].Reversed = 0
            FreeCADDOC.Objects[2 * k + 1].Midplane = 0
            FreeCADDOC.Objects[2 * k + 1].Length2 = wallsheight  # 100.000000
            FreeCADDOC.Objects[2 * k + 1].Type = 0
            FreeCADDOC.Objects[2 * k + 1].UpToFace = None


        FreeCADDOC.recompute() # maybe it should be in the loop.

        outputfileFull = os.path.abspath(outputfile)
        Mesh.export(FreeCADDOC.Objects, outputfileFull)

        self.addCacheDocument(type=self.doctype,
                              resource=outputfileFull,
                              dataFormat=datatypes.STRING)


        self.logger.info(f"toSTL: end. Maxheight {maxheight}")
        return maxheight





# class convert():
#
#     _projectName = None
#     _FilesDirectory = None
#     _projectMultiDB = None
#     _GISdatalayer = None
#     _manipulator = None
#
#     def __init__(self, projectName, FilesDirectory, users=[None], useAll=False):
#
#         self._FilesDirectory = FilesDirectory
#         self._projectName = projectName
#         self._GISdatalayer = GIS_datalayer(projectName=projectName, FilesDirectory=FilesDirectory, users=users, useAll=useAll)
#         self._projectMultiDB = hera.datalayer.project.ProjectMultiDB(projectName=projectName, databaseNameList=users, useAll=useAll)
#         self._manipulator = dataManipulations()
#
#     def addSTLtoDB(self, path, NewFileName, points, xMin, xMax, yMin, yMax, zMin, zMax, dxdy, **kwargs):
#         """
#         Adds a path to the dataframe under the type 'stlType'.
#
#         Parameters:
#
#             path: The path (string)
#             NewFileName: A name for the file (string)
#             kwargs: Additional parameters for the document.
#
#         Returns: A list that contains the stl string and a dict that holds information about it.
#         -------
#
#         """
#
#
#         self._projectMultiDB.addMeasurementsDocument(desc=dict(name = NewFileName, bounds = points, dxdy=dxdy,
#                                                                xMin=xMin, xMax=xMax, yMin=yMin, yMax=yMax, zMin=zMin, zMax=zMax, **kwargs),
#                                                       type="stlFile",
#                                                       resource=path,
#                                                       dataFormat="string")
#
#     def toSTL(self, data, NewFileName, dxdy=50, save=True, addtoDB=True, flat=None, path=None, **kwargs):
#         """
#         Converts a geopandas dataframe data to an stl file.
#
#         Parameters:
#         -----------
#
#             data: pandas.DataFrame
#                 The data that should be converted to stl. May be a dataframe or a name of a saved polygon in the database.
#             NewFileName: str
#                 A name for the new stl file, also used in the stl string. (string)
#             dxdy: float
#                 the dimention of each cell in the mesh in meters, the default is 50.
#             save: bool
#                 Default is True. If True, the new stl string is saved as a file and the path to the file is added to the database.
#             flat:
#                 Default is None. Else, it assumes that the area is flat and the value of flat is the height of the mesh cells.
#             path:
#                 Default is None. Then, the path in which the data is saved is the given self.FilesDirectory. Else, the path is path. (string)
#             kwargs:
#                 Any additional metadata to be added to the new document in the database.
#
#         Returns
#         -------
#
#         """
#
#         if type(data) == str:
#             polygon = self._GISdatalayer.getGeometry(data)
#             dataframe = self._GISdatalayer.getGISDocuments(geometry=data, geometryMode="contains")[0].getDocFromDB()
#             geodata = self._manipulator.PolygonDataFrameIntersection(polygon=polygon, dataframe=dataframe)
#         elif type(data) == geopandas.geodataframe.GeoDataFrame:
#             geodata = data
#         else:
#             raise KeyError("data should be geopandas dataframe or a polygon.")
#         xmin = geodata['geometry'].bounds['minx'].min()
#         xmax = geodata['geometry'].bounds['maxx'].max()
#
#         ymin = geodata['geometry'].bounds['miny'].min()
#         ymax = geodata['geometry'].bounds['maxy'].max()
#         points = [xmin, ymin, xmax, ymax]
#         if len(datalayer.Measurements.getDocuments(projectName=self._projectName, type="stlFile", bounds=points, dxdy=dxdy)) >0:
#             stlstr = datalayer.Measurements.getDocuments(projectName=self._projectName, type="stlFile", bounds=points, dxdy=dxdy)[0].getDocFromDB()
#             newdict = datalayer.Measurements.getDocuments(projectName=self._projectName, type="stlFile", bounds=points, dxdy=dxdy)[0].asDict()
#             newdata = pandas.DataFrame(dict(gridxMin=[newdict["desc"]["xMin"]], gridxMax=[newdict["desc"]["xMax"]],
#                                             gridyMin=[newdict["desc"]["yMin"]], gridyMax=[newdict["desc"]["yMax"]],
#                                             gridzMin=[newdict["desc"]["zMin"]], gridzMax=[newdict["desc"]["zMax"]]))
#         else:
#             stlstr, newdata = self.Convert_geopandas_to_stl(gpandas=geodata, points=points, flat=flat, NewFileName=NewFileName, dxdy=dxdy)
#
#         if save:
#             p = self._FilesDirectory if path is None else path
#             new_file_path = p + "/" + NewFileName + ".stl"
#             new_file = open(new_file_path, "w")
#             new_file.write(stlstr)
#             newdata = newdata.reset_index()
#             if addtoDB:
#                 self.addSTLtoDB(p, NewFileName, points=points, xMin=newdata["gridxMin"][0], xMax=newdata["gridxMax"][0],
#                                 yMin=newdata["gridyMin"][0], yMax=newdata["gridyMax"][0], zMin=newdata["gridzMin"][0], zMax=newdata["gridzMax"][0], dxdy=dxdy, **kwargs)
#
#         return stlstr, newdata
#
#     def _make_facet_str(self, n, v1, v2, v3):
#         facet_str = 'facet normal ' + ' '.join(map(str, n)) + '\n'
#         facet_str += '  outer loop\n'
#         facet_str += '      vertex ' + ' '.join(map(str, v1)) + '\n'
#         facet_str += '      vertex ' + ' '.join(map(str, v2)) + '\n'
#         facet_str += '      vertex ' + ' '.join(map(str, v3)) + '\n'
#         facet_str += '  endloop\n'
#         facet_str += 'endfacet\n'
#         return facet_str
#
#     def _makestl(self, X, Y, elev, NewFileName):
#         """
#             Takes a mesh of x,y and elev and convert it to stl file.
#
#             X - matrix of x coordinate. [[ like meshgrid ]]
#             Y - matrix of y coordinate. [[ like meshgrid ]]
#             elev - matrix of elevation.
#
#         """
#         base_elev = elev.min() - 10
#         stl_str = 'solid ' + NewFileName + '\n'
#         print(elev.shape[0] - 1)
#         for i in range(elev.shape[0] - 1):
#             print(i)
#             for j in range(elev.shape[1] - 1):
#
#                 x = X[i, j];
#                 y = Y[i, j]
#                 v1 = [x, y, elev[i, j]]
#
#                 x = X[i + 1, j];
#                 y = Y[i, j]
#                 v2 = [x, y, elev[i + 1, j]]
#
#                 x = X[i, j];
#                 y = Y[i, j + 1]
#                 v3 = [x, y, elev[i, j + 1]]
#
#                 x = X[i + 1, j + 1];
#                 y = Y[i + 1, j + 1]
#                 v4 = [x, y, elev[i + 1, j + 1]]
#
#                 # dem facet 1
#                 n = cross(array(v1) - array(v2), array(v1) - array(v3))
#                 n = n / sqrt(sum(n ** 2))
#                 stl_str += self._make_facet_str(n, v1, v2, v3)
#
#                 # dem facet 2
#                 n = cross(array(v2) - array(v3), array(v2) - array(v4))
#                 n = n / sqrt(sum(n ** 2))
#                 # stl_str += self._make_facet_str( n, v2, v3, v4 )
#                 stl_str += self._make_facet_str(n, v2, v4, v3)
#
#                 # base facets
#                 v1b = list(v1)
#                 v2b = list(v2)
#                 v3b = list(v3)
#                 v4b = list(v4)
#
#                 v1b[-1] = base_elev
#                 v2b[-1] = base_elev
#                 v3b[-1] = base_elev
#                 v4b[-1] = base_elev
#
#                 n = [0.0, 0.0, -1.0]
#
#                 stl_str += self._make_facet_str(n, v1b, v2b, v3b)
#                 stl_str += self._make_facet_str(n, v2b, v3b, v4b)
#
#                 vlist = [v1, v2, v3, v4]
#                 vblist = [v1b, v2b, v3b, v4b]
#
#                 # Now the walls.
#                 for k, l in [(0, 1), (0, 2), (1, 3), (2, 3)]:
#                     # check if v[i],v[j] are on boundaries.
#                     kboundary = False
#                     if vlist[k][0] == X.min() or vlist[k][0] == X.max():
#                         kboundary = True
#
#                     lboundary = False
#                     if vlist[l][1] == Y.min() or vlist[l][1] == Y.max():
#                         lboundary = True
#
#                     if (kboundary or lboundary):
#                         # Add i,j,j-base.
#                         n = cross(array(vlist[k]) - array(vlist[l]), array(vblist[l]) - array(vlist[l]))
#                         n = n / sqrt(sum(n ** 2))
#                         stl_str += self._make_facet_str(n, vlist[k], vblist[l], vlist[l])
#
#                         # add j-base,i-base,i
#                         n = cross(array(vlist[k]) - array(vblist[k]), array(vlist[k]) - array(vblist[l]))
#                         n = n / sqrt(sum(n ** 2))
#                         stl_str += self._make_facet_str(n, vlist[k], vblist[k], vblist[l])
#
#         stl_str += 'endsolid ' + NewFileName + '\n'
#         return stl_str
#
#     def Convert_geopandas_to_stl(self, gpandas, points, NewFileName, dxdy=50, flat=None):
#         """
#             Gets a shape file of topography.
#             each contour line has property 'height'.
#             Converts it to equigrid xy mesh and then build the STL.
#         """
#
#         # 1. Convert contour map to regular height map.
#         # 1.1 get boundaries
#         xmin = points[0]
#         xmax = points[2]
#
#         ymin = points[1]
#         ymax = points[3]
#
#         print("Mesh boundaries x=(%s,%s) ; y=(%s,%s)" % (xmin, xmax, ymin, ymax))
#         # 1.2 build the mesh.
#         grid_x, grid_y = numpy.mgrid[(xmin):(xmax):dxdy, (ymin):(ymax):dxdy]
#         # 2. Get the points from the geom
#         Height = []
#         XY = []
#         for i, line in enumerate(gpandas.iterrows()):
#             if isinstance(line[1]['geometry'], LineString):
#                 linecoords = [x for x in line[1]['geometry'].coords]
#                 lineheight = [line[1]['HEIGHT']] * len(linecoords)
#                 XY += linecoords
#                 Height += lineheight
#             else:
#                 for ll in line[1]['geometry']:
#                     linecoords = [x for x in ll.coords]
#                     lineheight = [line[1]['HEIGHT']] * len(linecoords)
#                     XY += linecoords
#                     Height += lineheight
#         if flat is not None:
#             for i in range(len(Height)):
#                 Height[i] = flat
#         grid_z2 = griddata(XY, Height, (grid_x, grid_y), method='cubic')
#         grid_z2 = self.organizeGrid(grid_z2)
#
#         stlstr = self._makestl(grid_x, grid_y, grid_z2, NewFileName)
#
#         data = pandas.DataFrame({"XY": XY, "Height": Height, "gridxMin":grid_x.min(), "gridxMax":grid_x.max(),
#                                  "gridyMin":grid_y.min(), "gridyMax":grid_y.max(), "gridzMin":grid_z2.min(), "gridzMax":grid_z2.max(),})
#
#         return stlstr, data
#
#     def organizeGrid(self, grid):
#
#         for row in grid:
#             for i in range(len(row)):
#                 if math.isnan(row[i]):
#                     pass
#                 else:
#                     break
#             for n in range(i):
#                 row[n] = row[i]
#             for i in reversed(range(len(row))):
#                 if math.isnan(row[i]):
#                     pass
#                 else:
#                     break
#             for n in range(len(row)-i):
#                 row[-n-1] = row[i]
#         return grid