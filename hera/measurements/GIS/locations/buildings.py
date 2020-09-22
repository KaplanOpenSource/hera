import os
import logging
import pandas
import geopandas
from .abstractLocation import datalayer as locationDatalayer
from ....datalayer import datatypes

import matplotlib.pyplot as plt

from shapely.geometry import Point,box,MultiLineString, LineString

try:
    from freecad import app as FreeCAD
    import Part
    import Mesh
except ImportError as e:
    logging.warning("FreeCAD not Found, cannot convert to STL")


class datalayer(locationDatalayer):

    _publicProjectName = None
    _analysis = None

    @property
    def doctype(self):
        return 'BuildingSTL'

    @property
    def analysis(self):
        return self._analysis

    def __init__(self, projectName, FilesDirectory="", databaseNameList=None, useAll=False,publicProjectName="Buildings",Source="BNTL"):

        self._publicProjectName = publicProjectName
        super().__init__(projectName=projectName,publicProjectName=self.publicProjectName,FilesDirectory=FilesDirectory,databaseNameList=databaseNameList,useAll=useAll,Source=Source)
        self._analysis = analysis(projectName=projectName, dataLayer=self)

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
        self.logger.info("toSTL: begin")

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


class analysis():

    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, projectName, dataLayer=None, FilesDirectory="", databaseNameList=None, useAll=False,
                 publicProjectName="Topography", Source="BNTL"):

        self._datalayer = datalayer(projectName=projectName, FilesDirectory=FilesDirectory, publicProjectName=publicProjectName,
                         databaseNameList=databaseNameList, useAll=useAll, Source=Source) if datalayer is None else dataLayer

    def ConvexPolygons(self, data, buffer=100):
        """
        Returns polygons of groups of buildings.
        """
        data = data.reset_index()
        d = data.buffer(buffer)
        indicelist = [[0]]
        for i in range(1, len(data)):
            found = False
            for g in range(len(indicelist)):
                for n in indicelist[g]:
                    if d[i].intersection(d[n]).is_empty:
                        continue
                    else:
                        indicelist[g].append(i)
                        found = True
                        break
                if found:
                    break
                if g == len(indicelist) - 1:
                    indicelist.append([i])

        geo = data.loc[indicelist[0]].unary_union.convex_hull
        gpd = geopandas.GeoDataFrame.from_dict([{"geometry": geo, "area": geo.area}])
        for indice in indicelist[1:]:
            geo = data.loc[indice].unary_union.convex_hull
            gpd = pandas.concat([gpd, geopandas.GeoDataFrame.from_dict([{"geometry": geo, "area": geo.area}])])

        gpd = gpd.sort_values(by="area", ascending=False).reset_index()
        found = False
        for i in range(len(gpd)):
            for j in range(i + 1, len(gpd)):
                if gpd.loc[i].geometry.intersection(gpd.loc[j].geometry).is_empty:
                    continue
                else:
                    found = True
                    break
            if found:
                break
        if found:
            gpd = self.ConvexPolygons(gpd, buffer=1)

        return gpd

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