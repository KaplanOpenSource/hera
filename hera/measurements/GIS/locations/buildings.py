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
        super().__init__(projectName=projectName,publicProjectName=self.publicProjectName,FilesDirectory=FilesDirectory,useAll=useAll,Source=Source)
        self._analysis = analysis(projectName=projectName, dataLayer=self)
        self.setConfig()

    def setConfig(self, Source="BNTL", dbName=None, **kwargs):
        config = dict(source=Source,**kwargs)
        super().setConfig(config, dbName=dbName)

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


class analysis():

    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, projectName, dataLayer=None, FilesDirectory="", databaseNameList=None, useAll=False,
                 publicProjectName="Buildings", Source="BNTL"):

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
