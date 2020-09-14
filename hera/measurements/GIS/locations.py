import logging

import hera.datalayer.project
from ... import datalayer
import getpass
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from shapely.geometry import Point,box

try:
    from freecad import app as FreeCAD
except ImportError as e:
    logging.warning("FreeCAD not Found, cannot convert to STL")



class locationImage(hera.datalayer.project.ProjectMultiDBPublic):
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

    def load(self, path, locationName, extents):
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


class topography(hera.datalayer.project.ProjectMultiDBPublic):
    """
        Holds a polygon with description. Allows querying on the location of the shapefile.

        The projectName in the public DB is 'locationGIS'

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
        super.__init__(projectName=projectName,publicProjectName='locationGIS')


    def load(self, points, CutName, mode="Contour", desc=None, useOwn=False):
        """
        Load a polygon document that holds the path of a GIS shapefile.

        Parameters:
        -----------
            points: list, dict
                Holds the ITM coordinates of a rectangle.
                If it is a list: [minimum x, maximum x, minimum y, maximum y]
                If it is a dic : dict(xmin=,xmax=,ymin=,ymax=).

            CutName: str
                Used as part of a new file's name. (string)\n
            mode:
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

        descDict = dict(CutName=CutName,points=points,mode=mode)

        if desc is not None:
            descDict.update(desc)

        documents = self.getMeasurementsDocumentsAsDict(points=points, mode=mode)
        if len(documents) == 0:

            FileName = "%s//%s%s-%s.shp" % (self._FilesDirectory, self._projectName, CutName, mode)

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





    def toSTL(self,doc):
        """

        Convert the data to STL.

        Parameters
        ----------

        doc: datalayer document or hera.datalayer.nonDBMetadata


        Returns
        -------
        str, the STL string file.
        """
        pass


    def query(self, imageName=None, point=None, **query):
        """
            query the existing topography.

        :param imageName:
        :param point:
        :param query:
        :return:
        """
        pass

class buildings(hera.datalayer.project.ProjectMultiDBPublic):
    """
        Holds the list of buildings.
    """



    def toSTL(self, doc11, flat=None):
        """
            Converts the document to the stl.

            Parameters
            ----------

            doc: hera.datalayer.document.MetadataFrame, hera.datalayer.
                The document with the data to convert.


            Returns
            -------
            str
            The string with the STL format.

        """

        maxheight = -500

        FreeCADDOC = FreeCAD.newDocument("Unnamed")
        shp =
        k = -1
        for j in range(len(shp)):  # converting al the buildings
            try:
                walls = shp['geometry'][j].exterior.xy
            except:
                continue
            print('bad j !!!', j)
        if j % 100 == 0:
            print('100j', j, k, len(shp))  # just to see that the function is still working

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
        FreeCADDOC.recompute()  # maybe it can go outside the for loop


        FreeCADDOC.recompute()
        Mesh.export(FreeCADDOC.Objects, "file-buildings.stl")

        return maxheight
