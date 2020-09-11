from ... import datalayer
import getpass
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from shapely.geometry import Point,box

class locationImage(datalayer.ProjectMultiDBPublic):
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


class topography(datalayer.ProjectMultiDBPublic):
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
            publicproject = datalayer.ProjectMultiDB(projectName="PublicData", databaseNameList=["public"])
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





    def toSTL(self,data,**query):
        """
            convert the data to STL.

            data can be either doc or str (the name of the location).
            query - a mongoDB query.

        :return:
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

class buildings(datalayer.ProjectMultiDBPublic):
    """
        Holds the list of buildings.
    """
    pass