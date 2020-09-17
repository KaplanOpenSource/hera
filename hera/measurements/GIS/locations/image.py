import os
import logging
import numpy
from .abstractLocation import datalayer as locationDatalayer
from ....datalayer import datatypes
from ....datalayer import project
import matplotlib.pyplot as plt

from shapely.geometry import Point,box,MultiLineString, LineString

class datalayer(project.ProjectMultiDBPublic):
    """
        A class to handle an image that represents a location.

        Looks up the location in the public database in project 'imageLocation'.

    """

    _projectName = None

    def __init__(self, projectName, databaseNameList=None, useAll=False,publicProjectName="Images"):

        self._projectName = projectName
        super().__init__(projectName=projectName, publicProjectName=publicProjectName,databaseNameList=databaseNameList,useAll=useAll)

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

    def load(self, path, imageName, extents):
        """
        Parameters:
        -----------
        projectName: str
                    The project name
        path:  str
                    The image path
        imageName: str
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
                   desc=dict(imageName=imageName,
                             xmin=extentList[0],
                             xmax=extentList[1],
                             ymin=extentList[2],
                             ymax=extentList[3]
                             )
                   )
        if self._databaseNameList[0] == "public" or self._databaseNameList[0] == "Public" and len(
                self._databaseNameList) > 1:
            userName = self._databaseNameList[1]
        else:
            userName = self._databaseNameList[0]
        self.addMeasurementsDocument(**doc,users=[userName])

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


