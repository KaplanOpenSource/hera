import os
import logging
import numpy
from ....datalayer import project
from ....datalayer import datatypes

import matplotlib.pyplot as plt

from shapely.geometry import Point,box,MultiLineString, LineString




class datalayer(project.ProjectMultiDBPublic):
    """
        A class to handle an image that represents a location.

        Looks up the location in the public database in project 'imageLocation'.

    """

    @property
    def publicProjectName(self):
        return 'imageLocation'


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

        super().__init__(projectName=projectName,publicProjectName=self.publicProjectName)

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


