import os
import logging
import numpy
from .abstractLocation import datalayer as locationDatalayer
from ....datalayer import datatypes
from ....datalayer import project
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from shapely.geometry import Point,box,MultiLineString, LineString

class datalayer(project.ProjectMultiDBPublic):
    """
        A class to handle an image that represents a location.

        Looks up the location in the public database in project 'imageLocation'.

    """

    _projectName = None
    _presentation = None

    docType = "GIS_Image"

    @property
    def presentation(self):
        return self._presentation

    def __init__(self, projectName, databaseNameList=None, useAll=False,publicProjectName="Images"):

        self._projectName = projectName
        super().__init__(projectName=projectName, publicProjectName=publicProjectName,databaseNameList=databaseNameList,useAll=useAll)
        self._presentation = presentation(projectName=projectName,dataLayer=self)

    def load(self, path, imageName, extents,**desc):
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

        desc: additional description of the figure.

        Returns
        -------
        """
        check = self.getMeasurementsDocuments(dataFormat='image',type=self.docType ,imageName=imageName)
        if len(check)>0:
            raise KeyError("The imageName is already used.")
        if isinstance(extents,dict):
            extentList = [extents['xmin'],extents['xmax'],extents['ymin'],extents['ymax']]
        elif isinstance(extents,list):
            extentList = extents
        else:
            raise ValueError("extents is either a list(xmin, xmax, ymin, ymax) or dict(xmin=, xmax=, ymin=, ymax=) ")

        imageparams = dict(imageName=imageName,
             xmin=extentList[0],
             xmax=extentList[1],
             ymin=extentList[2],
             ymax=extentList[3]
             )

        imageparams.update(desc)

        doc = dict(resource=path,
                   dataFormat='image',
                   type='GIS',
                   desc=imageparams)

        if self._databaseNameList[0] == "public" or self._databaseNameList[0] == "Public" and len(
                self._databaseNameList) > 1:
            userName = self._databaseNameList[1]
        else:
            userName = self._databaseNameList[0]
        self.addMeasurementsDocument(**doc,users=[userName])

    def listImages(self,**filters):
        return self.getMeasurementsDocuments(type=self.docType ,imageName=imageName)

class presentation():

    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, projectName, dataLayer=None, databaseNameList=None, useAll=False,
                 publicProjectName="Images"):

        self._datalayer = datalayer(projectName=projectName, publicProjectName=publicProjectName,
                         databaseNameList=databaseNameList, useAll=useAll) if datalayer is None else dataLayer

    def plot(self, imageName, ax=None, **query):
        """
        :param imageName: The location name
        :param query: Some more specific details to query on
        :return:
        """
        doc = self.datalayer.getMeasurementsDocuments(dataFormat='image',type='GIS',imageName=imageName,**query)
        if len(doc) > 1:
            raise ValueError('More than 1 documents fills those requirements')
        doc = doc[0]

        if ax is None:
            fig, ax = plt.subplots()
        else:
            plt.sca(ax)

        path = doc.resource
        extents = [doc.desc['xmin'], doc.desc['xmax'], doc.desc['ymin'], doc.desc['ymax']]
        image = mpimg.imread(path)
        ax = plt.imshow(image, extent=extents)
        return ax

