from . import abstractLocation
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from shapely.geometry import Point,box,MultiLineString, LineString

class ImageToolkit(abstractLocation.AbstractLocationToolkit):
    """
        A class to handle an image that represents a location.

        Looks up the location in the public database in project 'imageLocation'.

    """

    _presentation = None

    @property
    def presentation(self):
        return self._presentation

    def __init__(self, projectName,FilesDirectory=None):

        super().__init__(projectName=projectName,
                         toolkitName="Image",
                         FilesDirectory=FilesDirectory)

        self._presentation = presentation(dataLayer=self)

    def loadImageToDB(self, path, imageName, extents,**desc):
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
            return the document.
        """
        check = self.getLocationByRegion(imageName)

        if check is not None:
            raise KeyError(f"Image {imageName} already exists in project {self.projectName}")


        if isinstance(extents,dict):
            extentList = [extents['xmin'],extents['xmax'],extents['ymin'],extents['ymax']]
        elif isinstance(extents,list):
            extentList = extents
        else:
            raise ValueError("extents is either a list(xmin, xmax, ymin, ymax) or dict(xmin=, xmax=, ymin=, ymax=) ")

        imageparams = {abstractLocation.TOOLKIT_LOCATION_REGIONNAME : imageName,
                         "xmin":extentList[0],
                         "xmax":extentList[1],
                         "ymin":extentList[2],
                         "ymax":extentList[3]
                       }

        imageparams.update(desc)

        doc = self.addMeasurementsDocument(resource=path,
                                     dataFormat=abstractLocation.toolkit.datatypes.IMAGE,
                                     type=self.name,
                                     desc=imageparams)

        return doc


    def getImage(self,imageName,**filters):
        qry = {abstractLocation.TOOLKIT_LOCATION_REGIONNAME: imageName}
        qry.update(filters)
        docList = self.getMeasurementsDocuments(type=self.name,**qry)

        return None if len(docList) ==0 else docList[0]


    def listImages(self,**filters):
        return self.getMeasurementsDocuments(type=self.name,**filters)

class presentation:

    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self,dataLayer):
        self._datalayer = dataLayer

    def plot(self, imageNameOrData, ax=None):
        """
            Plot the image

        Parameters
        ----------

        imageName: str or image
            The image name from the

        Returns
        -------
            return the ax of the figure.
        """

        if isinstance(imageNameOrData,str):
            doc = self.datalayer.getImage(imageNameOrData)
            image = doc.getData()
        else:
            image = imageNameOrData


        if ax is None:
            fig, ax = plt.subplots()
        else:
            plt.sca(ax)

        path = doc.resource
        extents = [doc.desc['xmin'], doc.desc['xmax'], doc.desc['ymin'], doc.desc['ymax']]
        image = mpimg.imread(path)
        ax = plt.imshow(image, extent=extents)
        return ax

