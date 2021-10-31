
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import os
from hera.toolkit import TOOLKIT_SAVEMODE_NOSAVE,TOOLKIT_SAVEMODE_ONLYFILE,TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,TOOLKIT_SAVEMODE_FILEANDDB,TOOLKIT_SAVEMODE_FILEANDDB_REPLACE
from hera.datalayer import datatypes, nonDBMetadataFrame

class RasterToolkit():
    """
        A class to handle an image that represents a location.

        Looks up the location in the public database in project 'imageLocation'.

    """

    @property
    def doctype(self):
        return f"{self.toolkitName}_PNG"


    def __init__(self, projectName, filesDirectory=None):

        super().__init__(projectName=projectName,
                         toolkitName="Image",
                         filesDirectory=filesDirectory)

        self._presentation = presentation(dataLayer=self)


    def loadData(self, fileNameOrData, extents, saveMode=TOOLKIT_SAVEMODE_NOSAVE, regionName=None, additionalData=dict()):
        """
            Loading a data from file, and possibly store the region in the database.

        Parameters
        ----------

        fileNameOrData: str
                If str , the datafile to load
                If other objects - convert the

        extents: dict or list
            if list:
                [minX, maxX, minY, maxY]

            if dict:
                the dictionary with the keys:
                    - minX
                    - maxX
                    - minY
                    - maxY


        saveMode: str
                Can be either:

                    - TOOLKIT_SAVEMODE_NOSAVE   : Just load the data from file and return the datafile

                    - TOOLKIT_SAVEMODE_ONLYFILE : Loads the data from file and save to a file.
                                                  raise exception if file exists.

                    - TOOLKIT_SAVEMODE_ONLYFILE_REPLACE: Loads the data from file and save to a file.
                                                  Replace the file if it exists.

                    - TOOLKIT_SAVEMODE_FILEANDDB : Loads the data from file and save to a file and store to the DB as a source.
                                                    Raise exception if the entry exists.

                    - TOOLKIT_SAVEMODE_FILEANDDB_REPLACE: Loads the data from file and save to a file and store to the DB as a source.
                                                    Replace the entry in the DB if it exists.


            kwargs: Contains:
                regionName: If fileNameOrData is an object, required if saveMode is not NOSAVE.
                additionalData: additional metadata to add if adding to the DB.

        Returns
        -------

            The doc (DB document or nonDB document).
        """

        localSave = True
        if isinstance(fileNameOrData,str):
            if os.path.exists(os.path.abspath(fileNameOrData)):
                localSave = False
                regionName = os.path.basename(fileNameOrData).split(".")[0] if regionName is None else regionName
                data = mpimg.imread(os.path.abspath(fileNameOrData))
            else:
                raise FileNotFoundError(f"The {fileNameOrData} does not exist.")

        else:
            regionName = None
            data = fileNameOrData

        doc = None

        if saveMode in [TOOLKIT_SAVEMODE_ONLYFILE,
                        TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,
                        TOOLKIT_SAVEMODE_FILEANDDB,
                        TOOLKIT_SAVEMODE_FILEANDDB_REPLACE]:

            outputFileName = os.path.join(self.FilesDirectory, f"{regionName}.png")

            if saveMode in [TOOLKIT_SAVEMODE_FILEANDDB,TOOLKIT_SAVEMODE_ONLYFILE]:
                if os.path.exists(outputFileName):
                    raise FileExistsError(f"{outputFileName} exists in project {self.projectName}")

            if localSave:
                resource= outputFileName
                mpimg.imsave(outputFileName,data)
            else:
                resource = fileNameOrData

            if saveMode in [TOOLKIT_SAVEMODE_FILEANDDB,TOOLKIT_SAVEMODE_FILEANDDB_REPLACE]:

                doc = self.getImage(regionName)
                if doc is not None and saveMode==TOOLKIT_SAVEMODE_FILEANDDB:
                    raise ValueError(f"{regionName} exists in DB for project {self.projectName}")

                if isinstance(extents, dict):
                    extentList = [extents['minX'], extents['maxX'], extents['minY'], extents['maxY']]
                elif isinstance(extents, list):
                    extentList = extents
                else:
                    raise ValueError(
                        "extents is either a list(minX, maxX, minY, maxY) or dict(minX=, maxX=, minY=, maxY=) ")

                additionalData.update({abstractLocation.TOOLKIT_LOCATION_REGIONNAME: regionName,
                                       abstractLocation.toolkit.TOOLKIT_DATASOURCE_NAME: regionName,
                                       abstractLocation.toolkit.TOOLKIT_TOOLKITNAME_FIELD: self.toolkitName,
                                       "minX": extentList[0],
                                       "maxX": extentList[1],
                                       "minY": extentList[2],
                                       "maxY": extentList[3]
                                       })

                if doc is None:
                    self.addCacheDocument(
                        type = self.doctype,
                        resource=resource,
                        dataFormat=datatypes.IMAGE,
                        desc = additionalData
                    )

                else:
                    doc['resource'] = resource
                    doc.desc = additionalData
                    doc.save()

        return nonDBMetadataFrame(data) if doc is None else doc


    def getImage(self, regionName, **filters):
        qry = {abstractLocation.TOOLKIT_LOCATION_REGIONNAME: regionName,
               abstractLocation.toolkit.TOOLKIT_TOOLKITNAME_FIELD: self.toolkitName}
        qry.update(filters)
        docList = self.getCacheDocuments(type=self.doctype, **qry)

        return None if len(docList) ==0 else docList[0]

    def listImages(self,**filters):
        return self.getMeasurementsDocuments(type=self.doctype, **filters)

class presentation:

    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self,dataLayer):
        self._datalayer = dataLayer

    def plot(self, imageNameOrData,extents=None, ax=None,**filters):
        """
            Plot the image

        Parameters
        ----------

        imageName: str or image (numpy.array)
            str - the resource name in the DB .
            if its array then its the data

        Returns
        -------
            return the ax of the figure.
        """
        if isinstance(imageNameOrData,str):
            doc = self.datalayer.getImage(imageNameOrData,**filters)
            extents = [doc.desc['minX'], doc.desc['maxX'], doc.desc['minY'], doc.desc['maxY']]
            image = doc.getData()
        else:
            image = imageNameOrData
            if extents is None:
                    raise ValueError("extents must be supplied if imageNameOrData is image")

            if isinstance(extents, dict):
                extents = [extents['minX'], extents['maxX'], extents['minY'], extents['maxY']]
            elif isinstance(extents, list):
                extents = extents
            else:
                raise ValueError("extents is either a list(minX, maxX, minY, maxY) or dict(minX=, maxX=, minY=, maxY=)")


        if ax is None:
            fig, ax = plt.subplots()
        else:
            plt.sca(ax)


        ax = plt.imshow(image, extent=extents)
        return ax

