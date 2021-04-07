import geopandas
import io
from descartes import PolygonPatch
from matplotlib.collections import PatchCollection

import matplotlib.pyplot as plt
from ... import toolkit

from ...datalayer import datatypes, nonDBMetadataFrame
import os

class ShapesToolKit(toolkit.abstractToolkit):
    """
        Manages geoJSON shape files.
        Holds the geoJSON in the resource of the mongoDB
    """

    @property
    def doctype(self):
        return f"{self.toolkitName}_GeoJSON"



    def __init__(self, projectName):

        super().__init__(projectName=projectName,
                         toolkitName="Shapes")

        self._presentation = presentation(dataLayer=self)

    def loadData(self, fileNameOrData, saveMode=toolkit.TOOLKIT_SAVEMODE_NOSAVE, regionName=None, additionalData=dict()):
        """
            Loading shapes from a geojson file (or other file that geopandas read) and return
            as geoDataFrame.  Possibly store as a region in the database (as a geoJSON string in the resource).

        Parameters
        ----------
        fileNameOrData: str, geopandas.geoDataFrame
                If str , the datafile to load or the geoJSON str
                geopandas.geoDataFrame the geodata frame

        saveMode: str
                Can be either:

                    - TOOLKIT_SAVEMODE_NOSAVE   : Just load the data from file and return the datafile


                    - TOOLKIT_SAVEMODE_FILEANDDB : Loads the data from file and save to a file and store to the DB as a source.
                                                    Raise exception if the entry exists.

                    - TOOLKIT_SAVEMODE_FILEANDDB_REPLACE: Loads the data from file and save to a file and store to the DB as a source.
                                                    Replace the entry in the DB if it exists.


            kwargs: Contains:
                regionName: If fileNameOrData is an object, required if saveMode is not NOSAVE.
                additionalData: additional metadata to add if adding to the DB.

        Returns
        -------
            The document with a DB.
        """

        if isinstance(fileNameOrData,str):
            if os.path.exists(os.path.abspath(fileNameOrData)):
                regionName = os.path.basename(fileNameOrData).split(".")[0] if regionName is None else regionName
                data = geopandas.read_file(os.path.abspath(fileNameOrData))
            else:

                data = geopandas.read_file(io.StringIO(fileNameOrData))
        elif isinstance(fileNameOrData,geopandas.geodataframe.GeoDataFrame):
            data = fileNameOrData
        else:
            raise ValueError(f"fileNameOrData must be a filename, geoJSON str or geoPandas.geodataframa, got type {type(fileNameOrData)}")

        doc = None

        if saveMode in [toolkit.TOOLKIT_SAVEMODE_FILEANDDB,
                        toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE]:

            if regionName is None:
                raise ValueError("Must suply regionName when when saving to the DB.")

            doc = self.getShape(regionName)

            if doc is None:
                additionalData.update({
                                       toolkit.TOOLKIT_DATASOURCE_NAME: regionName,
                                       toolkit.TOOLKIT_TOOLKITNAME_FIELD: self.toolkitName
                                       })

                doc = self.addCacheDocument(
                            type = self.doctype,
                            resource=data.to_json(),
                            dataFormat=datatypes.JSON_GEOPANDAS,
                            desc = additionalData
                        )
            else:
                    doc['resource'] = data.to_json()
                    doc.desc = additionalData
                    doc.save()

        return nonDBMetadataFrame(data) if doc is None else doc


    def getShape(self, regionName):
        """
        Get the shape with the region name.

        Parameters:
        ----------
            regionName: str
                The name of the shape region

        """
        qry = {toolkit.TOOLKIT_DATASOURCE_NAME: regionName,
               toolkit.TOOLKIT_TOOLKITNAME_FIELD: self.toolkitName}

        docList = self.getCacheDocuments(**qry)
        return None if len(docList)==0 else docList[0]

class presentation():

    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self,dataLayer):

        self._datalayer = dataLayer


    def plot(self, regionNameOrData, ax=None,**patchParams):
        """
        Plots saved geometry shapes.

        Parameters:
        -----------
            regionNameOrData: str, geopandas.DataFrame
                    The name of the shape to plot, a geoJSON str or the geopandas.DataFrame.
            ax: matplotlib ax
                The ax of the plot.
        Returns
        --------
            The ax of the plot
        """
        if ax is None:
            fig, ax = plt.subplots(1,1)
        else:
            plt.sca(ax)

        if isinstance(regionNameOrData,str):
            shapesDoc = self.datalayer.getShape(regionNameOrData)

            if shapesDoc is None:
                shapes = geopandas.read_file(io.StringIO(regionNameOrData))
            else:
                shapes = shapesDoc.getData()
        elif isinstance(regionNameOrData,geopandas.geodataframe):
            shapes = regionNameOrData
        else:
            raise ValueError("must supply a region name, the geoJSON str or a geopandas")

        patches = []
        for pol in shapes:
            patches.append(PolygonPatch(pol,**patchParams))

        ax.add_collection(PatchCollection(patches, match_original=True))
        return ax
