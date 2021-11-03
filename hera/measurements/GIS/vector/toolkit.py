import io
from hera import toolkit
import geopandas as gp
from shapely.geometry import Polygon, box
from hera.datalayer.datahandler import datatypes

TOOLKIT_VECTOR_REGIONNAME = "regionName"


class VectorToolkit(toolkit.abstractToolkit):


    def __init__(self, projectName, toolkitName = 'VectorToolkit', filesDirectory=None):
        """
            Initializes vector data toolkit.

        Parameters
        ----------
        projectName: str
            The project Name that the toolkit is initialized on
        toolkitName: str
            the specific toolkit, getting from the child.

        FilesDirectory: str or None
                The path to save a regions files when they are created.

                if str then represents a path (relative or absolute) to save the files in. The directory is created automatically.

                if None, then tries to get the default path of the project from the config. if it does not
                exist, then use the current directory.

        """
        super().__init__(projectName=projectName,toolkitName=toolkitName,filesDirectory=filesDirectory)

        # FileDirectoryCheck in abstractToolkit

    # @staticmethod
    # def geopandasToGeoJson(geoData):
    #     features = []
    #     insert_features = lambda X: features.append(
    #         geojson.Feature(geometry=geojson.Point((X["long"],
    #                                                 X["lat"],
    #                                                 X["elev"])),
    #                         properties=dict(name=X["name"],
    #                                         description=X["description"])))
    #     geoData.apply(insert_features, axis=1)
    #     with open('map1.geojson', 'w', encoding='utf8') as fp:
    #         geojson.dump(geojson.FeatureCollection(features), fp, sort_keys=True, ensure_ascii=False)

    @staticmethod
    def geopandasToGeoJson(geoData):

        if isinstance(geoData,gp.GeoDataFrame):
            dataHandler = io.BytesIO()
            geoData.to_file(dataHandler,driver='GeoJSON')
            return dataHandler.getvalue().decode('ascii')

        else:
            raise ValueError("Function receives only GeoDataFrame")

    def _setGeoPandasFromRegionData(self,regionData, crs = None):
        """
            Converts a shapte to geopandas.

        Parameters
        ----------
        regionData: GeoDataFrame,dict,list,Polygon
            The shape data to use.

            list - a box with  [xmin,ymin,xmax,ymax]

        crs : int
            The CRS EPSG identification.

        Returns
        -------

        """

        if isinstance(regionData,gp.GeoDataFrame):

            data = regionData

        elif isinstance(regionData, dict):
            data = gp.GeoDataFrame.from_features(regionData['features'])

        elif isinstance(regionData, list):
            try:
                data =gp.GeoDataFrame({'geometry':[box(regionData[0], regionData[1], regionData[2], regionData[3])]})
            except Exception as e:
                raise ValueError(f"Can't create region shape from regionData list, The list should contain the following input: xmin,ymin,xmax,ymax ")

        elif isinstance(regionData, Polygon):
            data = gp.GeoDataFrame({'geometry':[regionData]})
        else:
            raise ValueError(f" regionData paremeter must be: GeoDataFrame/GeoJSON/list/Polygon")

        if crs is not None:
            try:
                data.crs = crs
            except Exception as e:
                raise ValueError(f" Please pass a valid crs number ")

        return data

    def addRegion(self,regionData, regionName, crs = None): #isBounds - save the shape as Bbox or not , use can add CRS for the point list or shape
        """
                Add region shape to DB.

                Parameters
                ----------
                regionData: geoPandas , geoJson, shaply geometry
                    The project Name that the toolkit is initialized on
                regionName: str
                   The document name to save with in DB.

                isBounds: Boolean
                    If the requested region

                crs: int
                    the EPSG number represent the region coordinate system.
                    * The crs is not saved to the db with the region when it is specified in the geoPandas, only when sendding
                      the crs to the function.

        """
        desc ={TOOLKIT_VECTOR_REGIONNAME: regionName,'crs':crs}

        data = self._setGeoPandasFromRegionData(regionData,crs = crs).to_json()
        # data = self.geopandasToGeoJson(data)

        self.addCacheDocument(resource = data, dataFormat=datatypes.GEOPANDAS, desc=desc)


    def cutRegionFromSource(self,shapeDataOrName,dataSourceName, isBounds = False, crs = None): # If  shapeDataOrName is data: if is Bounds = True: use the Bbox of shape as the region, else use the shpae as the region
        """
            Cuts a the shape from the requested datasource

        Parameters
        ----------
        shapeDataOrName: GeoDataFrame,dict,list,Polygon
            The shape data to use.

                        list - a box with the corners in [xmin,ymin,xmax,ymax]

        dataSourceName : str
            The name of the satasource to cur from.

        isBounds : bool
            If true, use the bounding box fo the polygon.

        crs : int
            The EPSG of the coordinate system of the shape (if it is a shape and not in the dtasource ative coordinates).

        Returns
        -------

        """
        self.logger.info("-- Start --")


        if isinstance(shapeDataOrName, str):
            shape= self.getRegionData(shapeDataOrName)
        else:
            shape = self._setGeoPandasFromRegionData(shapeDataOrName,crs = crs)

            shape.crs=2039

        self.logger.debug(f"The crs of the input is {shape.crs}")
        self.logger.debug(f"The shape is {shape.iloc[0]}")

        dct = dict(bbox=shape) if isBounds else dict(mask=shape)

        if isinstance(dataSourceName, str):
            doc = self.getDatasourceDocument(datasourceName=dataSourceName)
            self.logger.debug(f"The datasource {dataSourceName} is pointing to {doc.resource}")


            import pdb
            pdb.set_trace()

            if doc is None:
                sourceList = self.getDataSourceTable()['datasourceName'].str.cat()
                err = f"The Data sources available in the project are: " + sourceList
                self.logger.error(err)
                raise ValueError(err)
            else:
                return doc.getData(**dct)

    def getRegionNameList(self):
        """
            Return the list of region names.

        :return:
        """

        return [doc.desc[TOOLKIT_VECTOR_REGIONNAME] for doc in self.getMeasurementsDocuments()]

    def getRegionData(self, regionName):

        """

                Parameters
                ---------

                    regionName: str
                            The name of the region
                Returns
                -------
                    Return the vectorData with the region anme
                """
        qry ={TOOLKIT_VECTOR_REGIONNAME: regionName}
        shapeDoc = self.getCacheDocuments(**qry)
        return None if len(shapeDoc)==0 else shapeDoc[0].getData()

    def getRegionDocumentByName(self, regionName):
        """

        Parameters
        ---------

            regionName: str
                    The name of the region
        Returns
        -------
            Return the vectorData with the region anme
        """

        return self.getDatasourceDocument(regionName)

    def getRegionsTable(self):

        return self.getDataSourceTable()

    def deleteRegion(self,regionName):
        self.deleteCacheDocuments(desc ={TOOLKIT_VECTOR_REGIONNAME: regionName})




