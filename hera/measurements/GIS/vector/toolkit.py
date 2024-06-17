import io
from hera import toolkit
import geopandas 
from shapely.geometry import Polygon, box
from ..utils import ITM,ED50_ZONE36N,WSG84
from ....utils.logging import get_classMethod_logger


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

    @staticmethod
    def geopandasToGeoJson(geoData):
        if isinstance(geoData,geopandas.GeoDataFrame):
            dataHandler = io.BytesIO()
            geoData.to_file(dataHandler,driver='GeoJSON')
            return dataHandler.getvalue().decode('ascii')

        else:
            raise ValueError("Function receives only GeoDataFrame")

    def _RegionToGeopandas(self, regionData, crs = None):
        """
            Converts a shape to geopandas.

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
        if isinstance(regionData,geopandas.GeoDataFrame):
            data = regionData
        elif isinstance(regionData, dict):
            data = geopandas.GeoDataFrame.from_features(regionData['features'])

        elif isinstance(regionData, list):
            try:
                data =geopandas.GeoDataFrame({'geometry':[box(regionData[0], regionData[1], regionData[2], regionData[3])]})
            except Exception as e:
                raise ValueError(f"Can't create region shape from regionData list, The list should contain the following input: xmin,ymin,xmax,ymax ")

        elif isinstance(regionData, Polygon):
            data = geopandas.GeoDataFrame({'geometry':[regionData]})
        else:
            raise ValueError(f" regionData paremeter must be: GeoDataFrame/GeoJSON/list/Polygon")

        if crs is not None:
            try:
                data.crs = crs
            except Exception as e:
                raise ValueError(f" Please pass a valid crs number ")

        return data

    def cutRegionFromSource(self, datasourceDocument, shape, isBounds = False, inputCRS = WSG84):
        """
            Cuts a the shape from the requested datasource


        Parameters
        ----------
        datasourceDocument: hera document.
                The datasource.

        shape : list, dict, geopandas.
            list - a box with the corners in [xmin,ymin,xmax,ymax]
            dict - geoJSON

        datasourceName : str
            The name of the satasource to cur from.

        isBounds : bool
            If true, use the bounding box fo the polygon.

        inputCRS : int, default is WWSG84=4326.
            The EPSG of the coordinate system of the shape (if it is a shape and not in the dtasource ative coordinates).

        Returns
        -------
            geopandas of the right shape.
        """
        logger = get_classMethod_logger(self, "regionToSTL")

        regionWithCRS = self._RegionToGeopandas(shape, crs = inputCRS)

        logger.debug(f"The crs of the input is {regionWithCRS.crs}")
        logger.debug(f"The shape is {regionWithCRS.iloc[0]}")
        dct = dict(bbox=regionWithCRS) if isBounds else dict(mask=regionWithCRS)

        if regionWithCRS.crs is None:
            logger.execution("The region was defined without crs. Using the crs of the datasource.")
            regionWithCRS.crs = datasourceDocument.desc['crs']
        elif regionWithCRS.crs.to_epsg() != datasourceDocument.desc['crs']:
            logger.execution("shape and region crs mismatch. Converting the shape to the crs of the datasource.")
            regionWithCRS = regionWithCRS.to_crs(datasourceDocument.desc['crs'])
        else:
            logger.execution("shape and region crs match.")

        return datasourceDocument.getData(**dct)


