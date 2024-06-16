import matplotlib.pyplot as plt
import geopandas
import numpy
import math
from ....utils.logging import get_classMethod_logger
from .... import toolkit
from itertools import product
from PIL import Image
import requests
from io import BytesIO

WSG84 = 4326
ITM = 2039
ED50_ZONE36N = 23036

class TilesToolkit(toolkit.abstractToolkit):
    """
        A class to handle an image that represents a location.

        Looks up the location in the public database in project 'imageLocation'.

    """
    WSG84 = 4326
    ITM    = 2039
    ED50_ZONE36N = 23036

    Z0RES = 156543.03

    @property
    def doctype(self):
        return f"{self.toolkitName}_PNG"


    def __init__(self, projectName, filesDirectory=None):

        super().__init__(projectName=projectName,
                         toolkitName="Tiles",
                         filesDirectory=filesDirectory)

        self._presentation = presentation(dataLayer=self)

    def tileScaleAtLatLonZoom(self,latitude,longitude,zoomlevel):
        """
            Returns the scale of a tile im meters at the location and zoom level
        Parameters
        ----------
        latitude :  float
                The latitude in WGS84
        longitude : float
                The longitude in WGS
        zoomlevel

        Returns
        -------
            float,
            The scale in m of one tile
        """
        return self.Z0RES / (2 ** zoomlevel) * numpy.cos(numpy.deg2rad(latitude))

    def getImageFromCorners(self, left,right,bottom,top, zoomlevel,tileServer,inputCRS=WSG84,outputCRS=WSG84):
        """
            Gets the image from the lower left corner and upper right cornet.
            The lowerLeft,upperRight are given in WGS84 (degrees) if dgrees are True
            and in Israel 1993 / Israeli TM Grid

        Parameters
        ----------
        left,right,bottom,top : float
            The left, right bottom and up of the image.
        zoomlevel :
        tileServer
        inputCRS : int
                The ESPG of the input coordinates.
        outputCRS: int
                The ESPG of the output coordinates.

        Returns
        -------
            Tuple,
                Image and extent [left,right,bottom,top] in the
                coordinate system of the outputCRS.

        """
        logger = get_classMethod_logger(self,name="getImageFromTiles")
        logger.info("------- Start")

        lon = [top,bottom]
        lat = [left, right]

        gdf = geopandas.GeoDataFrame(
            None, geometry=geopandas.points_from_xy(lat, lon), crs=inputCRS  # "EPSG:4326"
        )

        logger.info(f"Converting the input coordinates from EPSG {inputCRS} to WGS84 (EPSG:4326)")
        gdf = gdf.to_crs(WSG84)

        tileULX, tileULY = self.deg2tile(gdf.iloc[0].geometry.y,gdf.iloc[0].geometry.x, zoomlevel)
        tileLRX, tileLRY = self.deg2tile(gdf.iloc[1].geometry.y,gdf.iloc[1].geometry.x, zoomlevel)

        tileurl = self.getDataSourceData(tileServer)
        tileurl = tileServer if tileurl is None else tileurl
        img,extent =  self._getImageFromTiles([tileULX, tileULY],[tileLRX, tileLRY], zoomlevel,tileurl)

        logger.info(f"Converting the output extent {extent} from WGS84 (EPSG:4326) to EPSG {outputCRS}")
        lat = [extent[0], extent[1]]
        lon = [extent[2], extent[3]]

        gdf = geopandas.GeoDataFrame(
            None, geometry=geopandas.points_from_xy(lat,lon), crs=WSG84
        )
        gdf = gdf.to_crs(outputCRS)

        extent = [gdf.iloc[0].geometry.x,gdf.iloc[1].geometry.x,gdf.iloc[1].geometry.y,gdf.iloc[0].geometry.y]
        return img,extent




    def _getImageFromTiles(self, ulTiles, lrTiles, zoomLevel, tileServer, square=True):
        """
            Creates an image with range of the tiles and compute their
            extent in degrees

            The pro

        Parameters
        ----------
        llTiles : list, tuple
                    The lower left corner tile. (latitude (NE),longitude (EW))
        urTiles : list, tuple
                    The upper right corner tile. (latitude (NE),longitude (EW))

        zoomLevel : integer
                    The zoom to retrieve. usually up to ~19. (highest).

        square : bool
                If true, return a square image by extening the minimal dimension.

        tileServer: string
                The url of the tile server. We assume that the zoom, x and y are give by {z},{x} and {y}.
        Returns
        -------
            Tuple. :
                The Image and a the extent of the image in degrees (WGS84): [left,right,bottom,top]
        """
        logger = get_classMethod_logger(self,name="getImageFromTiles")
        logger.info("------- Start")

        height = numpy.abs(ulTiles[1] - lrTiles[1])
        width = numpy.abs(ulTiles[0] - lrTiles[0])

        logger.debug(f"The number of tiles in x: {ulTiles[1]}-{lrTiles[1]}={width} and y: {ulTiles[0]}-{lrTiles[0]} {height}.")

        if square:
            sqrx = numpy.max([width, height])
            sqry = numpy.max([width, height])
            logger.debug(f"Squaring the image to be of equal height and width: {sqrx}")
        else:
            logger.debug(f"Get the image as requested (no squaring)")
            sqry = height
            sqrx = width

        lrTiles[0] = ulTiles[0] + sqrx
        lrTiles[1] = ulTiles[1]+ sqry

        TILESIZE = 256
        finalimg = Image.new('RGB', (TILESIZE * sqrx, TILESIZE * sqry))
        logger.info(f"Getting the Image  (tiles {sqrx}x{sqry}) from {tileServer}")
        for i, j in product(range(sqrx), range(sqry)):
            #response = requests.get(f"http://192.168.14.118/resat_tiles/{zoomlevel}/{tileLLX + i}/{}.png")
            logger.debug(f"Getting address {tileServer.format(z=zoomLevel,x=ulTiles[0] + i,y=ulTiles[1] + j)}")
            response = requests.get(tileServer.format(z=zoomLevel,x=ulTiles[0] + i,y=ulTiles[1] + j))
            img = Image.open(BytesIO(response.content))
            logger.debug(f"({i},{j}) --> {ulTiles[0] + i}, {ulTiles[1] + j}")
            finalimg.paste(img, (i * TILESIZE, j * TILESIZE))

        boundULY, boundULX = self.tile2deg(ulTiles[0], ulTiles[1] , zoomLevel)
        boundLRY, boundLRX = self.tile2deg(lrTiles[0], lrTiles[1], zoomLevel)
        logger.info(f"Got the image extent (upper left -> lower right): ({ulTiles[0]},{ulTiles[1]})={[boundULY, boundULX]} and ({lrTiles[0]},{lrTiles[1]})={[boundLRY, boundLRX]}")

        return finalimg,[boundULX,boundLRX,boundULY,boundLRY]


    def tile2deg(self, xtile, ytile, zoom):
        n = 2.0 ** zoom
        lon_deg = xtile / n * 360.0 - 180.0
        lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
        lat_deg = math.degrees(lat_rad)
        return (lat_deg, lon_deg)

    def deg2tile(self, lat_deg, lon_deg, zoom):
        lat_rad = math.radians(lat_deg)
        n = 2.0 ** zoom
        xtile = int((lon_deg + 180.0) / 360.0 * n)
        ytile = int((1.0 - math.asinh(math.tan(lat_rad)) / math.pi) / 2.0 * n)
        return (xtile, ytile)

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
        elif isinstance(imageNameOrData,tuple):
            if extents is not None:
                raise ValueError("extents must be None if imageNameOrData is the tuple (image,extents)")
            image = imageNameOrData[0]
            extents = imageNameOrData[1]
        else:
            image = imageNameOrData
            if extents is None:
                    raise ValueError("extents must be supplied if imageNameOrData is image")

            if isinstance(extents, dict):
                if 'minX' in extents:
                    extents = [extents['minX'], extents['maxX'], extents['minY'], extents['maxY']]
                else:
                    extents = [extents['left'], extents['right'], extents['bottom'], extents['top']]
            elif isinstance(extents, list):
                extents = extents
            else:
                raise ValueError("extents is either a list(minX, maxX, minY, maxY), dict(minX=, maxX=, minY=, maxY=), or dict(left=, right=, bottom=, top=)")

        if ax is None:
            fig, ax = plt.subplots()
        else:
            plt.sca(ax)


        ax = plt.imshow(image, extent=extents)
        return ax



    #
    # def getImageAndStore(self,regionName, center,zoomlevel):
    #     """
    #         Gets an image and stores it in the project.
    #
    #     Parameters
    #     ----------
    #     regionName : The name of the image.
    #     center : a tuple with te center is WSG84 coordinates.
    #     zoomlevel : The zoom level to get
    #
    #
    #     Returns
    #     -------
    #
    #     """
    #     qry = {abstractLocation.TOOLKIT_LOCATION_REGIONNAME: regionName,
    #            abstractLocation.toolkitExtension.TOOLKIT_TOOLKITNAME_FIELD: self.toolkitName}
    #     qry.update(filters)
    #     docList = self.getCacheDocuments(type=self.doctype, **qry)
    #
