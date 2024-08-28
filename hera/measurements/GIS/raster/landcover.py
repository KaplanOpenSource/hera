from .... import toolkit
from ..utils import stlFactory,convertCRS,ITM,WSG84,ED50_ZONE36N,BETA,KARMAN
from ....utils.logging import get_classMethod_logger
import numpy
import math
from osgeo import gdal
from pyproj import Transformer
from itertools import product
import pandas
import os
import xarray
import geopandas
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from hera import toolkitHome
import matplotlib.colors as mcolors


class LandCoverToolkit(toolkit.abstractToolkit):
    # """
    # We should use different roughness to buildings and to topo
    #
    #  MCD12Q1 - MODIS/Terra+Aqua Land Cover Type Yearly L3 Global 500m SIN Grid
    #
    #   data is from https://zenodo.org/records/8367523
    #   type are in https://ladsweb.modaps.eosdis.nasa.gov/filespec/MODIS/6/MCD12Q1
    #
    #   Land Cover Type 1: IGBP global vegetation classification scheme
    #   Land Cover Type 2: University of Maryland (UMD) scheme
    #
    #   DataField    Name                    Data      Dimensions
    #                                        Type
    #
    #   DataField_1  Land_Cover_Type_1       UINT8     Dimension_1
    #                                                  Dimension_2
    #
    #           Description:    land cover types (IGBP)
    #
    #           DataField_1 HDF Attributes:
    #
    #              long_name      STRING  1   PGE    "Land_Cover_Type_1"
    #              units          STRING  1   PGE     "class number"
    #              valid_range    uint8    2   PGE     0 16
    #              _FillValue     uint8    1   PGE     255
    #              Water          uint8    1   PGE     0
    #              Evergreen needleleaf forest
    #                             uint8    1   PGE     1
    #              Evergreen broadleaf forest
    #                             uint8    1   PGE     2
    #              Deciduous needleleaf forest
    #                             uint8    1   PGE     3
    #              Deciduous broadleaf forest
    #                             uint8    1   PGE     4
    #              Mixed forests  unit8    1   PGE     5
    #              Closed shrubland
    #                             unit8    1   PGE     6
    #              Open shrublands
    #                             unit8    1   PGE     7
    #              Woody savannas unit8    1   PGE     8
    #              Savannas       unit8    1   PGE     9
    #              Grasslands     unit8    1   PGE    10
    #              Permanent wetlands
    #                             unit8    1   PGE    11
    #              Croplands      unit8    1   PGE    12
    #              Urban and built-up
    #                             unit8    1   PGE    13
    #              Cropland/natural vegetation mosaic
    #                             unit8    1   PGE    14
    #              Snow and ice   unit8    1   PGE    15
    #              Barren or sparsely vegetated
    #                             unit8    1   PGE    16
    #
    #   DataField_2  Land_Cover_Type_2       UINT8     Dimension_1
    #                                                  Dimension_2
    #
    #           Description:    land cover types (UMD)
    #
    #           DataField_2 HDF Attributes:
    #              long_name      STRING   1   PGE    "Land_Cover_Type_2"
    #              units          STRING   1   PGE     "class number"
    #              valid_range    uint8    2   PGE     0 16
    #              _FillValue     uint8    1   PGE     255
    #              Water          uint8    1   PGE     0
    #              Evergreen needleleaf forest
    #                             uint8    1   PGE     1
    #              Evergreen broadleaf forest
    #                             uint8    1   PGE     2
    #              Deciduous needleleaf forest
    #                             uint8    1   PGE     3
    #              Decidous broadleaf forest
    #                             uint8    1   PGE     4
    #              Mixed forests  uint8    1   PGE     5
    #              Closed shrublands
    #                              unit8   1   PGE     6
    #              Open shrubland  unit8   1   PGE     7
    #              Woody savannas  unit8   1   PGE     8
    #              Savannas       unit8    1   PGE     9
    #              Grasslands     uint8    1   PGE    10
    #              Croplands      unit8    1   PGE    12
    #              Urban and built-up
    #                             unit8    1   PGE    13
    #              Barren or sparsely vegetated
    #                             unit8    1   PGE    16
    #
    #   Land cover may change between winter and summer
    #   urban area may change over the years
    # """

    def __init__(self, projectName, filesDirectory=None):
        """
            Initializes land cover data toolkit.

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
        super().__init__(projectName=projectName, toolkitName = 'LandCoverToolkit', filesDirectory=filesDirectory)
        self._presentation = presentation(dataLayer=self)


    def getLandCoverAtPoint(self,lon,lat,inputCRS=WSG84, dataSourceName=None):
        """
        Get the landcover type integer value in a specific point.

        Parameters
        ----------
        lon : float
            The longitude coodinates.

        lat : float
            The latitute coordinates.

        inputCRS : int , default=WSG84
            The ESPRG of the coordinates.

        dataSourceName : string , default=None
            The name of the data source to use.

        Returns
        -------
            int

        """
        dataSourceName = self.getConfig()['defaultLandCover'] if dataSourceName is None else dataSourceName
        datasourceDocument = self.getDataSourceDocument(dataSourceName)
        ds = self.getDataSourceData(dataSourceName)
        img = ds.GetRasterBand(1).ReadAsArray()
        gt = ds.GetGeoTransform()

        width = ds.RasterXSize
        height = ds.RasterYSize
        minx = gt[0]
        maxx = gt[0] + width * gt[1] + height * gt[2]
        miny = gt[3] + width * gt[4] + height * gt[5]
        maxy = gt[3]
        x = math.floor((lat - gt[3]) / gt[5])
        y = math.floor((lon - gt[0]) / gt[1])
        return img[x, y]

    def getLandCover(self,minlon,minlat,maxlon,maxlat,dxdy = 30, inputCRS=WSG84, dataSourceName=None):
        """
        Get Xarray LandCover map.

        Parameters
        ----------
        minlon : float
            Minimum value of longitude coodinates.

        minlat : float
            Minimum value of latitute coodinates.

        maxlon: float
            Maximum value of longitude coodinates.

        maxlat: float
            Maximum value of latitute coodinates.

        dxdy: int, default=30
            Spatial resolution of the output land cover map. Defines the step size for the grid points
            in both the x (longitude) and y (latitude) directions within the specified bounding box.
            Smaller values result in finer resolution.

        inputCRS : int , default=WSG84
            The ESPRG of the coordinates.

        dataSourceName : string , default=None
            The name of the data source to use.

        Returns
        -------
            xarray.DataArray

        """

        def vectorizeLandCoverCalc(lat, lon, img, lonUpperLeft, lonResolution, latUpperLeft, latResolution):
            ilat = math.floor((lat - latUpperLeft) / latResolution)
            ilon = math.floor((lon - lonUpperLeft) / lonResolution)
            return img[ilat, ilon]

        min_pp = convertCRS(points=[[minlon, minlat]], inputCRS=inputCRS, outputCRS=ITM)[0]
        max_pp = convertCRS(points=[[maxlon, maxlat]], inputCRS=inputCRS, outputCRS=ITM)[0]
        x = numpy.arange(min_pp.x, max_pp.x, dxdy)
        y = numpy.arange(min_pp.y, max_pp.y, dxdy)
        xx = numpy.zeros((len(x), len(y)))
        yy = numpy.zeros((len(x), len(y)))
        for ((i, vx), (j, vy)) in product([(i, vx) for (i, vx) in enumerate(x)], [(j, vy) for (j, vy) in enumerate(y)]):
            print((i, j), end="\r")
            newpp = convertCRS(points=[[vx, vy]], inputCRS=ITM, outputCRS=WSG84)[0]
            xx[i, j] = newpp.x
            yy[i, j] = newpp.y

        ds = self.getDataSourceData(dataSourceName)
        img = ds.GetRasterBand(1).ReadAsArray()
        lonUpperLeft, lonResolution, lonRotation, latUpperLeft, latRotation, latResolution = ds.GetGeoTransform()
        vectorizedLandCover = numpy.vectorize(lambda tlat, tlon: vectorizeLandCoverCalc(lat=tlat,
                                                                                        lon=tlon,
                                                                                        img=img,
                                                                                        lonUpperLeft=lonUpperLeft,
                                                                                        lonResolution=lonResolution,
                                                                                        latUpperLeft=latUpperLeft,
                                                                                        latResolution=latResolution))
        landcover = vectorizedLandCover(xx, yy)
        ### Transform to XArray
        i = np.arange(landcover.shape[0])
        j = np.arange(landcover.shape[1])
        xarray = xr.DataArray(
            landcover,
            coords={
                'i': i,
                'j': j,
                'lat': (['i', 'j'], xx),
                'lon': (['i', 'j'], yy),
                'landcover': (['i', 'j'], landcover),
                'dxdy': dxdy
                },
            dims=['i', 'j']
            )
        xarray.attrs['landcover_description'] = self.getCodingMap(dataSourceName)
        return xarray

    def getRoughnessAtPoint(self,lon,lat,inputCRS=WSG84, dataSourceName=None):
        """
        Get the roughness value of a specific point in the map.

        Parameters
        ----------
        lon : float
            The longitude coodinates.

        lat : float
            The latitute coordinates.

        inputCRS : int , default=WSG84
            The ESPRG of the coordinates.

        dataSourceName : string , default=None
            The name of the data source to use.

        Returns
        -------
            float

        """
        dataSourceName = self.getConfig()['defaultLandCover'] if dataSourceName is None else dataSourceName
        datasourceDocument = self.getDataSourceDocument(dataSourceName)
        landcover = self.getLandCoverAtPoint(lon=lon,lat=lat,inputCRS=inputCRS,dataSourceName=dataSourceName)

        handlerFunction = getattr(self, f"_handleType{datasourceDocument['desc']['type']}")
        return handlerFunction(landcover)


    def getRoughnessFromLandcover(self,landcover,windMeteorologicalDirection=None,resolution=None,isBuilding=False,dataSourceName=None,GIS_BUILDINGS_dataSourceName=None):
        """
        Adds Roughness field (z0) to landcover Xarray.

        Parameters
        ----------
        landcover : Xarray
            Landcover Xarray map in which the Roughness field will be added to.

        windMeteorologicalDirection: double,default=None
            The meteorological angle. Must be specified only if data includes urbanic area.

        resolution: double,default=None
            The size of the squares. Must be specified only if data includes urbanic area.

        isBuilding : bool, default=False
            Is the landcover containts urbanic area.

        dataSourceName : string , default=None
            The name of the data source to use.

        GIS_BUILDINGS_dataSourceName: string , default=None
            The name of the GIS Buildings datasource name. Relevant if landcover contains Urban areas.

        Returns
        -------
            xarray.DataArray

        """
        dataSourceName = self.getConfig()['defaultLandCover'] if dataSourceName is None else dataSourceName
        datasourceDocument = self.getDataSourceDocument(dataSourceName)
        if not isBuilding:
            handlerFunction = getattr(self, f"_handleType{datasourceDocument['desc']['type']}")
            roughness_values = np.vectorize(handlerFunction)(landcover['landcover'])
            landcover = landcover.assign_coords(z0=(['i', 'j'], roughness_values))
        else:
            if not windMeteorologicalDirection or not resolution:
                raise ValueError("windMeteorologicalDirection and reolution must be specified for calculating urban area.")
            landcover = self._getUrbanRoughnessFromLandCover(landcover,windMeteorologicalDirection,resolution,GIS_BUILDINGS_dataSourceName)
        return landcover

    def getRoughness(self,minlon,minlat,maxlon,maxlat,dxdy = 30, inputCRS=WSG84, dataSourceName=None,isBuilding=False,GIS_BUILDINGS_dataSourceName=None):
        """
        Returns Xarray LandCover map with Roughness (zo) field. Just as applying getLandCover and getRoughnessFromLandcover at the same time.

        Parameters
        ----------
        minlon : float
            Minimum value of longitude coodinates.

        minlat : float
            Minimum value of latitute coodinates.

        maxlon: float
            Maximum value of longitude coodinates.

        maxlat: float
            Maximum value of latitute coodinates.

        dxdy: int, default=30
            Spatial resolution of the output land cover map. Defines the step size for the grid points
            in both the x (longitude) and y (latitude) directions within the specified bounding box.
            Smaller values result in finer resolution.

        inputCRS : int , default=WSG84
            The ESPRG of the coordinates.

        dataSourceName : string , default=None
            The name of the data source to use.

        isBuilding : bool, default=False
            Is the landcover containts building area.

        Returns
        -------
            xarray.DataArray
        """
        landcover = self.getLandCover(minlon,minlat,maxlon,maxlat,dxdy = dxdy, inputCRS=inputCRS, dataSourceName=dataSourceName)
        landcover = self.getRoughnessFromLandcover(landcover,isBuilding=isBuilding,dataSourceName=dataSourceName,GIS_BUILDINGS_dataSourceName=GIS_BUILDINGS_dataSourceName)
        return landcover

    def _handleType1(self,landcover):
        """
        Converting land type of Type-1 to roughness.
        Based on the paper:
            * https://wes.copernicus.org/articles/6/1379/2021/ table a2
            * https://doi.org/10.5194/wes-6-1379-2021 Satellite-based estimation of roughness lengths and displacement heights for wind resource modelling, Rogier Floors, Merete Badger, Ib Troen, Kenneth Grogan, and Finn-Hendrik Permien

        Parameters
        ----------
        landcover : int
            Landcover type value.

        Returns
        -------
            float
        """
        roughnessDict = {
            0 : 0.0001, # Water,  depends on waves that depends on wind
            1 : 1,       # Evergreen needleleaf forest
            2:  1,         # Evergreen broadleaf forest
            3:  1,  # Deciduous needleleaf forest
            4:  1, # Deciduous broadleaf forest
            5:  1, # Mixed forests
            6:  0.05, # Closed shrubland
            7:  0.06, # Open shrublands
            8: 0.05, # Woody savannas
            9: 0.15,  # Savannas
            10: 0.12,  # Grasslands
            11: 0.3, # Permanent wetlands
            12: 0.15, # Croplands
            13: 0.8, # Urban and built-up
            14: 0.14, # Cropland/natural vegetation mosaic
            15: 0.001, # Snow and ice
            16: 0.01, # Barren or sparsely vegetated
        }
        return roughnessDict.get(landcover,0.05)


    def getCodingMap(self,datasourceName):
        """
        Returns dictionary that maps landcover int value to string of landcover.

        Parameters
        ----------
        datasourceName: str
            Datasource type (Type-1, Type-2 etc.)

        Returns
        -------
            dict
        """
        dict = {}
        if datasourceName=='Type-1':
            dict = {
                    0: "Water",
                    1: "Evergreen needleleaf forest",
                    2: "Evergreen broadleaf forest",
                    3: "Deciduous needleleaf forest",
                    4: "Deciduous broadleaf forest",
                    5: "Mixed forests",
                    6: "Closed shrubland",
                    7: "Open shrublands",
                    8: "Woody savannas",
                    9: "Savannas",
                    10: "Grasslands",
                    11: "Permanent wetlands",
                    12: "Croplands",
                    13: "Urban and built-up",
                    14: "Cropland/natural vegetation mosaic",
                    15: "Snow and ice",
                    16: "Barren or sparsely vegetated"
                }

        return dict

    def _getUrbanRoughnessFromLandCover(self,landcover,windMeteorologicalDirection,resolution,GIS_BUILDINGS_dataSourceName):
        """
        Add Roughness for Urban areas to landcover Xarray. z0 and dd fields are calculated from LambdaP, LambdaF and HC of each Urban area.
        """
        gis_building_tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_BUILDINGS, projectName=self.projectName)

        min_pp = convertCRS(points=[[float(landcover[0,0].lat.values), float(landcover[0,0].lon.values)]], inputCRS=WSG84, outputCRS=ITM)[0]

        max_i, max_j = int(max(landcover.i).values) , int(max(landcover.j).values)
        max_pp = convertCRS(points=[[float(landcover[max_i,max_j].lat.values), float(landcover[max_i,max_j].lon.values)]], inputCRS=WSG84, outputCRS=ITM)[0]


        buildings = gis_building_tk.getBuildingsFromRectangle(minx=min_pp.x,miny=min_pp.y,maxx=max_pp.x,maxy=max_pp.y,dataSourceName=GIS_BUILDINGS_dataSourceName,inputCRS=ITM)
        lambdaGrid = gis_building_tk.analysis.LambdaFromBuildingData(windMeteorologicalDirection, resolution, buildings)

        lambdaGrid.crs = ITM
        lambdaGrid = self._getRoughnessFromBuildingsDataFrame(lambdaGrid)
        square_size = float(landcover.dxdy.values)

        landcover['z0'] = (('i', 'j'), np.full((landcover.sizes['i'], landcover.sizes['j']), np.nan))
        landcover['dd'] = (('i', 'j'), np.full((landcover.sizes['i'], landcover.sizes['j']), np.nan))

        for i, arr in enumerate(landcover):
            for j, x in enumerate(arr):
                shape = convertCRS([[x.lat, x.lon]], inputCRS=WSG84, outputCRS=ITM)[0]
                lambdas = lambdaGrid.loc[shape.intersects(lambdaGrid['geometry'])]
                landcover['z0'].values[i, j] = lambdas['zz0'].values[0]
                landcover['dd'].values[i, j] = lambdas['dd'].values[0]

        return landcover

    def _getRoughnessFromBuildingsDataFrame(self,lambdaGrid):
        lambdaGrid.loc[(lambdaGrid.hc < 2), "lambdaF"] = 0.25
        lambdaGrid.loc[(lambdaGrid.hc < 2), "lambdaP"] = 0.25
        lambdaGrid.loc[(lambdaGrid.hc < 2), "hc"] = 2
        lambdaGrid.loc[(lambdaGrid.lambdaF > 0.4), "lambdaF"] = 0.4
        lambdaGrid.loc[(lambdaGrid.lambdaP > 0.6), "lambdaP"] = 0.6
        lambdaGrid["Lc"] = lambdaGrid["hc"] * (1 - lambdaGrid["lambdaP"]) / lambdaGrid["lambdaF"]
        lambdaGrid["ll"] = 2 * (BETA ** 3) * lambdaGrid["Lc"]
        lambdaGrid["zz0"] = lambdaGrid["ll"] / KARMAN * np.exp(-KARMAN / BETA)
        lambdaGrid["dd"] = lambdaGrid["ll"] / KARMAN
        return lambdaGrid


class presentation:
    """
    Presentation Layer of LandCover toolkit.
    """
    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, dataLayer):
        self._datalayer = dataLayer
        self.landcover_colors_map = {
            0:'red',
            1:'blue'
        }

    def plotRoughness(self,plot,landcover,alpha=0.5,figsize=(28,28)):
        """
        Plot Roughness upon Given Axes.

        Parameters
        ----------
        plot: matplotlib.image.AxesImage
            Satile Image to plot Polygons on.

        landcover: xarray
            Landcover xarray of plot area.

        alpha: float, default=0.2
            Opaqueness level.

        figsize: tuple, default=(28,28)
            Figure size.

        Returns
        -------
        """
        fig, ax = plt.subplots(figsize=figsize)
        rectangles = self._getRoughnessRectangles(landcover)
        self._plotWithRectangles(ax, plot, rectangles, alpha)
        plt.show()


    def plotLandcover(self,plot,landcover,alpha=0.2,figsize=(28,28)):
        """
        Plot LandCover upon Given Axes.

        Parameters
        ----------
        plot: matplotlib.image.AxesImage
            Satile Image to plot Polygons on.

        landcover: xarray
            Landcover xarray of plot area.

        alpha: float, default=0.2
            Opaqueness level.

        figsize: tuple, default=(28,28)
            Figure size.

        Returns
        -------
        """

        fig, ax = plt.subplots(figsize=figsize)
        rectangles = self._getLandcoverRectangles(landcover)
        self._plotWithRectangles(ax,plot,rectangles,alpha)
        plt.show()

    def _plotWithRectangles(self,ax,plot,rectangles,alpha):
        ax.imshow(plot.get_array(), extent=plot.get_extent())
        ax = self._adddRectanglesToPlot(ax, rectangles,alpha)

        ax.set_xlim(plot.get_extent()[0], plot.get_extent()[1])
        ax.set_ylim(plot.get_extent()[2], plot.get_extent()[3])

        return ax

    def _adddRectanglesToPlot(self,ax,rectangles,alpha):
        for rect in rectangles:
            x, y, width, height, color = rect
            rectangle = patches.Rectangle(
                (x, y),  # Bottom-left corner
                width,  # Width
                height,  # Height
                linewidth=1,
                edgecolor=color,
                facecolor=color,
                alpha=alpha
            )
            ax.add_patch(rectangle)

        return ax
    def _getLandcoverRectangles(self,landcover):
        rectangles = []
        for arr in landcover:
            for x in arr:
                shape = convertCRS([[x.lat, x.lon]], inputCRS=WSG84, outputCRS=ITM)[0]
                rectangles.append((shape.x, shape.y, float(landcover.dxdy.values), float(landcover.dxdy.values), self.landcover_colors_map.get(int(x.values))))

        return list(set(rectangles))

    def _getRoughnessRectangles(self,landcover):
        rectangles = []
        colormap = plt.cm.viridis
        norm = mcolors.Normalize(vmin=landcover.z0.min().values, vmax=landcover.z0.max().values)

        for arr in landcover:
            for x in arr:
                shape = convertCRS([[x.lat, x.lon]], inputCRS=WSG84, outputCRS=ITM)[0]
                color = colormap(norm(x['z0'].values))
                rectangles.append((shape.x, shape.y, float(landcover.dxdy.values), float(landcover.dxdy.values),
                                   color))

        return list(set(rectangles))