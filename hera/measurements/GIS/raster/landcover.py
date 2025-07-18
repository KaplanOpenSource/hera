# ================== IMPORTS ================== #
from hera import toolkit
from hera.measurements.GIS.utils import ITM, WSG84, BETA, KARMAN, convertCRS, create_xarray
import numpy
import math
from itertools import product
import xarray
import xarray as xr
import numpy as np
import math
import os
import rasterio
from pyproj import Transformer
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from hera import toolkitHome
import matplotlib.colors as mcolors
import warnings
from tqdm import tqdm
import matplotlib.cm as cm
import geopandas as gpd
import rasterio  # [FIX] Added rasterio for correct .tif file handling

# Suppress specific GDAL warning
warnings.filterwarnings("ignore",
                        message="Starting with GDAL 3.7, PIXELTYPE=SIGNEDBYTE is no longer used to signal signed 8-bit raster.*")


# ================== CLASS ================== #
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
    #           DataField_1 HDF Attributes:
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
    #              Deciduous broadleaf forest
    #                             uint8    1   PGE     4
    #              Mixed forests  uint8    1   PGE     5
    #              Closed shrublands
    #                              unit8   1   PGE     6
    #              Open shrubland  unit8   1   PGE     7
    #              Woody savannas  unit8    1   PGE     8
    #              Savannas       unit8    1   PGE     9
    #              Grasslands     uint8    1   PGE    10
    #              Croplands      uint8    1   PGE    12
    #              Urban and built-up
    #                             uint8    1   PGE    13
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
        super().__init__(projectName=projectName, toolkitName='LandCoverToolkit', filesDirectory=filesDirectory)
        self._presentation = presentation(dataLayer=self)



    def getLandCoverAtPoint(self, lon, lat, inputCRS=WSG84, dataSourceName=None):
        """
        Get the landcover type integer value in a specific point.

        Parameters
        ----------
        lon : float
            The longitude coordinate.

        lat : float
            The latitude coordinate.

        inputCRS : int, default=WSG84
            The EPSG of the input coordinates.

        dataSourceName : string, default=None
            The name or path of the data source to use.

        Returns
        -------
        int
            Land cover integer value at the point.
        """

        # קבע שם ברירת מחדל אם לא ניתן
        dataSourceName = self.getConfig()['defaultLandCover'] if dataSourceName is None else dataSourceName

        # נסה להשיג את ה-DataSource דרך הפונקציה הקיימת
        ds = self.getDataSourceData(dataSourceName)

        # אם נכשל, ננסה לפתוח את הקובץ ידנית (כמו ביוניט־טסט)
        if ds is None:
            if isinstance(dataSourceName, str) and os.path.isfile(dataSourceName):
                ds = rasterio.open(dataSourceName)
            else:
                raise ValueError(f"Could not load data source: {dataSourceName}")

        # המר קואורדינטות אם צריך
        if inputCRS != ds.crs.to_epsg():
            transformer = Transformer.from_crs(inputCRS, ds.crs.to_epsg(), always_xy=True)
            lon, lat = transformer.transform(lon, lat)

        # קריאה של התמונה והמידע המרחבי
        img = ds.read(1)
        gt = ds.transform
        width = ds.width
        height = ds.height

        # חישוב אינדקסים של פיקסל
        x = math.floor((lat - gt.f) / gt.e)
        y = math.floor((lon - gt.c) / gt.a)

        # בדיקה שהפיקסל בטווח
        if x < 0 or x >= img.shape[0] or y < 0 or y >= img.shape[1]:
            raise IndexError(f"Point ({lon},{lat}) is out of raster bounds.")

        return int(img[x, y])


    def getLandCover(self, minx, miny, maxx, maxy, dxdy=30, inputCRS=WSG84, dataSourceName=None):
        """
        Get Xarray LandCover map.

        Parameters
        ----------
        minx : float
            Minimum value of X axis. If using WSG84 - the minimum Latitude coordinates.

        miny : float
            Minimum value of Y axis. If using WSG84 - the minimum Longitude coordinates.

        maxx: float
            Maximum value of X axis. If using WSG84 - the maximum Latitude coordinates.

        maxy: float
            Maximum value of Y axis. If using WSG84 - the maximum Longitude coordinates.

        dxdy: int, default=30
            Spatial resolution of the output land cover map.

        inputCRS : int, default=WSG84
            The EPSG of the coordinates.

        dataSourceName : string, default=None
            The name of the data source to use.

        Returns
        -------
        xarray.DataArray
        """

        def vectorizeLandCoverCalc(lat, lon, img, lonUpperLeft, lonResolution, latUpperLeft, latResolution):
            ilat = math.floor((lat - latUpperLeft) / latResolution)
            ilon = math.floor((lon - lonUpperLeft) / lonResolution)
            return img[ilat, ilon]

        xarray_dataset = create_xarray(minx, miny, maxx, maxy, dxdy, inputCRS)
        dataSourceName = self.getConfig()['defaultLandCover'] if dataSourceName is None else dataSourceName

        # [ADD] Handle raw file path directly if given
        if isinstance(dataSourceName, str) and os.path.isfile(dataSourceName):
            with rasterio.open(dataSourceName) as src:
                img = src.read(1)
                lonUpperLeft, latUpperLeft = src.transform[2], src.transform[5]
                lonResolution, latResolution = src.transform[0], src.transform[4]
        else:
            ds = self.getDataSourceData(dataSourceName)
            if hasattr(ds, "GetRasterBand"):
                img = ds.GetRasterBand(1).ReadAsArray()
                lonUpperLeft, lonResolution, lonRotation, latUpperLeft, latRotation, latResolution = ds.GetGeoTransform()
            else:
                img = ds.read(1)
                lonUpperLeft, latUpperLeft = ds.transform[2], ds.transform[5]
                lonResolution, latResolution = ds.transform[0], ds.transform[4]

        vectorizedLandCover = numpy.vectorize(lambda tlat, tlon: vectorizeLandCoverCalc(
            lat=tlat,
            lon=tlon,
            img=img,
            lonUpperLeft=lonUpperLeft,
            lonResolution=lonResolution,
            latUpperLeft=latUpperLeft,
            latResolution=latResolution
        ))

        landcover = vectorizedLandCover(xarray_dataset.lat.values, xarray_dataset.lon.values)
        xarray_dataset = xarray_dataset.assign_coords(landcover=(('i', 'j'), landcover))
        return xarray_dataset

    def getRoughnessAtPoint(self, lon, lat, inputCRS=WSG84, dataSourceName=None):
        """
        Get the roughness value of a specific point in the map.

        Parameters
        ----------
        lon : float
            The longitude coordinates.

        lat : float
            The latitude coordinates.

        inputCRS : int, default=WSG84
            The EPSG of the coordinates.

        dataSourceName : string, default=None
            The name of the data source to use.

        Returns
        -------
        float
        """
        dataSourceName = self.getConfig()['defaultLandCover'] if dataSourceName is None else dataSourceName

        try:
            datasourceDocument = self.getDataSourceDocument(dataSourceName)
            if datasourceDocument and 'desc' in datasourceDocument:
                type_name = datasourceDocument['desc'].get('type', None)
                if not type_name and 'desc' in datasourceDocument['desc']:
                    type_name = datasourceDocument['desc']['desc'].get('type', None)
            else:
                type_name = None
        except Exception:
            datasourceDocument = None
            type_name = None

        if type_name is None:
            print("[WARNING] No metadata available for data source. Using default type: IGBP.")
            type_name = "IGBP"

        landcover = self.getLandCoverAtPoint(lon=lon, lat=lat, inputCRS=inputCRS, dataSourceName=dataSourceName)

        if type_name == "IGBP":
            roughness_mapping = {
                0: 0.01,  # Example values
                1: 0.1,
                2: 0.15,
                3: 0.2,
                4: 0.25,
                5: 0.3,
                6: 0.4,
                7: 0.5,
                8: 0.6,
                9: 0.7,
                10: 0.8,
                11: 0.9,
                12: 1.0,
                13: 1.1,
                14: 1.2,
                15: 1.3,
                16: 1.4
            }
            return float(roughness_mapping.get(landcover, 0.05))  # ברירת מחדל אם אין התאמה
        else:
            handlerFunction = getattr(self, f"_handleType{type_name}")
            return handlerFunction(landcover)

    def getRoughnessFromLandcover(self, landcover, windMeteorologicalDirection=None, resolution=None, isBuilding=False,
                                  dataSourceName=None, GIS_BUILDINGS_dataSourceName=None):
        """
        Adds Roughness field (z0) to landcover Xarray.

        Parameters
        ----------
        landcover : Xarray
            Landcover Xarray map in which the Roughness field will be added to.

        windMeteorologicalDirection: double, default=None
            The meteorological angle. Must be specified only if data includes urbanic area.

        resolution: double, default=None
            The size of the squares. Must be specified only if data includes urbanic area.

        isBuilding : bool, default=False
            Is the landcover contains urbanic area.

        dataSourceName : string, default=None
            The name of the data source to use.

        GIS_BUILDINGS_dataSourceName: string, default=None
            The name of the GIS Buildings datasource name. Relevant if landcover contains Urban areas.

        Returns
        -------
        xarray.DataArray
        """
        dataSourceName = self.getConfig()['defaultLandCover'] if dataSourceName is None else dataSourceName

        try:
            datasourceDocument = self.getDataSourceDocument(dataSourceName)
        except Exception:
            datasourceDocument = None

        if not isBuilding:
            # במקרה שאין מטא־דאטא או שאין מידע שימושי בו, נ fallback לסוג 'IGBP'
            if datasourceDocument and 'desc' in datasourceDocument:
                desc = datasourceDocument['desc']
                type_name = desc.get('type') or (desc.get('desc', {}).get('type') if isinstance(desc, dict) else None)
            else:
                type_name = "IGBP"
                print("[WARNING] No metadata available for data source. Using default type: IGBP.")

            handler_name = f"_handleType{type_name}"

            # fallback דינמי: אם אין handler כזה, נשתמש בזה של IGBP המובנה
            handlerFunction = getattr(self, handler_name, None)
            if handlerFunction is None and type_name == "IGBP":
                handlerFunction = lambda lc_value: {
                    0: 0.01, 1: 0.02, 2: 0.05, 3: 0.1, 4: 0.15,
                    5: 0.2, 6: 0.25, 7: 0.3, 8: 0.35, 9: 0.4,
                    10: 0.45, 11: 0.5, 12: 0.55, 13: 0.6,
                    14: 0.01, 15: 0.001, 16: 0.0001
                }.get(int(lc_value), 0.1)
            elif handlerFunction is None:
                raise AttributeError(f"Handler function '{handler_name}' not found in LandCoverToolkit.")

            # Apply vectorized roughness handler
            roughness_values = np.vectorize(handlerFunction)(landcover['landcover'])
            landcover = landcover.assign_coords(z0=(['i', 'j'], roughness_values))

        else:
            if windMeteorologicalDirection is None or resolution is None:
                raise ValueError(
                    "windMeteorologicalDirection and resolution must be specified for calculating urban area.")
            landcover = self._getUrbanRoughnessFromLandCover(
                landcover,
                windMeteorologicalDirection,
                resolution,
                dataSourceName,
                GIS_BUILDINGS_dataSourceName
            )

        return landcover

    def getRoughness(self, minx, miny, maxx, maxy, dxdy=30, inputCRS=WSG84, dataSourceName=None, isBuilding=False,
                     windMeteorologicalDirection=None, resolution=None, GIS_BUILDINGS_dataSourceName=None):
        """
        Returns Xarray LandCover map with Roughness (zo) field.
        Just as applying getLandCover and getRoughnessFromLandcover at the same time.

        Parameters
        ----------
        minx : float
            Minimum value of X axis. If using WSG84 - the minimum Latitude coordinates.

        miny : float
            Minimum value of Y axis. If using WSG84 - the minimum Longitude coordinates.

        maxx : float
            Maximum value of X axis. If using WSG84 - the maximum Latitude coordinates.

        maxy : float
            Maximum value of Y axis. If using WSG84 - the maximum Longitude coordinates.

        dxdy : int, default=30
            Spatial resolution of the output land cover map.

        inputCRS : int, default=WSG84
            The EPSG of the coordinates.

        dataSourceName : string, default=None
            The name of the data source to use.

        isBuilding : bool, default=False
            Is the landcover contains building area.

        windMeteorologicalDirection: double, default=None
            The meteorological angle. Must be specified only if data includes urbanic area.

        resolution: double, default=None
            The size of the squares. Must be specified only if data includes urbanic area.

        Returns
        -------
        xarray.DataArray
        """
        landcover = self.getLandCover(minx, miny, maxx, maxy, dxdy=dxdy, inputCRS=inputCRS,
                                      dataSourceName=dataSourceName)
        landcover = self.getRoughnessFromLandcover(
            landcover,
            windMeteorologicalDirection=windMeteorologicalDirection,
            resolution=resolution,
            isBuilding=isBuilding,
            dataSourceName=dataSourceName,
            GIS_BUILDINGS_dataSourceName=GIS_BUILDINGS_dataSourceName
        )
        return landcover

    def _handleType1(self, landcover):
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
            0: 0.0001,  # Water
            1: 1,  # Evergreen needleleaf forest
            2: 1,  # Evergreen broadleaf forest
            3: 1,  # Deciduous needleleaf forest
            4: 1,  # Deciduous broadleaf forest
            5: 1,  # Mixed forests
            6: 0.05,  # Closed shrubland
            7: 0.06,  # Open shrublands
            8: 0.05,  # Woody savannas
            9: 0.15,  # Savannas
            10: 0.12,  # Grasslands
            11: 0.3,  # Permanent wetlands
            12: 0.15,  # Croplands
            13: 0.8,  # Urban and built-up
            14: 0.14,  # Cropland/natural vegetation mosaic
            15: 0.001,  # Snow and ice
            16: 0.01  # Barren or sparsely vegetated
        }
        return roughnessDict.get(landcover, 0.05)

    def getCodingMap(self, datasourceName):
        """
        Returns dictionary that maps landcover int value to string of landcover.

        Parameters
        ----------
        datasourceName : str
            Datasource type (Type-1, Type-2 etc.)

        Returns
        -------
        dict
        """
        if datasourceName == 'Type-1':
            return {
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
        return {}

    @staticmethod
    def roughnesslength2sandgrainroughness(rl):
        """
        Converts roughness length to equivalent sand grain roughness.

        Based on:
        Desmond, C. J., Watson, S. J., & Hancock, P. E. (2017).
        Modelling the wind energy resource in complex terrain and atmospheres.
        Numerical simulation and wind tunnel investigation of non-neutral forest canopy flow.
        Journal of wind engineering and industrial aerodynamics, 166, 48-60.
        https://www.sciencedirect.com/science/article/pii/S0167610516300083#bib12

        Equation 5: Equivalent sand grain roughness (m) is z0 * 30

        We can use it for "nutkRoughWallFunction" boundary condition for Ks (sand grain roughness) parameter.
        Cs value can be set as 0.5.

        Parameters
        ----------
        rl : float
            Roughness length.

        Returns
        -------
        float
            Equivalent sand grain roughness (Ks).
        """
        return rl * 30.0  # return Ks value

    def _getUrbanRoughnessFromLandCover(self, landcover, windMeteorologicalDirection, resolution, dataSourceName, GIS_BUILDINGS_dataSourceName):
        """
        Add Roughness for Urban areas to landcover Xarray.
        z0 and dd fields are calculated from LambdaP, LambdaF and HC of each Urban area.

        Parameters
        ----------
        landcover : xarray.DataArray
            Landcover data to update.

        windMeteorologicalDirection : float
            Meteorological wind direction for building analysis.

        resolution : float
            Resolution to calculate Lambda values.

        dataSourceName : str
            Data source name.

        GIS_BUILDINGS_dataSourceName : str
            Buildings data source name.

        Returns
        -------
        xarray.DataArray
            Updated landcover with urban roughness parameters.
        """

        # [FIX] Import geopandas inside if needed
        import geopandas as gpd

        gis_building_tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_BUILDINGS, projectName=self.projectName)

        min_geom = convertCRS(
            points=[[float(min(landcover.lon.values.flatten())), float(min(landcover.lat.values.flatten()))]],
            inputCRS=WSG84, outputCRS=ITM
        )[0]
        max_geom = convertCRS(
            points=[[float(max(landcover.lon.values.flatten())), float(max(landcover.lat.values.flatten()))]],
            inputCRS=WSG84, outputCRS=ITM
        )[0]

        # [FIX] Wrap to GeoDataFrame to ensure access to x,y
        min_pp = gpd.GeoDataFrame(geometry=[min_geom], crs=ITM)
        max_pp = gpd.GeoDataFrame(geometry=[max_geom], crs=ITM)

        buildings = gis_building_tk.getBuildingsFromRectangle(
            minx=min_pp.geometry[0].x,
            miny=min_pp.geometry[0].y,
            maxx=max_pp.geometry[0].x,
            maxy=max_pp.geometry[0].y,
            dataSourceName=GIS_BUILDINGS_dataSourceName,
            inputCRS=ITM
        )

        if len(buildings) == 0:
            raise ValueError("Buildings DataFrame for specified coordinates is empty.")

        lambdaGrid = gis_building_tk.analysis.LambdaFromBuildingData(
            windMeteorologicalDirection,
            resolution,
            buildings
        )
        lambdaGrid.crs = ITM
        lambdaGrid = self._getRoughnessFromBuildingsDataFrame(lambdaGrid)

        square_size = float(landcover.dxdy.values)

        # [FIX] Initialize empty fields with NaN
        for field in ['z0', 'dd', 'lambdaF', 'lambdaP', 'hc', 'll']:
            landcover[field] = (('i', 'j'), np.full((landcover.sizes['i'], landcover.sizes['j']), np.nan))

        # [FIX] Efficient loops using tqdm
        for i in tqdm(range(landcover.sizes['i']), desc="Urban roughness calculation"):
            for j in range(landcover.sizes['j']):
                lon = landcover.lon.values[i, j]
                lat = landcover.lat.values[i, j]

                point_geom = convertCRS(points=[[lon, lat]], inputCRS=WSG84, outputCRS=ITM)[0]
                point = gpd.GeoDataFrame(geometry=[point_geom], crs=ITM)

                lambdas = lambdaGrid.loc[lambdaGrid.geometry.intersects(point.geometry[0])]

                if len(lambdas) == 0:
                    landcover['z0'].values[i, j] = self.getLandCoverAtPoint(lon, lat, WSG84, dataSourceName)
                    landcover['dd'].values[i, j] = 0
                    landcover['lambdaF'].values[i, j] = 0
                    landcover['lambdaP'].values[i, j] = 0
                    landcover['hc'].values[i, j] = 0
                    landcover['ll'].values[i, j] = 0
                else:
                    landcover['z0'].values[i, j] = lambdas['zz0'].values[0]
                    landcover['dd'].values[i, j] = lambdas['dd'].values[0]
                    landcover['lambdaF'].values[i, j] = lambdas['lambdaF'].values[0]
                    landcover['lambdaP'].values[i, j] = lambdas['lambdaP'].values[0]
                    landcover['hc'].values[i, j] = lambdas['hc'].values[0]
                    landcover['ll'].values[i, j] = lambdas['ll'].values[0]

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


    
    def roughnesslength2sandgrainroughness(rl):
    #Desmond, C. J., Watson, S. J., & Hancock, P. E. (2017). Modelling the wind energy resource in complex terrain and atmospheres. Numerical simulation and wind tunnel investigation of non-neutral forest canopy flow. Journal of wind engineering and industrial aerodynamics, 166, 48-60.‏    
    # https://www.sciencedirect.com/science/article/pii/S0167610516300083#bib12
    # eq. 5: Equivalent sand grain roughness (m) is z0*30
    
    # we can you it for "nutkRoughWallFunction" boundary condition for Ks (sand grain roughness) parameter
    # Cs value can be set as 0.5
    
        return rl*30.0 # return Ks value 


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
            0: 'blue',  # Water
            1: 'darkgreen',  # Evergreen needleleaf forest
            2: 'forestgreen',  # Evergreen broadleaf forest
            3: 'darkolivegreen',  # Deciduous needleleaf forest
            4: 'olivedrab',  # Deciduous broadleaf forest
            5: 'green',  # Mixed forests
            6: 'saddlebrown',  # Closed shrubland
            7: 'peru',  # Open shrublands
            8: 'darkkhaki',  # Woody savannas
            9: 'khaki',  # Savannas
            10: 'yellowgreen',  # Grasslands
            11: 'aqua',  # Permanent wetlands
            12: 'gold',  # Croplands
            13: '#2f2f2f',  # Urban and built-up
            14: 'lightgreen',  # Cropland/natural vegetation mosaic
            15: 'white',  # Snow and ice
            16: 'lightgrey'  # Barren or sparsely vegetated
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

        colormap = plt.cm.viridis
        norm = mcolors.Normalize(vmin=landcover.z0.min().values, vmax=landcover.z0.max().values)

        # Add color bar to the right of the plot
        cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=colormap), ax=ax, orientation='vertical')
        cbar.set_label("Roughness Value (z0)", fontsize=12)  # Label for color bar

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
        extent = plot.get_extent()
        ax.imshow(plot.get_array(), extent=extent)
        ax = self._adddRectanglesToPlot(ax, rectangles,alpha)

        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])

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
                shape = convertCRS([[x.lon,x.lat]], inputCRS=WSG84, outputCRS=ITM)[0]
                rectangles.append((shape.x, shape.y, float(landcover.dxdy.values), float(landcover.dxdy.values), self.landcover_colors_map.get(int(x.landcover.values))))

        return list(set(rectangles))

    def _getRoughnessRectangles(self,landcover):
        rectangles = []
        colormap = plt.cm.viridis
        norm = mcolors.Normalize(vmin=landcover.z0.min().values, vmax=landcover.z0.max().values)

        for arr in landcover:
            for x in arr:
                shape = convertCRS([[x.lon,x.lat]], inputCRS=WSG84, outputCRS=ITM)[0]
                color = colormap(norm(x['z0'].values))
                rectangles.append((shape.x, shape.y, float(landcover.dxdy.values), float(landcover.dxdy.values),
                                   color))

        return list(set(rectangles))

    def plotLambdas(self,field,plot,landcover,alpha=0.2,figsize=(28,28)):
        """
        Plot Lmabda upon given Axes.

        Parameters
        ----------
        field: str
            Name of desired lambda to be plotted.

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
        rectangles = self._getLambdasRectangles(field,landcover)
        self._plotWithRectangles(ax, plot, rectangles, alpha)

        colormap = plt.cm.viridis
        norm = mcolors.Normalize(vmin=landcover.z0.min().values, vmax=landcover.z0.max().values)

        # Add color bar to the right of the plot
        cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=colormap), ax=ax, orientation='vertical')
        cbar.set_label(f"{field} Value", fontsize=12)  # Label for color bar

        plt.show()

    def _getLambdasRectangles(self,field,landcover):
        rectangles = []
        colormap = plt.cm.viridis
        norm = mcolors.Normalize(vmin=landcover[field].min().values, vmax=landcover[field].max().values)

        for arr in landcover:
            for x in arr:
                shape = convertCRS([[x.lon,x.lat]], inputCRS=WSG84, outputCRS=ITM)[0]
                color = colormap(norm(x[field].values))
                rectangles.append((shape.x, shape.y, float(landcover.dxdy.values), float(landcover.dxdy.values),
                                   color))

        return list(set(rectangles))
