"""
    Just pieces of code from the vector.topography.

    Should be organized when needed.

"""
from .... import toolkit
from ..utils import stlFactory,convertCRS,ITM,WSG84,ED50_ZONE36N
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

class TopographyToolkit(toolkit.abstractToolkit):

    def __init__(self, projectName, filesDirectory=None):
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
        super().__init__(projectName=projectName, toolkitName = 'TopographyToolkit', filesDirectory=filesDirectory)
        self._analysis = topographyAnalysis(self)

    def getPointElevation(self,lat, long,dataSourceName=None):
        """
        Get the elevation above sea level in meters from DEM for a point SRTM30 means 30 meters resolution.

        Parameters
        ----------
        lat : float
            Latitiute in the WSG projection.

        long : float
            Longitude in the WSG projection.

        dataSourceName: string,default=None
            The name of the datasources loaded. If none, use the default datasource.

        Returns
        -------
            float
        """
        # How to use
        # ----------
        # from hera import toolkitHome
        # tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_RASTER)
        # tk.getPointElevation(lat = 33.459,long = 35.755)
        #
        # This should return ~  ~820 according to amud anan.

        logger = get_classMethod_logger(self,"getPointElevation")

        filename = 'N'+str(int(lat))+'E'+str(int(long)).zfill(3)+'.hgt'
        logger.info(f"Getting elevation from file {filename}")

        if dataSourceName is None:
            dataSourceName = self.getConfig()['defaultSRTM']

        if dataSourceName is None:
            raise ValueError(f"The datasource is not found!")

        logger.debug(f"Getting the data source {dataSourceName}")
        fheight = os.path.join(self.getDataSourceData(dataSourceName),filename)

        ds =  gdal.Open(fheight)
        myarray = numpy.array(ds.GetRasterBand(1).ReadAsArray())
        myarray[myarray < -1000] = 0  # deal with problematic points
        gt = ds.GetGeoTransform()
        rastery = (long - gt[0]) / gt[1]
        rasterx = (lat - gt[3]) / gt[5]
        height11 = myarray[int(rasterx), int(rastery)]
        if (int(rasterx)+1>=myarray.shape[0]) or (int(rastery)+1>=myarray.shape[1]):
            height12 = height11
            height21 = height11
            height22 = height11
        else:
            height12 = myarray[int(rasterx) + 1, int(rastery)]
            height21 = myarray[int(rasterx), int(rastery) + 1]
            height22 = myarray[int(rasterx) + 1, int(rastery) + 1]
        height1 = (1. - (rasterx - int(rasterx))) * height11 + (rasterx - int(rasterx)) * height12
        height2 = (1. - (rasterx - int(rasterx))) * height21 + (rasterx - int(rasterx)) * height22
        elevation = (1. - (rastery - int(rastery))) * height1 + (rastery - int(rastery)) * height2

        return elevation

    def getPointListElevation(self, pointList, dataSourceName=None,inputCRS=WSG84):
        """
        Return the elevation of the point list.
        For now we assume that pointList is in WSG84 coordinates.

        Parameters
        ----------
        pointList

        dataSourceName

        Returns
        -------

        """
        totalNumberOfPoints = pointList.shape[0]
        logger = get_classMethod_logger(self, "getPointListElevation")
        logger.info(f"Computing the elevation for a list of {totalNumberOfPoints} points")
        logger.debug("Computing the file name for each point")


        if isinstance(pointList,geopandas.geoseries.GeoSeries):
            pointList = pointList.to_frame("point").assign(filename=pointList.apply(lambda row: 'N' + str(int(row.y)) + 'E' + str(int(row.x)).zfill(3) + '.hgt' ),
                                                           lon=pointList.x,
                                                           lat=pointList.y,
                                                           elevation=0)

        elif isinstance(pointList,geopandas.geodataframe.GeoDataFrame):
            if 'point' not in pointList:
                error = "GeoDataFrame must contain a field 'point' that contain the points of interest"
                logger.error(error)
                raise ValueError(error)
            pointList = pointList.assign(filename=pointList.apply(lambda row: 'N' + str(int(row.point.y)) + 'E' + str(int(row.point.x)).zfill(3) + '.hgt', axis=1),
                                         lon=pointList.point.x,
                                         lat=pointList.point.y,
                                         elevation=0)
        else:
            pointList = pointList.assign(filename=pointList.apply(lambda x: 'N' + str(int(x.lat)) + 'E' + str(int(x.lon)).zfill(3) + '.hgt' ,axis=1),elevation=0)


        if dataSourceName is None:
            dataSourceName = self.getConfig()['defaultSRTM']
        logger.info(f"Getting the elevation for the points. Using datasource {dataSourceName}")

        if dataSourceName is None:
            err = "The datasource is not found!"
            logger.error(err)
            raise ValueError(err)

        processed = 0

        for grpid,grp in pointList.groupby("filename"):
            logger.info(f"\tProcessing the group {grpid} with {grp.shape[0]}. Processed so far {processed}/{totalNumberOfPoints}")
            fheight = os.path.join(self.getDataSourceData(dataSourceName), grpid)
            logger.debug(f"Getting data from the file {fheight}")
            ds = gdal.Open(fheight)
            myarray = numpy.array(ds.GetRasterBand(1).ReadAsArray())
            myarray[myarray < -1000] = 0  # deal with problematic points
            gt = ds.GetGeoTransform()
            for itemindx,item in grp.iterrows():
                rastery = (item.lon - gt[0]) / gt[1]
                rasterx = (item.lat - gt[3]) / gt[5]
                height11 = myarray[int(rasterx), int(rastery)]
                if (int(rasterx)+1>=myarray.shape[0]) or (int(rastery)+1>=myarray.shape[1]):
                    height12 = height11
                    height21 = height11
                    height22 = height11
                else:
                    height12 = myarray[int(rasterx) + 1, int(rastery)]
                    height21 = myarray[int(rasterx), int(rastery) + 1]
                    height22 = myarray[int(rasterx) + 1, int(rastery) + 1]
                height1 = (1. - (rasterx - int(rasterx))) * height11 + (rasterx - int(rasterx)) * height12
                height2 = (1. - (rasterx - int(rasterx))) * height21 + (rasterx - int(rasterx)) * height22
                elevation = (1. - (rastery - int(rastery))) * height1 + (rastery - int(rastery)) * height2
                pointList.loc[itemindx,'elevation'] = elevation
                processed += grp.shape[0]
        return pointList

    def getDomainElevation(self, xmin, ymin, xmax, ymax, dxdy = 30, inputCRS=WSG84, outputCRS=ITM, dataSourceName=None):
        """
            Returns the eleveation between the lowerleft point and the upper right points, at the requested resolution.

            Assumes the coordinates are given in the input CRS, and outputted in the output CRS.

            Note, that in order to get 30m resultion, we convert the coordinates to ITM and then convert them back to
            WSG84. This will work well in IL.

        Parameters
        ----------
        lowerleft_point : float
                The lower left corner
        upperright_point: float
                The upper right corner

        dxdy : float
                The resolution in m (default m).
        inputCRS : The ESPG code of the input projection.

        outputCRS : The ESPG code of the output projection.
                    [Default ITM]
        Returns
        -------
            xarrat with the fields
            - X : X[i,j] - the x coordinate
            - Y:  Y[i,j] - the y coordinate
            - Height: Height[i,j] - the height of the surface at point  (X[i,j],Y[i,j])
        """
        logger = get_classMethod_logger(self,"getDomainElevation")

        logger.info("Converting input coordinates to ITM (to have meters)")
        itm_points = convertCRS([[xmin, ymin], [xmax, ymax]], inputCRS=inputCRS, outputCRS=ITM)
        logger.debug(f"... Got {itm_points}")
        logger.info("Getting the coordinates")
        itm_lowerleft, itm_upperright = itm_points

        itm_gridx = numpy.arange(itm_lowerleft.x,itm_upperright.x,dxdy)
        itm_gridy = numpy.arange(itm_lowerleft.y,itm_upperright.y,dxdy)

        logger.info("Converting the coordinates to WSG84 for SRTM data")
        wsg_coords = convertCRS([x for x in product(itm_gridx,itm_gridy)], inputCRS=ITM, outputCRS=WSG84)

        logger.info("Getting output coords")
        if outputCRS != ITM:
            output_coords = convertCRS([x for x in product(itm_gridx,itm_gridy)], inputCRS=ITM, outputCRS=outputCRS)
        else:
            XX,YY = zip(*[x for x in product(itm_gridx,itm_gridy)])
            output_coords = geopandas.points_from_xy(XX,YY)

        grid_x = numpy.zeros([len(itm_gridx), len(itm_gridy)])
        grid_y = numpy.zeros_like(grid_x)
        grid_z = numpy.zeros_like(grid_x)

        pointList = pandas.DataFrame([dict(lon=pnt.x,lat=pnt.y,lon_output=pnt_output.x,lat_output=pnt_output.y) for pnt,pnt_output in zip(wsg_coords,output_coords)])
        pointHeights = self.getPointListElevation(pointList) #.set_index(["lon","lat"])

        logger.info("Converting the output to the desired coordinates.")
        grid_y = pointHeights['lat_output'].values.reshape(len(itm_gridx), len(itm_gridy))
        grid_x = pointHeights['lon_output'].values.reshape(len(itm_gridx), len(itm_gridy))
        grid_z = pointHeights['elevation'].values.reshape(len(itm_gridx), len(itm_gridy))

        retArray =  xarray.Dataset(data_vars=dict(X=(["i","j"],grid_x),
                                                  Y=(["i","j"],grid_y),
                                                  Elevation=(["i", "j"], grid_z)),coords=dict(i=("i",numpy.arange(grid_y.shape[0])),j=("j",numpy.arange(grid_y.shape[1]) )))

        return retArray

    def getDomainElevation_STL(self, minx, miny, maxx, maxy, dxdy = 30, inputCRS=WSG84, outputCRS=ITM, dataSourceName=None, solidName="Topography", shiftx=0, shifty=0):
        """
        Return the STL string from xarray dataset with the following fields.

        Parameters
        ----------
        lowerleft_point : float
                The lower left corner.

        upperright_point: float
                The upper right corner.

        dxdy : float,default=30
                The resolution in m (default m).

        inputCRS : int,default=4326
            The ESPG code of the input projection. Default is WSG84.

        outputCRS : int, default=2039
                The ESPG code of the output projection. Default is ITM.
                    
        shiftx : float,default=0
            Used when one wants to set another point as origin center
                    
        shifty : float,default=0
            Used when one wants to set another point as origin center

        Returns
        -------
            str
        """
        elevation = self.getDomainElevation(xmin=minx, xmax=maxx, ymax=maxy, ymin=miny, dxdy=dxdy, inputCRS=inputCRS, outputCRS=outputCRS, dataSourceName=dataSourceName)
        stlstr = stlFactory().rasterToSTL(elevation['X'].values - shiftx, elevation['Y'].values - shifty , elevation['Elevation'].values,solidName=solidName)
        return stlstr

    # def setDomainElevation(point1, point2, outputfile, amplitude=100., cosx=0., cosy=0., projectName = "default"):
    #     """
    #         create STL from equation
    #
    #     Parameters
    #     ----------
    #     point1 (lat, long)
    #     point2 (lat, long)
    #     outputfile - where to write the stl
    #     amplitude - of the wave [m]
    #     cosx - x axis frequency, optional
    #     cosy - y axis frequency, optional
    #
    #     Returns
    #     -------
    #
    #     How to use
    #     ----------
    #     from hera.measurements.GIS.raster.topography import TopographyToolkit
    #     TopographyToolkit.setDomainElevation([200000, 740000], [201000, 741000], 'test1.stl', amplitude=10, cosx=.5, cosy=.2)
    #
    #
    #     What to fix
    #     -----------
    #
    #     """
    #     miny = point1[1]  # y
    #     maxy = point2[1]
    #     minx = point1[0]  # x
    #     maxx = point2[0]
    #     dxdy = 30
    #
    #     NX = 1 + math.ceil((maxx - minx) / dxdy)
    #     NY = 1 + math.ceil((maxy - miny) / dxdy)
    #     grid_x = numpy.zeros([NX, NY])
    #     grid_y = numpy.zeros_like(grid_x)
    #     grid_z = numpy.zeros_like(grid_x)
    #     i = minx
    #     j = miny
    #     for i in range(NX):
    #         for j in range(NY):
    #             x = minx + i * dxdy
    #             y = miny + j * dxdy
    #             z = amplitude*(2.+math.cos(math.pi+i*cosx))*(2.+math.cos(math.pi+j*cosy))
    #             grid_x[i, j] = x
    #             grid_y[i, j] = y
    #             grid_z[i, j] = z
    #
    #     a = stlFactory()
    #     stlstr = a.rasterToSTL(grid_x, grid_y, grid_z, 'solidName')
    #
    #     with open(outputfile, "w") as stlfile:
    #         stlfile.write(stlstr)
    #
    #     return
    #
    #
    # def latlongtoITM(point1):
    #     """
    #         change the coordinate system, from lat long to x y
    #
    #     Parameters
    #     ----------
    #     lat
    #     long
    #
    #     Returns
    #     -------
    #     x, y in ITM
    #
    #     How to use
    #     ----------
    #     from hera.measurements.GIS.raster.topography import TopographyToolkit
    #     TopographyToolkit.latlongtoITM([35.159,32.755])
    #
    #     What to fix
    #     -----------
    #     """
    #     ITM = "epsg:2039"
    #     WGS84 = "epsg:4327"
    #     lat = point1[0]
    #     long = point1[1]
    #     coordTransformer = Transformer.from_crs(WGS84, ITM, always_xy=True)
    #     x, y = coordTransformer.transform(lat, long)
    #     return x,y




# srtmDocs = self.getMeasurementsDocumentsAsDict(name="SRTM")
# self._srtmBounds = []
# if type(srtmDocs) == dict:
#     for doc in srtmDocs["documents"]:
#         self._srtmBounds.append(doc["desc"]["Bounds"])

#
# def get_altitdue_ip(lat, lon):
#     """
#     returning the altitude of the point, uses free mapquest data that is limited in the amount of calls per month, it uses Nir BAMBA Benmoshe key
#
#     param lat - the path where we save the stl
#     param lon - the width of the domain in the x direction
#
#     return:
#     altitude - meters above sea level
#     """
#
#     resp = requests.get(
#         'http://open.mapquestapi.com/elevation/v1/profile?key=D5z9RSebQJLbUs4bohANIB4TzJdbvyvm&shapeFormat=raw&latLngCollection=' + str(
#             lat) + ',' + str(lon))
#
#     height = resp.json()['elevationProfile'][0]['height']
#
#     return height
# def get_altitdue_gdal(lat, lon):
#     #        if lat<=30:
#     # USGS EROS Archive - Digital Elevation - Global Multi-resolution Terrain Elevation Data 2010 (GMTED2010)
#     # fheight = r'/data3/nirb/10N030E_20101117_gmted_med075.tif'
#     # fheight = r'/data3/nirb/30N030E_20101117_gmted_med075.tif'
#     # https://dds.cr.usgs.gov/srtm/version2_1/SRTM3/Africa/   # 90m resolution
#     if lat > 29 and lat < 30 and lon > 34 and lon < 35:
#         fheight = r'/data3/nirb/N29E034.hgt'
#     elif lat > 29 and lat < 30 and lon > 35 and lon < 36:
#         fheight = r'/data3/nirb/N29E035.hgt'
#     elif lat > 30 and lat < 31 and lon > 34 and lon < 35:
#         fheight = r'/data3/nirb/N30E034.hgt'
#     elif lat > 30 and lat < 31 and lon > 35 and lon < 36:
#         fheight = r'/data3/nirb/N30E035.hgt'
#     elif lat > 31 and lat < 32 and lon > 34 and lon < 35:
#         fheight = r'/data3/nirb/N31E034.hgt'
#     elif lat > 31 and lat < 32 and lon > 35 and lon < 36:
#         fheight = r'/data3/nirb/N31E035.hgt'
#     elif lat > 32 and lat < 33 and lon > 34 and lon < 35:
#         fheight = r'/data3/nirb/N32E034.hgt'
#     elif lat > 32 and lat < 33 and lon > 35 and lon < 36:
#         fheight = r'/data3/nirb/N32E035.hgt'
#     elif lat > 33 and lat < 33 and lon > 35 and lon < 36:
#         fheight = r'/data3/nirb/N33E035.hgt'
#     else:
#         print('!!!!NOT in Israel !!!!!!!!')
#         # taken from https://earthexplorer.usgs.gov/
#         fheight = r'/ibdata2/nirb/gt30e020n40.tif'
#
#     ds = gdal.Open(fheight)
#     myarray = numpy.array(ds.GetRasterBand(1).ReadAsArray())
#     myarray[myarray < -1000] = 0
#     gt = ds.GetGeoTransform()
#     rastery = (lon - gt[0]) / gt[1]
#     rasterx = (lat - gt[3]) / gt[5]
#     height11 = myarray[int(rasterx), int(rastery)]
#     height12 = myarray[int(rasterx) + 1, int(rastery)]
#     height21 = myarray[int(rasterx), int(rastery) + 1]
#     height22 = myarray[int(rasterx) + 1, int(rastery) + 1]
#     height1 = (1. - (rasterx - int(rasterx))) * height11 + (rasterx - int(rasterx)) * height12
#     height2 = (1. - (rasterx - int(rasterx))) * height21 + (rasterx - int(rasterx)) * height22
#     height = (1. - (rastery - int(rastery))) * height1 + (rastery - int(rastery)) * height2
#
#     return height
#
#


# def makeRegion_tif(self, points, outputFileName, **kwargs):
#     availableBounds = {}
#     for bounds in self.srtmBounds:
#         if points[0] >= bounds[0] and points[0] <= bounds[2] and points[2] >= bounds[0] and points[2] <= bounds[2] and \
#                 points[1] >= bounds[1] and points[1] <= bounds[3] and points[3] >= bounds[1] and points[3] <= bounds[3]:
#             availableBounds["Bounds"] = bounds
#             break
#     if len(availableBounds) > 0:
#         FileName = f"{outputFileName}.parquet"
#         allData = self.getMeasurementsDocuments(name="SRTM", type=toolkit.TOOLKIT_DATASOURCE_TYPE, **availableBounds)[
#             0].getData()
#         dataArray = allData.read(1)
#         xs = numpy.linspace(allData.bounds.left, allData.bounds.right, dataArray.shape[1])
#         ys = numpy.linspace(allData.bounds.top, allData.bounds.bottom, dataArray.shape[0])
#         xmin = numpy.where(xs > points[0])[0].min()
#         ymin = numpy.where(ys < points[1])[0].min()
#         xmax = numpy.where(xs < points[2])[0].max()
#         ymax = numpy.where(ys > points[3])[0].max()
#         xColumn = self.getConfig()["xColumn"] if "xColumn" in self.getConfig().keys() else "x"
#         yColumn = self.getConfig()["yColumn"] if "yColumn" in self.getConfig().keys() else "y"
#         heightColumn = self.getConfig()["heightColumn"] if "heightColumn" in self.getConfig().keys() else "height"
#         data = xarray.DataArray(data=dataArray[ymax:ymin, xmin:xmax], dims=["yWGS84", "xWGS84"],
#                                 coords=[ys[ymax:ymin], xs[xmin:xmax]]).to_dataframe(heightColumn).reset_index()
#         gdf = geopandas.GeoDataFrame(data, geometry=geopandas.points_from_xy(data["xWGS84"], data["yWGS84"]))
#         gdf.crs = {"init": "epsg:4326"}
#         gdf = gdf.to_crs({"init": "epsg:2039"})
#         gdf[xColumn] = gdf.geometry.x
#         gdf[yColumn] = gdf.geometry.y
#         data = dask.dataframe.from_pandas(gdf.drop(columns="geometry"), npartitions=1)
#         data.to_parquet(FileName, compression="GZIP")
#     else:
#         raise KeyError("Couldn't find SRTM files which include the requested points!")
#     return
#
#
# def getHeight(self, latitude, longitude, source):
#     return getattr(self, f"getHeight_{source}")(longitude=longitude, latitude=latitude)
#
#
# def getHeight_mapquest(self, latitude, longitude):
#     resp = requests.get(
#         'http://open.mapquestapi.com/elevation/v1/profile?key=D5z9RSebQJLbUs4bohANIB4TzJdbvyvm&shapeFormat=raw&latLngCollection=' + str(
#             latitude) + ',' + str(longitude))
#
#     height = resp.json()['elevationProfile'][0]['height']
#
#     return height
#
#
# def getHeight_USGS(self, latitude, longitude):
#     path = self.getMeasurementsDocuments(name="USGS")[0].getData()
#
#     if latitude > 29 and latitude < 30 and longitude > 34 and longitude < 35:
#         fheight = r'%s/N29E034.hgt' % path
#     elif latitude > 29 and latitude < 30 and longitude > 35 and longitude < 36:
#         fheight = r'%s/N29E035.hgt' % path
#     elif latitude > 30 and latitude < 31 and longitude > 34 and longitude < 35:
#         fheight = r'%s/N30E034.hgt' % path
#     elif latitude > 30 and latitude < 31 and longitude > 35 and longitude < 36:
#         fheight = r'%s/N30E035.hgt' % path
#     elif latitude > 31 and latitude < 32 and longitude > 34 and longitude < 35:
#         fheight = r'%s/N31E034.hgt' % path
#     elif latitude > 31 and latitude < 32 and longitude > 35 and longitude < 36:
#         fheight = r'%s/N31E035.hgt' % path
#     elif latitude > 32 and latitude < 33 and longitude > 34 and longitude < 35:
#         fheight = r'%s/N32E034.hgt' % path
#     elif latitude > 32 and latitude < 33 and longitude > 35 and longitude < 36:
#         fheight = r'%s/N32E035.hgt' % path
#     elif latitude > 33 and latitude < 33 and longitude > 35 and longitude < 36:
#         fheight = r'%s/N33E035.hgt' % path
#     else:
#         raise KeyError("The point is outside Isreal.")
#
#     ds = gdal.Open(fheight)
#     myarray = numpy.array(ds.GetRasterBand(1).ReadAsArray())
#     myarray[myarray < -1000] = 0
#     gt = ds.GetGeoTransform()
#     rastery = (longitude - gt[0]) / gt[1]
#     rasterx = (latitude - gt[3]) / gt[5]
#     height11 = myarray[int(rasterx), int(rastery)]
#     height12 = myarray[int(rasterx) + 1, int(rastery)]
#     height21 = myarray[int(rasterx), int(rastery) + 1]
#     height22 = myarray[int(rasterx) + 1, int(rastery) + 1]
#     height1 = (1. - (rasterx - int(rasterx))) * height11 + (rasterx - int(rasterx)) * height12
#     height2 = (1. - (rasterx - int(rasterx))) * height21 + (rasterx - int(rasterx)) * height22
#     height = (1. - (rastery - int(rastery))) * height1 + (rastery - int(rastery)) * height2
#
#     return height



class topographyAnalysis:

    datalayer = None

    def __init__(self,datalayer):
        self.datalayer = datalayer

    def calculateStastics(self,elevation):
        """
        Calculate the domain elevation statistics. Returns domain size (in meters for ITM, assuming small domain (rectangle),
        max elevation + location, min elevation + location, mean elevation and std of elevation.

        Parameters
        ----------
        elevation: pd.DataFrame
            elevation data including coordinates

        Returns
        -------
            dict
        """
        # example
        # -------
        # after creating a project (hera-project project create project-name)
        # from hera import toolkitHome
        # tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_RASTER_TOPOGRAPHY,projectName="topotest")
        # elevation = tk.getDomainElevation(xmin=34.529,xmax=34.531,ymin=31.160,ymax=31.162)
        # print(tk.analysis.calculateStastics(elevation))

        xmin = elevation['X'].values.min()
        xmax = elevation['X'].values.max()
        ymin = elevation['Y'].values.min()
        ymax = elevation['Y'].values.max()
        domainsize = (xmax - xmin) * (ymax - ymin)
        domainmean = elevation['Elevation'].values.mean()
        domainstd = elevation['Elevation'].values.std()
        domainmax = elevation['Elevation'].values.max()
        maxpos = numpy.argmax(elevation['Elevation'].values.ravel())
        domainmaxlocation = (elevation['X'].values.ravel()[maxpos],elevation['Y'].values.ravel()[maxpos])
        domainmin = elevation['Elevation'].values.min()
        minpos = numpy.argmin(elevation['Elevation'].values.ravel())
        domainminlocation = (elevation['X'].values.ravel()[minpos],elevation['Y'].values.ravel()[minpos])
        elevation_statistics = {
                "xmin": xmin,
                "xmax": xmax,
                "ymin": ymin,
                "ymax": ymax,
                "size": domainsize,
                "mean": domainmean,
                "std": domainstd,
                "domainmax": domainmax,
                "domainmaxlocation": domainmaxlocation,
                "domainmin": domainmin,
                "domainminlocation": domainminlocation,
                }
        return elevation_statistics
