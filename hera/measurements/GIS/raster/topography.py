"""
    Just pieces of code from the vector.topography.

    Should be organized when needed.

"""
from .... import toolkit
#import geopandas
#from shapely.geometry import Polygon, box
#from hera.datalayer.datahandler import datatypes

from hera import toolkitHome
import numpy as np
import math
from osgeo import gdal

class TopographyToolkit(toolkit.abstractToolkit):


    def __init__(self, projectName, toolkitName = 'TopographyToolkit', filesDirectory=None):
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


    def getPointElevation(lat, long, projectName = "default"):
        """
            get the elevation above sea level in meters from DEM for a point
            SRTM30 means 30 meters resolution

        Parameters
        ----------
        lat
        long
        source - can be SRTM30, SRTM90

        Returns
        -------
        eleveation above sea level

        How to use
        ----------
        from hera.measurements.GIS.raster.topography import TopographyToolkit
        TopographyToolkit.getPointElevation(33.459,35.755)

        What to fix
        -----------
        We need to make sure that we have "default" projectName
        We need to give warning if outside Israel
        We need to be sure that we use SRTM30 and not SRTM90

        """
        # lat = 33.459 # elevation should be ~820 according to amud anan
        # long = 35.755 # elevation should be ~820 according to amud anan
        filename = 'N'+str(int(lat))+'E'+str(int(long)).zfill(3)+'.hgt'

        tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_RASTER_TOPOGRAPHY,projectName=projectName)
        doc = tk.getDatasourceDocumentsList()[0]
        fheight = doc.resource + '/' + filename

        ds = gdal.Open(fheight)
        myarray = np.array(ds.GetRasterBand(1).ReadAsArray())
        myarray[myarray < -1000] = 0  # deal with problematic points
        gt = ds.GetGeoTransform()
        rastery = (long - gt[0]) / gt[1]
        rasterx = (lat - gt[3]) / gt[5]
        height11 = myarray[int(rasterx), int(rastery)]
        height12 = myarray[int(rasterx) + 1, int(rastery)]
        height21 = myarray[int(rasterx), int(rastery) + 1]
        height22 = myarray[int(rasterx) + 1, int(rastery) + 1]
        height1 = (1. - (rasterx - int(rasterx))) * height11 + (rasterx - int(rasterx)) * height12
        height2 = (1. - (rasterx - int(rasterx))) * height21 + (rasterx - int(rasterx)) * height22
        elevation = (1. - (rastery - int(rastery))) * height1 + (rastery - int(rastery)) * height2

        return elevation



    def getDomainElevation(point1, point2, projectName = "default"):
        """
            get the elevation above sea level in meters from DEM for a domain
            SRTM30 means 30 meters resolution

        Parameters
        ----------
        point1 (lat, long)
        point2 (lat, long)
        source - can be SRTM30, SRTM90

        Returns
        -------
        eleveation above sea level

        How to use
        ----------
        from hera.measurements.GIS.raster.topography import TopographyToolkit
        TopographyToolkit.getPointElevation(33.459,35.755)

        What to fix
        -----------
        We need to make sure that we have "default" projectName
        We need to give warning if outside Israel
        We need to be sure that we use SRTM30 and not SRTM90

        """
        # lat = 33.459 # elevation should be ~820 according to amud anan
        # long = 35.755 # elevation should be ~820 according to amud anan

        minlat = min(point1[0], point2[0])
        maxlat = max(point1[0], point2[0])
        minlong = min(point1[1], point2[1])
        maxlong = max(point1[1], point2[1])
        if (int(maxlong)-int(minlong)>0) or (int(maxlat)-int(minlat)>0):
            assert ("This domain is too big for us, we have to change the code to deal with it")
        heightmap=[] # append([lat long height])
        print ('bounding',minlat, maxlat, minlong, maxlong)

        filename = 'N'+str(int(minlat))+'E'+str(int(minlong)).zfill(3)+'.hgt'

        tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_RASTER_TOPOGRAPHY,projectName=projectName)
        doc = tk.getDatasourceDocumentsList()[0]
        fheight = doc.resource + '/'+ filename

        ds = gdal.Open(fheight)

        if int(maxlong) - int(minlong) == 1:
            matlong = 2
        else:
            matlong = 1
        if int(maxlat) - int(minlat) == 1:
            matlat = 2
        else:
            matlat = 1
        mat=np.zeros([np.array(ds.GetRasterBand(1).ReadAsArray()).shape[0]*matlong,np.array(ds.GetRasterBand(1).ReadAsArray()).shape[0]*matlat])
        if matlong==1 and matlat==1:
            filename = 'N' + str(int(minlat)) + 'E' + str(int(minlong)).zfill(
                3) + '.hgt'  # filename=r'/data3/GIS_Data/Elevation/SRTM30m/N33E035.hgt'
            fheight = doc.resource + '/' + filename
            ds = gdal.Open(fheight)
            myarray = np.array(ds.GetRasterBand(1).ReadAsArray())
            myarray[myarray < -1000] = 0  # deal with problematic points
            gt = ds.GetGeoTransform()
            mat[:,:]=myarray
            gt0 = gt[0]
            gt3 = gt[3]
            gt1 = gt[1]
            gt5 = gt[5]
        gt = ds.GetGeoTransform()
        rasterymin = math.floor((minlong - gt0) / gt1)
        rasterymax = math.ceil((maxlong - gt0) / gt1)
        rasterxmax = math.floor((minlat - gt3) / gt5)
        rasterxmin = math.ceil((maxlat - gt3) / gt5)
        points = []
        print('gt',gt)
        print('rasterx',rasterxmin, rasterxmax,rasterymin, rasterymax)
        for y in range(rasterymin,rasterymax):
            for x in range(rasterxmin, rasterxmax):
                xlat = x * gt[5] + gt[3]
                ylong = y * gt[1] + gt[0]
                points.append([xlat, ylong, myarray[x, y]])
        return points


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
