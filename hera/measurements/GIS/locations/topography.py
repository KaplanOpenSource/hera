import logging
from .abstractLocation import datalayer as locationDatalayer
import geopandas
from scipy.interpolate import griddata
import numpy
import pandas
import random
import gdal
import numpy as np
from .shapes import datalayer as shapeDatalayer
from shapely.geometry import Point,box,MultiLineString, LineString
import requests
from unum.units import *
import scipy
from ..topography.STL import stlFactory
from ....simulations.utils import coordinateHandler
import os
import dask

class datalayer(locationDatalayer):

    _analysis = None

    @property
    def analysis(self):
        return self._analysis

    def __init__(self, projectName, FilesDirectory="", databaseNameList=None, useAll=False,publicProjectName="Topography",
                 source="BNTL", dxdy=50, heightSource="USGS", xColumn="x",yColumn="y",heightColumn="height", dbName=None,**kwargs):

        super().__init__(projectName=projectName,publicProjectName=publicProjectName,FilesDirectory=FilesDirectory,databaseNameList=databaseNameList,useAll=useAll)
        self._analysis = analysis(projectName=projectName, dataLayer=self)
        documents = self.getCacheDocuments(type="__config__")
        if len(documents)==0 or len(kwargs.keys())>0:
            self.setConfig(dbName=dbName,source=source,dxdy=dxdy,heightSource=heightSource,xColumn=xColumn,yColumn=yColumn,heightColumn=heightColumn,**kwargs)

    def getHeight(self, latitude, longitude):

        try:
            height = self.__getattribute__("getHeight_%s" % self.getConfig()["heightSource"])(longitude=longitude,latitude=latitude)
        except AttributeError:
            raise KeyError("The height source is not known.")

        return height

    def getHeight_mapquest(self, latitude, longitude):
        resp = requests.get(
            'http://open.mapquestapi.com/elevation/v1/profile?key=D5z9RSebQJLbUs4bohANIB4TzJdbvyvm&shapeFormat=raw&latLngCollection=' + str(
                latitude) + ',' + str(longitude))

        height = resp.json()['elevationProfile'][0]['height']

        return height

    def getHeight_USGS(self, latitude, longitude):

        path = self.getMeasurementsDocuments(type="Height",**self.getConfig())[0].getData()

        if latitude > 29 and latitude < 30 and longitude > 34 and longitude < 35:
            fheight = r'%s/N29E034.hgt' % path
        elif latitude > 29 and latitude < 30 and longitude > 35 and longitude < 36:
            fheight = r'%s/N29E035.hgt'% path
        elif latitude > 30 and latitude < 31 and longitude > 34 and longitude < 35:
            fheight = r'%s/N30E034.hgt'% path
        elif latitude > 30 and latitude < 31 and longitude > 35 and longitude < 36:
            fheight = r'%s/N30E035.hgt'% path
        elif latitude > 31 and latitude < 32 and longitude > 34 and longitude < 35:
            fheight = r'%s/N31E034.hgt'% path
        elif latitude > 31 and latitude < 32 and longitude > 35 and longitude < 36:
            fheight = r'%s/N31E035.hgt'% path
        elif latitude > 32 and latitude < 33 and longitude > 34 and longitude < 35:
            fheight = r'%s/N32E034.hgt'% path
        elif latitude > 32 and latitude < 33 and longitude > 35 and longitude < 36:
            fheight = r'%s/N32E035.hgt'% path
        elif latitude > 33 and latitude < 33 and longitude > 35 and longitude < 36:
            fheight = r'%s/N33E035.hgt'% path
        else:
            raise KeyError("The point is outside Isreal.")
            # print('!!!!NOT in Israel !!!!!!!!')
            # # taken from https://earthexplorer.usgs.gov/
            # fheight = r'/ibdata2/nirb/gt30e020n40.tif'

        ds = gdal.Open(fheight)
        myarray = np.array(ds.GetRasterBand(1).ReadAsArray())
        myarray[myarray < -1000] = 0
        gt = ds.GetGeoTransform()
        rastery = (longitude - gt[0]) / gt[1]
        rasterx = (latitude - gt[3]) / gt[5]
        height11 = myarray[int(rasterx), int(rastery)]
        height12 = myarray[int(rasterx) + 1, int(rastery)]
        height21 = myarray[int(rasterx), int(rastery) + 1]
        height22 = myarray[int(rasterx) + 1, int(rastery) + 1]
        height1 = (1. - (rasterx - int(rasterx))) * height11 + (rasterx - int(rasterx)) * height12
        height2 = (1. - (rasterx - int(rasterx))) * height21 + (rasterx - int(rasterx)) * height22
        height = (1. - (rastery - int(rastery))) * height1 + (rastery - int(rastery)) * height2

        return height

class analysis():

    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, projectName, dataLayer=None, FilesDirectory="", databaseNameList=None, useAll=False,
                 publicProjectName="Topography", Source="BNTL"):

        self._datalayer = datalayer(projectName=projectName, FilesDirectory=FilesDirectory, publicProjectName=publicProjectName,
                         databaseNameList=databaseNameList, useAll=useAll, Source=Source) if datalayer is None else dataLayer

    def PolygonDataFrameIntersection(self, dataframe, polygon):
        """
        Creates a new dataframe based on the intersection of a dataframe and a polygon.
        Parameters:
        ----------
        dataframe: A geopandas dataframe.
        polygon: A shapely polygon

        Returns: A new geopandas dataframe

        """

        newlines = []
        for line in dataframe["geometry"]:
            newline = polygon.intersection(line)
            newlines.append(newline)
        dataframe["geometry"] = newlines
        dataframe = dataframe[~dataframe["geometry"].is_empty]

        return dataframe

    def toSTL(self, data, NewFileName, save=True, addToDB=True, flat=None, path=None, **kwargs):

        """
        Converts a geopandas dataframe data to an stl file.

        Parameters:

            data: The data that should be converted to stl. May be a dataframe or a name of a saved polygon in the database.
            NewFileName: A name for the new stl file, also used in the stl string. (string)
            dxdy: the dimention of each cell in the mesh in meters, the default is 50.
            save: Default is True. If True, the new stl string is saved as a file and the path to the file is added to the database.
            flat: Default is None. Else, it assumes that the area is flat and the value of flat is the height of the mesh cells.
            path: Default is None. Then, the path in which the data is saved is the given self.FilesDirectory. Else, the path is path. (string)
            kwargs: Any additional metadata to be added to the new document in the database.

        Returns
        -------

        """

        if type(data) == str:
            polygon = shapeDatalayer(projectName=self.datalayer.projectName).getShape(data)
            dataframe = self.datalayer.getDocuments(Shape=data, ShapeMode="contains")[0].getData()
            geodata = self.PolygonDataFrameIntersection(polygon=polygon, dataframe=dataframe)
        elif type(data) == geopandas.geodataframe.GeoDataFrame:
            geodata = data
        else:
            raise KeyError("data should be geopandas dataframe or a polygon.")
        xmin = geodata['geometry'].bounds['minx'].min()
        xmax = geodata['geometry'].bounds['maxx'].max()

        ymin = geodata['geometry'].bounds['miny'].min()
        ymax = geodata['geometry'].bounds['maxy'].max()
        points = [xmin, ymin, xmax, ymax]
        documents = self.datalayer.getMeasurementsDocuments(type="stlFile", bounds=points, dxdy=self._datalayer.getConfig()["dxdy"])
        if len(documents) >0:
            stlstr = documents[0].getData()
            newdict = documents[0].asDict()
            newdata = pandas.DataFrame(dict(gridxMin=[newdict["desc"]["xMin"]], gridxMax=[newdict["desc"]["xMax"]],
                                            gridyMin=[newdict["desc"]["yMin"]], gridyMax=[newdict["desc"]["yMax"]],
                                            gridzMin=[newdict["desc"]["zMin"]], gridzMax=[newdict["desc"]["zMax"]]))
        else:
            stl = stlFactory()
            stlstr, newdata = stl.Convert_geopandas_to_stl(gpandas=geodata, points=points, flat=flat, NewFileName=NewFileName,dxdy=self._datalayer.getConfig()["dxdy"])
            if save:
                p = self.datalayer.FilesDirectory if path is None else path
                new_file_path = os.path.join(p,f"{NewFileName}.stl")
                with open(new_file_path, "w") as new_file:
                    new_file.write(stlstr)
                newdata = newdata.reset_index()
                if addToDB:
                    self.datalayer.addMeasurementsDocument(desc=dict(name=NewFileName, bounds=points, dxdy=self._datalayer.getConfig()["dxdy"],
                                                                           xMin=newdata["gridxMin"][0], xMax=newdata["gridxMax"][0], yMin=newdata["gridyMin"][0],
                                                                           yMax=newdata["gridyMax"][0], zMin=newdata["gridzMin"][0], zMax=newdata["gridzMax"][0], **kwargs),
                                                                 type="stlFile",
                                                                 resource=stlstr,
                                                                 dataFormat="string")
        return stlstr, newdata

    def getDEM(self,data,**kwargs):

        """
        Creats a dem format string of topography.
        data is either a CutName of a topography or a geopandas/pandas/dask dataframe.
        any additional kwargs may be used to locate the dataframe.
        """

        if type(data) == str:
            geodata = self.datalayer.getDocuments(CutName=data,**kwargs)[0].getData()
        elif type(data) == geopandas.geodataframe.GeoDataFrame or type(data) == pandas.core.frame.DataFrame or type(data) == dask.dataframe.core.DataFrame:
            geodata = data
        else:
            raise KeyError("data should be geopandas/pandas/dask dataframe or a string.")
        if type(data) == geopandas.geodataframe.GeoDataFrame:
            xmin = geodata['geometry'].bounds['minx'].min()
            xmax = geodata['geometry'].bounds['maxx'].max()

            ymin = geodata['geometry'].bounds['miny'].min()
            ymax = geodata['geometry'].bounds['maxy'].max()
            Nx = int(((xmax - xmin) / self._datalayer.getConfig()["dxdy"]))
            Ny = int(((ymax - ymin) / self._datalayer.getConfig()["dxdy"]))

            print("Mesh boundaries x=(%s,%s) ; y=(%s,%s); N=(%s,%s)" % (xmin, xmax, ymin, ymax, Nx, Ny))
            dx = (xmax - xmin) / (Nx)
            dy = (ymax - ymin) / (Ny)
            print("Mesh increments: D=(%s,%s); N=(%s,%s)" % (dx, dy, Nx, Ny))
            # 2.2 build the mesh.
            grid_x, grid_y = numpy.mgrid[xmin:xmax:self._datalayer.getConfig()["dxdy"], ymin:ymax:self._datalayer.getConfig()["dxdy"]]
            # 3. Get the points from the geom
            Height = []
            XY = []

            for i, line in enumerate(geodata.iterrows()):
                if isinstance(line[1]['geometry'], LineString):
                    linecoords = [x for x in line[1]['geometry'].coords]
                    lineheight = [line[1]['HEIGHT']] * len(linecoords)
                    XY += linecoords
                    Height += lineheight
                else:
                    for ll in line[1]['geometry']:
                        linecoords = [x for x in ll.coords]
                        lineheight = [line[1]['HEIGHT']] * len(linecoords)
                        XY += linecoords
                        Height += lineheight

            grid_z2 = griddata(XY, Height, (grid_x, grid_y), method='cubic')

            if np.ma.is_masked(grid_z2):
                nans = grid_z2.mask
            else:
                nans = np.isnan(grid_z2)
            notnans = np.logical_not(nans)
            grid_z2[nans] = scipy.interpolate.griddata((grid_x[notnans], grid_y[notnans]), grid_z2[notnans],
                                                 (grid_x[nans], grid_y[nans]), method='nearest').ravel()
        else:
            if type(geodata) == dask.dataframe.core.DataFrame:
                geodata = geodata.compute()
            xmin = geodata[self.datalayer.getConfig()["xColumn"]].min()
            xmax = geodata[self.datalayer.getConfig()["xColumn"]].max()

            ymin = geodata[self.datalayer.getConfig()["yColumn"]].min()
            ymax = geodata[self.datalayer.getConfig()["yColumn"]].max()
            Nx = int(((xmax - xmin) / self._datalayer.getConfig()["dxdy"]))
            Ny = int(((ymax - ymin) / self._datalayer.getConfig()["dxdy"]))
            grid_z2=coordinateHandler.regularizeTimeSteps(data=geodata, fieldList=[self.datalayer.getConfig()["heightColumn"]],
                                                          coord1=self.datalayer.getConfig()["xColumn"], coord2=self.datalayer.getConfig()["yColumn"],
                                                          n=(Nx, Ny), addSurface=False, toPandas=False)[0][self.datalayer.getConfig()["heightColumn"]]

        DEMstring = f"{xmin} {xmax} {ymin} {ymax}\n{int(Nx)} {int(Ny)}\n"
        for i in range(Ny):
            for j in range(Nx):
                DEMstring += f"{(float(grid_z2[j, i]))} "
            DEMstring += "\n"

        return DEMstring

    def addHeight(self,data,groundData,coord1="x",coord2="y",coord3="z",resolution=10,savePandas=False,addToDB=False,file=None,fillna=0,**kwargs):
        """
        adds a column of height from ground for a dataframe which describes a mesh.
        params:
        data = The dataframe of the mesh.
        groundData = A dataframe that holds vertical coordinates values of the ground.
        coord1 = An horizontal coordinate name, default is "x"
        coord2 = A second horizontal coordinate name, default is "y"
        coord3 = A vertical coordinate name, default is "z"
        resolution = The cell length used in the conversion of the ground dataframe to a regular grid
        file = The new file's name
        savePandas = Boolian, whether to save the dataframe
        addToDB = Boolian, whether to add the dataframe to the DB
        kwargs = Any additional parameters to use as descriptors of the new file
        """
        nx = int((groundData[coord1].max()-groundData[coord1].min())/resolution)
        ny = int((groundData[coord2].max() - groundData[coord2].min()) / resolution)
        xarrayGround = coordinateHandler.regularizeTimeSteps(data=groundData, fieldList=[coord3],coord1=coord1, coord2=coord2, n=(nx,ny), addSurface=False, toPandas=False)[0]
        ground = []
        for i in range(len(data)):
            x = data.loc[i][coord1]
            y = data.loc[i][coord2]
            ground.append(float(xarrayGround.interp(**{coord1:x,coord2:y}).fillna(fillna)[coord3]))
            if i > 9999 and i % 10000 == 0:
                print(i)
        data["ground"]=ground
        data["height"] = data[coord3] - data["ground"]
        data.loc[data.height < 0, "height"] = 0
        if savePandas:
            file = os.path.join(self.datalayer.FilesDirectory,"cellData.parquet") if file is None else file
            data.to_parquet(file, compression="gzip")
            if addToDB:
                self.datalayer.addCacheDocument(resource=file, dataFormat="parquet",
                                   type="cellData", desc=dict(resolution=resolution,**kwargs))

        return data

def get_altitdue_ip(lat, lon):
    """
    returning the altitude of the point, uses free mapquest data that is limited in the amount of calls per month, it uses Nir BAMBA Benmoshe key

    param lat - the path where we save the stl
    param lon - the width of the domain in the x direction

    return:
    altitude - meters above sea level
    """

    resp = requests.get(
        'http://open.mapquestapi.com/elevation/v1/profile?key=D5z9RSebQJLbUs4bohANIB4TzJdbvyvm&shapeFormat=raw&latLngCollection=' + str(
            lat) + ',' + str(lon))

    height = resp.json()['elevationProfile'][0]['height']

    return height


def get_altitdue_gdal(lat, lon):
    #        if lat<=30:
    # USGS EROS Archive - Digital Elevation - Global Multi-resolution Terrain Elevation Data 2010 (GMTED2010)
    # fheight = r'/data3/nirb/10N030E_20101117_gmted_med075.tif'
    # fheight = r'/data3/nirb/30N030E_20101117_gmted_med075.tif'
    # https://dds.cr.usgs.gov/srtm/version2_1/SRTM3/Africa/   # 90m resolution
    if lat > 29 and lat < 30 and lon > 34 and lon < 35:
        fheight = r'/data3/nirb/N29E034.hgt'
    elif lat > 29 and lat < 30 and lon > 35 and lon < 36:
        fheight = r'/data3/nirb/N29E035.hgt'
    elif lat > 30 and lat < 31 and lon > 34 and lon < 35:
        fheight = r'/data3/nirb/N30E034.hgt'
    elif lat > 30 and lat < 31 and lon > 35 and lon < 36:
        fheight = r'/data3/nirb/N30E035.hgt'
    elif lat > 31 and lat < 32 and lon > 34 and lon < 35:
        fheight = r'/data3/nirb/N31E034.hgt'
    elif lat > 31 and lat < 32 and lon > 35 and lon < 36:
        fheight = r'/data3/nirb/N31E035.hgt'
    elif lat > 32 and lat < 33 and lon > 34 and lon < 35:
        fheight = r'/data3/nirb/N32E034.hgt'
    elif lat > 32 and lat < 33 and lon > 35 and lon < 36:
        fheight = r'/data3/nirb/N32E035.hgt'
    elif lat > 33 and lat < 33 and lon > 35 and lon < 36:
        fheight = r'/data3/nirb/N33E035.hgt'
    else:
        print('!!!!NOT in Israel !!!!!!!!')
        # taken from https://earthexplorer.usgs.gov/
        fheight = r'/ibdata2/nirb/gt30e020n40.tif'

    ds = gdal.Open(fheight)
    myarray = numpy.array(ds.GetRasterBand(1).ReadAsArray())
    myarray[myarray < -1000] = 0
    gt = ds.GetGeoTransform()
    rastery = (lon - gt[0]) / gt[1]
    rasterx = (lat - gt[3]) / gt[5]
    height11 = myarray[int(rasterx), int(rastery)]
    height12 = myarray[int(rasterx) + 1, int(rastery)]
    height21 = myarray[int(rasterx), int(rastery) + 1]
    height22 = myarray[int(rasterx) + 1, int(rastery) + 1]
    height1 = (1. - (rasterx - int(rasterx))) * height11 + (rasterx - int(rasterx)) * height12
    height2 = (1. - (rasterx - int(rasterx))) * height21 + (rasterx - int(rasterx)) * height22
    height = (1. - (rastery - int(rastery))) * height1 + (rastery - int(rastery)) * height2

    return height



if __name__ == "__main__":
    longitude = random.randint(35750, 35800) / 1000.0  # Hermon
    latitude = random.randint(33250, 33800) / 1000.0
    #    lon = 34.986008  // elevation should be 273m according to amud anan
    #    lat = 32.808486  // elevation should be 273m according to amud anan
    #    lon = 35.755  // elevation should be ~820 according to amud anan
    #    lat = 33.459  // elevation should be ~820 according to amud anan
    longitude = 35.234987  # 744m
    latitude = 31.777978  # 744m
    #alt1 = get_altitdue_ip(lat, longitude)
    alt2 = get_altitdue_gdal(latitude, longitude)
    print("the altitude at position: ", latitude, longitude, " is ", alt1, alt2)

