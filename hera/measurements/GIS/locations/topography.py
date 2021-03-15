from .import abstractLocation
import geopandas
import random
import gdal
import numpy as np

from shapely.geometry import LineString
import math
import pandas
from numpy import array, cross, sqrt
import numpy
from scipy.interpolate import griddata

import requests

from ....simulations.utils import coordinateHandler


import os

class abstractLocationToolkit(abstractLocation.AbstractLocationToolkit):

    _analysis = None

    @property
    def analysis(self):
        return self._analysis

    def __init__(self, projectName, FilesDirectory=""):

        super().__init__(projectName=projectName,FilesDirectory=FilesDirectory,toolkitName="topography")
        self._analysis = analysis(projectName=projectName, dataLayer=self)

    #
    # def setConfig(self, Source="BNTL", dxdy=50, heightSource="USGS", dbName=None, **kwargs):
    #     config = dict(source=Source,dxdy=dxdy,heightSource=heightSource,**kwargs)
    #     super().setConfig(config, dbName=dbName)

    def getHeight(self, latitude, longitude):

        try:
            height = self.__getattribute__("getHeight_%s" % self.getConfig()["heightSource"])(longitude=longitude,latitude=latitude)
        except AttributeError:
            print("The height source is not known.")

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

    def toSTL(self, regionNameOrData, outputFileName, flat=None, saveMode=abstractLocation.toolkit.TOOLKIT_SAVEMODE_ONLYFILE, **kwargs):
        """
            Converts a topographic data to STL file.

            The topographic data can be given as a regionName in the DB, geopandas.Dataframe
            or a geoJSON string.

        Parameters:
        -----------

            regionNameOrData: str or geopandas.DataFrame
                The data that should be converted to stl.
                May be a dataframe or a name of a saved polygon in the database.


            outputFileName: str
                A name for the new stl file, also used in the stl string. (string)

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


            dxdy: float.
                the dimension of each cell in the mesh in meters, the default is 50.


            flat: float or None.
                If not None, it assumes that the area is flat and the value of flat is the height of the mesh cells.

            kwargs:
                Any additional metadata to be added to the new document in the database.

        Returns
        -------
            DB document object or the STL as a string (if no saving to DB).

        """

        if isinstance(regionNameOrData,str):
            regionDoc self.getRegionByName(regionNameOrData)



        documents = self.datalayer.getMeasurementsDocuments(type="stlFile", bounds=points, dxdy=self._datalayer.getConfig()["dxdy"])
        if len(documents) >0:
            stlstr = documents[0].getData()
            newdict = documents[0].asDict()
            newdata = pandas.DataFrame(dict(gridxMin=[newdict["desc"]["xMin"]], gridxMax=[newdict["desc"]["xMax"]],
                                            gridyMin=[newdict["desc"]["yMin"]], gridyMax=[newdict["desc"]["yMax"]],
                                            gridzMin=[newdict["desc"]["zMin"]], gridzMax=[newdict["desc"]["zMax"]]))
        else:
            stl = stlFactory()
            stlstr, newdata = stl.Convert_geopandas_to_stl(gpandas=geodata, points=points, flat=flat, NewFileName=outputFileName, dxdy=self._datalayer.getConfig()["dxdy"])
            if save:
                p = self.datalayer.FilesDirectory if path is None else path
                new_file_path = os.path.join(p, f"{outputFileName}.stl")
                with open(new_file_path, "w") as new_file:
                    new_file.write(stlstr)
                newdata = newdata.reset_index()
                if addToDB:
                    self.datalayer.addMeasurementsDocument(desc=dict(name=outputFileName, bounds=points, dxdy=self._datalayer.getConfig()["dxdy"],
                                                                     xMin=newdata["gridxMin"][0], xMax=newdata["gridxMax"][0], yMin=newdata["gridyMin"][0],
                                                                     yMax=newdata["gridyMax"][0], zMin=newdata["gridzMin"][0], zMax=newdata["gridzMax"][0], **kwargs),
                                                                 type="stlFile",
                                                                 resource=stlstr,
                                                                 dataFormat="string")
        return stlstr, newdata

    def getDEM(self,data,**kwargs):

        """
        Creats a dem format string of topography.
        data is either a CutName of a topography or a geopandas dataframe.
        any additional kwargs may be used to locate the dataframe.
        """

        if type(data) == str:
            geodata = self.datalayer.getDocuments(CutName=data,**kwargs)[0].getData()
        elif type(data) == geopandas.geodataframe.GeoDataFrame:
            geodata = data
        else:
            raise KeyError("data should be geopandas dataframe or a string.")
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

        DEMstring = f"{xmin} {xmax} {ymin} {ymax}\n{int(Nx)} {int(Ny)}\n"

        for i in range(Nx):
            for j in range(Ny):
                DEMstring += f"{(float(grid_z2[j, i]))} "
            DEMstring += "\n"

        return DEMstring


class analysis():

    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, projectName, dataLayer):
        self._datalayer = dataLayer


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




class stlFactory():

    def _make_facet_str(self, n, v1, v2, v3):
        facet_str = 'facet normal ' + ' '.join(map(str, n)) + '\n'
        facet_str += '  outer loop\n'
        facet_str += '      vertex ' + ' '.join(map(str, v1)) + '\n'
        facet_str += '      vertex ' + ' '.join(map(str, v2)) + '\n'
        facet_str += '      vertex ' + ' '.join(map(str, v3)) + '\n'
        facet_str += '  endloop\n'
        facet_str += 'endfacet\n'
        return facet_str

    def _makestl(self, X, Y, elev, NewFileName):
        """
            Takes a mesh of x,y and elev and convert it to stl file.

            X - matrix of x coordinate. [[ like meshgrid ]]
            Y - matrix of y coordinate. [[ like meshgrid ]]
            elev - matrix of elevation.

        """
        base_elev = elev.min() - 10
        stl_str = 'solid ' + NewFileName + '\n'
        for i in range(elev.shape[0] - 1):
            for j in range(elev.shape[1] - 1):

                x = X[i, j];
                y = Y[i, j]
                v1 = [x, y, elev[i, j]]

                x = X[i + 1, j];
                y = Y[i, j]
                v2 = [x, y, elev[i + 1, j]]

                x = X[i, j];
                y = Y[i, j + 1]
                v3 = [x, y, elev[i, j + 1]]

                x = X[i + 1, j + 1];
                y = Y[i + 1, j + 1]
                v4 = [x, y, elev[i + 1, j + 1]]

                # dem facet 1
                n = cross(array(v1) - array(v2), array(v1) - array(v3))
                n = n / sqrt(sum(n ** 2))
                stl_str += self._make_facet_str(n, v1, v2, v3)

                # dem facet 2
                n = cross(array(v2) - array(v3), array(v2) - array(v4))
                n = n / sqrt(sum(n ** 2))
                # stl_str += self._make_facet_str( n, v2, v3, v4 )
                stl_str += self._make_facet_str(n, v2, v4, v3)

                # base facets
                v1b = list(v1)
                v2b = list(v2)
                v3b = list(v3)
                v4b = list(v4)

                v1b[-1] = base_elev
                v2b[-1] = base_elev
                v3b[-1] = base_elev
                v4b[-1] = base_elev

                n = [0.0, 0.0, -1.0]

                stl_str += self._make_facet_str(n, v1b, v2b, v3b)
                stl_str += self._make_facet_str(n, v2b, v3b, v4b)

                vlist = [v1, v2, v3, v4]
                vblist = [v1b, v2b, v3b, v4b]

                # Now the walls.
                for k, l in [(0, 1), (0, 2), (1, 3), (2, 3)]:
                    # check if v[i],v[j] are on boundaries.
                    kboundary = False
                    if vlist[k][0] == X.min() or vlist[k][0] == X.max():
                        kboundary = True

                    lboundary = False
                    if vlist[l][1] == Y.min() or vlist[l][1] == Y.max():
                        lboundary = True

                    if (kboundary or lboundary):
                        # Add i,j,j-base.
                        n = cross(array(vlist[k]) - array(vlist[l]), array(vblist[l]) - array(vlist[l]))
                        n = n / sqrt(sum(n ** 2))
                        stl_str += self._make_facet_str(n, vlist[k], vblist[l], vlist[l])

                        # add j-base,i-base,i
                        n = cross(array(vlist[k]) - array(vblist[k]), array(vlist[k]) - array(vblist[l]))
                        n = n / sqrt(sum(n ** 2))
                        stl_str += self._make_facet_str(n, vlist[k], vblist[k], vblist[l])

        stl_str += 'endsolid ' + NewFileName + '\n'
        return stl_str

    def Convert_geopandas_to_stl(self, gpandas, points, NewFileName, dxdy=50, flat=None):
        """
            Gets a shape file of topography.
            each contour line has property 'height'.
            Converts it to equigrid xy mesh and then build the STL.
        """

        # 1. Convert contour map to regular height map.
        # 1.1 get boundaries
        xmin = points[0]
        xmax = points[2]

        ymin = points[1]
        ymax = points[3]

        print("Mesh boundaries x=(%s,%s) ; y=(%s,%s)" % (xmin, xmax, ymin, ymax))
        # 1.2 build the mesh.
        grid_x, grid_y = numpy.mgrid[(xmin):(xmax):dxdy, (ymin):(ymax):dxdy]
        # 2. Get the points from the geom
        Height = []
        XY = []
        for i, line in enumerate(gpandas.iterrows()):
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
        if flat is not None:
            for i in range(len(Height)):
                Height[i] = flat
        grid_z2 = griddata(XY, Height, (grid_x, grid_y), method='cubic')
        grid_z2 = self.organizeGrid(grid_z2)
        stlstr = self._makestl(grid_x, grid_y, grid_z2, NewFileName)

        data = pandas.DataFrame({"XY": XY, "Height": Height, "gridxMin":grid_x.min(), "gridxMax":grid_x.max(),
                                 "gridyMin":grid_y.min(), "gridyMax":grid_y.max(), "gridzMin":grid_z2[~numpy.isnan(grid_z2)].min(), "gridzMax":grid_z2[~numpy.isnan(grid_z2)].max(),})

        return stlstr, data

    def organizeGrid(self, grid):

        for row in grid:
            for i in range(len(row)):
                if math.isnan(row[i]):
                    pass
                else:
                    break
            for n in range(i):
                row[n] = row[i]
            for i in reversed(range(len(row))):
                if math.isnan(row[i]):
                    pass
                else:
                    break
            for n in range(len(row)-i):
                row[-n-1] = row[i]
        return grid

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

