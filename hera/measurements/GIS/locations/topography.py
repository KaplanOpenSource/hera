import geopandas
import random
import gdal
import io
import scipy
import numpy
import math
import os
import requests
from scipy.interpolate import griddata
from shapely.geometry import LineString
from .... import toolkit
import xarray
import dask

from .import abstractLocation
from ....simulations.utils import coordinateHandler

from ....datalayer import datatypes


class TopographyToolkit(abstractLocation.AbstractLocationToolkit):

    _stlFactory = None
    _srtmBounds = None
    _toolkitname = None

    @property
    def doctype(self):
        return f"{self.toolkitName}_STL"

    @property
    def stlFactory(self):
        return self._stlFactory

    @property
    def srtmBounds(self):
        return self._srtmBounds

    @property
    def toolkitname(self):
        return self._toolkitname

    def __init__(self, projectName, FilesDirectory=""):
        toolkitName = "Topography"
        super().__init__(projectName=projectName,FilesDirectory=FilesDirectory,toolkitName=toolkitName)
        self._analysis = analysis(projectName=projectName, dataLayer=self)
        self._toolkitname = toolkitName
        self._stlFactory = stlFactory()
        srtmDocs = self.getMeasurementsDocumentsAsDict(name="SRTM")
        self._srtmBounds = []
        if type(srtmDocs)==dict:
            for doc in srtmDocs["documents"]:
                self._srtmBounds.append(doc["desc"]["Bounds"])

    def makeRegion_SRTM(self, points,outputFileName,**kwargs):

        availableBounds = {}
        for bounds in self.srtmBounds:
            if points[0] >= bounds[0] and points[0] <= bounds[2] and points[2] >= bounds[0] and points[2] <= bounds[2] and points[1] >= bounds[1] and points[1] <= bounds[3] and points[3] >= bounds[1] and points[3] <= bounds[3]:
                availableBounds["Bounds"] = bounds
                break
        if len(availableBounds)>0:
            FileName = f"{outputFileName}.parquet"
            import pdb
            pdb.set_trace()
            allData = self.getMeasurementsDocuments(name="SRTM",type=toolkit.TOOLKIT_DATASOURCE_TYPE, **availableBounds)[0].getData()
            dataArray = allData.read(1)
            xs = numpy.linspace(allData.bounds.left,allData.bounds.right,dataArray.shape[1])
            ys = numpy.linspace(allData.bounds.top,allData.bounds.bottom,dataArray.shape[0])
            xmin = numpy.where(xs>points[0])[0].min()
            ymin = numpy.where(ys<points[1])[0].min()
            xmax = numpy.where(xs<points[2])[0].max()
            ymax = numpy.where(ys>points[3])[0].max()
            xColumn = self.getConfig()["xColumn"] if "xColumn" in self.getConfig().keys() else "x"
            yColumn = self.getConfig()["yColumn"] if "yColumn" in self.getConfig().keys() else "y"
            heightColumn = self.getConfig()["heightColumn"] if "heightColumn" in self.getConfig().keys() else "height"
            data = xarray.DataArray(data=dataArray[ymax:ymin,xmin:xmax],dims=["yWGS84","xWGS84"],coords=[ys[ymax:ymin],xs[xmin:xmax]]).to_dataframe(heightColumn).reset_index()
            gdf = geopandas.GeoDataFrame(data, geometry=geopandas.points_from_xy(data["xWGS84"], data["yWGS84"]))
            gdf.crs = {"init": "epsg:4326"}
            gdf = gdf.to_crs({"init": "epsg:2039"})
            gdf[xColumn] = gdf.geometry.x
            gdf[yColumn] = gdf.geometry.y
            data = dask.dataframe.from_pandas(gdf.drop(columns="geometry"),npartitions=1)
            data.to_parquet(FileName,compression="GZIP")
        else:
            raise KeyError("Couldn't find SRTM files which include the requested points!")
        return

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
        myarray = numpy.array(ds.GetRasterBand(1).ReadAsArray())
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

    def toSTL(self, regionNameOrData, outputFileName, dxdy=50, saveMode=abstractLocation.toolkit.TOOLKIT_SAVEMODE_ONLYFILE, additionalData=dict()):
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

            additional_data: dict
                Any additional metadata to be added to the new document in the database.

        Returns
        -------
            DB document object or the STL as a string (if no saving to DB).

        """

        if isinstance(regionNameOrData,str):
            region = self.getRegionByName(regionNameOrData)
            if region is None:
                region = geopandas.read_file(io.StringIO(regionNameOrData))
        elif isinstance(regionNameOrData,geopandas.geodataframe):
            region = regionNameOrData
        else:
            raise ValueError(f"The regionNameOrData must be region name, geoJSON str or geodataframe")

        stlstr = self.stlFactory.vectorToSTL(region,dxdy=dxdy)

        if saveMode in [abstractLocation.toolkit.TOOLKIT_SAVEMODE_FILEANDDB,
                        abstractLocation.toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE,
                        abstractLocation.toolkit.TOOLKIT_SAVEMODE_ONLYFILE,
                        abstractLocation.toolkit.TOOLKIT_SAVEMODE_ONLYFILE_REPLACE]:

            if "stl" not in outputFileName:
                outputFileName = f"{outputFileName.split('.')[0]}.stl"

            fullfileNames = os.path.join(self.FilesDirectory,outputFileName)

            if os.path.exists(fullfileNames) and saveMode in [abstractLocation.toolkit.TOOLKIT_SAVEMODE_ONLYFILE,
                                                              abstractLocation.toolkit.TOOLKIT_SAVEMODE_FILEANDDB]:
                raise FileExistsError(f"The output STL file {fullfileNames} exists")

            with open(fullfileNames,'w') as outfiles:
                outfiles.write(stlstr)


            if saveMode in [abstractLocation.toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE,
                            abstractLocation.toolkit.TOOLKIT_SAVEMODE_FILEANDDB]:

                regionDoc = self.getSTL(regionNameOrData)

                if regionDoc is not None and saveMode==abstractLocation.toolkit.TOOLKIT_SAVEMODE_FILEANDDB:
                    raise ValueError(f"{regionNameOrData} already exists in the DB")

                xmin, ymin, xmax, ymax = region.bounds

                additionalData.update({abstractLocation.TOOLKIT_LOCATION_REGIONNAME: outputFileName.split('.')[0],
                                       abstractLocation.toolkit.TOOLKIT_DATASOURCE_NAME: outputFileName.split('.')[0],
                                       abstractLocation.toolkit.TOOLKIT_TOOLKITNAME_FIELD: self.toolkitName,
                                       "xmin": xmin,
                                       "xmax": xmax,
                                       "ymin": ymin,
                                       "ymax": ymax
                                       })

                if regionDoc is None:
                    regionDoc = self.addCacheDocument(type=self.doctype,
                                                    resource=fullfileNames,
                                                    dataFormat=datatypes.STRING,
                                                    desc=additionalData)
                else:
                    regionDoc.resource = fullfileNames
                    regionDoc.desc = additionalData
                    regionDoc.save()


        return stlstr if regionDoc is None else regionDoc

    def getDEM(self,regionNameOrData,dxdy=50,**filters):

        """
            Creates a Data Elevation Model (DEM) format string of topography.

            Parameters
            ----------

            regionNameOrData: str, geoJSON string, or geodataframe.
                The name of the regions, geoJSON string or geodataframe.

            dxdy: float
                The resolution of the conversion .

            Returns
            -------

                The DEM string format.
        """

        if isinstance(regionNameOrData,str):
            data = self.getDatasourceData(datatypes=regionNameOrData,**filters)
            if data is None:
                data = geopandas.read_file(io.StringIO(regionNameOrData))
        elif isinstance(regionNameOrData,geopandas.geodataframe):
            data = regionNameOrData
        else:
            raise ValueError("regionNameOrData must be wither region name in the DB, geoJSON string or a geopandas.geodataframe")

        rasterized = self.stlFactory.rasterize(data,dxdy=dxdy)

        if numpy.ma.is_masked(rasterized['height']):
            nans = rasterized['height'].mask
        else:
            nans = numpy.isnan(rasterized['height'])

        notnans = numpy.logical_not(nans)

        grid_x = rasterized['x']
        grid_y = rasterized['y']

        xmin = grid_x.min()
        xmax = grid_x.max()
        ymin = grid_y.min()
        ymax = grid_y.max()
        Nx   = grid_x.shape[0]
        Ny   = grid_x.shape[1]


        rasterized['height'][nans] = scipy.interpolate.griddata((grid_x[notnans], grid_y[notnans]), rasterized['height'][notnans],
                                             (grid_x[nans], grid_y[nans]), method='nearest').ravel()

        DEMstring = f"{xmin} {xmax} {ymin} {ymax}\n{Nx} {Ny}\n"

        for i in range(Nx):
            for j in range(Ny):
                DEMstring += f"{(float(rasterized['height'][j, i]))} "
            DEMstring += "\n"

        return DEMstring


    def getSTL(self,regionNameSTL):
        """
            Return the topography STL

        Parameters
        -----------
        regionNameSTL: str
            The name of the region


        Returns
        -------
            Return the ATL document
        """
        desc = {
                abstractLocation.TOOLKIT_LOCATION_REGIONNAME: regionNameSTL,
                "type" : self.doctype
                }
        docList = self.getCacheDocuments(**desc)
        return None if len(docList)==0 else docList[0]

    def getSTLDocumentsList(self):
        """
            Return a list with all the STL documents list
        Returns
        -------
            list with STL documents
        """
        desc = {
                "type" : self.doctype
                }
        return self.getCacheDocuments(**desc)

    def getSTLList(self):
        """
            Return a list with all the STL documents names
        Returns
        -------
            list with STL documents
        """
        desc = {
                "type" : self.doctype
                }
        return [x.desc[abstractLocation.TOOLKIT_LOCATION_REGIONNAME] for x in self.getCacheDocuments(**desc)]

    def loadData(self):
        """
         BUILD
        :return:
        """
        pass

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


class stlFactory:
    """
        Helper class to convert geopandas to STL
    """

    _heightColumnsNames=None

    def __init__(self):
        self._heightColumnsNames= "HEIGHT"

    @property
    def heightColumnsNames(self):
        return self._heightColumnsNames

    @heightColumnsNames.setter
    def heightColumnsNames(self, value):
        self._heightColumnsNames = value

    def _make_facet_str(self, n, v1, v2, v3):
        facet_str = 'facet normal ' + ' '.join(map(str, n)) + '\n'
        facet_str += '  outer loop\n'
        facet_str += '      vertex ' + ' '.join(map(str, v1)) + '\n'
        facet_str += '      vertex ' + ' '.join(map(str, v2)) + '\n'
        facet_str += '      vertex ' + ' '.join(map(str, v3)) + '\n'
        facet_str += '  endloop\n'
        facet_str += 'endfacet\n'
        return facet_str

    def rasterize(self,gpandas, dxdy=50):
        """
            Convert the height-contour in geopandas to raster (regular mesh)


        Parameters
        -----------
        gpandas: geopandas.geodataframe
                The height contours as geopandas.

                The structure of the geodataframe is geometry
                which are height lines and a height column (its name given in self.heightColumnsNames, default is HEIGHT).

        dxdy: float
            The resolution of the regularization


        Returns
        -------
            A dict after interpolating to regular mesh.

            Return dict with the keys:
              x:    A 2D map of the x coordinate
              y :   A 2D map of the y coordinate
              height: A 2D map of the height


        """
        # 1. Convert contour map to regular height map.
        # 1.1 get boundaries
        xmin,ymin,xmax,ymax = gpandas.bounds

        #print("Mesh boundaries x=(%s,%s) ; y=(%s,%s)" % (xmin, xmax, ymin, ymax))
        # 1.2 build the mesh.
        grid_x, grid_y = numpy.mgrid[(xmin):(xmax):dxdy, (ymin):(ymax):dxdy]

        # 2. Get the points from the geom
        Height = []
        XY = []
        for i, line in enumerate(gpandas.iterrows()):
            if isinstance(line[1]['geometry'], LineString):
                linecoords = [x for x in line[1]['geometry'].coords]
                lineheight = [line[1][self.heightColumnsNames]] * len(linecoords)
                XY += linecoords
                Height += lineheight
            else:
                for ll in line[1]['geometry']:
                    linecoords = [x for x in ll.coords]
                    lineheight = [line[1][self.heightColumnsNames]] * len(linecoords)
                    XY += linecoords
                    Height += lineheight

        grid_z2 = griddata(XY, Height, (grid_x, grid_y), method='cubic')
        grid_z2 = self._organizeGrid(grid_z2)

        return dict(x=grid_x,y=grid_y,height=grid_z2)

    def rasterToSTL(self, grid_x, grid_y, grid_z, solidName):
        """

        Parameters
        ----------
        grid_x: numpy.array

            2D array with the x grid
        grid_y: numpy.array

            2D array with the y grid
        grid_z: numpy.array

            2D array with the z (height) grid

        solidName: str
            The name of the topography solid

        Returns
        -------
            The STL as string.
        """
        base_elev = grid_z.min() - 10
        stl_str = f"solid {solidName}\n"
        for i in range(grid_z.shape[0] - 1):
            for j in range(grid_z.shape[1] - 1):

                x = grid_x[i, j];
                y = grid_y[i, j]
                v1 = [x, y, grid_z[i, j]]

                x = grid_x[i + 1, j];
                y = grid_y[i, j]
                v2 = [x, y, grid_z[i + 1, j]]

                x = grid_x[i, j];
                y = grid_y[i, j + 1]
                v3 = [x, y, grid_z[i, j + 1]]

                x = grid_x[i + 1, j + 1];
                y = grid_y[i + 1, j + 1]
                v4 = [x, y, grid_z[i + 1, j + 1]]

                # dem facet 1
                n = numpy.cross(array(v1) - array(v2), array(v1) - array(v3))
                n = n / sqrt(sum(n ** 2))
                stl_str += self._make_facet_str(n, v1, v2, v3)

                # dem facet 2
                n = numpy.cross(array(v2) - array(v3), array(v2) - array(v4))
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
                    if vlist[k][0] == grid_x.min() or vlist[k][0] == grid_x.max():
                        kboundary = True

                    lboundary = False
                    if vlist[l][1] == grid_y.min() or vlist[l][1] == grid_y.max():
                        lboundary = True

                    if (kboundary or lboundary):
                        # Add i,j,j-base.
                        n = numpy.cross(array(vlist[k]) - array(vlist[l]), array(vblist[l]) - array(vlist[l]))
                        n = n / sqrt(sum(n ** 2))
                        stl_str += self._make_facet_str(n, vlist[k], vblist[l], vlist[l])

                        # add j-base,i-base,i
                        n = numpy.cross(array(vlist[k]) - array(vblist[k]), array(vlist[k]) - array(vblist[l]))
                        n = n / sqrt(sum(n ** 2))
                        stl_str += self._make_facet_str(n, vlist[k], vblist[k], vblist[l])

        stl_str += f"endsolid {solidName}\n"
        return stl_str


    def vectorToSTL(self,gpandas,dxdy=50,solidName="Topography"):
        """
            Convert the vector to topography

        Parameters
        ----------
        gpandas: geopandas.geodataframe
            The height contours as geopandas.

        dxdy: float
            The resolution of the conversion

        solidName: str
            The name of the solid. The default is 'Topogrphay'.

        Returns
        -------
            string of the STL
        """

        rasterMap = self.rasterize(gpandas,dxdy=dxdy)
        return self.rasterToSTL(grid_x=rasterMap['x'],
                                grid_y=rasterMap['y'],
                                grid_z=rasterMap['height'],
                                solidName=solidName)


    def _organizeGrid(self, grid):

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

