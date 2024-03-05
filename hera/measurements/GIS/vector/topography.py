import geopandas
import io
import scipy
import numpy
import math
import os
from scipy.interpolate import griddata
from shapely.geometry import LineString
from hera import toolkit
import dask
from numpy import array, sqrt
import pandas

from . import toolkit
from ....simulations.utils import coordinateHandler
from ....utils.logging import get_classMethod_logger

from hera.toolkit import TOOLKIT_SAVEMODE_ONLYFILE


class TopographyToolkit(toolkit.VectorToolkit):

    @property
    def stlFactory(self):
        return self._stlFactory

    def __init__(self, projectName, filesDirectory=""):
        toolkitName = "Topography"
        super().__init__(projectName=projectName, filesDirectory=filesDirectory, toolkitName=toolkitName)
        self._analysis = analysis(projectName=projectName, dataLayer=self)
        self._toolkitname = toolkitName
        self._stlFactory = stlFactory()


    def cutRegionFromSource(self, shapeDataOrName, datasourceName, isBounds = False, crs = None): # If  shapeDataOrName is data: if is Bounds = True: use the Bbox of shape as the region, else use the shpae as the region
        """
            Cuts a the shape from the requested datasource

            It overrides the parent procedure because we need to perform intersection operation on the results.

        Parameters
        ----------
        shapeDataOrName: GeoDataFrame,dict,list,Polygon
            The shape data to use.

                        list - a box with the corners in [xmin,ymin,xmax,ymax]

        datasourceName : str
            The name of the satasource to cur from.

        isBounds : bool
            If true, use the bounding box fo the polygon.

        crs : int
            The EPSG of the coordinate system of the shape (if it is a shape and not in the dtasource ative coordinates).

        Returns
        -------

        """
        topography = super().cutRegionFromSource(shapeDataOrName=shapeDataOrName, datasourceName=datasourceName, isBounds =isBounds, crs =crs)
        logger = get_classMethod_logger(self, "cutRegionFromSource")
        if isinstance(shapeDataOrName, str):
            shape = self.getRegionData(shapeDataOrName)
        else:
            shape = self._RegionToGeopandas(shapeDataOrName, crs=crs)
            doc = self.getDatasourceDocument(datasourceName=datasourceName)
            logger.debug(f"The datasource {datasourceName} is pointing to {doc.resource}")
            doc.desc['desc'].update({'crs': 2039})
            if 'crs' not in doc.desc['desc']:
                logger.error(f"The datasource {datasourceName} has no CRS defined in the metadata. please add it")
                raise ValueError(f"The datasource {datasourceName} has no CRS defined in the metadata. please add it")

            if shape.crs is None:
                logger.execution("The region was defined without crs. Using the crs of the datasource.")
                shape.crs = doc.desc['desc']['crs']
                shape = shape.to_crs(topography.crs)

        topography['geometry'] = topography.intersection(shape.iloc[0].geometry)
        return topography



    def geoPandasToSTL(self,gpandas, dxdy=50, solidName="Topography"):
        """
            Transforsm the gpandas to STL.

        Parameters
        ----------
        gpandas
        dxdy
        solidName

        Returns
        -------

        """
        return self.stlFactory.vectorToSTL(gpandas, dxdy=dxdy, solidName="Topography")

    def regionToSTL(self, shapeDataOrName, dxdy, datasourceName, crs=None):
        """
            Converts a region in a vector height map (contours) to STL at requested resolution

            We always use the bounding box of the input shape .

        Parameters:
        -----------

            shapeDataOrName: str or polygon
                The shape that we will use.

            dxdy: float.
                the dimension of each cell in the mesh in meters

            dataSourceName: str
                The datasource to use.

            crs : int [optional]
                Used if the CRS of the shapeData is different than the CRS of the datasource.

        Returns
        -------
            str
            The STL string.

        """
        logger = get_classMethod_logger(self, "regionToSTL")
        logger.info("-- Start --")

        topography = self.cutRegionFromSource(shapeDataOrName, datasourceName=datasourceName, isBounds=True, crs=crs)


        if len(topography) == 0:
            logger.warning("The requested region is empty. ")
            stlstr = None
        else:
            stlstr = self.stlFactory.vectorToSTL(topography,dxdy=dxdy)

        logger.info("-- End --")

        return stlstr

    def toDEM(self,regionNameOrData,dxdy=50,**filters):

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
            data = self.getDatasourceData(datasourceName=regionNameOrData,**filters)
            if data is None:
                data = geopandas.read_file(io.StringIO(regionNameOrData))
        elif isinstance(regionNameOrData,geopandas.geodataframe.GeoDataFrame):
            data = regionNameOrData
        else:
            raise ValueError("regionNameOrData must be wither region name in the DB, geoJSON string or a geopandas.geodataframe")

        rasterized = self.stlFactory.rasterizeGeopandas(data, dxdy=dxdy)

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



class analysis():

    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, projectName, dataLayer):
        self._datalayer = dataLayer


    def addHeight(self, data, groundData, coord1="x", coord2="y", coord3="z", resolution=10, saveMode=TOOLKIT_SAVEMODE_ONLYFILE, file=None, fillna=0, **kwargs):
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
        valuesDict = {}
        for i in range(len(data)):
            x = data.loc[i][coord1]
            y = data.loc[i][coord2]
            if x not in valuesDict.keys():
                valuesDict[x] = {}
            if y in valuesDict[x].keys():
                ground.append(valuesDict[x][y])
            else:
                val = float(xarrayGround.interp(**{coord1:x,coord2:y}).fillna(fillna)[coord3])
                ground.append(val)
                valuesDict[x][y] = val
            if i > 9999 and i % 10000 == 0:
                print(f"Finished calculation for {i} cells")

        data["ground"]=ground
        data["height"] = data[coord3] - data["ground"]
        data.loc[data.height < 0, "height"] = 0
        if saveMode in [toolkit.TOOLKIT_SAVEMODE_FILEANDDB,
                        toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE,
                        toolkit.TOOLKIT_SAVEMODE_ONLYFILE,
                        toolkit.TOOLKIT_SAVEMODE_ONLYFILE_REPLACE]:

            file = os.path.join(self.datalayer.filesDirectory, "cellData.parquet") if file is None else file

            if os.path.exists(file) and saveMode in [toolkit.TOOLKIT_SAVEMODE_ONLYFILE,
                                                     toolkit.TOOLKIT_SAVEMODE_FILEANDDB]:
                raise FileExistsError(f"The output file {file} exists")

            data.to_parquet(file, compression="gzip")


            if saveMode in [toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE,
                            toolkit.TOOLKIT_SAVEMODE_FILEANDDB]:

                regionDoc = self.datalayer.getCacheDcouments(resource=file, dataFormat="parquet",type="cellData", desc=dict(resolution=resolution,**kwargs))

                if len(regionDoc) >0 and saveMode== toolkit.TOOLKIT_SAVEMODE_FILEANDDB:
                    raise ValueError(f"{file} already exists in the DB")
                else:
                    self.datalayer.addCacheDocument(resource=file, dataFormat="parquet",
                                                    type="cellData", desc=dict(resolution=resolution, **kwargs))

        return data


class stlFactory:
    """
        Helper class to convert topography geopandas to STL
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

    def rasterizeGeopandas(self, gpandas, dxdy=50):
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
        xmin = gpandas.bounds.min()["minx"]
        ymin = gpandas.bounds.min()["miny"]
        xmax = gpandas.bounds.max()["maxx"]
        ymax = gpandas.bounds.max()["maxy"]
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

        grid_z2 = griddata(XY, Height, (grid_x, grid_y), method='linear')
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
        if isinstance(gpandas, geopandas.GeoDataFrame):
            rasterMap = self.rasterizeGeopandas(gpandas, dxdy=dxdy)

        elif isinstance(gpandas, pandas.DataFrame) or isinstance(gpandas, dask.dataframe.DataFrame):
            rasterMap = self.rasterizePandas(gpandas,dxdy=dxdy)


        return self.rasterToSTL(grid_x=rasterMap['x'],
                                    grid_y=rasterMap['y'],
                                    grid_z=rasterMap['height'],
                                    solidName=solidName)

    def rasterizePandas(self, gpandas, dxdy=50.,xColumn="x",yColumn="y", heightColumn="height"):
        """
            Gets a shape file of topography.
            each contour line has property 'height'.
            Converts it to equigrid xy mesh and then build the STL.
        """

        # 1. Convert contour map to regular height map.
        # 1.1 get boundaries
        xmin = gpandas[xColumn].min()
        xmax = gpandas[xColumn].max()

        ymin = gpandas[yColumn].min()
        ymax = gpandas[yColumn].max()

        # 1.2 build the mesh.
        grid_x, grid_y = numpy.mgrid[(xmin):(xmax):dxdy, (ymin):(ymax):dxdy]
        # 2. Get the points from the geom
        Nx = int(((xmax - xmin) / dxdy))
        Ny = int(((ymax - ymin) / dxdy))
        grid_z2 = coordinateHandler.regularizeTimeSteps(data=gpandas, fieldList=[heightColumn],
                                              coord1=xColumn,
                                              coord2=yColumn,
                                              n=(Nx, Ny), addSurface=False, toPandas=False)[0][heightColumn]
        grid_z2 = self._organizeGrid(grid_z2)
        return {"x":grid_x,"y":grid_y,"height":grid_z2.values}

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

