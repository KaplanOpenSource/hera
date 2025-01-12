import geopandas
import io
import scipy
import numpy
import math
import os
from scipy.interpolate import griddata
from shapely.geometry import LineString
import dask
from numpy import array, sqrt
import pandas

from . import toolkit
from ....simulations.utils import coordinateHandler
from ....utils.logging import get_classMethod_logger

from ....toolkit import TOOLKIT_SAVEMODE_ONLYFILE
from .toolkit import VectorToolkit
from ..utils import stlFactory

class TopographyToolkit(VectorToolkit):

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
                logger.debug("The region was defined without crs. Using the crs of the datasource.")
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


