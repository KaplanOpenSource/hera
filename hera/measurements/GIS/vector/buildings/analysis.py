from collections import OrderedDict
from itertools import product
import io
import geopandas
import geopandas as gpd
import matplotlib as mpl
import numpy
import numpy as np
import pandas
import pandas as pd
import shapely.wkt
from shapely.geometry import box, Polygon


class analysis():

    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, dataLayer):
        self._datalayer = dataLayer

    def ConvexPolygons(self, regionNameOrData, buffer=100):
        """
            Returns convex polygons of groups of buildings.

        Parameters
        ----------
        :param data: str or geopandas DataFrame.
                The data to get the convex of.

        :param buffer:

        Returns
        -------

        """
        if isinstance(regionNameOrData,str):
            data = self.getRegions(regionNameOrData).getData()
        else:
            data = regionNameOrData

        data = data.reset_index().buffer(buffer)
        indicelist = [[0]]
        for i in range(1, len(data)):
            found = False
            for g in range(len(indicelist)):
                for n in indicelist[g]:
                    if data[i].intersection(data[n]).is_empty:
                        continue
                    else:
                        indicelist[g].append(i)
                        found = True
                        break
                if found:
                    break
                if g == len(indicelist) - 1:
                    indicelist.append([i])

        geo = data.loc[indicelist[0]].unary_union.convex_hull
        gpd = geopandas.GeoDataFrame.from_dict([{"geometry": geo, "area": geo.area}])
        for indice in indicelist[1:]:
            geo = data.loc[indice].unary_union.convex_hull
            gpd = pandas.concat([gpd, geopandas.GeoDataFrame.from_dict([{"geometry": geo, "area": geo.area}])])

        gpd = gpd.sort_values(by="area", ascending=False).reset_index()
        found = False
        for i in range(len(gpd)):
            for j in range(i + 1, len(gpd)):
                if gpd.loc[i].geometry.intersection(gpd.loc[j].geometry).is_empty:
                    continue
                else:
                    found = True
                    break
            if found:
                break
        if found:
            gpd = self.ConvexPolygons(gpd, buffer=1)

        return gpd

    def LambdaOfDomain(self, wind,resolution,buildingsDataSourceNameOrData = None,exteriorBlockNameOrData = None, crs = None):
        """
        Calculate average λf and λp of a rectangular domain from BNTL buildings layer.

        Parameters
        ----------
        wind: float
              The meteorological wind direction.
        resolution: int
                    The size of the block in meter.
        newDataSourceName: string
                    If buildingsDataSourceNameOrData is Data, check if the newDataSourceName is not None and save buildingsDataSourceNameOrData
                    as dataSource in project by the name newDataSourceName

        buildingsDataSourceNameOrData: str or GeoDataframe
        exteriorBlockNameOrData: GeoDataFrame,str,list,polygon

        Returns
        -------
        GeoDataFrame of the grid geometry with the calculated parameters: lambdaP,lambdaF,hc
        """

        extent = None
        data = None

        if buildingsDataSourceNameOrData is not None:

            if isinstance(buildingsDataSourceNameOrData,str):

                data = self._datalayer.getDatasourceData(buildingsDataSourceNameOrData)
                if data is None:
                    data = geopandas.read_file(io.StringIO(buildingsDataSourceNameOrData))
            else:
                data = buildingsDataSourceNameOrData
                if isinstance(data, gpd.GeoDataFrame):
                    if data.crs is None:
                        if crs is None:
                            raise ValueError(f" You must provide the crs parameter of buildingsDataSourceNameOrData ")
                        else:
                            data.crs = crs
                else:
                    raise ValueError(f" buildingsDataSource must be GeoDataFrame ")

        if exteriorBlockNameOrData is not None:

            if  isinstance(exteriorBlockNameOrData,str):
                    extent = self._datalayer.getRegionData(exteriorBlockNameOrData)
            else:
                    extent = self._datalayer._setGeoPandasFromRegionData(exteriorBlockNameOrData)

            if isinstance(extent, gpd.GeoDataFrame):
                if extent.crs is not None:
                    crs = extent.crs
                elif crs is None:
                    raise ValueError(f" You must provide the crs parameter of exteriorBlockNameOrData ")

                extent = gpd.GeoDataFrame({'geometry': [box(*(extent.total_bounds))]},crs=crs)

            else:
                raise ValueError(f" exteriorBlockNameOrData must be: GeoDataFrame/GeoJSON/list/Polygon/string")

        if data is not None and extent is not None:
            data = self._datalayer.cutRegionFromSource(extent, dataSourceName=data, crs=crs)

        elif extent is not None:
                data = self._datalayer.cutRegionFromSource(extent,dataSourceName="BNTL")

        elif data is not None:
            extent = gpd.GeoDataFrame({'geometry': [box(*(data.total_bounds))]}, crs=data.crs)

        else:
            raise ValueError(f" You must provide one or more of the followng parameters: buildingsDataSourceNameOrData ,exteriorBlockNameOrData ")

        domainLambda = Blocks(level=0, df=extent, size=resolution).iterBlocks(size=resolution).Lambda(data, windDirection=wind)

        return domainLambda


class Blocks(object):
    """
        Splite the domaine into blocks.
        Use BNTL to calculate urban parameters for each Block.
    """
    _Level = None
    _ExteriorBlock = None
    _Division = None  # will be a map axis->{type: , parameters}
    _DivisionType = None  # a map of axis type
    _Df = None
    _Buildings = None # The BNTL map of the Block
    _listOfDicts = None
    _hc = 0 # Average buildings height of the Block
    _LambdaF = 0 # Lambda F of the Block
    _LambdaP = 0 # Lambda P of the Block
    _blockArea = 0 # Total erea of the Block

    def _SplitFunction(self, min, max, axis):
        funcDict = {"size": numpy.arange, "count": (lambda x, y, z: numpy.linspace(x, y, z)[:-1])}
        coords = funcDict[self._DivisionType[axis]](min, max, self._Division[axis])
        dx = (coords[1] - coords[0]) if len(coords) > 1 else (max - min)

        return coords, dx

    def __init__(self, level, df, exteriorBlock=None, **kwargs):
        """
            Initializes a block object that calculates the properties for a certain domain in the map.

        Parameters
        ----------
        level: int
                The level of the block (for nested blocks calulations). The first block is level 0.
        df:  geopandas
                The data to use
        exteriorBlock: Block
                Reference to the parent block. The exteriorBlock of the root (lowest level) is None.
        kwargs:
                parmeters that determine how to divide a block.

                - size: The size in meter of the block in the x,y direction.
                        Using this feature disallows using npxy.
                - npxy: The number of blocks in x,y directions.
                - width: the width of each block in meter
                        Using this feature disallows using npy.
                - height: the height of each block in meter
                -
        """
        self._Level = level
        self._ExteriorBlock = exteriorBlock
        self._Df = df
        self._Division = {}
        self._DivisionType = {}
        self._Buildings = gpd.GeoDataFrame()
        self._listOfDicts = None

        if "size" in kwargs:
            if "npxy" in kwargs:
                raise ValueError("Got size with npxy")

            self._Division["x"] = kwargs["size"]
            self._Division["y"] = kwargs["size"]
            self._DivisionType = {"x": "size", "y": "size"}

        elif "npxy" in kwargs:
            self._Division["x"] = kwargs["npxy"] + 1
            self._Division["y"] = kwargs["npxy"] + 1
            self._DivisionType = {"x": "count", "y": "count"}

        else:
            if "width" in kwargs:
                self._DivisionType["x"] = "size"
                self._Division["x"] = kwargs["width"]
                if "npx" in kwargs:
                    raise ValueError("Got width with npx")
            elif "npx" in kwargs:
                self._DivisionType["x"] = "count"
                self._Division["x"] = kwargs["npx"] + 1

            if "height" in kwargs:
                self._DivisionType["y"] = "size"
                self._Division["y"] = kwargs["height"]
                if "npy" in kwargs:
                    raise ValueError("Got height with npy")
            elif "npy" in kwargs:
                self._DivisionType["y"] = "count"
                self._Division["y"] = kwargs["npy"] + 1

        if any([x not in self._DivisionType for x in ["x", "y"]]):
            raise ValueError("Either x or y axis is not set.")

    def _GetBlocks(self):
        for x in self._listOfDicts:
            yield x

    # Creates a dictionary of
    def _BuildIndexList(self):
        self._listOfDicts = []
        enum = lambda L: [x for x in enumerate(L)]
        currentLevel = self._Level
        totalBounds = self._Df.total_bounds

        if self._ExteriorBlock is None:
            X, width = self._SplitFunction(totalBounds[0], totalBounds[2], "x")
            Y, height = self._SplitFunction(totalBounds[1], totalBounds[3], "y")

            for ((i, xMin), (j, yMin)) in product(enum(X), enum(Y)):
                xMax, yMax = min(xMin + width, totalBounds[2]), min(yMin + height, totalBounds[3])

                currenDict = {'i0': [i], 'j0': [j], 'xMin0': [xMin], 'yMin0': [yMin], 'xMax0': [xMax], 'yMax0': [yMax]}
                self._listOfDicts.append(currenDict)
        else:
            exteriorBlockListOfDicts = self._ExteriorBlock._BuildIndexList()._listOfDicts
            for extDict in exteriorBlockListOfDicts:

                exteriorXMax, exteriorYMax = extDict['xMax%s' % (currentLevel - 1)][0], \
                                             extDict['yMax%s' % (currentLevel - 1)][0]
                exteriorXMin, exteriorYMin = extDict['xMin%s' % (currentLevel - 1)][0], \
                                             extDict['yMin%s' % (currentLevel - 1)][0]

                X, width = self._SplitFunction(exteriorXMin, exteriorXMax, "x")
                Y, height = self._SplitFunction(exteriorYMin, exteriorYMax, "y")

                for ((i, xMin), (j, yMin)) in product(enum(X), enum(Y)):
                    currentInnerDict = dict(extDict)
                    xMax, yMax = min(xMin + width, exteriorXMax), min(yMin + height, exteriorYMax)

                    currentInnerDict.update({'i%s' % currentLevel: [i], 'j%s' % currentLevel: [j],
                                             'xMin%s' % currentLevel: [xMin], 'yMin%s' % currentLevel: [yMin],
                                             'xMax%s' % currentLevel: [xMax], 'yMax%s' % currentLevel: [yMax]})
                    self._listOfDicts.append(currentInnerDict)

        return self

    def iterBlocks(self, **kwargs):
        return Blocks(level=(self._Level + 1), exteriorBlock=self, df=self._Df, **kwargs)._BuildIndexList()

    def initBuildingsBlock(self, blockDict):

        extent = box(blockDict['xMin%s' % self._Level][0], blockDict['yMin%s' % self._Level][0],
                           blockDict['xMax%s' % self._Level][0], blockDict['yMax%s' % self._Level][0])
        boxToIntersect = gpd.GeoDataFrame([extent], columns=['geometry'])
        self._Buildings.crs = boxToIntersect.crs
        buildings = gpd.overlay(self._Buildings, boxToIntersect, how='intersection')
        OnebuildingsBlock = Blocks(level=0, df=boxToIntersect, size=200)  # buildings,box
        OnebuildingsBlock._Buildings = buildings
        OnebuildingsBlock._ExteriorBlock = extent
        OnebuildingsBlock._blockArea = extent.area

        return OnebuildingsBlock

    def _LambdaP(self):
        if self._Buildings.empty:
            return 0

        Map_A_p = 0
        i = 0
        area_mltp_h = 0
        sum_area_mltp_h = 0
        indexes = self._Buildings['BLDG_HT'].index
        numberOfBld = len(indexes)

        for i in indexes:
            area = 0
            if (self._Buildings['FTYPE'][i] == 16) or (self._Buildings['FTYPE'][i] == 14) or (
                    self._Buildings['BLDG_HT'][i] == 0):
                Map_A_p = Map_A_p
            else:
                if self._Buildings.geometry[i].geom_type == 'MultiPolygon':
                    for poly in self._Buildings.geometry[i]:
                        Map_A_p = Map_A_p + poly.area
                        area = area + poly.area
                else:
                    Map_A_p = Map_A_p + self._Buildings.geometry[i].area
                    area = self._Buildings.geometry[i].area

            area_mltp_h = area * self._Buildings['BLDG_HT'][i]
            sum_area_mltp_h = sum_area_mltp_h + area_mltp_h

        lambdaP = Map_A_p / self._blockArea
        if (Map_A_p != 0):
            self._hc = sum_area_mltp_h / Map_A_p

        return lambdaP

    def _A_f(self, buildingGeometry, height, windDirection):
        A_f = 0
        g = gpd.GeoSeries([buildingGeometry])

        gdf = gpd.GeoDataFrame(geometry=g)
        gdf['angle'] = [windDirection]

        for index, row in gdf.iterrows():
            rotated = Polygon(shapely.affinity.rotate(row['geometry'], row['angle']))
            gdf.loc[index, 'geometry'] = rotated
            bounds = rotated.bounds
            xMin = bounds[0]
            xMax = bounds[2]
            A_f = (xMax - xMin) * height
        return A_f

    def _LambdaF(self, windDirection=270):
        """
        # Calculate average lambda F of a block
        """
        bldHeightToReduce = 0
        errorBuildings = 0
        if self._Buildings.empty:
            return 0
        Map_A_f = 0
        indexes = self._Buildings['BLDG_HT'].index
        numberOfBld = len(indexes)

        for i in indexes:
            if (self._Buildings['FTYPE'][i] == 16) or (self._Buildings['FTYPE'][i] == 14):
                bldHeightToReduce = bldHeightToReduce + self._Buildings['BLDG_HT'][i]
                numberOfBld = numberOfBld - 1
            else:
                bldHeight = self._Buildings['BLDG_HT'][i]
                if bldHeight < 2:
                    bldHeight = self._Buildings['HI_PNT_Z'][i] - self._Buildings['HT_LAND'][i]
                    if bldHeight < 2:
                        self._Buildings.at[i, 'BLDG_HT'] = 0
                        errorBuildings = errorBuildings + 1
                        continue
                    else:
                        self._Buildings.at[i, 'BLDG_HT'] = bldHeight

                if self._Buildings.geometry[i].geom_type == 'MultiPolygon':
                    for poly in self._Buildings.geometry[i]:
                        Map_A_f = Map_A_f + self._A_f(Polygon(poly), bldHeight, windDirection)
                else:

                    Map_A_f = Map_A_f + self._A_f(Polygon(self._Buildings.geometry[i]), bldHeight, windDirection)

        lambda_f = Map_A_f / self._blockArea
        # totalBldHeight = self._blockBuildings['BLDG_HT'].sum() - bldHeightToReduce
        return lambda_f

    def getHc(self):

        if self._LambdaP is not None:

            return self._hc
        else:

            return None

    def Lambda(self, buildings, windDirection=270):
        """
        This method calculates average Lambda P and F of each block in the domain.

        :return: geopandas with the calculated average lambda.
        """
        listOfBuildingsBlock = []
        self._Buildings = buildings
        currentDict = {'lambdaP': [], 'lambdaF':[], 'hc': [], 'geometry': []}


        for i,blockDict in enumerate( self._GetBlocks()):
            # currentPandas = pandas.DataFrame.from_dict(blockDict)
            listOfBuildingsBlock.append(self.initBuildingsBlock(blockDict))
            currentDict['lambdaF'].append(listOfBuildingsBlock[i]._LambdaF(windDirection=windDirection))
            currentDict['lambdaP'].append(listOfBuildingsBlock[i]._LambdaP())
            currentDict['hc'].append(listOfBuildingsBlock[i].getHc()) # The calculation of HC is done at LambdaP function
            currentDict['geometry'].append(listOfBuildingsBlock[i]._ExteriorBlock)


        df = pd.DataFrame.from_dict(currentDict, orient='columns')
        return gpd.GeoDataFrame(df, geometry=df['geometry'])




#divide blocks
