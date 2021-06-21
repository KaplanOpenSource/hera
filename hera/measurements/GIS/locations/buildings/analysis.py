from collections import OrderedDict
from itertools import product

import geopandas
import geopandas as gpd
import matplotlib as mpl
import numpy
import numpy as np
import pandas
import pandas as pd
import shapely.wkt
# from .building import Building as bld
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

    def LambdaGrid(self, regionNameOrData,wind,rez):
        """
         Calculate average λf and λp of a rectangular domain from BNTL buildings layer.

         Parameter:
                   buildings: str or GeoDataframe of rectangular area of the buildings layer of BNTL.
                   wind: The meteorological wind direction.
                   rez: The cell's resolution of the average lambda grid(output).
         Return:  geoPandas of lambda grid
                  Save the geoPandas to excel file('output_lambda.xlsx')

        """
        if isinstance(regionNameOrData,str):
            data = self.datalayer.getRegions(regionNameOrData).getData()
        else:
            data = regionNameOrData

        nf = field(data, rez, wind, self.datalayer.FilesDirectory)
        lambdaGrid = nf.getLambdaGrid()
        #nf.dataToExcel('output_lambda')
        return lambdaGrid


class field(object):
    """
    A class to calculate lambda F+P of a domain by creating grid cells from
    the domain and calculate the average lambda F+P and average buildings height of each cell.
    """
    _buildings = None
    _buildingsGrid = list()
    _lambdaGrid = list()
    _rez = None
    _wind = None
    _GridMedianHeight = list()
    _GridMeanHeight = list()
    _GridBuildingHeightArray = list()
    _ErrorBLD  = list()

    def __init__(self, buildings, rez, wind, directory, Frame=None):
        self._buildingsGrid = []
        self._GridMedianHeight = []
        self._GridMeanHeight = []
        self._GridBuildingHeightArray = []
        self._ErrorBLD = []
        self._buildings = buildings
        self._rez = rez
        self._wind = wind
        self._path = directory
        self.CreateNewField(Frame)

    @classmethod
    def setFieldFromFile(cls, FileName, buildings, rez, wind):
        cls._rez = rez
        cls._wind = wind
        cls._buildings = buildings
        cls._lambdaGrid = []
        cls.dataFromExcel(cls, FileName)
        return cls

    def validate(self):
        return True

    # Splite the GeoDataFrame buildings map to blockes in rezolution (=rez) and calculate lambdaF+P for each block
    def CreateNewField(self, Frame=None):  # Use Fram for resolution test

        # Parameters
        Grid = []
        lambda_f = []
        lambda_p = []
        built = []
        unbuilt = []
        errorBLDArea = []
        numberOfBld = []
        totalBldHeight = []
        errorBLDNum = []
        totalBldAreaLessThanTwo = []
        numberOfBldLessThanTwo = []
        hc = []

        # Create the Grid
        f = self.generateBlocks(Frame)
        grid = f._BuildIndexList()

        # Calculate lambda for each block in grid

        for i in range(len(grid)):
            # initize block by grid frame
            block, ibox = self.initBlock(Grid, grid, i)

            # lambda caculations for each block
            self.processBlock(Grid, block, built, i, lambda_f, lambda_p, unbuilt)

            # orgenize all calculated parameters
            self.buildDataBox(block, errorBLDArea, errorBLDNum, ibox, numberOfBld, numberOfBldLessThanTwo,
                              totalBldAreaLessThanTwo, totalBldHeight, hc)

        # building with heihgt  = 0
        ZERO_BLD = {'geometry': self._ErrorBLD}
        self._ErrorBLD = gpd.GeoDataFrame(ZERO_BLD)

        # Insert calculated data to GeoDataFrame
        d = {'lambda_f': lambda_f, 'lambda_p': lambda_p, 'built': built, 'unbuilt': unbuilt, 'errorArea': errorBLDArea,
             'MedianHeight': self._GridMedianHeight, 'MeanHeight': self._GridMeanHeight, 'numberOfBld': numberOfBld,
             'totalBldHeight': totalBldHeight, 'errorBLDNum': errorBLDNum,
             'totalAreaLessThanTwo': totalBldAreaLessThanTwo,
             'numberOfBldLessThanTwo': numberOfBldLessThanTwo, 'hc': hc, 'geometry': Grid}

        self._lambdaGrid = gpd.GeoDataFrame(d, crs="EPSG:2039")

    # use :def CreateNewField - orgenize all calculated one block parameters
    def buildDataBox(self, block, errorBLDArea, errorBLDNum, ibox, numberOfBld, numberOfBldLessThanTwo,
                     totalBldAreaLessThanTwo, totalBldHeight, hc):

        errors = gpd.GeoDataFrame({'geometry': block.ErrorBuildings})
        errorBLDNum.append(len(errors.index))
        errors['area'] = errors.area
        errorArea = errors['area'].sum()
        errorBLDArea.append(errorArea)
        numberOfBld.append(block._numberOfBld)
        totalBldHeight.append(block._totalBldHeight)
        hc.append(block.getHc())
        self._ErrorBLD.append(block.ErrorBuildings)
        totalBldAreaLessThanTwo.append(block._totalBldAreaLessThanTwo)
        numberOfBldLessThanTwo.append(block._numberOfBldLessThanTwo)
        self._buildingsGrid.append(ibox)
        self.BuildingsMedianHeight(ibox)  ## correct to mean without errors
        self.BuildingsMeanHeight(ibox)  ## correct to mean without errors

    # use :def CreateNewField - one block lambda calculation
    def processBlock(self, Grid, block, built, i, lambda_f, lambda_p, unbuilt):

        lambda_f.append(block.LambdaF(self._wind, Grid[i].area))
        lamP, A_p, area = block.LambdaP(Grid[i].area)
        if (area - A_p) < 0:  ## print error message
            built.append(area)
            unbuilt.append(0)
        else:
            built.append(A_p)
            unbuilt.append(area - A_p)
        lambda_p.append(lamP)

    # use :def CreateNewField - Set the block building map by grid frame
    def initBlock(self, Grid, grid, i):

        iFrame = grid[i]
        Grid.append(box(iFrame['xMin0'][0], iFrame['yMin0'][0], iFrame['xMax0'][0], iFrame['yMax0'][0], 1))
        iGeoFrame = gpd.GeoDataFrame([Grid[i]], columns=['geometry'])
        self._buildings.crs = iGeoFrame.crs
        ibox = gpd.overlay(self._buildings, iGeoFrame, how='intersection')
        self._GridBuildingHeightArray.append(ibox['BLDG_HT'])
        block = blk(ibox)
        return block, ibox

    # use :def CreateNewField - create field grid
    def generateBlocks(self, Frame):

        if Frame != None:
            GeoFrame = gpd.GeoDataFrame([Frame], columns=['geometry'])
            f = blks(level=0, df=GeoFrame, size=self._rez)
        else:
            f = blks(level=0, df=self._buildings, size=self._rez)
        return f

    def plotBuildings(self):
        fig, ax = plt.subplots(figsize=(1, 1))
        sns.set()
        l = sns.color_palette("coolwarm", 6)
        cmap = mpl.colors.ListedColormap(l)
        # cmap = mpl.colors.ListedColormap(['cyan','blue','yellowgreen', 'orange','red'])
        bounds = [0, 5, 10, 20, 40, 60, 200]
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        fig.colorbar(
            mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
        )
        self._buildings.plot(ax=ax, alpha=0.8, cmap=cmap, edgecolor='gray')
        ax.set_title("תל-אביב: גובה מבנים", fontsize=20)
        self._lambdaGrid.plot(ax=ax, legend=True, alpha=0.1, edgecolor='black')
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
        plt.show()
        return fig

    def plotMap(self, data):

        cmaps = OrderedDict()
        fig, ax = plt.subplots(1, 1)
        self._buildings.plot(ax=ax)
        self._lambdaGrid.plot(column=data, ax=ax, legend=True, alpha=0.5, edgecolor='black', vmin=0., vmax=0.6,
                              cmap='viridis')
        # self._ErrorBLD.plot(ax=ax,alpha=0.7,edgecolor='red')
        # for x, y, label in zip(self._lambdaGrid.geometry.centroid.x, self._lambdaGrid.geometry.centroid.y, self._lambdaGrid[data]):
        # ax.annotate(label, xy=(x, y), xytext=(3, 3), textcoords="offset points")
        m = plt.cm.ScalarMappable(cmap='viridis')
        m.set_clim(0., 0.6)
        ax.set_title(data + "/" + str(self._rez) + "Rez/Wind" + str(self._wind), fontsize=20)
        return fig
        # plt.show()

    def BuildingsMedianHeight(self, box):

        if box.empty:
            self._GridMedianHeight.append(0)
        else:
            medianHeight = (np.median(box['BLDG_HT']))
            self._GridMedianHeight.append(medianHeight)

    def BuildingsMeanHeight(self, box):

        if box.empty:
            self._GridMeanHeight.append(0)
        else:
            meanHeight = (np.mean(box['BLDG_HT']))
            self._GridMeanHeight.append(meanHeight)

    def dataToExcel(self, XlName):
        with pd.ExcelWriter(self._path + XlName + '.xlsx') as writer:
            self._lambdaGrid.to_excel(writer, sheet_name='Sheet_name_1')

    def dataFromExcel(self, XlName):
        geometryT = []
        FILE_NAME = self._path + XlName
        df = pd.read_excel(FILE_NAME, index_col=0, header=0)
        c = df['geometry']
        for poly in c:
            geometryT.append(shapely.wkt.loads(poly))
        df.drop('geometry', axis='columns')
        self._lambdaGrid = gpd.GeoDataFrame(df, geometry=geometryT)
        return self

    def getLambdaGrid(self):

        return self._lambdaGrid


class Blocks(object):
    """
        Splite the domaine into blocks.
    """
    _Level = None
    _ExteriorBlock = None
    _Division = None # will be a map axis->{type: , parameters}
    _DivisionType = None # a map of axis type
    _Df = None

    def _SplitFunction(self, min, max, axis):
        funcDict = { "size" : numpy.arange, "count" : (lambda x, y, z : numpy.linspace(x, y, z)[:-1])}
        coords = funcDict[self._DivisionType[axis]](min, max, self._Division[axis])
        dx = (coords[1] - coords[0]) if len(coords) > 1 else (max - min)

        return coords, dx

    def __init__(self, level, df, exteriorBlock = None, **kwargs):
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

        if "size" in kwargs:
            if "npxy" in kwargs:
                raise ValueError("Got size with npxy")

            self._Division["x"] = kwargs["size"]
            self._Division["y"] = kwargs["size"]
            self._DivisionType = {"x" : "size", "y" : "size"}

        elif "npxy" in kwargs:
            self._Division["x"] = kwargs["npxy"] + 1
            self._Division["y"] = kwargs["npxy"] + 1
            self._DivisionType = {"x" : "count", "y" : "count"}

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
        for x in self._BuildIndexList():
            yield x

    # Creates a dictionary of
    def _BuildIndexList(self):
        listOfDicts = []
        enum = lambda L: [x for x in enumerate(L)]
        currentLevel = self._Level
        totalBounds = self._Df.total_bounds

        if self._ExteriorBlock is None:
            X, width = self._SplitFunction(totalBounds[0], totalBounds[2], "x")
            Y, height = self._SplitFunction(totalBounds[1], totalBounds[3], "y")

            for ((i, xMin), (j, yMin)) in product(enum(X), enum(Y)):
                xMax, yMax = min(xMin + width, totalBounds[2]), min(yMin + height, totalBounds[3])
                currenDict = {'i0': [i], 'j0': [j], 'xMin0': [xMin], 'yMin0': [yMin], 'xMax0' : [xMax], 'yMax0' : [yMax]}
                listOfDicts.append(currenDict)
        else:
            exteriorBlockListOfDicts = self._ExteriorBlock._BuildIndexList()
            for extDict in exteriorBlockListOfDicts:

                exteriorXMax, exteriorYMax = extDict['xMax%s' % (currentLevel - 1)][0], extDict['yMax%s' % (currentLevel - 1)][0]
                exteriorXMin, exteriorYMin = extDict['xMin%s' % (currentLevel - 1)][0], extDict['yMin%s' % (currentLevel - 1)][0]

                X, width  = self._SplitFunction(exteriorXMin, exteriorXMax, "x")
                Y, height  = self._SplitFunction(exteriorYMin, exteriorYMax, "y")

                for ((i, xMin), (j, yMin)) in product(enum(X), enum(Y)):
                    currentInnerDict = dict(extDict)
                    xMax, yMax = min(xMin + width, exteriorXMax), min(yMin + height, exteriorYMax)

                    currentInnerDict.update({'i%s' % currentLevel: [i], 'j%s' % currentLevel: [j],
                                             'xMin%s' % currentLevel: [xMin], 'yMin%s' % currentLevel: [yMin],
                                             'xMax%s' % currentLevel: [xMax], 'yMax%s' % currentLevel: [yMax]})
                    listOfDicts.append(currentInnerDict)

        return listOfDicts

    def iterBlocks(self, **kwargs):
        return Blocks(level = (self._Level + 1), exteriorBlock = self, df = self._Df, **kwargs)


class Block:

    _block = None
    lambda_f = None
    lambda_p = None
    ErrorBuildings = []
    _numberOfBld = 0
    _totalBldHeight = 0
    _totalBldAreaLessThanTwo = 0
    _numberOfBldLessThanTwo = 0
    _hc = 0


    def __init__(self, box):
        self._hc = 0
        self.BLDG_HT = []
        self.ErrorBuildings = []
        self._block = box
        self.BLDG_HT = box['BLDG_HT']

    #
    # def LambdaF(self, MeteoAngle ,TotalArea=None):
    #     """
    #     # Calculate average lambda F of a block
    #     """
    #     bldHeightToReduce = 0
    #     if self._block.empty:
    #         return 0
    #     Map_A_f = 0
    #     i = 0
    #     j = 0
    #     indexes = self._block['BLDG_HT'].index
    #     self._numberOfBld = len(indexes)
    #
    #     for i in indexes:
    #         if (self._block['FTYPE'][i] == 16) or (self._block['FTYPE'][i] == 14):
    #             bldHeightToReduce = self.updateNAbuilding(bldHeightToReduce, i)
    #         else:
    #                 bldHeight = self._block['BLDG_HT'][i]
    #                 if bldHeight < 2:
    #                     bldHeight = self._block['HI_PNT_Z'][i]-self._block['HT_LAND'][i]
    #                     if bldHeight < 2:
    #                         self.updateErrorBuilding(i)
    #                         j = j + 1
    #                         continue
    #                     else:
    #                         self._block.at[i,'BLDG_HT'] = bldHeight
    #
    #                 if self._block.geometry[i].geom_type == 'MultiPolygon':
    #                     for poly in self._block.geometry[i]:
    #                         building = bld(Polygon(poly),bldHeight)
    #                         Map_A_f = Map_A_f + building.A_f2(MeteoAngle)
    #                 else:
    #                     building = bld(self._block.geometry[i], bldHeight)
    #                     Map_A_f = Map_A_f + building.A_f2(MeteoAngle)
    #
    #         j = j+1
    #
    #     if TotalArea == None:
    #         bounds = self._block.total_bounds
    #         TotalArea = ((bounds[2] - bounds[0])) * ((bounds[3] - bounds[1]))
    #
    #     self.lambda_f = Map_A_f / TotalArea
    #     self._totalBldHeight = self._block['BLDG_HT'].sum() - bldHeightToReduce
    #     return self.lambda_f

    def updateErrorBuilding(self, index):
        """
        # Update data when building has no height or is less than 2 meter height
        """
        self.ErrorBuildings.append(self._block.geometry[index])
        self._block.at[index, 'BLDG_HT'] = 0
        self._totalBldAreaLessThanTwo = self._totalBldAreaLessThanTwo + self._block.geometry[index].area
        self._numberOfBldLessThanTwo = self._numberOfBldLessThanTwo + 1

    def updateNAbuilding(self, bldHeightToReduce, index):
        """
        # Update data when object is not a building or does not meet the criterion of a wind obstacle(see BNTL types of structures)
        """
        self._numberOfBld = self._numberOfBld - 1
        bldHeightToReduce = bldHeightToReduce + self._block['BLDG_HT'][index]
        return bldHeightToReduce

    def LambdaP(self, TotalArea=None):
        """

        """
        if self._block.empty:
            return 0 , 0 , TotalArea
        Map_A_p = 0
        sum_area_mltp_h = 0
        indexes = self._block['geometry'].index

        for i in indexes:
            area = 0
            if (self._block['FTYPE'][i] == 16) or (self._block['FTYPE'][i] == 14) or (self._block['BLDG_HT'][i]==0):
                    Map_A_p = Map_A_p
            else:
                if self._block.geometry[i].geom_type == 'MultiPolygon':
                    for poly in self._block.geometry[i]:
                            Map_A_p = Map_A_p + poly.area
                            area = area+poly.area
                else:
                    Map_A_p = Map_A_p + self._block.geometry[i].area
                    area = self._block.geometry[i].area

            area_mltp_h = area * self._block['BLDG_HT'][i]
            sum_area_mltp_h = sum_area_mltp_h + area_mltp_h

        if TotalArea == None:
            bounds = self._block.total_bounds
            TotalArea = ((bounds[2] - bounds[0])) * ((bounds[3] - bounds[1]))

        self.lambda_p = Map_A_p / TotalArea
        if(Map_A_p != 0 ):
            self._hc = sum_area_mltp_h/Map_A_p
        return self.lambda_p , Map_A_p , TotalArea

    def getHc(self):
        """
        # Return block average buildings height (Normalize to relative of buildings area to total block area)
        """
        return self._hc