from .blocks import Blocks as blks
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon,Point
from shapely.geometry import box
from .block import Block as blk
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import matplotlib as mpl
from collections import OrderedDict
import seaborn as sns
import shapely.wkt



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


def setPolygon(xmin , ymin,sizex,sizey):
    xmax = xmin + sizex
    ymax = ymin + sizey
    polygon = Polygon([(xmin,ymin),(xmax,ymin),(xmax,ymax),(xmin,ymax),(xmin,ymin)])
    return polygon

def cutAreafromBNTL(polygon,map): #### Try new method to reduce runtime

    BNTLBldMap = gpd.read_file(map)
    geoPoly = gpd.GeoDataFrame([polygon],columns=['geometry'])
    BNTLBldMap.crs = geoPoly.crs
    polygonMap = gpd.overlay(BNTLBldMap, geoPoly, how='intersection')
    return polygonMap

if __name__ == "__main__":
    # Parameters:
    xmin = 177787
    ymin = 663894
    sizeX = 1000
    sizeY = 1000
    polygon = setPolygon(xmin, ymin, sizeX, sizeY)
    dataFile = "/data3/hera-data/GIS_Data/BNTL_MALE_ARZI/BNTL_MALE_ARZI/BUILDINGS/BLDG.shp"
    buildings =cutAreafromBNTL(polygon,dataFile)

    path = '/home/project/lambda_Output/'
    rez = 200
    wind = 240

    # Run:
    nf = field(buildings, rez, wind,path)
    nf.dataToExcel('output_test')
    nf.plotMap('lambda_p')
    nf.plotMap('lambda_f')




