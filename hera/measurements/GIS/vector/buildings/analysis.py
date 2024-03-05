from itertools import product
import io
import geopandas
import numpy
import os
import pandas
import pandas as pd
import shapely.wkt
from shapely.geometry import box, Polygon
from hera.datalayer import datatypes,nonDBMetadataFrame
from hera.measurements.GIS.vector import topography


import logging
BUILDINGS_LAMBDA_WIND_DIRECTION = 'wind'
BUILDINGS_LAMBDA_RESOLUTION = 'resolution'

class analysis():

    _datalayer = None

    @property
    def logger(self):
        return self._datalayer.logger

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



    def LambdaFromBuildingData(self, windMeteorologicalDirection, resolution, buildingsData,externalShape=None, overwrite=False):
        """

        Parameters
        ----------
        windMeteorologicalDirection : double
                The meteorological angle

        resolution: double
                The size of the squares.

        buildingsData: geopandas.Geopandas
                Geopandas of the buildings
                The columns are:
                    -
                    -

                The crs of the building data must be set.

        externalShape : str, shape [optional]
                will be used in the crs of the building data.

                If None, then use the boundaries of the buidingData.

        overwrite : bool
            If true, write over data and update the database

        Returns
        -------

        """
#nir        self.logger.info("--- Start ---")

        if isinstance(buildingsData,str):
            #nir            self.logger.debug("Reading the bounds from StringIO")
            data = geopandas.read_file(io.StringIO(buildingsData))
        elif isinstance(buildingsData, geopandas.GeoDataFrame):
            #nir            self.logger.debug("Using the user input")
            data = buildingsData
        else:
            err = f"buildingsData must be str or geopandas.GeoDataFrame. Got {type(buildingsData)}"
            #nir self.logger.error(err)
            raise ValueError(err)

        if data.crs is None:
            err = "The buildingData crs must be set"
            #nir            self.logger.error(err)
            raise ValueError(err)

        if isinstance(externalShape, str):
            bounds= self.datalayer.getRegionData(externalShape)
        else:
            bounds = self.datalayer._RegionToGeopandas(externalShape, crs = data.crs)

        desc = {
                 "bounds" : bounds.total_bounds,
                 BUILDINGS_LAMBDA_WIND_DIRECTION: windMeteorologicalDirection,
                 BUILDINGS_LAMBDA_RESOLUTION : resolution,
                 "crs":data.crs.to_epsg()
               }

        #nir        self.logger.info(f"Check if cached data exists for data {desc}")
        dataDoc = self.datalayer.getCacheDocuments(type="Lambda_Buildings",**desc)

        #nir        if self.logger.isEnabledFor(logging.DEBUG):
        dbgstr = "Cached data Not found"  if len(dataDoc) == 0 else "Found data in the cache"
            #nir            self.logger.debug(dbgstr)

        if len(dataDoc)==0 or overwrite:
            #nir            self.logger.info(f"Calculatings lambda data with bounds {bounds.total_bounds}")
            domainLambda = Blocks(level=0, df=bounds, size=resolution).iterBlocks(size=resolution).Lambda(data, windDirection=windMeteorologicalDirection)

            if len(dataDoc) == 0:
                #nir                self.logger.debug("Adding new record to the DB")
                doc = self.datalayer.addCacheDocument(resource="",
                                                      dataFormat=datatypes.GEOPANDAS,
                                                      type="Lambda_Buildings",
                                                      desc=desc)
                filename = f"{str(doc.id)}.geojson"
                outputfileFull = os.path.abspath(os.path.join(self._datalayer.filesDirectory, filename))
                #nir                self.logger.debug(f"Writing Lambda data to {outputfileFull}")

            else:
                #nir                self.logger.info("Updating old record.")
                doc = dataDoc[0]
                outputfileFull = doc.resource

            domainLambda.to_file(outputfileFull,driver='GeoJSON')
            doc.resource = outputfileFull
            doc.save()

        else:
            #nir            self.logger.debug("Return found data in DB")
            try:
                domainLambda = dataDoc[0].getData()
            except fiona.errors.DriverError:
                errmsg = f"The cached data in location {dataDoc[0].resource} is not found on the disk. Maybe it was removed?. Use overwrite=True to recalculate and update the cache."
                #nir                self.logger.error(errmsg)
                raise FileNotFoundError(errmsg)

        #nir        self.logger.info("---- End ----")
        return domainLambda



    def LambdaFromDatasource(self, windMeteorologicalDirection, resolution, shapeDataOrName,datasourceName, crs=None, overwrite=False):
        """
        Calculate average λf and λp of a rectangular domain from BNTL buildings layer.

        Parameters
        ----------
        windMeteorologicalDirection: float
              The meteorological wind direction.
        resolution: int
                    The size of the block in meter.

        dataSourceName: str
                the name of the buildings datasource.

        shapeDataOrName:  str, polygon or [xmin,ymin,xmax,ymax]
                The definition of the shape to use.
        crs: [optional], int
            if None, use the crs of the datasource.
            if not None, set the


        overwrite: bool
            Calculate the results even if data in the db.
            if it exists in db, update.

        Returns
        -------
        GeoDataFrame of the grid geometry with the calculated parameters: lambdaP,lambdaF,hc
        """
#nir        self.logger.info("--- Start ---")

        buildingsData = self.datalayer.cutRegionFromSource(shapeDataOrName=shapeDataOrName,
                                                           datasourceName=datasourceName, crs=crs)

        if isinstance(shapeDataOrName, str):
            bounds= self.datalayer.getRegionData(shapeDataOrName)
        else:
            bounds = self.datalayer._RegionToGeopandas(shapeDataOrName, crs = crs)

        domainLambda = self.LambdaFromBuildingData(windMeteorologicalDirection=windMeteorologicalDirection,
                                                   resolution=resolution,
                                                   buildingsData=buildingsData,
                                                   externalShape=bounds,
                                                   overwrite=overwrite)

        if crs is not None:
            domainLambda.crs = crs

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
        self._Buildings = geopandas.GeoDataFrame()
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
        """
            Initializes the ... .

        Parameters
        ----------
        blockDict : dict
                The dictionary that ....

        Returns
        -------

        """

        extent = box(blockDict['xMin%s' % self._Level][0], blockDict['yMin%s' % self._Level][0],
                           blockDict['xMax%s' % self._Level][0], blockDict['yMax%s' % self._Level][0])
        boxToIntersect = geopandas.GeoDataFrame([extent], columns=['geometry'])
        self._Buildings.crs = boxToIntersect.crs
        buildings = geopandas.overlay(self._Buildings, boxToIntersect, how='intersection',keep_geom_type=True)
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
#        numberOfBld = len(indexes)

        for i in indexes:
            area = 0
            # if there is no land height data, than the building height will be incorrect so we will take the data from the nearby building
            if ((self._Buildings['HT_LAND'][i]==0.0) and (self._Buildings['BLDG_HT'][i]>0.0)):
                farest = 99999999999
                farheight=0
                for j in indexes:
                    try:
                        walls = self._Buildings['geometry'][j].exterior.xy
                    except:
                        continue
                    if (self._Buildings['HT_LAND'][j]!=0.0):
                        far = ((self._Buildings['geometry'][i].exterior.xy[0][0]-self._Buildings['geometry'][j].exterior.xy[0][0])**2+
                               (self._Buildings['geometry'][i].exterior.xy[1][0]-self._Buildings['geometry'][j].exterior.xy[1][0])**2)
                        if (farest**2.>far):
                            farest = far
                            farheight= self._Buildings['HT_LAND'][j]
                building_height = max((self._Buildings['BLDG_HT'][i]-farheight),0.0) # don't calculate underground building
            else:
                building_height = self._Buildings['BLDG_HT'][i]


            if (self._Buildings['FTYPE'][i] == 16) or (self._Buildings['FTYPE'][i] == 14) or (
                    self._Buildings['BLDG_HT'][i] == 0) or (building_height==0):
                Map_A_p = Map_A_p
            else:
                if self._Buildings.geometry[i].geom_type == 'MultiPolygon':
                    for poly in self._Buildings.geometry[i].geoms:
                        Map_A_p = Map_A_p + poly.area
                        area = area + poly.area
                else:
                    Map_A_p = Map_A_p + self._Buildings.geometry[i].area
                    area = self._Buildings.geometry[i].area



            area_mltp_h = area * building_height
            sum_area_mltp_h = sum_area_mltp_h + area_mltp_h

        lambdaP = Map_A_p / self._blockArea
        if (Map_A_p != 0):
            self._hc = sum_area_mltp_h / Map_A_p

        return lambdaP

    def _A_f(self, buildingGeometry, height, windDirection):
        A_f = 0
        g = geopandas.GeoSeries([buildingGeometry])

        gdf = geopandas.GeoDataFrame(geometry=g)
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
        # bldHeightToReduce = 0
#        errorBuildings = 0
        if self._Buildings.empty:
            return 0
        Map_A_f = 0
        indexes = self._Buildings['BLDG_HT'].index
#        numberOfBld = len(indexes)

        # topography.analysis.addHeight(self,'data', 'groundData')

        for i in indexes:

            # if there is no land height data, than the building height will be incorrect so we will take the data from the nearby building
            if ((self._Buildings['HT_LAND'][i]==0.0) and (self._Buildings['BLDG_HT'][i]>0.0)):
                farest = 99999999999
                farheight=0
                for j in indexes:
                    try:
                        walls = self._Buildings['geometry'][j].exterior.xy
                    except:
                        continue
                    if (self._Buildings['HT_LAND'][j]!=0.0):
                        far = ((self._Buildings['geometry'][i].exterior.xy[0][0]-self._Buildings['geometry'][j].exterior.xy[0][0])**2+
                               (self._Buildings['geometry'][i].exterior.xy[1][0]-self._Buildings['geometry'][j].exterior.xy[1][0])**2)
                        if (farest**2.>far):
                            farest = far
                            farheight= self._Buildings['HT_LAND'][j]
                building_height = max((self._Buildings['BLDG_HT'][i]-farheight),0.0) # don't calculate underground building
            else:
                building_height = self._Buildings['BLDG_HT'][i]


            if (self._Buildings['FTYPE'][i] == 16) or (self._Buildings['FTYPE'][i] == 14) or (building_height==0):
 #               bldHeightToReduce = bldHeightToReduce + self._Buildings['BLDG_HT'][i]
 #               numberOfBld = numberOfBld - 1
                pass
            else:
                bldHeight = self._Buildings['BLDG_HT'][i]
                if bldHeight < 2:
                    bldHeight = self._Buildings['HI_PNT_Z'][i] - self._Buildings['HT_LAND'][i]
                    if bldHeight < 2:
                        self._Buildings.at[i, 'BLDG_HT'] = 0
#                        errorBuildings = errorBuildings + 1
                        continue
                    else:
                        self._Buildings.at[i, 'BLDG_HT'] = bldHeight

                if self._Buildings.geometry[i].geom_type == 'MultiPolygon':
                    for poly in self._Buildings.geometry[i].geoms:
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

    def Lambda(self, buildings, windDirection):
        """
        This method calculates average Lambda P and F of each block in the domain.

        :return: geopandas with the calculated average lambda.
        """
        listOfBuildingsBlock = []
        self._Buildings = buildings
        currentDict = {'lambdaP': [], 'lambdaF':[], 'hc': [], 'geometry': [],'i0':[],'j0':[]}

        for i,blockDict in enumerate( self._GetBlocks()):
            # currentPandas = pandas.DataFrame.from_dict(blockDict)
            BuildingsBlock = self.initBuildingsBlock(blockDict)
            currentDict['lambdaF'].append(BuildingsBlock._LambdaF(windDirection=windDirection))
            currentDict['lambdaP'].append(BuildingsBlock._LambdaP())
            currentDict['hc'].append(BuildingsBlock.getHc()) # The calculation of HC is done at LambdaP function
            currentDict['geometry'].append(BuildingsBlock._ExteriorBlock)
            currentDict['i0'].append(blockDict['i0'][0])
            currentDict['j0'].append(blockDict['j0'][0])


        df = pd.DataFrame.from_dict(currentDict, orient='columns')
        return geopandas.GeoDataFrame(df, geometry=df['geometry'])


if __name__ == "__main__":
    from hera import toolkitHome
    bt = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_BUILDINGS,projectName="testbamba") # tlvbig
    tlvbounding = [175000, 658000, 185000, 668000]  #
    tlvbounding = [175500, 658000, 185000, 668000]
    bsbounding = [175000, 569000, 189000, 579000]
    bsboundingsmall = [181000, 577000, 182000, 578000]
    ashkelonbounding = [156000, 616000, 164000, 625000]
    natanyabounding = [184000, 689000, 192000, 697000]

    bounding = tlvbounding
    cityname = 'tlv1'

    bounding = bsboundingsmall
    cityname = 'bssm'

    bt.addRegion(bounding, cityname, crs=2039)
    if 5 == 5:
        reg = bt.cutRegionFromSource(cityname, datasourceName='BNTL', isBounds=True, crs=2039)
        #	    bt.regionToSTL(cityname,cityname+'-buildings.stl','BNTL')
        print('dddeeebbb')
        lm = bt.analysis.LambdaFromDatasource(270, 250, reg, 'BNTL', crs=2039, overwrite=True)
        print(lm)
        file = open(cityname + '-lambda1.csv', 'w')
        file.writelines('[')
        for i in range(len(lm)):
            lines = ['[', str(min(lm.iloc[i]['geometry'].exterior.coords.xy[0])), ', ',
                     str(max(lm.iloc[i]['geometry'].exterior.coords.xy[0])), ', ',
                     str(min(lm.iloc[i]['geometry'].exterior.coords.xy[1])), ', ',
                     str(max(lm.iloc[i]['geometry'].exterior.coords.xy[1])), ', ',
                     str(lm.iloc[i]['lambdaF']), ', ',
                     str(lm.iloc[i]['lambdaP']), ', ',
                     str(lm.iloc[i]['hc']), '],\n']
            file.writelines(lines)
        file.writelines(']')
        file.close()

    if 5 == 6:
        lm = bt._analysis.LambdaFromDatasource(270, 250, reg, 'BNTL', crs=2039)
        print(lm)
        file = open(cityname + '-lambda.csv', 'w')
        file.writelines('[')
        for i in range(len(lm)):
            lines = ['[', str(min(lm.iloc[i]['geometry'].exterior.coords.xy[0])), ', ',
                     str(max(lm.iloc[i]['geometry'].exterior.coords.xy[0])), ', ',
                     str(min(lm.iloc[i]['geometry'].exterior.coords.xy[1])), ', ',
                     str(max(lm.iloc[i]['geometry'].exterior.coords.xy[1])), ', ',
                     str(round(lm.iloc[i]['lambdaF'], 3)), ', ',
                     str(round(lm.iloc[i]['lambdaP'], 3)), ', ',
                     str(round(lm.iloc[i]['hc'], 3)), '],\n']
            file.writelines(lines)
        file.writelines(']')
        file.close()

    if 5 == 6:
        bt2 = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_TOPOGRAPHY, projectName="testbamba")
        bt2.addRegion(bounding, cityname, crs=2039)
        # reg = bt2.cutRegionFromSource('bs',datasourceName='BNTL',isBounds = True, crs = 2039)
        topo = bt2.regionToSTL(bounding, 50, 'BNTL')
        file1 = open(cityname + '-topo.stl', 'w')
        file1.write(topo)
        file1.close()
