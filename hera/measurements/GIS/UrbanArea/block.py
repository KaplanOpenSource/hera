from shapely.geometry import Polygon
from .building import Building as bld


class Block():

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


    def LambdaF(self, MeteoAngle ,TotalArea=None):
        """
        # Calculate average lambda F of a block
        """
        bldHeightToReduce = 0
        if self._block.empty:
            return 0
        Map_A_f = 0
        i = 0
        j = 0
        indexes = self._block['BLDG_HT'].index
        self._numberOfBld = len(indexes)

        for i in indexes:
            if (self._block['FTYPE'][i] == 16) or (self._block['FTYPE'][i] == 14):
                bldHeightToReduce = self.updateNAbuilding(bldHeightToReduce, i)
            else:
                    bldHeight = self._block['BLDG_HT'][i]
                    if bldHeight < 2:
                        bldHeight = self._block['HI_PNT_Z'][i]-self._block['HT_LAND'][i]
                        if bldHeight < 2:
                            self.updateErrorBuilding(i)
                            j = j + 1
                            continue
                        else:
                            self._block.at[i,'BLDG_HT'] = bldHeight

                    if self._block.geometry[i].geom_type == 'MultiPolygon':
                        for poly in self._block.geometry[i]:
                            building = bld(Polygon(poly),bldHeight)
                            Map_A_f = Map_A_f + building.A_f2(MeteoAngle)
                    else:
                        building = bld(self._block.geometry[i], bldHeight)
                        Map_A_f = Map_A_f + building.A_f2(MeteoAngle)

            j = j+1

        if TotalArea == None:
            bounds = self._block.total_bounds
            TotalArea = ((bounds[2] - bounds[0])) * ((bounds[3] - bounds[1]))

        self.lambda_f = Map_A_f / TotalArea
        self._totalBldHeight = self._block['BLDG_HT'].sum() - bldHeightToReduce
        return self.lambda_f

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
