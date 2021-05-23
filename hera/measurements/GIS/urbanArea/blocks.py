import numpy
from itertools import product


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

