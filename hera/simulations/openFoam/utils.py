# import pandas
import numpy
import os
from . import HERAMETADATA
from ...utils import loadJSON
# from ..utils.coordinateHandler import coordinateHandler
# handler = coordinateHandler()



def getCellDataAndGroundData(casePath,ground="ground"):
    f = open(os.path.join(casePath, "0", "cellCenters"), "r")
    lines = f.readlines()
    f.close()
    fboundary = open(os.path.join(casePath, "0", "Hmix"), "r")
    boundarylines = fboundary.readlines()
    fboundary.close()

    for i in range(len(lines)):
        if "internalField" in lines[i]:
            cellStart = i + 3
            nCells = int(lines[i + 1])
        if ground in lines[i]:
            break
    for boundaryLine in range(len(boundarylines)):
        if "boundaryField" in boundarylines[boundaryLine]:
            break

    nGroundValues = int(lines[i + 4])
    groundData = pandas.read_csv(os.path.join(casePath, "0", "cellCenters"), skiprows=i + 6,
                                 skipfooter=len(lines) - (i + 6 + nGroundValues),
                                 engine='python',
                                 header=None,
                                 delim_whitespace=True, names=['x', 'y', 'z'])

    groundData['x'] = groundData['x'].str[1:]
    groundData['z'] = groundData['z'].str[:-1]
    groundData = groundData.astype(float)

    cellData = pandas.read_csv(os.path.join(casePath, "0", "cellCenters"), skiprows=cellStart,
                               skipfooter=len(lines) - (cellStart + nCells),
                               engine='python',
                               header=None,
                               delim_whitespace=True, names=['x', 'y', 'z'])
    cellData['x'] = cellData['x'].str[1:]
    cellData['z'] = cellData['z'].str[:-1]
    cellData = cellData.astype(float)
    return cellData,groundData