import pandas
import numpy
import os
from ..utils.coordinateHandler import coordinateHandler
handler = coordinateHandler()
from ... datalayer import project
p = project.Project("openFoamdata")
def centersToPandas(skipend, filepath='C', skiphead = 22, saveToTxt=False, fileName="CellCenters.txt"):
    """
        Extract pandas from openfoam cell centers file.

        It is also possible to save it to  a new txt file

    Parameters
    -----------

    filepath: str
        The path to the cell centers file.

        Default: 'C'

    saveToTxt:  boolean
        whether to save the centers coords to txt file
        Default: False

    fileName: str
        The file name to save.
        Default: 'CellCenters.txt'

    Returns
    --------

    cellData: pandas DF

    """
    cellData = pandas.read_csv(filepath, skiprows=skiphead,
                               skipfooter=skipend,
                               engine='python',
                               header=None,
                               delim_whitespace=True, names=['x', 'y', 'z'])

    cellData['x'] = cellData['x'].str[1:]
    cellData['z'] = cellData['z'].str[:-1]
    cellData = cellData.astype(float)

    if saveToTxt:

        L = numpy.array(cellData[['x', 'y', 'z']])
        numberOfCells = L.shape[0]
        out=[]
        for i in range(numberOfCells):
            out.append(f"({L[i, 0]} {L[i, 1]} {L[i, 2]})\n")

        with open(fileName, "w") as outfile:
            outfile.write("".join(out))

    return cellData

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