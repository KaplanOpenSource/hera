import pandas
import numpy
import os

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

def groundToPandas(casePath,times, ground="ground",fillna=True,fillVal=0):
    f = open(os.path.join(casePath,"0","cellCenters"), "r")
    lines = f.readlines()
    f.close()
    fboundary = open(os.path.join(casePath,"0","Hmix"), "r")
    boundarylines = fboundary.readlines()
    fboundary.close()

    for i in range(len(lines)):
        if "internalField" in lines[i]:
            cellStart = i+3
            nCells = int(lines[i+1])
        if ground in lines[i]:
            break
    for boundaryLine in range(len(boundarylines)):
        if "boundaryField" in boundarylines[boundaryLine]:
            break
    nGroundValues = int(lines[i+4])
    groundData = pandas.read_csv(os.path.join(casePath,"0","cellCenters"), skiprows=i+6,
                               skipfooter=len(lines)-(i+6+nGroundValues),
                               engine='python',
                              header=None,
                                delim_whitespace=True, names=['x', 'y', 'z'])

    groundData['x'] = groundData['x'].str[1:]
    groundData['z'] = groundData['z'].str[:-1]
    groundData = groundData.astype(float)
    xarrayGround = groundData.set_index(["x", "y"]).to_xarray()
    cellData = pandas.read_csv(os.path.join(casePath,"0","cellCenters"), skiprows=cellStart,
                                skipfooter=len(lines)-(cellStart+nCells),
                              engine='python',
                              header=None,
                              delim_whitespace=True, names=['x', 'y', 'z'])
    cellData['x'] = cellData['x'].str[1:]
    cellData['z'] = cellData['z'].str[:-1]
    cellData = cellData.astype(float)
    interpetedGroundValues = xarrayGround.interp(x=cellData['x'],y=cellData['y']).to_dataframe()
    cellData = cellData.set_index(["x", "y"]).join(interpetedGroundValues.rename(columns={"z":"ground"}).reset_index().drop_duplicates(["x","y"]).set_index(["x","y"]),on=["x","y"])
    if fillna:
        cellData=cellData.fillna(fillVal)
    cellData["height"]=cellData["z"]-cellData["ground"]
    newFileString = ""
    cellData = cellData.reset_index()
    for i in range(cellStart):
         newFileString += lines[i]
    for i in range(nCells):
         newFileString += f"({cellData['x'][i]} {cellData['y'][i]} {cellData['height'][i]})\n"
    newFileString += ")\n;\n\n"
    for i in range(boundaryLine,len(boundarylines)):
         newFileString += boundarylines[i]
    for time in times:
         with open(os.path.join(casePath,str(time),"cellHeights"),"w") as newFile:
            newFile.write(newFileString)


