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

def makeCellHeights(casePath,times, ground="ground",fileName="cellHeights",fillna=True,fillVal=0, savePandas=False,addToDB=False):
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
    xarrayGround = handler.regularizeTimeSteps(data=groundData,fieldList=["z"],coord2="y",addSurface=False,toPandas=False)[0]
    cellData = pandas.read_csv(os.path.join(casePath,"0","cellCenters"), skiprows=cellStart,
                                skipfooter=len(lines)-(cellStart+nCells),
                              engine='python',
                              header=None,
                              delim_whitespace=True, names=['x', 'y', 'z'])
    cellData['x'] = cellData['x'].str[1:]
    cellData['z'] = cellData['z'].str[:-1]
    cellData = cellData.astype(float)
    nsteps = int(nCells/1000)
    interpList = []
    concatedList = []

    for i in range(1,nsteps):
        partition = cellData.loc[i*1000:(i+1)*1000]
        newInterp = xarrayGround.interp(x=partition['x'],y=partition['y']).to_dataframe()
        interpList.append(newInterp.drop_duplicates())
        if i>=100 and i%100==0:
            concatedList.append(pandas.concat(interpList))
            interpList = []
            print(f"Interpolated ground heights for another step")

    partition = cellData.loc[(i+1) * 1000:]
    newInterp = xarrayGround.interp(x=partition['x'], y=partition['y']).to_dataframe()
    interpList.append(newInterp)
    concatedList.append(pandas.concat(interpList))
    print("finished interpolations")
    interpolatedGroundValues = pandas.concat(concatedList)
    cellData = cellData.set_index(["x", "y"]).join(interpolatedGroundValues.rename(columns={"z":"ground"}).reset_index().drop_duplicates(["x","y"]).set_index(["x","y"]),on=["x","y"])
    if fillna:
        cellData=cellData.fillna(fillVal)
    cellData["height"]=cellData["z"]-cellData["ground"]
    cellData.loc[cellData.height < 0, "height"] = 0
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
         with open(os.path.join(casePath,str(time),fileName),"w") as newFile:
            newFile.write(newFileString)
    if savePandas:
        cellData.to_parquet(os.path.join(casePath,f"{fileName}.parquet"),compression="gzip")
        if addToDB:
            p.addCacheDocument(resource=os.path.join(casePath,f"{fileName}.parquet"),dataFormat="parquet",type="cellData",desc={"casePath":casePath})

def makeUstar(casePath,times,fileName="ustar", savePandas=False, addToDB=False):
    for time in times:
        f = open(os.path.join(casePath,str(time),"U"), "r")
        lines = f.readlines()
        f.close()
        fboundary = open(os.path.join(casePath,str(time),"Hmix"), "r")
        boundarylines = fboundary.readlines()
        fboundary.close()

        for i in range(len(lines)):
            if "internalField" in lines[i]:
                cellStart = i+3
                nCells = int(lines[i+1])
        for boundaryLine in range(len(boundarylines)):
            if "boundaryField" in boundarylines[boundaryLine]:
                break
        cellData = p.getCacheDocuments(type="cellData",casePath=casePath)[0].getData(usePandas=True)
        Ufield = pandas.read_csv(os.path.join(casePath,str(time),"U"), skiprows=cellStart,
                                    skipfooter=len(lines)-(cellStart+nCells),
                                  engine='python',
                                  header=None,
                                  delim_whitespace=True, names=['u', 'v', 'w'])
        Ufield['u'] = Ufield['u'].str[1:]
        Ufield['w'] = Ufield['w'].str[:-1]
        Ufield = Ufield.astype(float)
        Ufield["U"] = numpy.sqrt(Ufield['u']**2+Ufield['v']**2+Ufield['w']**2)
        data = cellData.join(Ufield)
        xarrayU = handler.regularizeTimeSteps(data=data.loc[data.U<5].loc[data.U>4].drop_duplicates(["x","y"]),fieldList=["U"],coord2="y",addSurface=False,toPandas=False)[0]
        nsteps = int(nCells/1000)
        interpList = []
        concatedList = []

        for i in range(1,nsteps):
            partition = data.loc[i*1000:(i+1)*1000]
            newInterp = xarrayU.interp(x=partition['x'],y=partition['y']).to_dataframe()
            interpList.append(newInterp.drop_duplicates())
            if i>=100 and i%100==0:
                concatedList.append(pandas.concat(interpList))
                interpList = []
                print(f"Interpolated velocity near ground for another step")

        partition = data.loc[(i+1) * 1000:]
        newInterp = xarrayU.interp(x=partition['x'], y=partition['y']).to_dataframe()
        interpList.append(newInterp)
        concatedList.append(pandas.concat(interpList))
        print("finished interpolations")
        interpolatedUValues = pandas.concat(concatedList)
        data = data.set_index(["x", "y"]).join(interpolatedUValues.rename(columns={"U":"UnearGround"}).reset_index().drop_duplicates(["x","y"]).set_index(["x","y"]),on=["x","y"])
        data.loc[data["UnearGround"]<1.5,"UnearGround"]=1.5
        data["ustar"]=data["UnearGround"]*0.4/numpy.log(data["height"]/0.15)
        data.loc[data["ustar"] < 0.2, "ustar"] = 0.2
        data = data.fillna(data.ustar.mean())
        newFileString = ""
        data = data.reset_index()
        for i in range(cellStart):
             newFileString += lines[i]
        newFileString = newFileString.replace("vector","scalar").replace("Vector","Scalar")
        for i in range(nCells):
             newFileString += f"{data['ustar'][i]}\n"
        newFileString += ")\n;\n\n"
        for i in range(boundaryLine,len(boundarylines)):
             newFileString += boundarylines[i]
        for time in times:
             with open(os.path.join(casePath,str(time),fileName),"w") as newFile:
                newFile.write(newFileString)
        print("wrote time ",time)
        if savePandas:
            data.to_parquet(os.path.join(casePath,f"{fileName}_{time}.parquet"),compression="gzip")
            if addToDB:
                p.addCacheDocument(resource=os.path.join(casePath,f"{fileName}.parquet"),dataFormat="parquet",type="ustar",desc={"casePath":casePath,"time":time})