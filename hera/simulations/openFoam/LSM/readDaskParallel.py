import pandas
import numpy
import os
import glob
from dask.delayed import delayed
from dask import dataframe
from itertools import product



def extractFile(path, columnNames, vector=True):
    """
        Extracts data from a csv file.

    Parameters
    ----------
    path: str
        The path of the file
    time: str
        The files' time step
    columnNames: list of str
        The names of the columns
    skiphead: int
        Number of lines to skip from the beginning of the file
    skipend: int
        Number of lines to skip from the ending of the file

    Returns
    -------
        Pandas with the data.
    """
    skiphead = 20
    skipend = 4

    cnvrt = lambda x: float(x.replace("(", "").replace(")", ""))
    cnvrtDict = dict([(x, cnvrt) for x in columnNames])

    try:
        newData = pandas.read_csv(path,
                                  skiprows=skiphead,
                                  skipfooter=skipend,
                                  engine='python',
                                  header=None,
                                  delim_whitespace=True,
                                  converters=cnvrtDict,
                                  names=columnNames)
    except ValueError:
        newData = []

    if len(newData) == 0:
        file = open(path, "r")
        lines = file.readlines()
        file.close()
        vals = lines[17]
        data = []

        if vector:
            if "{" in vals:
                inputs = vals.split("{")
                repeat = int(inputs[0])
                valuesList = inputs[1][inputs[1].find("(") + 1:inputs[1].find(")")]
                data = dict(
                    [(colname, [float(x)] * repeat) for colname, x in zip(columnNames, valuesList.split(" "))])
            else:
                for rcrdListTuple in vals.split("(")[2:]:
                    record = dict(
                        [(name, float(y)) for name, y in zip(columnNames, rcrdListTuple.split(")")[0].split(" "))])
                    data.append(record)
        else:

            if "{" in vals:
                inputs = vals.split("{")
                repeat = int(inputs[0])
                value = float(inputs[1].split("}")[0])
                data = [{columnNames[0]: value} for x in range(repeat)]

            else:
                valuesList = vals.split("(")[1]
                for rcrdListItem in valuesList.split(" "):
                    record = {columnNames[0]: float(rcrdListItem.split(")")[0])}
                    data.append(record)

        newData = pandas.DataFrame(data)

    return newData.astype(float)


def readRecord(timeName, casePath, withVelocity=False, withReleaseTimes=False, withMass=False,
               cloudName="kinematicCloud"):

    print(f"Processing {timeName}")

    columnsDict = dict(x=[], y=[], z=[], id=[], procId=[], globalID=[], globalX=[], globalY=[], globalZ=[])
    if withMass:
        columnsDict['mass'] = []
    if withReleaseTimes:
        columnsDict['age'] = []
    if withVelocity:
        columnsDict['U_x'] = []
        columnsDict['U_y'] = []
        columnsDict['U_z'] = []

    newData = pandas.DataFrame(columnsDict, dtype=numpy.float64)

    try:
        newData = extractFile(
            os.path.join(casePath, timeName, "lagrangian", cloudName, "globalSigmaPositions"),
            ['x', 'y', 'z'])
        for fld in ['x', 'y', 'z']:
            newData[fld] = newData[fld].astype(numpy.float64)

        dataID = extractFile(os.path.join(casePath, timeName, "lagrangian", cloudName, "origId"),
                             ['id'], vector=False)
        newData['id'] = dataID['id'].astype(numpy.float64)

        dataprocID = extractFile(
            os.path.join(casePath, timeName, "lagrangian", cloudName, "origProcId"), ['procId'],
            vector=False)
        newData['procId'] = dataprocID['procId'].astype(numpy.float64)

        newData = newData.ffill().assign(globalID=1000000000 * newData.procId + newData.id)

        dataGlobal = extractFile(
            os.path.join(casePath, timeName, "lagrangian", cloudName, "globalPositions"),
            ['globalX', 'globalY', 'globalZ'])

        for col in ['globalX', 'globalY', 'globalZ']:
            newData[col] = dataGlobal[col].astype(numpy.float64)

        if withVelocity:
            dataU = extractFile(os.path.join(casePath, timeName, "lagrangian", cloudName, "U"),
                                ['U_x', 'U_y', 'U_z'])
            for col in ['U_x', 'U_y', 'U_z']:
                newData[col] = dataU[col]

        if withReleaseTimes:
            dataM = extractFile(os.path.join(casePath, timeName, "lagrangian", cloudName, "age"),
                                ['age'], vector=False)
            # newData["releaseTime"] = dataM["time"] - dataM["age"] + releaseTime

        if withMass:
            dataM = extractFile(os.path.join(casePath, timeName, "lagrangian", cloudName, "mass"),
                                ['mass'], vector=False)
            try:
                newData["mass"] = dataM["mass"]
            except:
                newData = newData.compute()
                newData["mass"] = dataM["mass"]

    except:
        pass

    theTime = os.path.split(timeName)[-1]
    newData['time'] = float(theTime)
    return newData


def loadDataParallel(casePath,
             times=None,
             parallelCase=True,
             withVelocity=False,
             withReleaseTimes=False,
             withMass=False,
             cloudName="kinematicCloud"):
    """
        Extracts results of an LSM run.

    Parameters
    -----------

        times: list
             a list of time steps to extract.
             If None, it extracts all time steps in the casePath.
        file: str
            The name for a file, in which the data is saved. Default is "<cloud name>_data.parquet", in the current working directory.

        withVelocity: bool
                True: extract the particles' velocities in addition to their positions.
                Default is False.

        **kwargs:
            Any additional parameters to add to the description in the DB.
    Returns
    --------
        document of the data.
    """
    finalCasePath = os.path.abspath(casePath)
    loader = lambda timeName: readRecord(timeName,
                                         casePath=finalCasePath,
                                         withVelocity=withVelocity,
                                         withReleaseTimes=withReleaseTimes,
                                         withMass=withMass,
                                         cloudName=cloudName)

    if parallelCase:

        processorList = [os.path.basename(proc) for proc in glob.glob(os.path.join(finalCasePath, "processor*"))]
        if len(processorList) == 0:
            raise ValueError(f"There are no processor* directories in the case {finalCasePath}. Is it parallel?")

        timeList = sorted([x for x in os.listdir(os.path.join(finalCasePath, processorList[0])) if (
                os.path.isdir(os.path.join(finalCasePath, processorList[0], x)) and
                x.isdigit() and
                (not x.startswith("processor") and x not in ["constant", "system", "rootCase", 'VTK']))],
                          key=lambda x: int(x))

        data = dataframe.from_delayed(
            [delayed(loader)(os.path.join(processorName, timeName)) for processorName, timeName in
             product(processorList, timeList[1:])])
    else:

        timeList = sorted([x for x in os.listdir(finalCasePath) if (
                os.path.isdir(x) and
                x.isdigit() and
                (not x.startswith("processor") and x not in ["constant", "system", "rootCase", 'VTK']))],
                          key=lambda x: int(x))

        data = dataframe.from_delayed([delayed(loader)(timeName) for timeName in timeList])

    return data
