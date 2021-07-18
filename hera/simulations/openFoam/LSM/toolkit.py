import glob
import pandas
import os
import xarray
import numpy
from dask.delayed import delayed
from dask import dataframe

from unum.units import *
from ....utils import tounit, tonumber

from ....measurements.GIS.locations.topography import TopographyToolkit

from ..utils import getCellDataAndGroundData
from ...utils import coordinateHandler
from ....datalayer import nonDBMetadataFrame
from .... import toolkit

from .sourcesFactoryTool import sourcesFactoryTool
from itertools import product


class OFLSMToolkit(toolkit.abstractToolkit):
    _casePath = None
    _cloudName = None
    _sources = None
    _topography = None

    _parallelCase = None

    @property
    def analysis(self):
        return self._analysis

    @property
    def sourcesFactory(self):
        return self._sourcesFactory

    @property
    def casePath(self):
        return self._casePath

    @casePath.setter
    def casePath(self, newPath):
        self._casePath = newPath

    @property
    def topography(self):
        return self._topography

    @property
    def cloudName(self):
        return self._cloudName

    @cloudName.setter
    def cloudName(self, value):
        self._cloudName = str(value)

    @property
    def parallelCase(self):
        return self._parallelCase

    @parallelCase.setter
    def parallelCase(self, value):
        self._parallelCase = value

    def __init__(self, projectName, casePath=None, cloudName="kinematicCloud", FilesDirectory=None, parallelCase=False):
        """
        Parameters
        ----------
        casePath: str
            The path of the case

        cloudName: str,
            The name of the cloud, in which the particles' properties are saved.
        """
        super().__init__(projectName=projectName, toolkitName="OF_LSM", FilesDirectory=FilesDirectory)
        self._casePath = os.getcwd() if casePath is None else os.path.abspath(casePath)
        self._sourcesFactory = sourcesFactoryTool()
        self._cloudName = cloudName
        self._parallelCase = parallelCase
        self._topography = TopographyToolkit(projectName)
        self._analysis = Analysis(self)

    def _extractFile(self, path, columnNames, vector=True):
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
            self.logger.execute(f"{path} is not a cvs, going to a specialized parser")
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

    def _readRecord(self, timeName,casePath, withVelocity=False, withReleaseTimes=False, withMass=False):
        self.logger.debug(f"Starting the read record with timeName {timeName}")



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
            newData = self._extractFile(
                os.path.join(casePath, timeName, "lagrangian", self._cloudName, "globalSigmaPositions"),
                ['x', 'y', 'z'])
            for fld in ['x', 'y', 'z']:
                newData[fld] = newData[fld].astype(numpy.float64)

            dataID = self._extractFile(os.path.join(casePath, timeName, "lagrangian", self._cloudName, "origId"),
                                       ['id'], vector=False)
            newData['id'] = dataID['id'].astype(numpy.float64)

            dataprocID = self._extractFile(
                os.path.join(casePath, timeName, "lagrangian", self._cloudName, "origProcId"), ['procId'],
                vector=False)
            newData['procId'] = dataprocID['procId'].astype(numpy.float64)

            newData = newData.ffill().assign(globalID=1000000000 * newData.procId + newData.id)

            dataGlobal = self._extractFile(
                os.path.join(casePath, timeName, "lagrangian", self._cloudName, "globalPositions"),
                ['globalX', 'globalY', 'globalZ'])

            for col in ['globalX', 'globalY', 'globalZ']:
                newData[col] = dataGlobal[col].astype(numpy.float64)

            if withVelocity:
                dataU = self._extractFile(os.path.join(casePath, timeName, "lagrangian", self._cloudName, "U"),
                                          ['U_x', 'U_y', 'U_z'])
                for col in ['U_x', 'U_y', 'U_z']:
                    newData[col] = dataU[col]

            if withReleaseTimes:
                dataM = self._extractFile(os.path.join(self._casePath, timeName, "lagrangian", self._cloudName, "age"),
                                          ['age'], vector=False)
                # newData["releaseTime"] = dataM["time"] - dataM["age"] + releaseTime

            if withMass:
                dataM = self._extractFile(os.path.join(casePath, timeName, "lagrangian", self._cloudName, "mass"),
                                          ['mass'], vector=False)
                try:
                    newData["mass"] = dataM["mass"]
                except:
                    newData = newData.compute()
                    newData["mass"] = dataM["mass"]

        except:
            self.logger.debug(f"No data at time {timeName}")

        theTime = os.path.split(timeName)[-1]
        newData['time'] = float(theTime)
        self.logger.debug(f"The output is {newData}")
        return newData

    @property
    def doctype(self):
        return "LSMRuns"

    def to_paraview_CSV(self, data, outputdirectory, filename, timeFactor=1):
        """
            Writes the globalPositions (globalX,globalY,globalZ) as  CSV for visualization in paraview.
            In paraview, each timestep is a different file.

        Parameters
        -----------
        data: dask.dataframe or pandas.dataframe
            The data to present

        outputdirectory: str
            The directory to write the files in

        timeFactor : int
            Multiply the time by a factro to make the time step round (so that paraview will recognize it).

        filename: str
            The filename to write.

        Returns
        -------
            None
        """
        for times, timedata in data.groupby("time"):
            with open(os.path.join(outputdirectory, f"{filename}_{str(int(timeFactor * times)).replace('.', '_')}.csv"),
                      "w") as outputfile:
                outputfile.writelines(timedata[['globalX', 'globalY', 'globalZ']].to_csv(index=False))

    def loadData(self,
                 times=None,
                 saveMode=toolkit.TOOLKIT_SAVEMODE_NOSAVE,
                 withVelocity=False,
                 withReleaseTimes=False,
                 withMass=False,
                 releaseTime=0,
                 cloudName="kinematicCloud",
                 casePath=None,
                 **kwargs):
        """
            Extracts results of an LSM run.

        Parameters
        -----------

            times: list
                 a list of time steps to extract.
                 If None, it extracts all time steps in the casePath.
            file: str
                The name for a file, in which the data is saved. Default is "<cloud name>_data.parquet", in the current working directory.

            saveMode: str

                    - toolkit.TOOLKIT_SAVEMODE_NOSAVE   : Just load the data from file and return the datafile

                    - toolkit.TOOLKIT_SAVEMODE_ONLYFILE : Loads the data from file and save to a file.
                                                  raise exception if file exists.

                    - toolkit.TOOLKIT_SAVEMODE_ONLYFILE_REPLACE: Loads the data from file and save to a file.
                                                  Replace the file if it exists.

                    - toolkit.TOOLKIT_SAVEMODE_FILEANDDB : Loads the data from file and save to a file and store to the DB as a source.
                                                    Raise exception if the entry exists.

                    - toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE: Loads the data from file and save to a file and store to the DB as a source.
                                                    Replace the entry in the DB if it exists.


            withVelocity: bool
                    True: extract the particles' velocities in addition to their positions.
                    Default is False.

            **kwargs:
                Any additional parameters to add to the description in the DB.
        Returns
        --------
            document of the data.
        """

        if saveMode not in [toolkit.TOOLKIT_SAVEMODE_NOSAVE,
                            toolkit.TOOLKIT_SAVEMODE_FILEANDDB,
                            toolkit.TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,
                            toolkit.TOOLKIT_SAVEMODE_FILEANDDB,
                            toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE
                            ]:
            optins = ",".join([toolkit.TOOLKIT_SAVEMODE_NOSAVE,
                               toolkit.TOOLKIT_SAVEMODE_FILEANDDB,
                               toolkit.TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,
                               toolkit.TOOLKIT_SAVEMODE_FILEANDDB,
                               toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE
                               ])

            raise ValueError(f"saveMode must be one of [{optins}]. Got {saveMode}")

        finalCasePath = self.casePath if casePath is None else os.path.abspath(casePath)

        finalFileName = os.path.abspath(os.path.join(self.FilesDirectory, f"{self.cloudName}Data.parquet"))

        if saveMode in [toolkit.TOOLKIT_SAVEMODE_ONLYFILE, toolkit.TOOLKIT_SAVEMODE_FILEANDDB] and os.path.exists(
                finalFileName):
            raise FileExistsError(f"saveMode was set to no replace and file {finalFileName} already exists")

        loader = lambda timeName: self._readRecord(timeName,
                                                   casePath=finalCasePath,
                                                   withVelocity=withVelocity,
                                                   withReleaseTimes=withReleaseTimes,
                                                   withMass=withMass)

        if self.parallelCase:

            processorList = [os.path.basename(proc) for proc in glob.glob(os.path.join(finalCasePath, "processor*"))]
            if len(processorList) == 0:
                raise ValueError(f"There are no processor* directories in the case {finalCasePath}. Is it parallel?")

            timeList = sorted([x for x in os.listdir(os.path.join(finalCasePath, processorList[0])) if (
                    os.path.isdir(os.path.join(finalCasePath, processorList[0],x)) and
                    x.isdigit() and
                    (not x.startswith("processor") and x not in ["constant", "system", "rootCase", 'VTK']))],
                              key=lambda x: int(x))

            data = dataframe.from_delayed(
                [delayed(loader)(os.path.join(processorName, timeName)) for processorName, timeName in
                 product(processorList, timeList)])
        else:

            timeList = sorted([x for x in os.listdir(finalCasePath) if (
                    os.path.isdir(x) and
                    x.isdigit() and
                    (not x.startswith("processor") and x not in ["constant", "system", "rootCase", 'VTK']))],
                              key=lambda x: int(x))

            data = dataframe.from_delayed([delayed(loader)(timeName) for timeName in timeList])

        if saveMode in [toolkit.TOOLKIT_SAVEMODE_ONLYFILE,
                        toolkit.TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,
                        toolkit.TOOLKIT_SAVEMODE_FILEANDDB,
                        toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE]:

            data.to_parquet(finalFileName, compression="GZIP")

            if saveMode in [toolkit.TOOLKIT_SAVEMODE_FILEANDDB, toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE]:
                doc = self.getSimulationsDocuments(type=self.doctype, casePath=self.casePath, cloudName=self.cloudName)

                if doc is not None and saveMode == toolkit.TOOLKIT_SAVEMODE_FILEANDDB:
                    raise FileExistsError(f"Data already in the DB. save mode is set to no replace")
                elif doc is None:
                    self.addSimulationsDocument(resource=finalFileName,
                                                dataFormat=toolkit.datatypes.PARQUET,
                                                type=self.doctype,
                                                desc=dict(casePath=self.casePath, cloudName=self.cloudName, **kwargs))
                else:
                    doc.resource = finalFileName
                    doc.save()
        else:
            doc = nonDBMetadataFrame(data=data, type=self.doctype, casePath=self.casePath, cloudName=self.cloudName)

        return doc

    def makeSource(self, x, y, z, nParticles, type="Point", fileName="kinematicCloudPositions", **kwargs):
        """
            Writes a source position files.
            Saves the file to the constant directory of the case.
        Parameters:
        -----------
            x: float
                The x coordinate of the source.
            y: float
                The y coordinate of the source.
            z: float
             The z coordinate of the source.
            nParticles: int
                The number of particles.
            type: str
                The type of the source. Must be one of the sources.sourcesTypeList
            fileName: str
                The output filename (without the directory).
            kwargs:
                additional parameters for the source creation.

        Returns
        -------
            None


        """
        string = "/*--------------------------------*- C++ -*----------------------------------*\n " \
                 "| =========                 |                                                 |\n" \
                 "| \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n" \
                 "|  \    /   O peration     | Version:  dev                                   |\n" \
                 "|   \  /    A nd           | Web:      www.OpenFOAM.org                      |\n" \
                 "|    \/     M anipulation  |                                                 |\n" \
                 "\*---------------------------------------------------------------------------*/\n" \
                 "FoamFile\n{    version     2.0;\n    format      ascii;\n    class       vectorField;\n" \
                 "    object      kinematicCloudPositions;\n}\n" \
                 "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n" \
            f"{nParticles}\n(\n"
        source = self.sourcesFactory.makeSource(x=x, y=y, z=z, nParticles=nParticles, type=type, **kwargs)
        for i in range(nParticles):
            string += f"({source.loc[i].x} {source.loc[i].y} {source.loc[i].z})\n"
        string += ")\n"
        with open(os.path.join(self.casePath, "constant", fileName), "w") as writeFile:
            writeFile.write(string)

    def makeCellHeights(self, times, ground="ground", fileName="cellHeights", resolution=10,
                        saveMode=toolkit.TOOLKIT_SAVEMODE_ONLYFILE, fillna=0):
        """
        makes a file with the height of each cell.
        params:
        times = A list of time directories in which to save the new file
        ground = The name of the ground patch, default is "ground"
        fileName = The new file's name
        resolution = The cell length used in the conversion of the ground dataframe to a regular grid
        savePandas = Boolean, whether to save the dataframe
        addToDB = Boolean, whether to add the dataframe to the DB
        """
        documents = self.topography.getCacheDocuments(type="cellData", resolution=resolution, casePath=self.casePath)
        if len(documents) == 0:
            cellData, groundData = getCellDataAndGroundData(casePath=self.casePath, ground=ground)
            cellData = self.topography.analysis.addHeight(data=cellData, groundData=groundData, resolution=resolution,
                                                          file=os.path.join(self.casePath, f"{fileName}.parquet"),
                                                          casePath=self.casePath,
                                                          saveMode=saveMode, fillna=fillna)
        else:
            cellData = documents[0].getData(usePandas=True)
        f = open(os.path.join(self.casePath, "0", "cellCenters"), "r")
        lines = f.readlines()
        f.close()
        fboundary = open(os.path.join(self.casePath, "0", "Hmix"), "r")
        boundarylines = fboundary.readlines()
        fboundary.close()
        for i in range(len(lines)):
            if "internalField" in lines[i]:
                cellStart = i + 3
                nCells = int(lines[i + 1])
        for boundaryLine in range(len(boundarylines)):
            if "boundaryField" in boundarylines[boundaryLine]:
                break

        newFileString = ""
        for i in range(cellStart):
            newFileString += lines[i]
        for i in range(nCells):
            newFileString += f"({cellData['x'][i]} {cellData['y'][i]} {cellData['height'][i]})\n"
        newFileString += ")\n;\n\n"
        for i in range(boundaryLine, len(boundarylines)):
            newFileString += boundarylines[i]
        for time in times:
            with open(os.path.join(self.casePath, str(time), fileName), "w") as newFile:
                newFile.write(newFileString)

    def makeUstar(self, times, fileName="ustar", ground="ground", heightLimits=[1, 2], dField=None, dColumn="D",
                  saveMode=toolkit.TOOLKIT_SAVEMODE_ONLYFILE, resolution=10):
        """
        makes a file with the shear velocity in each cell.
        params:
        times = A list of time directories in which to save the new file
        fileName = The new file's name
        ground = The name of the ground patch, default is "ground"
        resolution = The cell length used in the conversion of the ground dataframe to a regular grid
        savePandas = Boolian, whether to save the dataframe
        addToDB = Boolian, whether to add the dataframe to the DB
        """

        documents = self.topography.getCacheDocuments(type="cellData", resolution=resolution, casePath=self.casePath)

        if len(documents) == 0:
            cellData, groundData = getCellDataAndGroundData(casePath=self.casePath, ground=ground)
            cellData = self.topography.analysis.addHeight(data=cellData, groundData=groundData, resolution=resolution,
                                                          file=os.path.join(self.casePath, f"{fileName}.parquet"),
                                                          casePath=self.casePath,
                                                          saveMode=saveMode)
        else:
            cellData = documents[0].getData(usePandas=True)

        for time in times:
            documents = self.getCacheDocuments(type="ustar", casePath=self.casePath, time=time)

            f = open(os.path.join(self.casePath, str(time), "U"), "r")
            lines = f.readlines()
            f.close()
            fboundary = open(os.path.join(self.casePath, str(time), "Hmix"), "r")
            boundarylines = fboundary.readlines()
            fboundary.close()

            for i in range(len(lines)):
                if "internalField" in lines[i]:
                    cellStart = i + 3
                    nCells = int(lines[i + 1])
            for boundaryLine in range(len(boundarylines)):
                if "boundaryField" in boundarylines[boundaryLine]:
                    break

            if len(documents) == 0:
                Ufield = pandas.read_csv(os.path.join(self.casePath, str(time), "U"), skiprows=cellStart,
                                         skipfooter=len(lines) - (cellStart + nCells),
                                         engine='python',
                                         header=None,
                                         delim_whitespace=True, names=['u', 'v', 'w'])

                Ufield['u'] = Ufield['u'].str[1:]
                Ufield['w'] = Ufield['w'].str[:-1]
                Ufield = Ufield.astype(float)
                Ufield["U"] = numpy.sqrt(Ufield['u'] ** 2 + Ufield['v'] ** 2 + Ufield['w'] ** 2)
                data = cellData.join(Ufield)
                nx = int((data["x"].max() - data["x"].min()) / resolution)
                ny = int((data["y"].max() - data["y"].min()) / resolution)
                xarrayU = coordinateHandler.regularizeTimeSteps(
                    data=data.loc[data.height < heightLimits[1]].loc[data.height > heightLimits[0]]
                        .drop_duplicates(["x", "y"]), n=(nx, ny),
                    fieldList=["U"], coord2="y", addSurface=False, toPandas=False)[0]
                if dField is None:
                    xarrayD = xarray.zeros_like(xarrayU).rename({"U": dColumn})
                else:
                    xarrayD = coordinateHandler.regularizeTimeSteps(data=dField, n=(nx, ny),
                                                                    coord1Lim=(data["x"].min(), data["x"].max()),
                                                                    coord2Lim=(data["y"].min(), data["y"].max()),
                                                                    fieldList=[dColumn], coord2="y", addSurface=False,
                                                                    toPandas=False)[0]
                fillU = xarrayU.U.mean()
                fillD = xarrayD[dColumn].mean()
                unearground = []
                ds = []
                valuesDict = {}
                for i in range(len(data)):
                    x = data.loc[i].x
                    y = data.loc[i].y
                    if x not in valuesDict.keys():
                        valuesDict[x] = {}
                    if y in valuesDict[x].keys():
                        unearground.append(valuesDict[x][y]["u"])
                        ds.append(valuesDict[x][y]["d"])
                    else:
                        valuesDict[x][y] = {}
                        uVal = float(xarrayU.interp(x=x, y=y).fillna(fillU).U)
                        dVal = float(xarrayD.interp(x=x, y=y).fillna(fillD)[dColumn])
                        valuesDict[x][y]["u"] = uVal
                        valuesDict[x][y]["d"] = dVal
                    if i > 9999 and i % 10000 == 0:
                        print(i)
                data["UnearGround"] = unearground
                data[dColumn] = ds
                data.loc[data["UnearGround"] < 1.5, "UnearGround"] = 1.5
                data["ustar"] = data["UnearGround"] * 0.4 / numpy.log(
                    (0.5 * (heightLimits[1] + heightLimits[0]) - data[dColumn]) / (0.15))
                data.loc[data["ustar"] < 0.2, "ustar"] = 0.2
                newFileString = ""
                data = data.reset_index()
                # if savePandas:
                #     data.to_parquet(os.path.join(self.casePath, f"{fileName}_{time}.parquet"), compression="gzip")
                #     if addToDB:
                #         self.addCacheDocument(resource=os.path.join(self.casePath, f"{fileName}.parquet"), dataFormat="parquet",
                #                            type="ustar", desc={"casePath": self.casePath, "time": time})
            else:
                data = documents[0].getData(usePandas=True)
            for i in range(cellStart):
                newFileString += lines[i]
            newFileString = newFileString.replace("vector", "scalar").replace("Vector", "Scalar")
            for i in range(nCells):
                newFileString += f"{data['ustar'][i]}\n"
            newFileString += ")\n;\n\n"
            for i in range(boundaryLine, len(boundarylines)):
                newFileString += boundarylines[i]

            with open(os.path.join(self.casePath, str(time), fileName), "w") as newFile:
                newFile.write(newFileString)
            print("wrote time ", time)

    def createRootCaseMeshLink(self, rootCase):
        """
            Creates the directories for run (currently only parallel).

            For each processorXX in the rootCase:

                    1. Copy the timestep

            If parallel, create all the processor** and link it.

        :param rootCase:
        :param parallel:
        :return:
        """
        for fl in glob.glob(os.path.join(rootCase, "processor*")):
            print(fl)
            fullpath = os.path.join(os.path.abspath(fl), lastTS)

            proc = os.path.split(fl)[-1]
            destination = os.path.join(os.path.abspath(proc), "3600")
            os.makedirs(os.path.dirname(destination), exist_ok=True)
            os.system(f"cp {fullpath} {destination} -rT")

            fullpath = os.path.abspath(os.path.join(fl, "constant", "polyMesh"))
            destination = os.path.join(os.path.abspath(proc), "constant", "polyMesh")
            os.makedirs(os.path.dirname(destination), exist_ok=True)
            os.system(f"ln -s {fullpath} {destination}")

            # link the root dir .
            curdir = os.path.abspath(os.path.join("rootCase", os.path.basename(fl)))
            targetdir = os.path.abspath(os.path.join(fl, "rootCase"))
            os.system(f"ln -s {curdir} {targetdir} ")


class Analysis:
    _datalayer = None

    def __init__(self, datalayer):
        self._datalayer = datalayer

    def getConcentration(self, endTime, startTime=1, Q=1 * kg, dx=10 * m, dy=10 * m, dz=10 * m, dt=10 * s, loc=None,
                         sigmaCoordinates=True,
                         Qunits=mg, lengthUnits=m, timeUnits=s, OFmass=False, save=False, addToDB=True, file=None,
                         releaseTime=0, nParticles=None, withID=False, **kwargs):
        """
        Calculates the concentration of dispersed particles.
        parmas:
        nParticles = The total number of particles induced to the system at all times.
        endTime = The final time step to read.
        startTime = The first time step to read. The default is 1.
        Q = The total mass of particles induced to the system. The default is 1 kg.
        dx = The distance between grid points of the concentration field in the x direction. The default is 10 meters.
        dy = The distance between grid points of the concentration field in the y direction. The default is 10 meters.
        dz = The distance between grid points of the concentration field in the z direction. The default is 10 meters.
        dt = The time measure between grid points of the concentration field. The default is 10 seconds.
        Qunits = The mass units of the concentration in the final dataframe. The default is mg.
        lengthUnits = The length units of the concentration in the final dataframe. The default is meter.
        timeUnits = The time units of the time steps in the final dataframe. The default is second.
        file = a name for a file, in which the data is saved. Default is "Concentration.parquet", in the current working directory.
        save = a boolian parameter, to choose whether to save the data. default is False.
        addToDB = a boolian parameter, to choose whether to save the data. default is True. It is used only if save is True.
        **kwargs = any additional parameters to add to the description in the DB.
        """
        if nParticles is None:
            with open(os.path.join(self._datalayer.casePath, "constant", "kinematicCloudPositions"), "r") as readFile:
                Lines = readFile.readlines()
            try:
                nParticles = int(Lines[15])
            except:
                raise KeyError("Couldn't find number of particles; please deliver it as nParticles")
        dx = tonumber(tounit(dx, lengthUnits), lengthUnits)
        dy = tonumber(tounit(dy, lengthUnits), lengthUnits)
        dz = tonumber(tounit(dz, lengthUnits), lengthUnits)
        dt = int(tonumber(tounit(dt, timeUnits), timeUnits))
        withReleaseTimes = False
        if type(Q) == list:
            if len(Q) != endTime - startTime + 1:
                raise KeyError("Number of values in Q must be equal to the number of time steps!")
            try:
                Q = [tonumber(tounit(q, Qunits), Qunits) for q in Q]
            except:
                Q = [tonumber(tounit(q, Qunits), Qunits / timeUnits) for q in Q]
            releaseTimes = [releaseTime + i for i in range(int(endTime - startTime + 1))]
            dataQ = pandas.DataFrame({"releaseTime": releaseTimes, "Q": Q})
            withReleaseTimes = True
        else:
            try:
                Q = tonumber(tounit(Q, Qunits), Qunits)
            except:
                Q = tonumber(tounit(Q, Qunits), Qunits / timeUnits)
        documents = self._datalayer.getSimulationsDocuments(type="openFoamLSMrun",
                                                            casePath=self._datalayer.casePath,
                                                            cloudName=self._datalayer.cloudName,
                                                            startTime=startTime, endTime=endTime, Q=Q, dx=dx,
                                                            dy=dy, dz=dz, dt=dt, nParticles=nParticles, **kwargs)

        if len(documents) == 0:
            datalist = []
            for time in [startTime + dt * i for i in range(int((endTime - startTime) / dt))]:
                data = self._datalayer._readRecord(time, withMass=OFmass, withID=withID,
                                                   sigmaCoordinates=sigmaCoordinates)
                for t in range(time + 1, time + dt):
                    data = data.append(self._datalayer._readRecord(t, withMass=OFmass, withID=withID,
                                                                   sigmaCoordinates=sigmaCoordinates))
                data["x"] = (data["x"] / dx).astype(int) * dx + dx / 2
                data["y"] = (data["y"] / dy).astype(int) * dy + dy / 2
                data["z"] = (data["z"] / dz).astype(int) * dz + dz / 2
                data["time"] = ((data["time"] - 1) / dt).astype(int) * dt + dt
                if OFmass:
                    data = data.groupby(["x", "y", "z", "time"]).sum().reset_index()

                    data["mass"] = data["mass"] * tonumber(tounit(1 * kg, Qunits), Qunits)
                    if time == startTime:
                        Qinsim = data.mass.sum() / dt
                        Qratio = Q / Qinsim
                    data["mass"] = data["mass"] * Qratio
                    data["Dosage"] = data["mass"] / (dx * dy * dz)
                else:
                    if type(Q) == list:
                        data = data.set_index("releaseTime").join(dataQ.set_index("releaseTime")).reset_index()
                    else:
                        data["Q"] = Q
                    data["Dosage"] = data["Q"] / (nParticles * dx * dy * dz)
                    data = data.drop(columns="Q")
                    data = data.groupby(["x", "y", "z", "time"]).sum().reset_index()
                if loc is not None:
                    data = data.loc[data[loc[0]] == loc[1]]
                data["Concentration"] = data["Dosage"] / dt
                datalist.append(data)
                print(f"Finished times {time} to {time + dt}")
            data = pandas.concat(datalist).reset_index()
            if save:
                cur = os.getcwd()
                file = os.path.join(cur, f"{self._datalayer.cloudName}Concentration.parquet") if file is None else file
                data.to_parquet(file, compression="GZIP")
                if addToDB:
                    self._datalayer.addSimulationsDocument(resource=file, dataFormat="parquet", type="openFoamLSMrun",
                                                           desc=dict(casePath=self._datalayer.casePath,
                                                                     cloudName=self._datalayer.cloudName,
                                                                     startTime=startTime, endTime=endTime, Q=Q, dx=dx,
                                                                     dy=dy, dz=dz, dt=dt, nParticles=nParticles,
                                                                     **kwargs))
        else:
            data = documents[0].getData(usePandas=True)
        return data
