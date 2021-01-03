import pandas
import dask
import os
from .... datalayer import project
from unum.units import *
from ...utils import toUnum, toNumber
from ....measurements.GIS import topography
import numpy
from ..utils import getCellDataAndGroundData
from ...utils import coordinateHandler

class preProcess(project.ProjectMultiDBPublic):

    _publicProjectName = None
    _casePath = None
    _cloudName = None

    @property
    def casePath(self):
        return self._casePath

    @casePath.setter
    def casePath(self,newPath):
        self._casePath = newPath

    @property
    def cloudName(self):
        return self._cloudName

    @cloudName.setter
    def cloudName(self,newName):
        self._cloudName = newName

    def __init__(self, projectName, casePath, cloudName, databaseNameList=None, useAll=False,publicProjectName="OpenFOAMLSM"):
        """
        params:
        casePath = the path of the case
        cloudName = the name of the cloud, in which the particles' properties are saved.
        """
        self._publicProjectName = publicProjectName
        super().__init__(projectName=projectName,publicProjectName=publicProjectName,databaseNameList=databaseNameList,useAll=useAll)
        self._casePath = casePath
        self._cloudName = cloudName


    def extractFile(self,path,time,names,skiphead=20,skipend=4,vector=True):
        """
        Extracts data from a csv file.
        params:
        path: the path of the file
        time: the files' time step
        names: the names of the columns
        skiphead: number of lines to skip from the beginning of the file
        skipend: number of lines to skip from the ending of the file
        """
        newData = dask.dataframe.read_csv(path, skiprows=skiphead,
                                  skipfooter=skipend, engine='python', header=None, delim_whitespace=True, names=names)
        if vector:
            newData[names[0]] = newData[names[0]].str[1:]
            newData[names[2]] = newData[names[2]].str[:-1]
        newData["time"] = time
        return newData.astype(float)

    def extractRunResult(self, times = None, file=None, save=False, addToDB=True, withVelocity=False, withReleaseTimes=False,releaseTime=0, **kwargs):
        """
        Extracts results of an LSM run.
        parmas:
        times = a list of time steps to extract. If None, it extracts all time steps in the casePath.
        file = a name for a file, in which the data is saved. Default is "Positions.parquet", in the current working directory.
        save = a boolian parameter, to choose whether to save the data. default is False.
        addToDB = a boolian parameter, to choose whether to save the data. default is True. It is used only if save is True.
        withVelocity = set to True in order to extract the particles' velocities in addition to their positions. Default is False.
        **kwargs = any additional parameters to add to the description in the DB.
        """
        data = dask.dataframe.from_pandas(pandas.DataFrame(),npartitions=2)

        times = os.listdir(self.casePath) if times is None else times
        for filename in times:
            try:
                newData = self.extractFile(f"{self.casePath}/{filename}/lagrangian/{self.cloudName}/globalSigmaPositions",filename,['x', 'y', 'height'])
                if withVelocity:
                    dataU = self.extractFile(f"{self.casePath}/{filename}/lagrangian/{self.cloudName}/U",filename,['U_x', 'U_y', 'U_z'])
                    for col in ['U_x', 'U_y', 'U_z']:
                        newData[col] = dataU[col]
                if withReleaseTimes:
                    dataM = self.extractFile(f"{self.casePath}/{filename}/lagrangian/{self.cloudName}/age",filename,['age'],vector=False)
                    newData["releaseTime"] = dataM["time"]-dataM["age"]+releaseTime
                data=data.append(newData)
            except:
                pass
        if save:
            cur = os.getcwd()
            file = os.path.join(cur,f"{self.cloudName}Positions.parquet") if file is None else file
            data.to_parquet(file, compression="GZIP")
            if addToDB:
                self.addSimulationsDocument(resource=file,dataFormat="parquet",type="LSMPositions",
                                            desc=dict(casePath=self.casePath,cloudName=self.cloudName,**kwargs))
        return data

    def getConcentration(self, endTime, startTime=1, Q=1*kg, dx=10 * m, dy=10 * m, dz =10 * m, dt =10 * s,
                         Qunits=mg, lengthUnits=m, timeUnits=s, save=False, addToDB=True, file=None,releaseTime=0, nParticles=None, **kwargs):
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
            with open(os.path.join(self.casePath,"constant","kinematicCloudPositions"),"r") as readFile:
                Lines = readFile.readlines()
            try:
                nParticles=int(Lines[15])
            except:
                raise KeyError("Couldn't find number of particles; please deliver it as nParticles")
        dx = toNumber(toUnum(dx, lengthUnits), lengthUnits)
        dy = toNumber(toUnum(dy, lengthUnits), lengthUnits)
        dz = toNumber(toUnum(dz, lengthUnits), lengthUnits)
        dt = int(toNumber(toUnum(dt, timeUnits), timeUnits))
        withReleaseTimes = False
        if type(Q) == list:
            if len(Q) != endTime - startTime + 1:
                raise KeyError("Number of values in Q must be equal to the number of time steps!")
            try:
                Q = [toNumber(toUnum(q, Qunits), Qunits) for q in Q]
            except:
                Q = [toNumber(toUnum(q, Qunits), Qunits / timeUnits) for q in Q]
            releaseTimes = [releaseTime + i for i in range(int(endTime - startTime + 1))]
            dataQ = pandas.DataFrame({"releaseTime": releaseTimes, "Q": Q})
            withReleaseTimes = True
        else:
            try:
                Q = toNumber(toUnum(Q, Qunits), Qunits)
            except:
                Q = toNumber(toUnum(Q, Qunits), Qunits / timeUnits)
        documents = self.getSimulationsDocuments(type="openFoamLSMrun",
                                                casePath=self.casePath, cloudName=self.cloudName,
                                                      startTime=startTime, endTime=endTime, Q=Q, dx=dx,
                                                      dy=dy, dz=dz, dt=dt,nParticles=nParticles, **kwargs)
        if len(documents)==0:
            datalist = []
            for time in [startTime + dt * i for i in range(int((endTime-startTime) / dt))]:
                data = self.extractRunResult(times=[t for t in range(time, time + dt)], withReleaseTimes=withReleaseTimes,releaseTime=releaseTime)
                data["x"] = (data["x"] / dx).astype(int) * dx + dx / 2
                data["y"] = (data["y"] / dy).astype(int) * dy + dy / 2
                data["height"] = (data["height"] / dz).astype(int) * dz + dz / 2
                data["time"] = ((data["time"]-1) / dt).astype(int) * dt + dt
                if type(Q)==list:
                    data = data.set_index("releaseTime").join(dataQ.set_index("releaseTime")).reset_index()
                else:
                    data["Q"] = Q
                data["Dosage"] = data["Q"] / (nParticles * dx * dy * dz)
                data = data.drop(columns="Q")
                data = data.groupby(["x","y","height","time"]).sum()
                data["Concentration"] = data["Dosage"] / dt
                datalist.append(data.compute())
                print(f"Finished times {time} to {time+dt}")
            data = pandas.concat(datalist).reset_index()

            if save:
                cur = os.getcwd()
                file = os.path.join(cur,f"{self.cloudName}Concentration.parquet") if file is None else file
                data.to_parquet(file, compression="GZIP")
                if addToDB:
                    self.addSimulationsDocument(resource=file, dataFormat="parquet", type="openFoamLSMrun",
                                                desc=dict(casePath=self.casePath, cloudName=self.cloudName,
                                                          startTime=startTime, endTime=endTime, Q=Q, dx=dx,
                                                          dy=dy, dz=dz, dt=dt,nParticles=nParticles, **kwargs))
        else:
            data = documents[0].getData(usePandas=True)
        return data

    def makeSource(self, x, y, z, nParticles,type="Point",**kwargs):
        """
        Writes an instantaneous point source for the run.
        params:
        x = The x coordinate of the source.
        y = The y coordinate of the source.
        z = The z coordinate of the source.
        nParticles = The number of particles.
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
        string = getattr(self,f"makeSource_{type}")(string,x,y,z,nParticles,**kwargs)
        with open(os.path.join(self.casePath,"constant","kinematicCloudPositions"),"w") as writeFile:
            writeFile.write(string)

    def makeSource_Point(self,string,x,y,z,nParticles,**kwargs):
        for i in range(nParticles):
            string += f"({x} {y} {z})\n"
        string += ")\n"
        return string

    def makeCellHeights(self,times, ground="ground", fileName="cellHeights", resolution=10,savePandas=False, addToDB=False):
        """
        makes a file with the height of each cell.
        params:
        times = A list of time directories in which to save the new file
        ground = The name of the ground patch, default is "ground"
        fileName = The new file's name
        resolution = The cell length used in the conversion of the ground dataframe to a regular grid
        savePandas = Boolian, whether to save the dataframe
        addToDB = Boolian, whether to add the dataframe to the DB
        """
        documents = topography.getCacheDocuments(type="cellData", resolution=resolution,casePath=self.casePath)
        if len(documents)==0:
            cellData, groundData = getCellDataAndGroundData(casePath=self.casePath,ground=ground)
            cellData = topography.analysis.addHeight(data=cellData,groundData=groundData,resolution=resolution,
                                                     file=os.path.join(self.casePath, f"{fileName}.parquet"),casePath=self.casePath,
                                                     savePandas=savePandas,addToDB=addToDB)
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

    def makeUstar(self, times, fileName="ustar", ground="ground",heightLimits=[1,2], dField=None, dColumn="D",savePandas=False, addToDB=False,resolution=10):
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

        documents = topography.getCacheDocuments(type="cellData", resolution=resolution,casePath=self.casePath)
        if len(documents)==0:
            cellData, groundData = getCellDataAndGroundData(casePath=self.casePath,ground=ground)
            cellData = topography.analysis.addHeight(data=cellData,groundData=groundData,resolution=resolution,
                                                     file=os.path.join(self.casePath, f"{fileName}.parquet"),casePath=self.casePath,
                                                     savePandas=savePandas,addToDB=addToDB)
        else:
            cellData = documents[0].getData(usePandas=True)

        for time in times:
            documents = self.getCacheDocuments(type="ustar",casePath=self.casePath,time=time)

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
                xarrayU = coordinateHandler.regularizeTimeSteps(data=data.loc[data.height < heightLimits[1]].loc[data.height > heightLimits[0]]
                                                                .drop_duplicates(["x", "y"]), n=(nx,ny),
                                                                fieldList=["U"], coord2="y", addSurface=False, toPandas=False)[0]
                nsteps = int(len(data) / 10000)
                interpList = []
                if dField is not None:
                    xarrayD = coordinateHandler.regularizeTimeSteps(data=dField, n=(nx,ny),coord1Lim=(data["x"].min(), data["x"].max()),
                                                                   coord2Lim=(data["y"].min(), data["y"].max()),
                                                                   fieldList=[dColumn], coord2="y", addSurface=False, toPandas=False)[0]
                    for i in range(nsteps):
                        partition = data.loc[i * 10000:(i + 1) * 10000]
                        newInterpU = xarrayU.interp(x=partition["x"],y=partition["y"]).to_dataframe().reset_index().drop_duplicates(["x", "y"])
                        nans = newInterpU.loc[newInterpU.z.isnull()]
                        interpNans = []
                        for l, line in enumerate(nans.iterrows()):
                            interpNans.append(float(xarrayU.sel(x=line[1]["x"], y=line[1]["y"], method="nearest")["U"]))
                        nans["U"] = interpNans
                        newInterpU = pandas.concat([newInterpU.dropna(), nans])
                        newInterpD = xarrayD.interp(x=partition['x'], y=partition['y']).to_dataframe()
                        newInterp = newInterpD.reset_index().drop_duplicates(["x", "y"]).set_index(["x", "y"]).join(
                                    newInterpU.rename(columns={"U": "UnearGround"}).reset_index().drop_duplicates(
                                    ["x", "y"]).set_index(["x", "y"]), on=["x", "y"])
                        cellData = partition.set_index(["x", "y"]).join(newInterp, on=["x", "y"])
                        interpList.append(cellData)
                        print("finished interpolating for another 10000 cells")
                    partition = data.loc[(i + 1) * 10000:]
                    # newInterpU = xarrayU.interp(x=partition['x'], y=partition['y']).to_dataframe()
                    newInterpU = xarrayU.interp(x=partition["x"],y=partition["y"]).to_dataframe().reset_index().drop_duplicates(["x", "y"])
                    nans = newInterpU.loc[newInterpU.z.isnull()]
                    interpNans = []
                    for l, line in enumerate(nans.iterrows()):
                        interpNans.append(float(xarrayU.sel(x=line[1]["x"], y=line[1]["y"], method="nearest")["U"]))
                    nans["U"] = interpNans
                    newInterpU = pandas.concat([newInterpU.dropna(), nans])
                    newInterpD = xarrayD.interp(x=partition['x'], y=partition['y']).to_dataframe()
                    newInterp = newInterpD.reset_index().drop_duplicates(["x", "y"]).set_index(["x", "y"]).join(
                        newInterpU.rename(columns={"U": "UnearGround"}).reset_index().drop_duplicates(
                            ["x", "y"]).set_index(["x", "y"]), on=["x", "y"])
                    cellData = partition.set_index(["x", "y"]).join(newInterp, on=["x", "y"])
                    interpList.append(cellData)
                else:
                    for i in range(nsteps):
                        partition = data.loc[i * 10000:(i + 1) * 10000]
                        newInterp = xarrayU.interp(x=partition['x'], y=partition['y']).to_dataframe()
                        cellData = partition.set_index(["x", "y"]).join(
                            newInterp.rename(columns={"U": "UnearGround"}).reset_index().drop_duplicates(
                                ["x", "y"]).set_index(["x", "y"]), on=["x", "y"])
                        interpList.append(cellData)
                        print("finished interpolating for another 10000 cells")
                    partition = data.loc[(i + 1) * 10000:]
                    newInterp = xarrayU.interp(x=partition['x'], y=partition['y']).to_dataframe()
                    cellData = partition.set_index(['x', 'y']).join(
                        newInterp.rename(columns={"U": "UnearGround"}).reset_index().drop_duplicates(
                            ["x", "y"]).set_index(["x", "y"]), on=["x", "y"])
                    interpList.append(cellData)

                data = pandas.concat(interpList)
                data.loc[data["UnearGround"] < 1.5, "UnearGround"] = 1.5
                data["ustar"] = data["UnearGround"] * 0.4 / numpy.log((heightLimits[1]-heightLimits[0])/(2*0.15))
                data.loc[data["ustar"] < 0.2, "ustar"] = 0.2
                newFileString = ""
                data = data.reset_index()
                if savePandas:
                    data.to_parquet(os.path.join(self.casePath, f"{fileName}_{time}.parquet"), compression="gzip")
                    if addToDB:
                        self.addCacheDocument(resource=os.path.join(self.casePath, f"{fileName}.parquet"), dataFormat="parquet",
                                           type="ustar", desc={"casePath": self.casePath, "time": time})
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
            for time in times:
                with open(os.path.join(self.casePath, str(time), fileName), "w") as newFile:
                    newFile.write(newFileString)
            print("wrote time ", time)


if __name__ == "__main__":

    xs = [0.01*h for h in range(200)]
    vels = []
    for h in xs:
        if h <= 0.39:
            vels.append((3.84*h+1.5))
        elif h >0.39 and h<1.18:
            vels.append((1.27901266*h+2.50064835))
        else:
            vels.append(4)