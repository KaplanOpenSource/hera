import pandas
import dask
import os
from .... datalayer import project
from ... utils import toNumber,toUnum
from unum.units import *

class preProcess(project.ProjectMultiDBPublic):

    _publicProjectName = None

    def __init__(self, projectName,databaseNameList=None,useAll=False,publicProjectName="OpenFOAMLSM"):
        self._publicProjectName = publicProjectName
        super().__init__(projectName=projectName,publicProjectName=publicProjectName,useAll=useAll)

    def extractFile(self,path,time,names,skiphead=20,skipend=4):
        newData = dask.dataframe.read_csv(path,skiprows=skiphead,skipfooter=skipend,engine="python",header=None,delim_whitespace=True,names=names)
        newData[names[0]] = newData[names[0]].str[1:]
        newData[names[2]] = newData[names[2]].str[:-1]
        newData["time"] = time
        return newData.astype(float)

    def extractRunResult(self,casePath,cloudName,times = None, file=None,save=False,addToDB=True,**kwargs):

        data = self.extractFile(f"{casePath}/constant/kinematicCloudPositions",0,["x","y","z"],skiphead=18)
        for col in ["U_x","U_y","U_z"]:
            data[col]=0
        times = os.listdir(casePath) if times is None else times
        for filename in times:
            try:
                newData = self.extractFile(f"{casePath}/{filename}/lagrangian/{cloudName}/globalPositions",filename,["x","y","z"])
                dataU = self.extractFile(f"{casePath}/{filename}/lagrangian/{cloudName}/U",filename,["U_x","U_y","U_z"])
                for col in ["U_x","U_y","U_z"]:
                    newData[col] = dataU[col]
                data=data.append(newData)
            except:
                pass
        if save:
            cur = os.getcwd()
            file = os.path.join(cur,f"{cloudName}Positions.parquet") if file is None else file
            data.to_parquet(file,compression="GZIP")
            if addToDB:
                self.addSimulationsDocument(resource=file,dataFormat="parquet",type="LSMPositions",desc=dict(casePath=casePath,cloudName=cloudName,**kwargs))
        return data

    def getConcentration(self,data,Q=1*kg, measureX=10*m, measureY=10*m,measureZ=10*m,measureTime=5*s,Qunits=mg,lengthUnits=m,timeUnits=s,nParticles=None):

        data = data.copy()
        measureX = toNumber(toUnum(measureX,lengthUnits),lengthUnits)
        measureY = toNumber(toUnum(measureY, lengthUnits), lengthUnits)
        measureZ = toNumber(toUnum(measureZ, lengthUnits), lengthUnits)
        measureTime = toNumber(toUnum(measureTime, timeUnits), timeUnits)
        Q = toNumber(toUnum(Q, Qunits), Qunits)
        nParticles = len(data.loc[data.time==data.time.min()]) if nParticles is None else nParticles
        data["x"] = (data["x"]/measureX).astype(int)*measureX+measureX/2
        data["y"] = (data["y"] / measureY).astype(int) * measureY + measureY / 2
        data["z"] = (data["z"] / measureZ).astype(int) * measureZ + measureZ / 2
        data["time"] = ((data["time"]-1) / measureTime).astype(int) * measureTime + measureTime
        data["Dosage"] = Q/(nParticles*measureX*measureY*measureZ)
        data = data.drop(columns=["U_x","U_y","U_z"]).groupby(["x","y","z","time"]).sum()
        data["Concentration"] = data["Dosage"]/measureTime
        return data.reset_index()

    def makeOnePointSource(self,casePath,x,y,z,nParticles):

        str = "/*--------------------------------*- C++ -*----------------------------------*\\ \n| =========                 |                                                 |" \
              "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n" \
              "|  \\    /   O peration     | Version:  dev                                   |\n" \
              "|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n" \
              "|    \\/     M anipulation  |                                                 |\n" \
              "\*---------------------------------------------------------------------------*/\n" \
              "FoamFile\n" \
              "{\n" \
              "    version     2.0;\n" \
              "    format      ascii;\n" \
              "    class       vectorField;\n" \
              "    object      kinematicCloudPositions;\n" \
              "}\n" \
              "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" \
              "\n" \
              f"{nParticles}\n" \
              "(\n"
        for i in range(nParticles):
            str += f"({x} {y} {z})\n"
        str += ")\n" \
               "\n" \
               "// ************************************************************************* //"
        with open(os.path.join(casePath,"constant","kinematicCloudPositions"),"w") as file:
            file.write(str)

if __name__=="__main__":
    from hera.simulations.openfoam.lsm.preProcess import preProcess
    pre = preProcess(projectName="debuggingLSM")
    directory = "/home/eyal/Projects/2020/LSMstudy/openFOAMruns/stepByStep"
    data = pre.extractRunResult(directory,"kinematicCloud", file="/home/eyal/Projects/2020/LSMstudy/openFOAMruns/stepByStep/Gauss.parquet",
                                save=False, addToDB=False,nParticles=200000, test="noRandom")