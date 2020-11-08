import pandas
import dask
import os
from .... datalayer import project
from unum.units import *
from ...utils import toUnum, toNumber

class preProcess(project.ProjectMultiDBPublic):

    _publicProjectName = None


    def __init__(self, projectName, databaseNameList=None, useAll=False,publicProjectName="OpenFOAMLSM"):

        self._publicProjectName = publicProjectName
        super().__init__(projectName=projectName,publicProjectName=publicProjectName,databaseNameList=databaseNameList,useAll=useAll)

    def extractFile(self,path,time,names,skiphead=20,skipend=4):
        newData = dask.dataframe.read_csv(path, skiprows=skiphead,
                                  skipfooter=skipend, engine='python', header=None, delim_whitespace=True, names=names)
        newData[names[0]] = newData[names[0]].str[1:]
        newData[names[2]] = newData[names[2]].str[:-1]
        newData["time"] = time
        return newData.astype(float)

    def extractRunResult(self, casePath, cloudName,file=None,save=False,addToDB=True, **kwargs):

        data = self.extractFile(f"{casePath}/constant/kinematicCloudPositions",0,['x','y','z'])
        for col in ['U_x', 'U_y', 'U_z']:
            data[col] = 0
        for filename in os.listdir(casePath):
            try:
                newData = self.extractFile(f"{casePath}/{filename}/lagrangian/{cloudName}/globalPositions",filename,['x', 'y', 'z'])
                dataU = self.extractFile(f"{casePath}/{filename}/lagrangian/{cloudName}/U",filename,['U_x', 'U_y', 'U_z'])
                for col in ['U_x', 'U_y', 'U_z']:
                    newData[col] = dataU[col]
                data=data.append(newData)
            except:
                pass
        #data = pandas.concat(datalist)
        if save:
            cur = os.getcwd()
            file = os.path.join(cur,f"{cloudName}Positions.parquet") if file is None else file
            data.to_parquet(file, compression="GZIP")
            if addToDB:
                self.addSimulationsDocument(resource=file,dataFormat="parquet",type="LSMPositions",desc=dict(casePath=casePath,cloudName=cloudName,**kwargs))
        return data

    def getConcentration(self, data,Q=1*kg, measureX=10*m, measureY=10*m, measureZ = 10*m, measureTime = 10*s,Qunits=mg,lengthUnits=m,timeUnits=s,nParticles=None):

        data = data.copy()
        measureX = toNumber(toUnum(measureX, lengthUnits),lengthUnits)
        measureY = toNumber(toUnum(measureY, lengthUnits),lengthUnits)
        measureZ = toNumber(toUnum(measureZ, lengthUnits),lengthUnits)
        measureTime = toNumber(toUnum(measureTime, timeUnits),timeUnits)
        Q = toNumber(toUnum(Q, Qunits),Qunits)
        nParticles = len(data.loc[data.time == data.time.min()]) if nParticles is None else nParticles
        data["x"] = (data["x"]/measureX).astype(int)*measureX+measureX/2
        data["y"] = (data["y"]/ measureY).astype(int)*measureY + measureY / 2
        data["z"] = (data["z"]/ measureZ).astype(int)*measureZ + measureZ / 2
        data["time"] = ((data["time"]-1) / measureTime).astype(int)*measureTime + measureTime
        data["Dosage"] = Q/(nParticles*measureX*measureY*measureZ)
        data = data.drop(columns=["U_x","U_y","U_z"]).groupby(["x","y","z","time"]).sum()
        data["Concentration"] = data["Dosage"]/measureTime

        return data.reset_index()


# def getA(ustar,u,v,w,sigmau,sigmav,sigmaw,dsu,dsv,dsw,Tl):
#     B = [numpy.sqrt(2*sigmau**2/Tl),numpy.sqrt(2*sigmav**2/Tl),numpy.sqrt(2*sigmaw**2/Tl)]
#     return [0.5*(-B[0]**2+dsu*w)*(sigmaw**2*u+ustar**2*w)/(sigmau**2*sigmaw**2-ustar**4),0.5*(-B[1]**2+dsv*w)*v/sigmav**2,
#             0.5*dsw+0.5*(-B[2]**2+dsw*w)*(sigmau**2*w+ustar**2*u)/(sigmau**2*sigmaw**2-ustar**4)]
#
# def getTurbVel(ustar,u,v,w,sigmau,sigmav,sigmaw,dsu,dsv,dsw, dt,rand,Tl):
#     B = [numpy.sqrt(2 * sigmau ** 2 / Tl), numpy.sqrt(2 * sigmav ** 2 / Tl), numpy.sqrt(2 * sigmaw ** 2 / Tl)]
#     A = [0.5*(-B[0]**2+dsu*w)*(sigmaw**2*u+ustar**2*w)/(sigmau**2*sigmaw**2-ustar**4),0.5*(-B[1]**2+dsv*w)*v/sigmav**2,
#             0.5*dsw+0.5*(-B[2]**2+dsw*w)*(sigmau**2*w+ustar**2*u)/(sigmau**2*sigmaw**2-ustar**4)]
#     vel = [u,v,w]
#     return [vel[i]+A[i]*dt+B[i]*numpy.sqrt(dt)*rand[i] for i in range(3)]
#
# def getTurbVelFortran(ustar,u,v,w,sigmau,sigmaw, dt,rand, tlm1h,tlm1v,dsw):
#     B = [sigmau*numpy.sqrt(1-numpy.exp(-2*dt*tlm1h)),sigmau*numpy.sqrt(1-numpy.exp(-2*dt*tlm1h)),sigmaw*numpy.sqrt(1-numpy.exp(-2*dt*tlm1v))]
#     return [u-sigmau**2*tlm1h*dt*(sigmaw**2*u+ustar**2*w)/(sigmau**2*sigmaw**2-ustar**4)+B[0]*rand[0],
#             v*numpy.exp(-dt*tlm1h)+B[1]*rand[1],
#             w+0.5*dsw*dt+(0.5*dsw*w-sigmaw**2*tlm1v)*dt*(sigmau**2*w+ustar**2*u)/(sigmau**2*sigmaw**2-ustar**4)+B[2]*rand[2]]

# from hera.simulations.LSM import LagrangianReader
# import os
# import xarray
# saveDir = os.getcwd()
# results_full_path = os.path.join(saveDir, "tozaot", "machsan", "OUTD2d03_3_")
# L = []
# i = 0
# for xray in LagrangianReader.toNetcdf(basefiles=results_full_path):
#     L.append(xray)
# data = xarray.concat(L, dim="datetime")
# data.to_netcdf(os.path.join("/home/ofir/Projects/2020/testLSM/checkGauss/netcdf", "data%s.nc" % i))
if __name__=="__main__":
    from hera.simulations.openfoam.lsm.preProcess import preProcess
    pre = preProcess(projectName="LSMOpenFOAM")
    directory = "/home/ofir/Projects/2020/LSM/horizontalInPlace"
    data = pre.extractRunResult(directory,"kinematicCloud",save=True,addToDB=True,nparticles=100000)