import pandas
import os
class preProcess():

    def extractFile(self,casePath,filename,cloudName,file,names,skiphead=20,skipend=4):
        newData = pandas.read_csv(f"{casePath}/{filename}/lagrangian/{cloudName}/{file}", skiprows=skiphead,
                                  skipfooter=skipend, engine='python', header=None, delim_whitespace=True, names=names)
        newData[names[0]] = newData[names[0]].str[1:]
        newData[names[2]] = newData[names[2]].str[:-1]
        newData["time"] = filename
        return newData.astype(float)

    def extractRunResult(self, casePath, cloudName):
        datalist = []
        for filename in os.listdir(casePath):
            try:
                print(float(filename))
                newData = self.extractFile(casePath,filename,cloudName,"globalPositions",['x', 'y', 'z'])
                dataU = self.extractFile(casePath,filename,cloudName,"U",['U_x', 'U_y', 'U_z'])
                for col in ['U_x', 'U_y', 'U_z']:
                    newData[col] = dataU[col]
                datalist.append(newData)
            except:
                pass
        return pandas.concat(datalist)