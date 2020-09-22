import numpy
import pandas
from .... import datalayer
import os

class dataManipulations():

    _projectName = None

    def __init__(self, projectName):

        self._projectName = projectName

    def arrangeSlice(self, data, twoD=False, xdir=True, ydir=True, save=False, path=None, key="Slice", addToDB=True, hdfName="ArrangedSlice", **kwargs):
        """
        Arranging data of a slice: adding distance downwind, velocity and height over terrain.
        Params:
        data: The data of the slice (pandas DataFrame)
        xdir: If true, the x component of the velocity is positive.
        ydir: If true, the y component of the velocity is positive.
        Returns:
            The arranged data
        """
        data["terrain"] = [numpy.nan for x in range(len(data))]
        if twoD:
            data = data.drop(columns="U_y")
            data["Velocity"] = numpy.sqrt(data["U_x"] * data["U_x"] + data["U_z"] * data["U_z"])
        else:
            data["Velocity"] = numpy.sqrt(data["U_x"] * data["U_x"] + data["U_y"] * data["U_y"] + data["U_z"] * data["U_z"])

        if xdir:
            x2 = data["x"] - data['x'].min()
        else:
            x2 = data["x"].max() - data['x']

        if ydir:
            y2 = data["y"] - data['y'].min()
        else:
            y2 = data["y"].max() - data['y']
        data["distance"] = numpy.sqrt(x2 * x2 + y2 * y2)

        base_data = data.query("Velocity==0").reset_index()
        for dist in base_data.distance.drop_duplicates():
            base = base_data.loc[base_data.distance==dist]["z"].max()
            index = list(data.loc[data.distance==dist].index)[0]
            data.at[index, "terrain"] = base

        data = data.sort_values(by="distance").set_index("distance").interpolate(method='index').reset_index().dropna()
        data["heightOverTerrain"] = data["z"] - data["terrain"]
        data = data.query("heightOverTerrain>=0")

        self.saveAndAddtoDB(save=save, addToDB=addToDB, data=data, path=path, key=key, projectName=self._projectName,
                            filter=hdfName, **kwargs)
        return data

    def findDetaiedLocations(self, data, n=10):

        optional = []
        for d in data.distance.drop_duplicates():
           if len(data.query("distance==@d")) > n:
               optional.append(d)
        optional.sort()
        return optional

    def arrangeClip(self, data, save=False, path=None, key="Clip", addToDB=True, hdfName="ArrangedClip", **kwargs):

        data["Velocity"] = numpy.sqrt(data["U_x"] * data["U_x"] + data["U_y"] * data["U_y"] + data["U_z"] * data["U_z"])
        data["terrain"] = [numpy.nan for x in range(len(data))]
        base_data = data.query("Velocity==0").reset_index()
        for i in range(len(base_data)):
            base = base_data.loc[i]["z"].max()
            x = base_data.loc[i]["x"]
            y = base_data.loc[i]["y"]

            indice = list(data.loc[data.x==x].loc[data.y==y].index)
            for index in indice:
                data.at[index, "terrain"] = base

        data = data.dropna()
        data["heightOverTerrain"] = data["z"] - data["terrain"]

        self.saveAndAddtoDB(save=save, addToDB=addToDB, data=data, path=path, key=key, projectName=self._projectName,
                            filter=hdfName, **kwargs)

        return data

    def makeSliceHeightData(self, data, height, variable, limit=100, save=False, path=None, key="Slice", addToDB=True, hdfName="SliceHeightData", **kwargs):

        stop = False
        delta = 1
        selecteddist = []
        finer = []
        newdata = data
        while stop==False:
            found=False
            while found==False:
                d = newdata.query("heightOverTerrain>@height-@delta and heightOverTerrain<@height+@delta").sort_values(by="distance")
                if d.empty:
                    delta += 1
                else:
                    found=True
            dif = [list(d.distance)[i + 1] - list(d.distance)[i] for i in range(len(d) - 1)]
            for i in range(len(dif)):
                if dif[i]<limit:
                    selecteddist.append(list(d.distance)[i])
                else:
                    finer.append(i)
            selecteddist.append(list(d.distance)[-1])
            if len(finer)==0:
                stop=True
            else:
                newdata=pandas.DataFrame()
                for i in finer:
                    mindist = list(d.distance)[i]
                    maxdist = list(d.distance)[i+1]
                    newd = data.query("distance>@mindist and distance<@maxdist")
                    newdata = newdata.append(newd)
                delta += 1
                finer = []
        values = []
        for dist in selecteddist:
            d = data.loc[data.distance == dist]
            d2 = pandas.DataFrame([dict(distance=dist, heightOverTerrain=height)])
            value = d.append(d2).sort_values(by="heightOverTerrain").set_index("heightOverTerrain").interpolate(
                method='index').loc[height][variable]
            values.append(value)

        returndata = pandas.DataFrame({"distance": selecteddist, variable: values}).dropna().sort_values(by="distance")
        self.saveAndAddtoDB(save=save, addToDB=addToDB, data=returndata, path=path, key=key, projectName=self._projectName,
                            filter=hdfName, height=height, **kwargs)
        return returndata

    def makeClipHeightData(self, data, height, variable, limit=100, skips=5, save=False, path=None, key="Clip", addToDB=True, hdfName="ClipHeightData", **kwargs):

        returndata = pandas.DataFrame()
        dx = list(data.x.drop_duplicates().sort_values())
        for i in range(0, len(dx),skips):
            print(i)
            x = dx[i]
            d = data.loc[data.x==x]
            if d.heightOverTerrain.max()<height or d.heightOverTerrain.min()>height:
                pass
            else:
                d["distance"] = d["y"]
                newdata = self.makeSliceHeightData(data=d, height=height, variable=variable, limit=limit)
                newdata = newdata.rename(columns={"distance":"y"})
                newdata["x"] = x
                returndata = returndata.append(newdata)
        self.saveAndAddtoDB(save=save, addToDB=addToDB, data=returndata, path=path, key=key, projectName=self._projectName,
                            filter=hdfName, height=height, **kwargs)

        return returndata

    def saveAndAddtoDB(self, save, addToDB, data, path, key, projectName, filter, **kwargs):

        path = path if path is not None else os.path.join("%s.hdf" % filter)
        if save:
            data.to_hdf(path, key=key, format="table")
            if addToDB:
                datalayer.Measurements.addDocument(projectName=projectName, desc=(dict(filter=filter, **kwargs)),
                                                   resource=dict(path=path, key=key), type="OFsimulation", dataFormat="HDF")

    def toSigmaCoordinate(self, data, sigmas=None):
        zmax = data["z"].max()
        data["sigma"] = data["z"]/(zmax-data["terrain"])
        return data

    def meanVelocity(self, data, sigmas, method="interpolate"):

        datalist = []
        if method=="interpolate":
            for dist in list(data.distance.drop_duplicates()):
                for sigma in sigmas:
                    if sigma > data.loc[data.distance == dist].sigma.min() and sigma < data.loc[
                        data.distance == dist].sigma.max():
                        data = data.append({"distance": dist, "sigma": sigma}, ignore_index=True)
                datalist.append(data.loc[data.distance == dist].sort_values(by=["distance", "sigma"]).set_index("sigma").interpolate(method="index").reset_index())
        else:
            for sigma in sigmas:
                data["difference"] = numpy.square(data["sigma"] - sigma)
                newData = data.loc[data.difference == data.difference.min()]
                newData["sigma"] = sigma
                datalist.append(newData)
        data = pandas.concat(datalist)
        datalist = []
        for sigma in sigmas:
            newdata = data.loc[data.sigma==sigma]
            newdata["disturbance"] = 100*(newdata.Velocity - newdata.Velocity.mean())/newdata.Velocity
            datalist.append(newdata)
        return pandas.concat(datalist)