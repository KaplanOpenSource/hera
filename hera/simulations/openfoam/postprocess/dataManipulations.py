import numpy
import pandas
from .... import datalayer
import os

class dataManipulations():

    _projectName = None

    def __init__(self, projectName):

        self._projectName = projectName

    def arrangeSlice(self, data, xdir=True, ydir=True, save=False, path=None, key="Clip", addToDB=True, **kwargs):
        """
        Arranging data of a slice: adding distance downwind, velocity and height over terrain.
        Params:
        data: The data of the slice (pandas DataFrame)
        xdir: If true, the x component of the velocity is positive.
        ydir: If true, the y component of the velocity is positive.
        Returns:
            The arranged data
        """
        path = path if path is not None else os.path.join("ArrangedClip.hdf")
        data["terrain"] = [numpy.nan for x in range(len(data))]
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
            # x = base_data.loc[i]["x"]
            # y = base_data.loc[i]["y"]
            # if data.loc[data["x"] == x].loc[data["y"] == y].empty:
            #     pass
            # else:
            index = list(data.loc[data.distance==dist].index)[0]
            data.at[index, "terrain"] = base

        data = data.sort_values(by="distance").set_index("distance").interpolate(method='index').reset_index().dropna()
        data["heightOverTerrain"] = data["z"] - data["terrain"]
        data = data.query("heightOverTerrain>=0")

        if save:
            data.to_hdf(path, key=key, format="table")
            if addToDB:
                datalayer.Measurements.addDocument(projectName=self._projectName, desc=(dict(filter="OrganizedSlice", **kwargs)),
                                                   resource=dict(path=path, key=key), type="OFsimulation", dataFormat="HDF")

        return data

    def findDetaiedLocations(self, data, n=10):

        optional = []
        for d in data.distance.drop_duplicates():
           if len(data.query("distance==@d")) > n:
               optional.append(d)
        optional.sort()
        return optional

    def arrangeClip(self, data, save=False, path=None, key="Clip", addToDB=True, **kwargs):

        path = path if path is not None else os.path.join("ArrangedClip.hdf")
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
        if save:
            data.to_hdf(path, key=key, format="table")
            if addToDB:
                datalayer.Measurements.addDocument(projectName=self._projectName, desc=(dict(filter="OrganizedClip", **kwargs)),
                                                   resource=dict(path=path, key=key), type="OFsimulation", dataFormat="HDF")

        return data

    def makeSliceHeightData(self, data, height, variable, limit=100):

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
        return returndata

    def makeClipHeightData(self, data, height, variable, limit=100, skips=5, save=False, path=None, key="Clip", addToDB=True, **kwargs):

        path = path if path is not None else os.path.join("ArrangedClip.hdf")
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
        if save:
            returndata.to_hdf(path, key=key, format="table")
            if addToDB:
                datalayer.Measurements.addDocument(projectName=self._projectName, desc=(dict(filter="OrganizedClip", height=height, **kwargs)),
                                                   resource=dict(path=path, key=key), type="OFsimulation", dataFormat="HDF")
        return returndata
