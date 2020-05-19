import numpy
import pandas

class dataManipulations():

    def arrangeSlice(self, data, xdir=True, ydir=True):
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
        return data

    def findDetaiedLocations(self, data):

        optional = []
        for d in data.distance.drop_duplicates():
           if len(data.query("distance==@d")) > 10:
               optional.append(d)
        optional.sort()
        return optional