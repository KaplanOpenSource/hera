import numpy

def addBase(data, base_data):
    """
    Adding the base height for points on a slice.

    Parameters
    ----------
    data: the data of the inletMesh patch of the slice
    base_data: the data of the terrain patch of the slice

    Returns the data with the base
    -------

    """

    data["base"] = [numpy.nan for x in range(len(data))]

    for i in range(len(base_data)):
        base = base_data.loc[i]["z"]
        x = base_data.loc[i]["x"]
        y = base_data.loc[i]["y"]
        if data.loc[data["x"] == x].loc[data["y"] == y].empty:
            pass
        else:
            index = list(data.loc[data["x"] == x].loc[data["y"] == y].index)[0]
            data.at[index, "base"] = base
    data["distance"] = numpy.sqrt(data["x"] * data["x"] + data["y"] * data["y"])
    data = data.set_index("distance").interpolate(method='index').reset_index().dropna()
    return data