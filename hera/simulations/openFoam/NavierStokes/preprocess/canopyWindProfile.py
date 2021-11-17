import numpy
import pandas
import geopandas
import os

def urbanLogNormalProfile(cellCenters, lambdaGrid, stations):
    """
        Creates a log-normal wind profile in the mesh.


    Parameters
    ----------
    cellCenters: pandas.DataFrame
            The cell centers in the grid in which the results should be calculated.
            Has columns:
                - x
                - y
                - z

    lambdaGrid: pandas.DataFrame, geopandas.GeoDataFrame.

            The lambda of the urban grid.
            Has the following fields:
                - lambda_f
                - lambda_p
                - hx

            if GeoDataFrame: then take the x,y of the lambda from the centroid of the geometry.
            Otherwise, they have to be supplied.

            Calculates the Lc,ll,zz0, beta, and dd (displacenemnt).

    stations: pandas.DataFrame
            The wind speed and direction above building height

            Has the columns:
                x
                y
                speed10m
                direction (for now  is mathematical degrees).

            If has 1 station - define it along the entire domain.
            If has more than 1 station - interpolate each cell according to the inverse distance
            (takes the height differences into account).

    Returns
    -------
        pandas.DataFrame
            columns U_x, U_y hold the velocity in each field.

    """
    print(stations)
    stations = stations.reset_index()
    lambdaGrid = lambdaGrid.copy()
    aa_c = 0.25
    ccc = 1.0

    beta = 0.2
    norm_fac = 0.000000000000000000000000000001

    lambdaGrid.loc[(lambdaGrid.hc < 2), "lambdaF"]=0.25
    lambdaGrid.loc[(lambdaGrid.hc < 2), "lambdaP"] = 0.25
    lambdaGrid.loc[(lambdaGrid.hc < 2), "hc"] = 2
    lambdaGrid.loc[(lambdaGrid.lambdaF > 0.4), "lambdaF"] = 0.4
    lambdaGrid.loc[(lambdaGrid.lambdaP > 0.6), "lambdaF"] = 0.6
    lambdaGrid["Lc"] = lambdaGrid["hc"] * (1 - lambdaGrid["lambdaP"]) / lambdaGrid["lambdaF"]
    lambdaGrid["ll"] = 2 * (beta ** 3) * lambdaGrid["Lc"]
    lambdaGrid["zz0"] = lambdaGrid["ll"] / 0.4 * numpy.exp(-0.4 / beta)
    lambdaGrid["dd"] = lambdaGrid["ll"] / 0.4
    if isinstance(lambdaGrid, geopandas.geodataframe.GeoDataFrame):
        lambdaGrid["x"] = lambdaGrid.geometry.centroid.x
        lambdaGrid["y"] = lambdaGrid.geometry.centroid.y
    dx = lambdaGrid.reset_index()["x"][1] - lambdaGrid.reset_index()["x"][0]
    dy = lambdaGrid.reset_index()["y"][1] - lambdaGrid.reset_index()["y"][0]

    pi = numpy.pi
    tgr = 300.0
    aa2 = aa_c * (dx ** 2 + dy ** 2)

    stations["ust1"] = stations["speed10m"]*(numpy.cos(stations["MeteorologicalDirection"]))
    stations["vst1"] = stations["speed10m"]*(numpy.sin(stations["MeteorologicalDirection"]))
    stations["ust2"] = 1.5*stations["ust1"]
    stations["vst2"] = 1.5 * stations["vst1"]
    stations["temperature"] = tgr


    if len(stations)==1:
        lambdaGrid["U2h"] = float(stations["speed10m"].max())
        lambdaGrid["theta"] = 270. - float(stations["MeteorologicalDirection"].max())
        lambdaGrid["theta"] = lambdaGrid["theta"] * pi / 180.0
    else:
        lambdaGrid["sum_u"] = 0.0
        lambdaGrid["sum_v"] = 0.0
        lambdaGrid["total_sum"] = 0.0
        for i in range(len(stations)):
            lambdaGrid["r2"] = (lambdaGrid["x"] - stations["x"][i]) ** 2 + (lambdaGrid["y"] - stations["y"][i]) ** 2
            lambdaGrid["z_diff_sqr"] = (2 * lambdaGrid["h_c"] - stations["height"][i]) ** 2
            lambdaGrid["weight_h"] = 1. / (1. + lambdaGrid["r2"] / aa2)
            lambdaGrid["weight_v"] = 1. / (1. + ccc * lambdaGrid["z_diff_sqr"] / aa2)
            lambdaGrid["total_weight"] = lambdaGrid["weight_h"] * lambdaGrid["weight_v"]
            lambdaGrid["sum_u"] += lambdaGrid["total_weight"] * stations["ust1"][i]
            lambdaGrid["sum_v"] += lambdaGrid["total_weight"] * stations["vst1"][i]
            lambdaGrid["total_sum"] += lambdaGrid["total_weight"]

        lambdaGrid["uu"] = 1.0 * lambdaGrid["sum_u"] / (lambdaGrid["total_sum"] + norm_fac)
        lambdaGrid["vv"] = 1.0 * lambdaGrid["sum_v"] / (lambdaGrid["total_sum"] + norm_fac)
        lambdaGrid["U2h"] = numpy.sqrt((lambdaGrid["uu"]) ** 2 + (lambdaGrid["vv"]) ** 2)
        lambdaGrid["theta"] = numpy.atan2(lambdaGrid["vv"], lambdaGrid["uu"])

    lambdaGrid["Uh"] = (lambdaGrid["U2h"] * 0.4) / (beta * numpy.log((lambdaGrid["hc"] + lambdaGrid["dd"]) / lambdaGrid["zz0"]))
    partLambdaGrid = pandas.DataFrame({"x":lambdaGrid["x"],
                                       "y":lambdaGrid["y"],
                                       "Uh":lambdaGrid["Uh"],
                                       "hc":lambdaGrid["hc"],
                                       "ll":lambdaGrid["ll"],
                                       "dd":lambdaGrid["dd"],
                                       "zz0":lambdaGrid["zz0"],
                                       "U2h":lambdaGrid["U2h"],
                                       "theta":lambdaGrid["theta"]})

    lambdaXarray = partLambdaGrid.set_index(["x","y"]).to_xarray()
    Uh = []
    hc = []
    ll = []
    dd = []
    zz0 = []
    U2h = []
    theta = []
    valuesDict = {}

    for i in range(len(cellCenters)):
        x = cellCenters.loc[i]["x"]
        y = cellCenters.loc[i]["y"]
        if x not in valuesDict.keys():
            valuesDict[x] = {}
        if y in valuesDict[x].keys():
            for variableList, name in zip([Uh, hc, ll, dd, zz0, U2h, theta],
                                          ["Uh", "hc", "ll", "dd", "zz0", "U2h", "theta"]):
                variableList.append(valuesDict[x][y][name])
        else:
            valuesDict[x][y] = {}
            for variableList, name in zip([Uh,hc,ll,dd,zz0,U2h,theta],["Uh","hc","ll","dd","zz0","U2h","theta"]):
                newVal = float(lambdaXarray.interp(**{"x": x, "y": y}).fillna(0)[name])
                variableList.append(newVal)
                valuesDict[x][y][name] = newVal
        if i%5000==0 and i>=5000:
            print(f"finished {i} steps")
    for variableList, name in zip([Uh,hc,ll,dd,zz0,U2h,theta],["Uh","hc","ll","dd","zz0","U2h","theta"]):
        cellCenters[name] = variableList
    cellCenters["Umag"] = 0

    cellCenters.loc[(cellCenters.z < 2 * cellCenters.hc), "Umag"]= (cellCenters.loc[cellCenters.z > cellCenters.hc].loc[cellCenters.z < 2 * cellCenters.hc].Uh * beta / 0.4) * \
                                                                   numpy.log(
                                                                       1.0 * (cellCenters.loc[cellCenters.z > cellCenters.hc].loc[cellCenters.z < 2 * cellCenters.hc].z -
                                                                              cellCenters.loc[cellCenters.z > cellCenters.hc].loc[cellCenters.z < 2 * cellCenters.hc].hc +
                                                                              cellCenters.loc[cellCenters.z > cellCenters.hc].loc[cellCenters.z < 2 * cellCenters.hc].dd) /
                                                                       cellCenters.loc[cellCenters.z > cellCenters.hc].loc[cellCenters.z < 2 * cellCenters.hc].zz0)
    cellCenters.loc[(cellCenters.z < cellCenters.hc), "Umag"]= cellCenters.loc[cellCenters.z < cellCenters.hc].Uh \
                                                               * numpy.exp(beta * (cellCenters.loc[cellCenters.z < cellCenters.hc].z - cellCenters.loc[cellCenters.z < cellCenters.hc].hc) / cellCenters.loc[cellCenters.z < cellCenters.hc].ll)
    cellCenters.loc[(cellCenters.z > 2 * cellCenters.hc), "Umag"] = cellCenters.loc[cellCenters.z > 2 * cellCenters.hc].U2h

    cellCenters["U_x"] = numpy.cos(cellCenters.theta) * cellCenters.Umag
    cellCenters["U_y"] = numpy.sin(cellCenters.theta) * cellCenters.Umag

    return cellCenters