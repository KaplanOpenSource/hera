import numpy
import pandas
import geopandas

def windProfile(cellData, LambdaGrid,stations):

    stations = stations.reset_index()
    cellData = cellData.copy()
    LambdaGrid = LambdaGrid.copy()
    aa_c = 0.25
    ccc = 1.0

    beta = 0.2
    norm_fac = 0.000000000000000000000000000001

    LambdaGrid.loc[(LambdaGrid.hc<2),"lambda_f"]=0.25
    LambdaGrid.loc[(LambdaGrid.hc < 2), "lambda_p"] = 0.25
    LambdaGrid.loc[(LambdaGrid.hc < 2), "hc"] = 2
    LambdaGrid.loc[(LambdaGrid.lambda_f > 0.4), "lambda_f"] = 0.4
    LambdaGrid.loc[(LambdaGrid.lambda_p > 0.6), "lambda_f"] = 0.6
    LambdaGrid["Lc"] = LambdaGrid["hc"] * (1 - LambdaGrid["lambda_p"]) / LambdaGrid["lambda_f"]
    LambdaGrid["ll"] = 2 * (beta ** 3) * LambdaGrid["Lc"]
    LambdaGrid["zz0"] = LambdaGrid["ll"] / 0.4 * numpy.exp(-0.4 / beta)
    LambdaGrid["dd"] = LambdaGrid["ll"] / 0.4
    if isinstance(LambdaGrid,geopandas.geodataframe.GeoDataFrame):
        LambdaGrid["x"] = LambdaGrid.geometry.centroid.x
        LambdaGrid["y"] = LambdaGrid.geometry.centroid.y
    dx = LambdaGrid.reset_index()["x"][1]-LambdaGrid.reset_index()["x"][0]
    dy = LambdaGrid.reset_index()["y"][1] - LambdaGrid.reset_index()["y"][0]

    pi = numpy.pi
    tgr = 300.0
    aa2 = aa_c * (dx ** 2 + dy ** 2)

    stations["ust1"] = stations["velocity"]*(numpy.cos(stations["direction"]))
    stations["vst1"] = stations["velocity"]*(numpy.sin(stations["direction"]))
    stations["ust2"] = 1.5*stations["ust1"]
    stations["vst2"] = 1.5 * stations["vst1"]
    stations["temperature"] = tgr

    if len(stations)==1:
        LambdaGrid["U2h"] = float(stations["velocity"].max())
        LambdaGrid["theta"] = 270. - float(stations["direction"].max())
        LambdaGrid["theta"] = LambdaGrid["theta"] * pi / 180.0
    else:
        LambdaGrid["sum_u"] = 0.0
        LambdaGrid["sum_v"] = 0.0
        LambdaGrid["total_sum"] = 0.0
        for i in range(len(stations)):
            LambdaGrid["r2"] = (LambdaGrid["x"] - stations["x"][i]) ** 2 + (LambdaGrid["y"] - stations["y"][i]) ** 2
            LambdaGrid["z_diff_sqr"] = (2 * LambdaGrid["h_c"] - stations["height"][i]) ** 2
            LambdaGrid["weight_h"] = 1. / (1. + LambdaGrid["r2"] / aa2)
            LambdaGrid["weight_v"] = 1. / (1. + ccc * LambdaGrid["z_diff_sqr"] / aa2)
            LambdaGrid["total_weight"] = LambdaGrid["weight_h"] * LambdaGrid["weight_v"]
            LambdaGrid["sum_u"] += LambdaGrid["total_weight"] * stations["ust1"][i]
            LambdaGrid["sum_v"] += LambdaGrid["total_weight"] * stations["vst1"][i]
            LambdaGrid["total_sum"] += LambdaGrid["total_weight"]

        LambdaGrid["uu"] = 1.0 * LambdaGrid["sum_u"] / (LambdaGrid["total_sum"] + norm_fac)
        LambdaGrid["vv"] = 1.0 * LambdaGrid["sum_v"] / (LambdaGrid["total_sum"] + norm_fac)
        LambdaGrid["U2h"] = numpy.sqrt((LambdaGrid["uu"]) ** 2 + (LambdaGrid["vv"]) ** 2)
        LambdaGrid["theta"] = numpy.atan2(LambdaGrid["vv"], LambdaGrid["uu"])

    LambdaGrid["Uh"] = (LambdaGrid["U2h"] * 0.4) / (beta * numpy.log((LambdaGrid["hc"] + LambdaGrid["dd"]) / LambdaGrid["zz0"]))
    partLambdaGrid = pandas.DataFrame({"x":LambdaGrid["x"],"y":LambdaGrid["y"],"Uh":LambdaGrid["Uh"],"hc":LambdaGrid["hc"],
                                       "ll":LambdaGrid["ll"],"dd":LambdaGrid["dd"],"zz0":LambdaGrid["zz0"],
                                       "U2h":LambdaGrid["U2h"],"theta":LambdaGrid["theta"]})
    lambdaXarray = partLambdaGrid.set_index(["x","y"]).to_xarray()
    Uh = []
    hc = []
    ll = []
    dd = []
    zz0 = []
    U2h = []
    theta = []
    valuesDict = {}
    for i in range(len(cellData)):
        x = cellData.loc[i]["x"]
        y = cellData.loc[i]["y"]
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
    for variableList, name in zip([Uh,hc,ll,dd,zz0,U2h,theta],["Uh","hc","ll","dd","zz0","U2h","theta"]):
        cellData[name] = variableList
    cellData["Umag"] = 0

    cellData.loc[(cellData.z<2*cellData.hc),"Umag"]= (1.0 * cellData.loc[cellData.z>cellData.hc].loc[cellData.z<2*cellData.hc].Uh * beta / 0.4) * \
                                                                                numpy.log(
                                                                   1.0 * (cellData.loc[cellData.z>cellData.hc].loc[cellData.z<2*cellData.hc].z -
                                                                          cellData.loc[cellData.z>cellData.hc].loc[cellData.z<2*cellData.hc].hc +
                                                                          cellData.loc[cellData.z>cellData.hc].loc[cellData.z<2*cellData.hc].dd) /
                                                                          cellData.loc[cellData.z>cellData.hc].loc[cellData.z<2*cellData.hc].zz0)
    cellData.loc[(cellData.z<cellData.hc),"Umag"]= cellData.loc[cellData.z<cellData.hc].Uh * numpy.exp(
                                                            beta * (cellData.loc[cellData.z<cellData.hc].z -
                                                            cellData.loc[cellData.z<cellData.hc].hc) /cellData.loc[cellData.z<cellData.hc].ll)
    cellData.loc[(cellData.z > 2*cellData.hc), "Umag"] = cellData.loc[cellData.z>2*cellData.hc].U2h

    cellData["U_x"] = numpy.cos(cellData.theta) * cellData.Umag
    cellData["U_y"] = numpy.sin(cellData.theta) * cellData.Umag

    return cellData