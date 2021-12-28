import dask.dataframe
import numpy
import pandas
import geopandas
import os
from .....utils.angle import toMathematicalAngle

def urbanLogNormalProfile(cellCenters, lambdaGrid, stations):
    """
        Creates a log-normal wind profile in the mesh.

        TODO: seemd wrong - check equations.


    Parameters
    ----------
    cellCenters: pandas.DataFrame
            The cell centers in the grid in which the results should be calculated.
            Has columns:
                - x
                - y
                - z

            TODO: We currently assume that the topography is flat!!.
                  When there is topography, we have to use the height above the ground and not z (which is absolute).


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

            We currently assumethat there is 1 station. It will be extended in the future.

    Returns
    -------
        pandas.DataFrame
            columns U_x, U_y hold the velocity in each field.

    """
    def calcU(row):
        z = row.z
        hc = row.hc
        dd = row.dd
        ll = row.ll
        Uh = row.Uh
        zz0 = row.zz0
        beta = 0.2

        surfaceLayerHeight = 200

        if (z<hc):
            ret = Uh*numpy.exp(beta * (z - hc) / ll)
        elif (z>hc) and (z<surfaceLayerHeight):
            ret = (Uh * beta / 0.4) * numpy.log((z - hc + dd) /zz0)
        else: # Above the surface layer
            z = surfaceLayerHeight
            ret = (Uh * beta / 0.4) * numpy.log((z - hc + dd) /zz0)
        return ret


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

    minx,miny,maxx,maxy = lambdaGrid.iloc[0].geometry.bounds
    dx = maxx-minx
    dy = maxy-miny

    ## Here we assume that there is only 1 station.
    U = float(stations["speed10m"].max())
    theta = numpy.deg2rad(toMathematicalAngle(stations["MeteorologicalDirection"].max()))
    lambdaGrid["Uh"] = (U * 0.4) / (beta * numpy.log((lambdaGrid["hc"] + lambdaGrid["dd"]) / lambdaGrid["zz0"]))


    cellCenters = cellCenters.assign(i0=(cellCenters.x-cellCenters.x.min())//dx,
                                     j0=(cellCenters.y-cellCenters.y.min())//dy)

    if isinstance(cellCenters,dask.dataframe.DataFrame):
        cellCenters = cellCenters.compute()

    # Sort by old index to maintain openfoam numbering
    data = cellCenters.reset_index().merge(lambdaGrid,on=['i0','j0']).sort_values('index')

    data["Umag"] = data.apply(calcU,axis=1)
    data["Ux"] =  numpy.cos(theta) * data.Umag
    data["Uy"] = numpy.sin(theta) * data.Umag
    return data
