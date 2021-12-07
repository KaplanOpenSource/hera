import dask.dataframe
import numpy
import pandas
import geopandas
import os
from .....utils.angle import toMathematicalAngle

karman = 0.41
def urbanLogExponentProfile(cellCenters, lambdaGrid, stations):
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
                direction (for now is mathematical degrees).

            We currently assume that there is 1 station. It will be extended in the future.

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

        # these variables help us to fine the center of the lambda matrix cell
        i01 = row.i01
        i01max = row.i01max
        j01 = row.j01
        j01max = row.j01max

        surfaceLayerHeight = 2000  # previously it was 200m, but it should be high above the buildings height

        if i01 != i01max//2 or j01 != j01max//2:
            # we want to put the profile only in the lambda matrix cells center
            ret = numpy.nan
        else:
            if z < hc:
                ret = Uh*numpy.exp(beta * (z - hc) / ll)
            elif (z > hc) and (z < surfaceLayerHeight):
                ret =  (Uh * beta / karman) * numpy.log((z - hc + dd) / zz0)
            else:  # Above the surface layer
                z = surfaceLayerHeight
                ret = (Uh * beta / karman) * numpy.log((z - hc + dd) / zz0)

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
    lambdaGrid.loc[(lambdaGrid.lambdaP > 0.6), "lambdaP"] = 0.6
    lambdaGrid["Lc"] = lambdaGrid["hc"] * (1 - lambdaGrid["lambdaP"]) / lambdaGrid["lambdaF"]
    lambdaGrid["ll"] = 2 * (beta ** 3) * lambdaGrid["Lc"]
    lambdaGrid["zz0"] = lambdaGrid["ll"] / karman * numpy.exp(-karman / beta)
    lambdaGrid["dd"] = lambdaGrid["ll"] / karman

    minx,miny,maxx,maxy = lambdaGrid.iloc[0].geometry.bounds
    # dx and dy of the lambda matrix cells
    dx = maxx-minx
    dy = maxy-miny

    ## Here we assume that there is only 1 station.
    #U = float(stations["speed10m"].max())
    U = float(stations["speedgeo"].max())
    geoheight=800.
    theta = numpy.deg2rad(toMathematicalAngle(stations["MeteorologicalDirection"].max()))
    #lambdaGrid["Uh"] = (U * karman) / (beta * numpy.log((lambdaGrid["hc"] + lambdaGrid["dd"]) / lambdaGrid["zz0"]))
    lambdaGrid["Uh"] = (U * karman) / (beta * numpy.log((geoheight - lambdaGrid["hc"] + lambdaGrid["dd"]) / lambdaGrid["zz0"]))

    cellCenters = cellCenters.assign(i0=(cellCenters.x-cellCenters.x.min())//dx,
                                     j0=(cellCenters.y-cellCenters.y.min())//dy)

    if isinstance(cellCenters,dask.dataframe.DataFrame):
        cellCenters = cellCenters.compute()

    # dx and dy of the domain cells
    dx2=cellCenters.groupby('x').x.mean().iloc[1]-cellCenters.groupby('x').x.mean().iloc[0] # dx of the domain
    dy2=cellCenters.groupby('y').y.mean().iloc[1]-cellCenters.groupby('y').y.mean().iloc[0] # dy of the domain

    # calculate the cells index inside a lambda cell
    cellCenters["i01"] = ((cellCenters.x-cellCenters.x.min())%dx)//dx2 #index of domain cell in lambda cell
    cellCenters["j01"] = ((cellCenters.y-cellCenters.y.min())%dy)//dy2 #index of domain cell in lambda cell

    # find the maximum index, it help us find the center (max divide into two)
    cellCenters["i01max"] = cellCenters.i01.max()
    cellCenters["j01max"] = cellCenters.j01.max()

    # Sort by old index to maintain openfoam numbering
    data = cellCenters.reset_index().merge(lambdaGrid,on=['i0','j0']).sort_values('index')

    data["Umag"] = data.apply(calcU,axis=1)
    data["Ux"] =  numpy.cos(theta) * data.Umag
    data["Uy"] = numpy.sin(theta) * data.Umag

#    data.z = data.z.round(4) # we want to round the elevation info of each cell
    interpolation = 'linear'
    # we want to interpolate and smooth the wind profile between the centers of cells, it will make easy life for the mass consistency code
    if interpolation=='linear':
        xdata = data.to_xarray()
        xdata = xdata.set_coords(['x','y','z'])
        xdata['Ux'] = xdata.Ux.interpolate_na(dim='index', method="linear", fill_value="extrapolate")
        xdata['Uy'] = xdata.Uy.interpolate_na(dim='index', method="linear", fill_value="extrapolate")
        data['Ux'] = xdata.Ux.data
        data['Uy'] = xdata.Uy.data
    if interpolation=='Elevation':
        print('I didnt finish write the code of interpolation using edw although it is better due to perfomance issues')
        a=100.
        c=5.
        pointswithdata = data['Ux'].isnull()
        k=0
        m=0
        for i in range(len(data)):
            w=0
            ix = data.iloc[i].x
            iy = data.iloc[i].y
            iz = data.iloc[i].z
            for j in range(len(data)):
                dx = ix - data.iloc[j].x
                dy = iy - data.iloc[j].y
                dz = iz - data.iloc[j].z
                r2= dx*dx+dy*dy
                if pointswithdata[i]:
                    if r2>0:
                        w += (1./(1.+c*dz*dz/a/a))*(1./r2)
            if i%1==0:
                print(i,w)
        print('fin')
    return data
