import numpy
import os
import pandas
import xarray
import dask.dataframe as dd
from scipy.interpolate import griddata
import glob

class coordinateHandler(object):
    """
    used to convert pandas dataframe to new coordinates.
    """

    def regularizeTimeSteps(self, data, fieldList=["U_x","U_z"], timelist=None, coordinateType="Regular", toPandas=True, **kwargs):
        """
        Converts a pandas dataframe of a slice into new grid points.

        data = The pandas dataframe.
        fieldList = A list of the columns being interpulated (default is U_x and U_z)
        timelist = A list of the time steps to be interpulated. The default is None, in which case the function assumes that all values belong to the same time step.
        coordinateType = The kind of new grid, supports "Regular" and "Sigma".
        toPandas = If true, returns a pandas dataframe, else, an xarray.
        **kwargs = Other parameters that are unique for differend coordinate types are listed in other functions.
        """

        method = getattr(self, "regularizeTimeStep_%s" % coordinateType)
        if timelist is None:
            if "time" in data.columns:
                timelist = data.time.drop_duplicates()
            else:
                timelist = [None]

        retlist = []
        for time in timelist:
            ret = method(data, fieldList=fieldList, time=time, toPandas=toPandas,**kwargs)
            retlist.append(ret)
        if toPandas:
            retlist = pandas.concat(retlist)
        return retlist

    def regularizeTimeStep_Regular(self, data, fieldList, time, coord1Lim=(None, None), coord2Lim=(None, None), coord1='x',
                                   coord2='z', n=(600, 600), toPandas=True):
        """
        Converts a pandas dataframe of a slice into regular grid points.

        data = The pandas dataframe.
        fieldList = A list of the columns being interpulated (default is U_x and U_z)
        time = The time step
        coord1Lim = a tuple with limits for the first coordinate.
        coord2Lim = a tuple with limits for the second coordinate.
        coord1 = the name of the first coordinate, default is "x".
        coord2 = the name of the second coordinate, default is "z".
        n = Number of points in each coordinate, a tuple.
        toPandas = If true, returns a pandas dataframe, else, an xarray.
        """

        if coord1Lim is None:
            c1min = data[coord1].min()
            c1max = data[coord1].max()
        else:
            c1min = data[coord1].min() if coord1Lim[0] is None else coord1Lim[0]
            c1max = data[coord1].max() if coord1Lim[1] is None else coord1Lim[1]

        if coord2Lim is None:
            c2min = data[coord2].min()
            c2max = data[coord2].max()
        else:
            c2min = data[coord2].min() if coord2Lim[0] is None else coord2Lim[0]
            c2max = data[coord2].max() if coord2Lim[1] is None else coord2Lim[1]

        grid_1, grid_2 = numpy.mgrid[numpy.min(c1min):numpy.max(c1max):complex(n[0]), \
                         numpy.min(c2min):numpy.max(c2max):complex(n[1])]
        points = numpy.vstack([data[coord1], data[coord2]]).T

        C1 = grid_1[:, 0]
        C2 = grid_2[0, :]

        # create the xarray
        ret = {}

        interpField_U = griddata(points, data['U_z'], (grid_1, grid_2), method='linear', fill_value=0)
        SurfaceIndex = []
        for i in range(len(C2)):
            SurfaceIndex.append(numpy.argmax(interpField_U[:, i] > 0) - 1)

        surface = numpy.take(grid_1[:, 0], SurfaceIndex)
        distance = grid_1 - surface

        for field in numpy.atleast_1d(fieldList):
            try:
                interpField = griddata(points, data[field], (grid_1, grid_2), method='linear', fill_value=0)
            except KeyError:
                raise KeyError("field %s not found. Available keys are %s" % (field, ",".join(data.columns)))
            ret[field] = ([coord1, coord2], interpField)

        ret['surfaceIndex'] = ([coord2], SurfaceIndex)
        ret['surfaceX'] = ([coord2], surface)
        ret['distance'] = ([coord1, 'z'], distance)
        if time is None:
            ret = xarray.Dataset(ret, coords={coord1: C1, coord2: C2})
        else:
            ret = xarray.Dataset(ret, coords={coord1: C1, coord2: C2, 'time': time})
        if toPandas:
            ret = ret.to_dataframe()
        return ret

    def regularizeTimeStep_Sigma(self,data, fieldList, coord1='x',coord2='z', n=10, time=0, terrain=None,toPandas=True):
        """
        Converts a pandas dataframe of a slice into sigma terrain-following grid points.

        data = The pandas dataframe.
        fieldList = A list of the columns being interpulated (default is U_x and U_z)
        time = The time step
        coord1 = the name of the first coordinate, default is "x".
        coord2 = the name of the second coordinate, default is "z".
        n = Number of sigma values between 0 to 1.
        toPandas = If true, returns a pandas dataframe, else, an xarray.
        terrain = a dataframe with x and z values of the terrain.
        """
        if terrain is None:
            if "terrain" not in data.columns:
                raise KeyError("Can't find a terrain. Deliver a name of the terrain column or a dataframe with the terrain properties.")
            else:
                terrain = pandas.DataFrame({coord1:data[coord1].drop_duplicates(),coord2:data.terrain.drop_duplicates()})
        elif type(terrain)==str:
            if terrain not in data.columns:
                raise KeyError("Can't find a terrain. Deliver a name of the terrain column or a dataframe with the terrain properties.")
            else:
                terrain = pandas.DataFrame({coord1:data[coord1].drop_duplicates(),coord2:data[terrain].drop_duplicates()})
        c2max = data[coord2].max()
        terrain = terrain.sort_values(coord1)
        grid_x = numpy.array([numpy.array([x for i in range(n)]) for x in terrain[coord1]])
        points = numpy.vstack([data[coord1], data[coord2]]).T

        grid_z = numpy.array([numpy.array([(sigma/n)*(c2max-z)+z for sigma in range(n)]) for z in terrain[coord2]])

        retmap = {}
        for field in numpy.atleast_1d(fieldList):
            try:
                interpField = griddata(points, data[field], (grid_x, grid_z), method='linear', fill_value=0)
            except KeyError:
                raise KeyError("field %s not found. Available keys are %s" % (field, ",".join(data.columns)))
            retmap[field] = (['x', "sigma",'time'], numpy.atleast_3d(interpField))

        retmap['z'] = (['x', "sigma",'time'], numpy.atleast_3d(grid_z))
        if time is None:
            ret = xarray.Dataset(retmap, coords={"sigma":[sigma/n for sigma in range(n)], 'x':grid_x[:, 0]})
        else:
            ret = xarray.Dataset(retmap, coords={"sigma": [sigma / n for sigma in range(n)], 'x': grid_x[:, 0], 'time': [time]})
        if toPandas:
            ret = ret.to_dataframe()
        return ret