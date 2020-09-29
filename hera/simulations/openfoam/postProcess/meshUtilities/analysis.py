from .presentation import presentation
import numpy
import pandas
import xarray
from scipy.interpolate import griddata

class analysis():

    _datalayer = None
    _presentation = None

    @property
    def presentation(self):
        return self._presentation
    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, dataLayer):

        self._datalayer = dataLayer
        self._presentation = presentation(dataLayer=self.datalayer, Analysis=self)

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
        data = data.copy()
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
            data.loc[data.distance==dist,"terrain"]=base
        data = data.sort_values(by="distance").set_index("distance").interpolate(method='index').reset_index().dropna()
        data["heightOverTerrain"] = data["z"] - data["terrain"]
        data = data.query("heightOverTerrain>=0")

        self.datalayer.saveAndAddtoDB(save=save, addToDB=addToDB, data=data, path=path, key=key, filter=hdfName, **kwargs)
        return data

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
        ret['Distance'] = ([coord1, 'z'], distance)
        if time is None:
            ret = xarray.Dataset(ret, coords={coord1: C1, coord2: C2})
        else:
            ret = xarray.Dataset(ret, coords={coord1: C1, coord2: C2, 'time': time})
        if toPandas:
            ret = ret.to_dataframe()
        return ret

    def regularizeTimeStep_Sigma(self,data, fieldList, coord1='x',coord2='z', n=10, time=0, terrain="terrain",toPandas=True):
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

        if type(terrain)==str:
            if terrain not in data.columns:
                raise KeyError("Can't find a terrain. Deliver a name of the terrain column or a dataframe with the terrain properties.")
            else:
                terraindata = data.drop_duplicates(coord1)
                terrain = pandas.DataFrame({coord1:terraindata[coord1],coord2:terraindata[terrain]})
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
            retmap[field] = ([coord1, "sigma",'time'], numpy.atleast_3d(interpField))

        retmap['z'] = ([coord1, "sigma",'time'], numpy.atleast_3d(grid_z))
        if time is None:
            ret = xarray.Dataset(retmap, coords={"sigma":[sigma/n for sigma in range(n)], coord1:grid_x[:, 0]})
        else:
            ret = xarray.Dataset(retmap, coords={"sigma": [sigma / n for sigma in range(n)], coord1: grid_x[:, 0], 'time': [time]})
        if toPandas:
            ret = ret.to_dataframe()
        return ret

    def regularizeTimeStep_constantHeight(self,data, fieldList, coord1='x',coord2='z', time=0, terrain="terrain",heightsList=None,toPandas=True):
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
        heightsList = [i*10 for i in range(50)] if heightsList is None else heightsList
        if type(terrain)==str:
            if terrain not in data.columns:
                raise KeyError("Can't find a terrain. Deliver a name of the terrain column or a dataframe with the terrain properties.")
            else:
                terraindata = data.drop_duplicates(coord1)
                terrain = pandas.DataFrame({coord1:terraindata[coord1],coord2:terraindata[terrain]})
        terrain = terrain.sort_values(coord1)
        grid_x = numpy.array([numpy.array([x for i in range(len(heightsList))]) for x in terrain[coord1]])
        points = numpy.vstack([data[coord1], data[coord2]]).T

        grid_z = numpy.array([numpy.array([height+z for height in heightsList]) for z in terrain[coord2]])

        retmap = {}
        for field in numpy.atleast_1d(fieldList):
            try:
                interpField = griddata(points, data[field], (grid_x, grid_z), method='linear', fill_value=0)
            except KeyError:
                raise KeyError("field %s not found. Available keys are %s" % (field, ",".join(data.columns)))
            retmap[field] = ([coord1, "sigma",'time'], numpy.atleast_3d(interpField))

        retmap['z'] = ([coord1, "sigma",'time'], numpy.atleast_3d(grid_z))
        if time is None:
            ret = xarray.Dataset(retmap, coords={"heightOverTerrain":heightsList, coord1:grid_x[:, 0]})
        else:
            ret = xarray.Dataset(retmap, coords={"heightOverTerrain": heightsList, coord1: grid_x[:, 0], 'time': [time]})
        if toPandas:
            ret = ret.to_dataframe()
        return ret