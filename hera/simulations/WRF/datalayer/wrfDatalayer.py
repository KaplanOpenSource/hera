import numpy
import wrf
import pandas
import geopandas
from netCDF4 import Dataset
import xarray

class wrfDatalayer():

    def find_i(self, request, xdata, dimension, coordinate):
        """
        Finds the i value of a dimension that may interpolate a requested value of a coordinate,
        for example, dimension west_east for coordinate XLONG

        Parameters
        ----------
        request: The requested value
        xdata: An xarray with the data
        dimension: The dimension, for example, west_east
        coordinate: The coordinate, for example, XLONG

        Returns: The value that may be used for an interpolation.
        -------

        """

        for i in range(len(xdata.__getattr__(dimension))):
            compare_dict = {dimension: i}
            compare_dictplus = {dimension: i + 1}
            if request > float(xdata.isel(Time=0, bottom_top=0, **compare_dict).__getattr__(coordinate).min()) - 0.01 and\
               request < float(xdata.isel(Time=0, bottom_top=0, **compare_dict).__getattr__(coordinate).max()) + 0.01:
                delta = float(xdata.isel(Time=0, bottom_top=0, **compare_dictplus).__getattr__(coordinate).mean()) - \
                        float(xdata.isel(Time=0, bottom_top=0, **compare_dict).__getattr__(coordinate).mean())
                request_delta = request - float(xdata.isel(Time=0, bottom_top=0, **compare_dict).__getattr__(coordinate).mean())
                break

        return i + request_delta / delta

    def getPandas(self, datapath, Time=None, lat=None, lon=None, heightLimit=None, compare_lat=650000, compare_lon=200000):
        """
        Makes a pandas dataframe that holds data of wrf simulation.
        Parameters
        ----------
        datapath: The path of the wrf results.
        Time: A datetime.time object with time in hours and minutes, or None for all times.
        lat: An ITM latitude coordinate. Either one of lat or lan must be given.
        lon: An ITM latitude coordinate. Either one of lat or lan must be given.
        heightLimit: An optional height limit
        compare_lat: A coordinate that is used to get the wpf-84 coordinates for the function's operation.
        compare_lon: A coordinate that is used to get the wpf-84 coordinates for the function's operation.

        Returns: A pandas dataframe.
        -------

        """

        d = pandas.DataFrame()
        data = Dataset(datapath)
        xdata = xarray.open_dataset(datapath)

        if Time is not None:
            pStart = pandas.Timestamp(wrf.getvar(data, "times").values)
            pTime = pandas.Timestamp(year=pStart.year, month=pStart.month, day=pStart.day, hour=Time.hour, minute=Time.minute)
            for i in range(len(xdata.Time)):
                if pTime >= wrf.getvar(data, "times", timeidx=i).values and  pTime <= wrf.getvar(data, "times", timeidx=i+1).values:
                    delta = (wrf.getvar(data, "times", timeidx=i+1).values-wrf.getvar(data, "times", timeidx=i).values)/numpy.timedelta64(1, 's')
                    request_delta = (pTime - wrf.getvar(data, "times", timeidx=i).values).total_seconds()
                    request_i_time = [i + request_delta/delta]
                    break
        else:
            request_i_time = [x for x in range(len(xdata.Time))]

        if lat is not None:
            geo = geopandas.GeoDataFrame(dict(geometry=geopandas.points_from_xy([compare_lon],[lat])),index=[0])
            geo.crs = 2039
            geo = geo.to_crs(epsg=4326)
            lat = geo.geometry[0].y
            changes = ["south_north", "south_north", "south_north_stag"]
            request_i_u = request_i = self.find_i(lat, xdata, "south_north", "XLAT")
            request_i_v = self.find_i(lat, xdata, "south_north_stag", "XLAT_V")

        elif lon is not None:
            geo = geopandas.GeoDataFrame(dict(geometry=geopandas.points_from_xy([lon],[compare_lat])),index=[0])
            geo.crs = 2039
            geo = geo.to_crs(epsg=4326)
            lon = geo.geometry[0].x
            changes = ["west_east", "west_east_stag", "west_east"]
            request_i_v = request_i = self.find_i(lon, xdata, "west_east", "XLONG")
            request_i_u = self.find_i(lon, xdata, "west_east_stag", "XLONG_U")
        else:
            raise KeyError("Choose longtitude or latitude.")

        H_vals = []
        if heightLimit is not None:
            query_dict = {changes[0]: request_i}
            ht = wrf.getvar(data, "z")
            for i in range(len(xdata.bottom_top)):
                if ht.isel(bottom_top=i).interp(**query_dict).min()<heightLimit:
                    H_vals.append(i)
        else:
            H_vals.append(len(xdata.bottom_top))


        ter = wrf.getvar(data, "ter")

        for i in request_i_time:
            if Time is None:
                height = wrf.getvar(data, "z", timeidx=i)
            else:
                height = wrf.getvar(data, "z", timeidx=wrf.ALL_TIMES).interp(Time=pTime)

            for j in range(max(H_vals)):

                query_dict = {changes[0]: request_i}
                query_dictU = {changes[1]: request_i_u}
                query_dictV = {changes[2]: request_i_v}
                xnewdata = xdata.isel(bottom_top=j, bottom_top_stag=j).interp(Time=i, **query_dict)
                if lat is not None:
                    new_u = xdata.isel(bottom_top=j).interp(Time=i, **query_dictU).U[:-1]
                    new_v = xdata.isel(bottom_top=j).interp(Time=i, **query_dictV).V
                elif lon is not None:
                    new_u = xdata.isel(bottom_top=j).interp(Time=i, **query_dictU).U
                    new_v = xdata.isel(bottom_top=j).interp(Time=i, **query_dictV).V[:-1]
                if Time is None:
                    new_t = wrf.getvar(data, "time", timeidx=i).values
                else:
                    new_t = pTime

                new_d = pandas.DataFrame(dict(height = height.isel(bottom_top=j).interp(**query_dict),
                                              U=new_u,
                                              V=new_v,
                                              W=xnewdata.W,
                                              Time=new_t,
                                              terrain=ter.interp(**query_dict),
                                              LAT=xnewdata.XLAT,
                                              LONG=xnewdata.XLONG))
                d = pandas.concat([d, new_d])

        gdf = geopandas.GeoDataFrame(d, geometry=geopandas.points_from_xy(d.LONG, d.LAT))
        gdf.crs = 4326
        gdf = gdf.to_crs(epsg=2039)
        gdf["LAT"] = gdf.geometry.y
        gdf["LONG"] = gdf.geometry.x
        gdf["height_over_terrain"] = gdf.height - gdf.terrain

        return gdf
