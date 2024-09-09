from hera import toolkit
import pandas as pd
import xarray
import numpy as np
from hera.utils.angle import toMathematicalAngle
from hera.utils import BETA,KARMAN

class WindProfileToolkit(toolkit.abstractToolkit):
    def __init__(self, projectName, filesDirectory=None):
        super().__init__(projectName=projectName, toolkitName='WindProfileToolkit', filesDirectory=filesDirectory)

    def getWindProfile(self, stations:pd.DataFrame, xarray:xarray.DataArray, height:float, dz:float):
        """
        Returns dataframe with u and v velocities in different specified heights.

        Parameters
        ----------
        stations: pd.DataFrame
            Stations dataframe, with lon, lat, wind and direction columns. Direction should be in meteorological direction.

        xarray: xarray.DataArray
            Xarray of stations region.

        height: float
            The maximum height to calculate velocities.

        dz: float
            Step size in meters.

        Returns
        -------
            pd.DataFrame
        """
        results = []
        fields = xarray.coords
        for _, station in stations.iterrows():
            lat = station['lat']
            lon = station['lon']
            wind_speed = station['wind']
            wind_direction = station['direction']
            i,j = self._find_lat_lon_index_in_xarray(lat,lon,xarray)
            z0 = float(xarray[i,j].z0.values)

            for z in np.arange(0, height + dz, dz):
                if 'dd' not in fields:
                    U = self.compute_log_wind_profile(z0, z, wind_speed)
                elif float(xarray[i,j].dd.values)==0:
                    U = self.compute_log_wind_profile(z0, z, wind_speed)
                else:
                    dd = float(xarray[i,j].dd.values)
                    hc = float(xarray[i,j].hc.values)
                    U = self.compute_exp_log_wind_profile(z,z0,hc,dd,wind_speed)

                u = U * np.cos(np.radians(toMathematicalAngle(wind_direction)))
                v = U * np.sin(np.radians(toMathematicalAngle(wind_direction)))

                speed = np.sqrt(u ** 2 + v ** 2)

                results.append({
                    'lat': lat,
                    'lon': lon,
                    'height': z,
                    'u': u,
                    'v': v,
                    'speed': speed,
                    'direction': wind_direction
                })

        return pd.DataFrame(results)
    def compute_exp_log_wind_profile(self,z,z0,hc,dd,wind_speed):
        return (wind_speed * KARMAN) / (BETA * np.log((z - hc + dd) / z0))

    def compute_log_wind_profile(self, z0, z, wind_speed):
        return (wind_speed * BETA / KARMAN) * np.log((z) / z0)

    def _find_lat_lon_index_in_xarray(self,lat,lon,xarray):
        latitudes = xarray['lat'].values
        longitudes = xarray['lon'].values

        # Find the index where lat=32 and lon=34
        lat_diff = np.abs(latitudes - lat)
        lon_diff = np.abs(longitudes - lon)

        # Find the combined minimum distance
        total_diff = lat_diff + lon_diff
        min_index = np.unravel_index(np.argmin(total_diff), latitudes.shape)

        # Extract i, j from the index
        i, j = min_index
        return i,j
