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
                U_star = (wind_speed * KARMAN) / np.log(z / z0)
                if 'hc' in fields:
                    hc = float(xarray[i,j].hc.values)
                    if hc > 2.0:          #Urban Area
                        if z > hc:
                            U_z = (U_star / KARMAN) * np.log(z / z0)
                        else:
                            U_hc = (U_star/KARMAN) * np.log(hc/z0)
                            U_z = U_hc * np.exp(BETA * (z - hc) / xarray[i,j].ll.values)
                    else:
                        continue
                else:                       #Non-Urban Area
                    U_z = (U_star / KARMAN) * np.log(z / z0)

                u = U_z * np.cos(np.radians(toMathematicalAngle(wind_direction)))
                v = U_z * np.sin(np.radians(toMathematicalAngle(wind_direction)))
                # speed = np.sqrt(u ** 2 + v ** 2)
                results.append({
                    'lat': lat,
                    'lon': lon,
                    'height': z,
                    'u': u,
                    'v': v,
                    'U_z': U_z,
                    'direction': wind_direction
                })

        return pd.DataFrame(results)

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
