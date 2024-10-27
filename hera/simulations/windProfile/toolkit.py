from hera import toolkit
import xarray as xr
import pandas as pd
import xarray
import numpy as np
from hera.utils.angle import toMathematicalAngle
from hera.utils import BETA,KARMAN,convertCRS,ITM,WSG84
from hera import toolkitHome
import json
from tqdm import tqdm
import requests
from hera.simulations.utils.interpolations import spatialInterpolate


class WindProfileToolkit(toolkit.abstractToolkit):
    def __init__(self, projectName, filesDirectory=None):
        super().__init__(projectName=projectName, toolkitName='WindProfileToolkit', filesDirectory=filesDirectory)

    def getWindProfile(self, df:pd.DataFrame, xarray:xarray.DataArray, height:float, dz:float):
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
        for _, row in df.iterrows():
            lat = row['lat']
            lon = row['lon']
            wind_speed = row['ws']
            wind_direction = row['wd']
            # z0 = df['z0']
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

    def getSpatialWind(self,minlon,minlat,maxlon,maxlat,IMS_TOKEN,dxdy=30,inputCRS=WSG84,landcover_DataSource=None):
        landcover_tk = toolkitHome.getToolkit(toolkitHome.GIS_LANDCOVER,projectName=self.projectName)
        xarray = landcover_tk.getLandCover(minlon, minlat, maxlon, maxlat, dxdy,inputCRS=inputCRS,dataSourceName=landcover_DataSource)
        xarray = landcover_tk.getRoughnessFromLandcover(xarray,dxdy)       #For now is simple - only with dictionary
        # stations = self._getStationsInRegion(minlon,minlat,maxlon,maxlat,inputCRS=inputCRS)
        # if len(stations)==0:
        #     raise ValueError(f"No Stations in Specified Region.")
        with open('wind_stations.json', 'r') as json_file:
            wind_stations = json.load(json_file)
        stations = [station for station in wind_stations]
        stations_with_data = self._getWindSpeedDirection(stations,IMS_TOKEN)
        return xarray,stations_with_data

    def _getStationsInRegion(self,minlon,minlat,maxlon,maxlat,inputCRS):
        min_pp = convertCRS(points=[[minlon, minlat]], inputCRS=inputCRS, outputCRS=ITM)[0]
        max_pp = convertCRS(points=[[maxlon, maxlat]], inputCRS=inputCRS, outputCRS=ITM)[0]
        with open('wind_stations.json', 'r') as json_file:
            wind_stations = json.load(json_file)
        stations_in_region = []
        for station in wind_stations:
            lat = station['attributes'][2]['value']['latitude']
            lon = station['attributes'][2]['value']['longitude']
            point_ITM = convertCRS(points=[[lon, lat]], inputCRS=WSG84, outputCRS=ITM)[0]
            if min_pp.x <= point_ITM.x and max_pp.x >= point_ITM.x and min_pp.y <= point_ITM.y and max_pp.y >= point_ITM.y:
                stations_in_region.append(station)

        return stations_in_region

    def _getWindSpeedDirection(self,stations,IMS_TOKEN):
        stations_with_data = []
        headers = {'Authorization': IMS_TOKEN['Authorization']}
        for station_dict in tqdm(stations):
            lat = station_dict['attributes'][2]['value']['latitude']
            lon = station_dict['attributes'][2]['value']['longitude']
            station_id = station_dict['attributes'][0]['value']
            height = station_dict['height_above_sea']
            trials = 0
            data = None
            while (trials < 20):
                try:
                    url = f'https://api.ims.gov.il/v1/envista/stations/{station_id}/data/daily'
                    response = requests.request('GET', url, headers=headers)
                    data = json.loads(response.text.encode('utf8'))
                    break
                except:
                    trials += 1
                    print(f"Trial {trials} for Station {station_id}")
            if data:
                for channel in data['data'][0]['channels']:
                    if channel['name'] == 'WS':
                        ws = channel['value']
                    if channel['name'] == 'WD':
                        wd = channel['value']

                stations_with_data.append([lat, lon, height, [ws, wd]])

        return stations_with_data

    def add_interpolated_ws_wd(self,xarray,stations):
        ws_wd = xr.apply_ufunc(
            self._interpolate_wd_ws,
            xarray['lon'], xarray['lat'],
            vectorize=True,
            kwargs={'stations_with_data': stations},
            output_core_dims=[[], []],
        )
        xarray['ws'] = ws_wd[0]
        xarray['wd'] = ws_wd[1]
        return xarray


    def _interpolate_wd_ws(self,lon,lat,stations_with_data):
        result = spatialInterpolate().interp([lon, lat, 10.0], stations_with_data)
        ws, wd = result[0], result[1]
        return float(ws), float(wd)