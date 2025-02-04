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
import os
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from hera.measurements.GIS.utils import create_xarray

class WindProfileToolkit(toolkit.abstractToolkit):
    def __init__(self, projectName, filesDirectory=None):
        super().__init__(projectName=projectName, toolkitName='WindProfileToolkit', filesDirectory=filesDirectory)
        self._presentation = presentation(dataLayer=self)

    def getWindProfile(self,xarray:xarray.DataArray, height:float, dz:float):
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

        df = pd.DataFrame()
        df['lon'] = xarray['lon'].values.flatten()
        df['lat'] = xarray['lat'].values.flatten()
        df['ws'] = xarray['ws'].values.flatten()
        df['wd'] = xarray['wd'].values.flatten()

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

        lat_diff = np.abs(latitudes - lat)
        lon_diff = np.abs(longitudes - lon)

        total_diff = lat_diff + lon_diff
        min_index = np.unravel_index(np.argmin(total_diff), latitudes.shape)

        i, j = min_index
        return i,j

    def getSpatialWind(self,minlat,minlon,maxlat,maxlon,IMS_TOKEN,dxdy=30,inputCRS=WSG84,landcover_DataSource=None):
        landcover_tk = toolkitHome.getToolkit(toolkitHome.GIS_LANDCOVER,projectName=self.projectName)
        xarray = landcover_tk.getLandCover(minlat,minlon,maxlat,maxlon, dxdy,inputCRS=inputCRS,dataSourceName=landcover_DataSource)
        xarray = landcover_tk.getRoughnessFromLandcover(xarray,dxdy)       #For now is simple - only with dictionary
        with open(f'{os.path.dirname(os.path.abspath(__file__))}/wind_stations.json', 'r') as json_file:
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
        stations_datetimes = []
        for station_dict in tqdm(stations):
            lat = station_dict['attributes'][2]['value']['latitude']
            lon = station_dict['attributes'][2]['value']['longitude']
            station_id = station_dict['attributes'][0]['value']
            height = station_dict['height_above_sea']
            trials = 0
            data = None
            datetime_str = None
            while (trials < 20):
                try:
                    url = f'https://api.ims.gov.il/v1/envista/stations/{station_id}/data/latest'
                    response = requests.request('GET', url, headers=headers)
                    data = json.loads(response.text.encode('utf8'))
                    datetime_str = data['data'][0]['datetime']
                    print(f"{data['data'][0]['datetime']}, station: {station_dict['name']}")
                    break
                except:
                    trials += 1
                    # print(f"Trial {trials} for Station {station_id}")
            if data:
                for channel in data['data'][0]['channels']:
                    if channel['name'] == 'WS':
                        ws = channel['value']
                    if channel['name'] == 'WD':
                        wd = channel['value']

                if abs(ws) < 20 and abs(wd) <= 360:
                    stations_with_data.append([lat, lon, height, [ws, wd]])
                    stations_datetimes.append(datetime_str)

        datetime_objects = [datetime.fromisoformat(dt) for dt in stations_datetimes]
        max_datetime = max(datetime_objects)
        threshold = timedelta(minutes=15)
        filtered_stations = [
            stations_with_data[i] for i, dt in enumerate(datetime_objects) if (max_datetime - dt) < threshold
        ]
        return filtered_stations

    def add_interpolated_ws_wd(self,xarray,stations):
        ws_wd = xr.apply_ufunc(
            self._interpolate_wd_ws,
            xarray['lon'], xarray['lat'],xarray['elevation'],
            vectorize=True,
            kwargs={'stations_with_data': stations},
            output_core_dims=[[], []],
        )
        xarray['ws'] = ws_wd[0]
        xarray['wd'] = ws_wd[1]
        return xarray


    def _interpolate_wd_ws(self,lon,lat,elevation,stations_with_data):
        result = spatialInterpolate().interp([lat, lon, elevation + 10.0], stations_with_data)
        ws, wd = result[0], result[1]
        return float(ws), float(wd)

class presentation:
    """
    Presentation Layer of WindProfile1D toolkit.
    """
    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, dataLayer):
        self._datalayer = dataLayer

    def plotWindDirections(self,plot,landcover,figsize=(28,28),color='black',arrow_length=10,width=10,head_width=50,head_length=70,show_photo=True):
        fig, ax = plt.subplots(figsize=figsize)
        ax.imshow(np.array(plot.get_array()), extent=plot.get_extent(), origin='upper')
        for lat, lon, wd in zip(landcover.lat.values.flatten(),
                                landcover.lon.values.flatten(),
                                landcover.wd.values.flatten()):
            itm_pp = convertCRS([[lon, lat]], inputCRS=WSG84, outputCRS=ITM)[0]
            arrow_x_index, arrow_y_index = itm_pp.x, itm_pp.y
            arrow_direction = wd
            arrow_dx = np.cos(np.radians(arrow_direction)) * arrow_length
            arrow_dy = np.sin(np.radians(arrow_direction)) * arrow_length
            ax.arrow(arrow_x_index, arrow_y_index, arrow_dx * arrow_length, arrow_dy * arrow_length,
                     width=width, head_width=head_width, head_length=head_length, fc=color, ec=color)

        ax.set_xlim(plot.get_extent()[0], plot.get_extent()[1])
        ax.set_ylim(plot.get_extent()[2], plot.get_extent()[3])

        if show_photo:
            plt.show()
        else:
            plt.close(fig)
            return fig
