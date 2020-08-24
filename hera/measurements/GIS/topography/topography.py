import requests
import random
from osgeo import gdal
import numpy as np

def get_altitdue_ip(lat, lon):
	"""
	returning the altitude of the point, uses free mapquest data that is limited in the amount of calls per month, it uses Nir BAMBA Benmoshe key

	param lat - the path where we save the stl
	param lon - the width of the domain in the x direction

	return:
	altitude - meters above sea level
	"""

	resp = requests.get('http://open.mapquestapi.com/elevation/v1/profile?key=D5z9RSebQJLbUs4bohANIB4TzJdbvyvm&shapeFormat=raw&latLngCollection='+str(lat)+','+str(lon))


	height = resp.json()['elevationProfile'][0]['height']

	return height

def get_altitdue_gdal(lat, lon):

        # taken from https://earthexplorer.usgs.gov/
        fheight = r'/ibdata2/nirb/gt30e020n40.tif'
        ds = gdal.Open(fheight)
        myarray = np.array(ds.GetRasterBand(1).ReadAsArray())
        myarray[myarray<-1000]=0
        gt = ds.GetGeoTransform()
        rastery = int((lon - gt[0]) / gt[1])
        rasterx = int((lat - gt[3]) / gt[5])
        height = myarray[rasterx,rastery]

if __name__ == "__main__":
    lon = random.randint(35750, 35800) / 1000.0 # Hermon
    lat = random.randint(33250, 33800) / 1000.0
    alt = get_altitdue_ip(lat,lon)
    alt = get_altitdue_gdal(lat,lon)
    print("the altitude at position: ",lat,lon," is ", alt)


