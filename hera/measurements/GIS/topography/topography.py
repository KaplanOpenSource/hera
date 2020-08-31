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

#        if lat<=30: 
        # USGS EROS Archive - Digital Elevation - Global Multi-resolution Terrain Elevation Data 2010 (GMTED2010)
        fheight = r'/data3/nirb/10N030E_20101117_gmted_med075.tif'    
        fheight = r'/data3/nirb/30N030E_20101117_gmted_med075.tif'    
        #https://dds.cr.usgs.gov/srtm/version2_1/SRTM3/Africa/   # 90m resolution     
        if lat>29 and lat<30 and lon>34 and lon<35:
            fheight = r'/data3/nirb/N29E034.hgt'
        elif lat>29 and lat<30 and lon>35 and lon<36:
            fheight = r'/data3/nirb/N29E035.hgt'
        elif lat>30 and lat<31 and lon>34 and lon<35:
            fheight = r'/data3/nirb/N30E034.hgt'
        elif lat>30 and lat<31 and lon>35 and lon<36:
            fheight = r'/data3/nirb/N30E035.hgt'
        elif lat>31 and lat<32 and lon>34 and lon<35:
            fheight = r'/data3/nirb/N31E034.hgt'
        elif lat>31 and lat<32 and lon>35 and lon<36:
            fheight = r'/data3/nirb/N31E035.hgt'
        elif lat>32 and lat<33 and lon>34 and lon<35:
            fheight = r'/data3/nirb/N32E034.hgt'
        elif lat>32 and lat<33 and lon>35 and lon<36:
            fheight = r'/data3/nirb/N32E035.hgt'
        elif lat>33 and lat<33 and lon>35 and lon<36:
            fheight = r'/data3/nirb/N33E035.hgt'
        else:
            print ('!!!!NOT in Israel !!!!!!!!')            
            # taken from https://earthexplorer.usgs.gov/
            fheight = r'/ibdata2/nirb/gt30e020n40.tif'

        ds = gdal.Open(fheight)
        myarray = np.array(ds.GetRasterBand(1).ReadAsArray())
        myarray[myarray<-1000]=0
        gt = ds.GetGeoTransform()
        rastery = (lon - gt[0]) / gt[1]
        rasterx = (lat - gt[3]) / gt[5]
        height11 = myarray[int(rasterx),int(rastery)]
        height12 = myarray[int(rasterx)+1,int(rastery)]
        height21 = myarray[int(rasterx),int(rastery)+1]
        height22 = myarray[int(rasterx)+1,int(rastery)+1]
        height1 = (1. - (rasterx - int(rasterx))) * height11 + (rasterx - int(rasterx)) * height12
        height2 = (1. - (rasterx - int(rasterx))) * height21 + (rasterx - int(rasterx)) * height22
        height = (1. - (rastery - int(rastery))) * height1 + (rastery - int(rastery)) * height2
        
        return height

if __name__ == "__main__":
    lon = random.randint(35750, 35800) / 1000.0 # Hermon
    lat = random.randint(33250, 33800) / 1000.0
#    lon = 34.986008  // elevation should be 273m according to amud anan
#    lat = 32.808486  // elevation should be 273m according to amud anan
#    lon = 35.755  // elevation should be ~820 according to amud anan
#    lat = 33.459  // elevation should be ~820 according to amud anan  
    lon = 35.234987 # 744m
    lat = 31.777978 # 744m
    alt1 = get_altitdue_ip(lat,lon)
    alt2 = get_altitdue_gdal(lat,lon)
    print("the altitude at position: ",lat,lon," is ", alt1, alt2) 

