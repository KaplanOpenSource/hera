import math
from osgeo import gdal

#        The download function is not ready yet, manually download is needed nowdays
#        Get the GFS field list from a file
#        The file will download from here https://www.nco.ncep.noaa.gov/pmb/products/gfs
#        We will choose the 0.25 degree resolution
#        the downloaded file can be '/ibdata2/nirb/Projects/test.grib2'


class GFS:

    def get_gfs_list(gfs_file):
        """
	Get the list of all the fields in the grib file 
       
        param gfs_file - the file we will analyze
        
        return:
        list of the fields
        """

        dataset = gdal.Open(gfs_file, gdal.GA_ReadOnly)
        number_of_bands = dataset.RasterCount

        for i in range(1,number_of_bands+1):
            band = dataset.GetRasterBand(i)
            metadata = band.GetMetadata()   
            print(i,metadata['GRIB_ELEMENT'], metadata['GRIB_COMMENT'],metadata['GRIB_SHORT_NAME'])
        # close the file
        del dataset, band

    def get_gfs_data(gfs_file, lat=31., lon=33., band_num = 1):
        """
	Get the list of all the fields in the grib file 
       
        param gfs_file - the file we will analyze
        
        return:
        list of the fields
        """

        dataset = gdal.Open(gfs_file, gdal.GA_ReadOnly)

        band = dataset.GetRasterBand(band_num)  # 442 and 443 are u and v at 10m
        arr = band.ReadAsArray()
    
	# get origin's coordinates, pixel width and pixel height
	# the GetGeoTransform method returns the skew in the x and y axis but you
	# can ignore these values
        ox, pw, xskew, oy, yskew, ph = dataset.GetGeoTransform()
	# calculate the indices (row and column)
        i = int(math.floor(-(oy - lat) / ph))
        j = int(math.floor((lon - ox) / pw))

	# close the file
        del dataset, band

        # index the array to return the correspondent value
        return arr[i, j]


if __name__ == "__main__":
    mygfs = gfs()
    filename2 = '/ibdata2/nirb/testme.grib'
    gfs.get_gfs_list(filename2)
    lat = 31.76
    lon = 35.21
    band_num = 442
    data = gfs.get_gfs_data(filename2, lat=lat, lon=lon, band_num=band_num)
    print('at lat:',lat, 'lon:', lon, 'the data of band ', band_num, 'is ', data)

