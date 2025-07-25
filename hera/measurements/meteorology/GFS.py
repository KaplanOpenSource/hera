import datetime
import sklearn.metrics
import math
import numpy as np    
from osgeo import gdal


#        The download function is not ready yet, manually download is needed nowdays
#        Get the GFS field list from a file
#        The file will download from here https://www.nco.ncep.noaa.gov/pmb/products/gfs
#        We will choose the 0.25 degree resolution
#        the downloaded file can be '/ibdata2/nirb/Projects/test.grib2'


class GFS:

    def get_gfs_list(self,gfs_file):
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

    def get_gfs_data(self, gfs_file, lat=31., lon=33., band_num = 1):
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
    
    
# https://stackoverflow.com/questions/67963199/xarray-from-grib-file-to-dataset    
        import cfgrib
        
        
        s = cfgrib.open_file(filename2)
        
        grib_data = cfgrib.open_datasets(filename2,backend_kwargs={'filter_by_keys': {'typeOfLevel': 'heightAboveGroundLayer'}})
        grib_data = cfgrib.open_datasets(filename2,backend_kwargs={'filter_by_keys': {'typeOfLevel': 'surface'}})
        sorted(grib_data.dimensions.items())
        grib_data[32]['v']
        
        
        ds = xr.load_dataset(filename2, engine="cfgrib",backend_kwargs={'filter_by_keys': {'typeOfLevel': 'surface'}})

        d = xr.open_dataset     (filename2,backend_kwargs={'filter_by_keys':{'typeOfLevel': 'heightAboveGroundLayer'}})


def getProbe(file, jumps = 1):
    # get openfoam probe data and convert it to python data
    f = open(file+"U","r")
    lines = f.readlines()
    times = []
    vectors = []
    probes=0
    while (lines[probes][:5]=="# Pro"):
        probes+=1
    # if len(lines)>100000:
    #     jumps = len(lines)//100000
    for i in range(len(lines)):
        if i % jumps ==0:
            if(i>probes):
                #first 20 chars of line are times
                times.append(float(lines[i][0:14]))
                #remaining chars are vector
                vector = lines[i][14:]
                #remove ( and )
                vector = vector.replace("\n", "")
                vector = vector.replace("(", "")
                vector = vector.replace(")", "")
                #split by space
                vectors.append(vector.split(" "))
    return times, vectors


if __name__ == "__main__":
    mygfs = GFS()
    filename2 = '/ibdata2/nirb/testme.grib'
    filename2 = '/ibdata2/nirb/gfs.t06z.pgrb2.0p25.f002'
    gfs_file = filename2 
    mygfs.get_gfs_list(filename2)
    lat = 31.76
    lon = 35.21
    band_num = 442
    data = mygfs.get_gfs_data(filename2, lat=lat, lon=lon, band_num=band_num)
    print('at lat:',lat, 'lon:', lon, 'the data of band ', band_num, 'is ', data)
    
    import netCDF4 as nc
    filenc='/ibdata2/nirb/gdas.t06z.atmf003.nc'
    filenc='/data4bk/nirb/Simulations/Haifa/gdas.t12z.atmf003.nc'
    filenc='/data4bk/nirb/Simulations/Haifa/gdas.t18z.atmf000.nc' # 2023061318gmt
    ds = nc.Dataset(filenc)
    print(ds.variables.keys())
    print(ds['ugrd'].shape)
    print(ds['pfull'])
    lat = 32.84 # Haifa
    lon = 34.94 # Haifa
    latindex=487
    longindex=298 
    pfull = ds.variables['pfull'][:].data
    delz = ds.variables['delz'][0,:,latindex,longindex].data
    latspan=ds.variables['lat'][:,0].data
    lonspan=ds.variables['lon'][0,:].data # zero is in UK
    
    print('delz')
    z=np.zeros(len(delz))
    ztotal=0
    for i in range(len(delz)):
        ztotal-=delz[len(delz)-1-i]
        z[len(delz)-1-i]=ztotal
        print(i,len(delz)-1-i,z[i],delz[i],ztotal)
    print(np.flipud(z[83:127]))
    print(np.flipud(ds['ugrd'][0,:,latindex,longindex].data[83:127]))
    print(np.flipud(ds['vgrd'][0,:,latindex,longindex].data[83:127]))
    plt.figure()
    plt.imshow(ds['ugrd'][0,126,:,:])
    plt.title('ugrd')
    plt.colorbar()
    plt.show()
    plt.figure()
    plt.imshow(ds['ugrd'][0,126,:,:])
    plt.title('vgrd')
    plt.colorbar()
    plt.show()
    
    uu = np.asarray([0.24, 1.62, 1.48, 0.64, 2.78, 2.98])
    vv = np.asarray([-3.98, -2.36, -1.89, -0.37, -1.43, -0.35])
    uu = np.asarray([.35,1.21, 1.29, 2.96, 1.65, 1.90])
    vv = np.asarray([-3.34, -1.36, -1.02, -1.20, -1.22, -0.42])

    uu = np.asarray([-1.7,1.29,.88, 1.46, 3.22,2.77])
    vv = np.asarray([-4.74, -2.72, -3.0,-.69,-2.34, -3.25])
    
    uu = np.asarray([-1.31,-.69, -1.41, .71, 1.20, 1.74])
    vv = np.asarray([-3.04, -2.47, -3.10, -2.99, -1.43, -1.32])

    uu = np.asarray([-1.52, .04,-1.55, 0.47, 0.93, 0.40]) #p2
    vv = np.asarray([-2.6, -.66, -2.66, -2.93, -3.53, -3.16])


    uu = np.asarray([-1.62, -0.73, -1.53, 0.74, 1.05, 0.91]) #p 30525
    vv = np.asarray([-3.25, -1.99, 2.8, -2.8, -2.27, -2.53])
    
    uu = np.asarray([-1.66, -0.74, -1.57, 0.72, .99, 0.72]) #p 34273
    vv = np.asarray([-3.38, -1.9, -2.72, -2.77, -2.44, -2.88])
    
    uu = np.asarray([-1.50, -0.04, -1.540, 0.47, .93, .39]) #p2wrf zero time
    vv = np.asarray([-2.58, -.77, -2.65, -2.92, -3.53, -3.15])

    uu = np.asarray([-1.29, -0.19, -.0, 0.50, 1.37, -.58]) #p2wrf
    vv = np.asarray([-2.18, -.8, -2.67, -2.92, -3.9, -2.96])
    
    uu = np.asarray([-1.29, -0.19, -.01, 0.49, 1.27, -.55]) #p2wrf convergence
    vv = np.asarray([-2.18, -.8, -2.67, -2.92, -3.9, -2.96])
    
    uu = np.asarray([-1.45, -1.23, -.01, 0.49, 1.27, -.61]) #p2wrf bisli
    vv = np.asarray([-2.46, -3.31, -2.67, -2.91, -3.9, -3.45])

    uu = np.asarray([-1.83, -1.08, -.04, 0.73, 1.32, -.66]) #p2wrf bisli 3m
    vv = np.asarray([-4.13, -3.62, -3.02, -3.55, -4.08, -3.95])

    uu = np.asarray([-1.45, -1.19, -.19, 0.36, 1.17, -.53]) #p2wrf bisli 5m
    vv = np.asarray([-2.46, -2.61, -1.94, -2.28, -3.56, -2.98])
    
    uu = np.asarray([-1.55, -1.19, -.19, 0.36, 0.15, -.43]) #p2wrf bisli 3m
    vv = np.asarray([-2.59, -2.61, -1.94, -2.28, -1.68, -2.38])
        
    uu = np.asarray([-1.57, 1.15, .63, 3.22, 3.18, 3.92]) #p4wrf bisli 3m
    vv = np.asarray([-2.64, -1.85, -1.71, -1.68, -2.1, -2.91])

    uu = np.asarray([0.86, 1.16, .90, 1.03, 2.07, 2.13]) #p5wrf WRF
    vv = np.asarray([-2.42, -1.81, -1.67, -1.83, -2.31, -2.25])

    uu = np.asarray([-.04  ,  -.35,  -.62,  -.15,  .92, 2.51]) #p5wrf of 10m
    vv = np.asarray([-1.95, -2.33, -1.31, -2.34, -2., -2.05])
    
    uu = np.asarray([-1.25  ,  1.36,  1.16,  4.34, 4.45, 6.30]) #p5wrf of 10m
    vv = np.asarray([-2.54, -1.24, -.41, -.49, -1.43, -3.23])
    
    uu = np.asarray([-1.22  ,  1.14,  1.05,  3.69, 3.03, 3.51]) #p5wrf of 3m
    vv = np.asarray([-1.52, -1.22, -.39, -.04, -1.09, -2.73])
    
    uu = np.asarray([0.23, 0.21, 0.05, 0.80, 0.71, 1.45]) #p5all2 of 10m
    vv = np.asarray([-1.19, -1.29, -0.31, -1.12, -1.61, -1.38])
    
    uu = np.asarray([0.42, 0.35, 0.27, 1.03, 0.94, 1.78]) #p5all2 of 20m
    vv = np.asarray([-1.69, -1.45, -0.75, -1.41, -2.16, -2.30])
    
    uu = np.asarray([0.73, 0.53, 0.18, 0.96, 1.57, 1.04]) #p5all3 of 10m
    vv = np.asarray([-1.73, -1.36, -0.70, -1.41, -1.50, -1.29])
    
    uu = np.asarray([0.53, 0.38, 0.25, 0.84, 0.71, 0.91]) #p5all3 of 10m
    vv = np.asarray([-1.51, -1.14, -0.47, -0.81, -1.30, -1.47])
    
    uu = np.asarray([0.52, 0.51, 0.40, 0.98, 0.82, 1.37]) #p5all3 of 20m
    vv = np.asarray([-1.92, -1.33, -0.90, -1.04, -1.65, -2.11])

    uu = np.asarray([0.56 ,  0.57,  0.45,  1.02,  0.87,  1.55]) #p5all3 of 30m
    vv = np.asarray([-2.20, -1.42, -1.01, -1.12, -1.81, -2.37])

    uu = np.asarray([0.53, 0.38, 0.26, .84, .70, .91]) #p5all3 interpolate 10m
    vv = np.asarray([-1.51, -1.15, -0.48, -.9, -1.27, -1.47])
    
    uu = np.asarray([0.17, 0.59, 0.14, 0.54, 0.92, .55]) #comfort0 (10m)
    vv = np.asarray([-.68, -.96 ,-.29, -.78, -0.95, -.56])
    
    uu = np.asarray([ 0.71,  1.66,  0.76,  1.54,  2.33,  2.22]) #comfort0 (10m new)
    vv = np.asarray([-2.83, -2.76 ,-1.53, -2.26, -2.44, -2.34])
    
    uu = np.asarray([ 0.62,  1.55,  0.46,  1.99,  2.66,  4.55]) #comfort0 (10m new)
    vv = np.asarray([-2.03, -2.34 ,-1.36, -1.81, -2.47, -2.22])
    
    uu = np.asarray([-0.34,  1.22,  1.15,  2.39,  2.57,  1.47]) #comfort0simple (10m)
    vv = np.asarray([-2.20, -1.41, -0.85, -1.03, -1.46, -1.69])
    
    uu = np.asarray([ 0.93,  1.92,  1.47,  1.60,  2.35,  2.13]) #comfortpressure (10m new)
    vv = np.asarray([-3.42, -3.18 ,-2.32, -2.97, -2.47, -2.17])
    vv = np.asarray([-1.51, -3.18 ,-0.39, -2.97, -2.47, -2.17])
    
    uu = np.asarray([ 0.70,  1.63 , 0.91,  1.34,  2.39,  2.16]) #comfortpressure (10m newnew)
    vv = np.asarray([-3.42, -2.84 ,-2.32, -2.97, -2.47, -2.17])
    vv = np.asarray([-1.30, -2.84 ,-0.10, -2.61, -2.49, -1.72])
    
    uu = np.asarray([-0.21,  0.31 ,-0.33, -1.14,  2.36,  0.97]) #comfortpressure4 (10m newnew)
    vv = np.asarray([-2.50, -2.35 ,-1.87,  0.30, -2.64, -2.31])

    uu = np.asarray([-0.01,  0.02 ,-0.69, -1.07,  0.39,  0.49]) #comfortpressure4 (10m newnew)
    vv = np.asarray([-2.25, -2.47 ,-1.57,  0.32, -1.86, -2.13])

    uu = np.asarray([ 0.06, -0.25 ,-0.46, -1.79,  0.14,  0.48]) #comfortpressure4 (7870)
    vv = np.asarray([-2.44, -2.57 ,-1.45, -0.67, -1.49, -2.12])
    uu = np.asarray([ 0.04, -0.24 ,-0.46, -1.77,  0.13,  0.50]) #comfortpressure4 (7899)
    vv = np.asarray([-2.43, -2.56 ,-1.46, -0.67, -1.48, -2.12])
    uu = np.asarray([ 0.02, -0.23 ,-0.45, -1.75,  0.11,  0.52]) #comfortpressure4 (7944)
    vv = np.asarray([-2.42, -2.54 ,-1.46, -0.68, -1.47, -2.12])
    uu = np.asarray([-0.11, -0.23 ,-0.48, -1.64, -0.03,  0.47]) #comfortpressure3 (8305)
    vv = np.asarray([-2.42, -2.44 ,-1.52, -0.76, -1.35, -2.12])
    uu = np.asarray([ 2.02,  2.43 , 2.34,  0.56,  2.69,  2.78]) #comfortpressure3 (1490m top=outlet)
    vv = np.asarray([-3.30, -2.60 ,-1.62, -0.58,  2.32, -2.42])
    uu = np.asarray([ 1.18,  1.68 , 2.17,  0.17,  0.08,  0.37]) #comfortpressure3 (2490m top=outlet)
    vv = np.asarray([-3.44, -3.88 ,-2.92, -3.48,  3.42, -3.11])
    uu = np.asarray([ 0.72,  1.00 , 0.81,  0.47,  2.35,  1.20]) #windarroundbuildingwrf 40000
    vv = np.asarray([-2.62, -1.95 ,-1.94, -1.20, -2.10, -2.14])
    uu = np.asarray([ 0.74,  0.82 , 0.75,  0.46,  2.35,  1.16]) #windarroundbuildingwrf 50000
    vv = np.asarray([-2.62, -1.85 ,-1.85, -1.08, -2.17, -2.17])
    uu = np.asarray([ 0.47,  0.56 , 0.99,  0.02,  1.91,  1.04]) #windarroundbuildingwrf 684000
    vv = np.asarray([-1.73, -0.86 ,-1.17, -1.15, -0.90, -1.39])
    uu = np.asarray([ 0.28 ,  0.99,  1.24,  1.47,  2.08,  2.24]) #windarroundbuildingwrf2 154000
    vv = np.asarray([-0.81 , -1.08, -0.98, -1.03, -1.49, -2.32])
    uu = np.asarray([ 0.52 ,  0.81,  1.11,  1.55,  1.82,  1.57]) #windarroundbuildingwrf2 471000
    vv = np.asarray([-1.41 , -0.74, -0.60, -0.42, -0.74, -1.75])
    uu = np.asarray([ 0.46 ,  0.78,  1.11,  1.50,  1.87,  1.42]) #windarroundbuildingwrf2 581000
    vv = np.asarray([-1.33 , -0.68, -0.54, -0.43, -0.75, -1.59])
    uu = np.asarray([ 0.16 , -0.19,  0.27, -0.57,  0.10, -0.762]) #windarroundbuildingwrf2 new 113000 (0.03 fvsolusion)
    vv = np.asarray([-1.07 , -1.42, -0.86, -1.61, -1.74, -3.18])
    uu = np.asarray([-0.02 , -0.05,  0.40, -0.66,  0.66,  1.20]) #windarroundbuildingwrf2 new 227000
    vv = np.asarray([-1.08 , -1.28, -0.98, -1.99, -1.31, -2.21])
    uu = np.asarray([-0.02 , -0.10,  0.40, -0.38,  0.93,  1.78]) #windarroundbuildingwrf2 new 256000 fin
    vv = np.asarray([-1.07 , -1.26, -0.94, -1.97, -1.31, -2.09])
    uu = np.asarray([ 0.26 ,  0.71,  0.94,  0.59,  0.75,  0.51]) #windarroundbuildingwrf2 new 113000 0.03 + temperature
    vv = np.asarray([-0.89 , -1.00, -0.95, -1.21, -1.58, -2.69])
    uu = np.asarray([ 0.18 ,  0.65,  0.92,  0.63,  0.84,  0.65]) #windarroundbuildingwrf2 new 113000 0.03 + temperature
    vv = np.asarray([-0.82 , -0.80, -0.81, -1.29, -1.50, -2.36])

    uu = np.asarray([-0.01 ,  0.43,  0.65, -0.06,  1.49,  1.04]) #windarroundbuildingwrf2 154000b
    vv = np.asarray([-1.21 , -2.20, -1.29, -1.08, -2.29, -3.15])
    uu = np.asarray([-0.20 ,  1.29, -0.54, -1.08,  0.20,  0.26]) #windarroundbuildingwrf2 260000b
    vv = np.asarray([-1.35 , -2.45, -0.77, -2.36, -2.06, -2.99])
    uu = np.asarray([ 0.46 ,  1.15,  0.21,  1.13,  1.88,  2.08]) #windarroundbuildingwrf3 40000
    vv = np.asarray([-2.35 , -2.26, -0.70, -1.64, -2.38, -1.72])
    uu = np.asarray([ 0.51 ,  0.92, -0.06,  1.01,  1.81,  1.89]) #windarroundbuildingwrf3 51000
    vv = np.asarray([-2.22 , -2.19, -0.42, -1.55, -2.35, -1.69])

    
    uu = np.asarray([-1.22  ,  1.14,  1.05,  3.69, 3.03, 3.51]) #wrf probe
    vv = np.asarray([-1.52, -1.22, -.39, -.04, -1.09, -2.73])
    
#################################################
    uu = np.asarray([1.08  ,  2.07,  1.12,  2.04,  2.09,  2.14]) #wrf2 202306131800 probe
    vv = np.asarray([-3.1  , -2.70, -2.41, -2.56, -2.72, -2.25])
    uu = np.asarray([ 0.47,  0.56 , 0.99,  0.02,  1.91,  1.04]) #windarroundbuildingwrf 684000
    vv = np.asarray([-1.73, -0.86 ,-1.17, -1.15, -0.90, -1.39])
    uu = np.asarray([ 0.46 ,  0.78,  1.11,  1.50,  1.87,  1.42]) #windarroundbuildingwrf2 581000
    vv = np.asarray([-1.33 , -0.68, -0.54, -0.43, -0.75, -1.59])
    uu = np.asarray([ 0.51 ,  0.92, -0.06,  1.01,  1.81,  1.89]) #windarroundbuildingwrf3 51000
    vv = np.asarray([-2.22 , -2.19, -0.42, -1.55, -2.35, -1.69])




    uu = np.asarray([ 2.50,  2.69,  2.80,  2.16,  1.74,  1.81]) #wrf 202306121800 probe *************
    vv = np.asarray([-0.88, -0.98, -1.03, -0.9 , -0.83, -0.80])
    uu = np.asarray([ 2.47,  2.65,  2.76,  2.21,  1.74,  1.81]) #wrf 202306121800 probe b2 *************
    vv = np.asarray([-0.87, -0.97, -1.02, -0.91, -0.83, -0.80])
    uu = np.asarray([ 1.01,  1.66,  1.61,  2.14,  0.95,  0.24]) #wrf 202306121800 c2s 227k
    vv = np.asarray([-0.11,  0.45, -0.06,  1.19, -1.11, -0.04])
    uu = np.asarray([ 1.37,  2.04,  1.89,  3.08,  1.29,  1.28]) #wrf 202306121800 c2s 227k + height
    vv = np.asarray([-0.27,  0.69, -0.67,  1.23, -1.42, -0.16])
    uu = np.asarray([ 2.61,  4.99,  4.95,  2.76,  2.30,  3.07]) #wrf 202306121800 a2az0 500k
    vv = np.asarray([ 0.16, -0.03,  0.07,  0.84, -2.24, -0.86])
    uu = np.asarray([ 1.28,  1.52,  2.36,  0.66,  0.61,  1.61]) #wrf 202306121800 a2aimpleb1 500k
    vv = np.asarray([-0.27, -0.54,  0.18, -0.54, -0.43, -0.90])
    uu = np.asarray([ 1.31,  4.64,  3.88,  4.87,  4.19,  4.75]) #wrf 202306121800 b2az1 500k
    vv = np.asarray([ 0.79, -0.11, -0.69,  0.69, -0.45, -0.50])


    uu = np.asarray([ 2.61 ,  5.19, 5.21, 2.96, 2.54, 3.33]) #wrf 202306121800 a2 43k fin
    vv = np.asarray([ 0.18 ,  0.04, 0.04, 1.33, -2.58, -0.95])
    uu = np.asarray([ 0.43 ,  4.45, 4.67, 4.62, 2.05, 2.57]) #wrf 202306121800 b2 533k
    vv = np.asarray([-0.52 , -0.13,-0.87, 0.69, -1.45, -1.25])
    uu = np.asarray([ 0.93 ,  5.73, 6.38, 6.07, 2.99, 3.51]) #wrf 202306121800 b2 941k
    vv = np.asarray([-0.15 ,  0.25,-0.62, 1.78, -1.29, -1.10])
    uu = np.asarray([ 1.82 ,  2.62, 2.59, 2.71, 1.70, 1.86]) #wrf 202306121800 c2 82k
    vv = np.asarray([-0.43 ,  0.07,-0.12,-0.19, -0.80, -0.51])
    uu = np.asarray([ 0.79 ,  1.97, 1.93, 3.61, 1.04, 1.00]) #wrf 202306121800 b2az0s 291k
    vv = np.asarray([ 0.12 ,  0.21,-0.65, 1.25, -1.50, -0.28])
    
    
    uu = np.asarray([ 1.64,  2.03,  1.36,  1.80,  1.83,  1.62]) #wrf 202306122200 probe
    vv = np.asarray([-0.66, -1.24, -0.70, -0.67, -0.99, -0.59])
    uu = np.asarray([ 1.98,  2.96,  2.25,  1.42,  1.18,  2.02]) #wrf 202306122200 1932k
    vv = np.asarray([-1.01, -0.91, -1.19, -0.36, -1.06, -1.18])
    uu = np.asarray([ 2.26,  4.42,  3.49,  3.17,  2.60,  2.65]) #wrf 202306122200b 878k
    vv = np.asarray([-0.51,  0.60,  1.13,  0.91, -2.52, -1.61])
  
    uu = np.asarray([ 3.43,  4.93,  4.31,  4.45,  2.90,  3.51]) #wrf 202306131200 probe a1z0*************
    vv = np.asarray([ 0.95,  1.96,  1.58,  1.94,  0.75, -2.70])
    uu = np.asarray([ 3.40,  4.83,  4.31,  4.42,  3.28,  3.95]) #wrf 202306131200 probe from b2z0s
    vv = np.asarray([ 0.91,  1.94,  1.76,  1.86,  0.42, -2.97])
    uu = np.asarray([ 1.29,  3.46,  3.67, -0.37,  0.88,  0.91]) #wrf 202306131200 b2az0simple 1000K
    vv = np.asarray([ 0.95,  0.14,  0.19, -0.12, -0.80, -1.51])
    uu = np.asarray([ 4.17, 10.04,  5.37,  8.10,  5.89,  6.74]) #wrf 202306131200 b2az0 900K
    vv = np.asarray([ 6.17,  8.12,  3.72,  9.90,  5.26,  3.77])
    uu = np.asarray([ 3.31,  6.97,  6.49,  3.64,  1.31,  4.80]) #wrf 202306131200 a4z0 800K
    vv = np.asarray([ 5.27,  5.41,  5.40,  4.67, -1.20, 3.65])
    uu = np.asarray([ 2.43,  4.41,  3.54,  0.61,  0.45, 4.64]) #wrf 202306131200 a2z0s 590K
    vv = np.asarray([ 1.30,  1.70,  0.96,  0.53, -0.37,-0.93])
    uu = np.asarray([ 0.73,  3.29,  2.52,  2.16,  2.08, 5.41]) #wrf 202306131200 b2z0s 174K
    vv = np.asarray([ 0.56,  1.26,  0.50,  1.68, -0.62,-0.03])
    uu = np.asarray([ 1.21,  3.61,  2.15,  1.43,  2.51, 4.65]) #wrf 202306131200 c2z0s 36K
    vv = np.asarray([ 0.62,  1.91,  0.47,  1.28, -0.51,-0.29])
    uu = np.asarray([ 1.77,  4.29,  2.80,  2.20,  2.87, 5.22]) #wrf 202306131200 c2z0s 36K + height
    vv = np.asarray([ 0.79,  2.56,  0.24,  1.60, -0.63,-0.22])


    
    
    
    # uu = np.asarray([]) #p5wrf probe10
    # vv = np.asarray([])
    
    # ws0=np.asarray([1.5,1,1.9,0.2,2.2,1.8]) #day1
    # wd0=np.asarray([58,35,73,258,289,283])
    
    import sklearn.metrics
    import math
    import numpy as np    
    
    # https://envihaifa.net/

    ws0=np.asarray([2.,3.9,3.2,3.8,5.5,4.4]) # 13/06/2023 12GMT
    wd0=np.asarray([243,226,256,232,240,313])
    uu0=np.round(-ws0*np.sin(wd0/180*math.pi),2) #[ 1.78, 2.81, 3.10, 2.99, 4.76, 3.22]
    vv0=np.round(-ws0*np.cos(wd0/180*math.pi),2) #[ 0.91, 2.71, 0.77, 2.34, 2.75,-3.00  ]



    ws0=np.asarray([1.5,2.,1.8,1.9,3.8,2.9]) # 13/06/2023 18GMT
    wd0=np.asarray([334,286,262,294,308,305])
    uu0=np.round(-ws0*np.sin(wd0/180*math.pi),2) #[ 0.66  1.92  1.78  1.74  2.99  2.38]
    vv0=np.round(-ws0*np.cos(wd0/180*math.pi),2) #[-1.35 -0.55  0.25 -0.77 -2.34 -1.66]



    ws0=np.asarray([1.8,2.1,2.6,2.2,3.9,3.6]) # 12/06/2023 18GMT
    wd0=np.asarray([286,267,270,296,294,292])
    ws0=np.asarray([1.8,1.1,1.8,2.6,0.2,2.9]) # 12/06/2023 18GMT new reading
    wd0=np.asarray([308,297,299,335,298,288])
    uu0=np.round(-ws0*np.sin(wd0/180*math.pi),2) #[ 1.73, 2.1 , 2.6 , 1.98, 3.56, 3.34]
    vv0=np.round(-ws0*np.cos(wd0/180*math.pi),2) #[-0.5 , 0.11, 0.  ,-0.96,-1.59,-1.35]
    
    ws0=np.asarray([0.7,1.5,1.3,1.4,3.2,2.4]) # 12/06/2023 22GMT
    wd0=np.asarray([299,282,292,304,307,296])
    uu0=np.round(-ws0*np.sin(wd0/180*math.pi),2) #[ 0.61,  1.47,  1.21,  1.16,  2.56,  2.16]
    vv0=np.round(-ws0*np.cos(wd0/180*math.pi),2) #[-0.34, -0.31, -0.49, -0.78, -1.93, -1.05]
    
    ws0=np.asarray([1.7,4.4,1.8,3.4,1.3,2.0]) # test
    wd0=np.asarray([271,252,305,238,279,325])
    uu0=np.round(-ws0*np.sin(wd0/180*math.pi),2) #[1.7 , 4.18, 1.47, 2.88, 1.28, 1.15]
    vv0=np.round(-ws0*np.cos(wd0/180*math.pi),2) #[-0.03,  1.36, -1.03,  1.8 , -0.2 , -1.64]
    
    # universiry 202240/740730 (475m) 12: 1.7/271 13: 4.4/252 ==> 
    # technion   202390/742120 (245m) 12: 1.8/305 13: 3.4/238 ==> 
    # Byalik     207632/746631 (6m) 12: 1.3/279 13: 2.0/325 ==> 

    ws01218=np.asarray([1.7,1.8,1.3,1.8,1.1,1.8,2.6,0.2,2.9]) # 12/06/2023 18GMT new reading
    wd01218=np.asarray([271,305,279,308,297,299,335,298,288])
    uu01218=np.round(-ws01218*np.sin(wd01218/180*math.pi),2) #[ 1.73, 2.1 , 2.6 , 1.98, 3.56, 3.34]
    vv01218=np.round(-ws01218*np.cos(wd01218/180*math.pi),2) #[-0.5 , 0.11, 0.  ,-0.96,-1.59,-1.35]

    uu =    np.asarray([ 3.60, 2.44, 1.73 , 2.50, 2.61, 2.76, 2.21, 1.74, 1.81]) #wrf 202306121800 c2z0s2 probe
    vv =    np.asarray([-1.28,-1.00,-0.85 ,-0.88,-0.96,-1.02,-0.91,-0.83,-0.80])

    uu =    np.asarray([ 2.54, 1.29, 0.12, 0.86,  1.80,  1.40,  2.14,  1.59, 0.53]) #wrf 202306121800 c2z0s2 4k
    vv =    np.asarray([ 0.73,-0.17, 0.02,-0.11,  0.50, -0.06,  1.31, -1.67, 0.07])
    uu = uu+np.asarray([ 3.11, 1.57, 0.54, 1.31,  2.14,  1.75,  3.11,  1.88, 1.64]) #wrf 202306121800 c2z0s2 4k + height
    vv = vv+np.asarray([ 0.46,-0.18,-0.26,-0.24,  0.76, -0.59,  1.36, -1.85,-0.20])    
    uu=uu/2
    vv=vv/2

    
    ws01312=np.asarray([4.4,3.4,2.0,2.0,3.9,3.2,3.8,5.5,4.4]) # 13/06/2023 12GMT
    wd01312=np.asarray([252,238,325,243,226,256,232,240,313])
    uu01312=np.round(-ws01312*np.sin(wd01312/180*math.pi),2) #[ 1.78, 2.81, 3.10, 2.99, 4.76, 3.22]
    vv01312=np.round(-ws01312*np.cos(wd01312/180*math.pi),2) #[ 0.91, 2.71, 0.77, 2.34, 2.75,-3.00  ]

    ws01316=np.asarray([1.7,1.0,1.4,2.0,3.9,3.2,3.8,5.5,4.4]) # 13/06/2023 16GMT
    wd01316=np.asarray([242,276,297,243,226,256,232,240,313])
    uu01316=np.round(-ws01312*np.sin(wd01312/180*math.pi),2) #[ 1.78, 2.81, 3.10, 2.99, 4.76, 3.22]
    vv01316=np.round(-ws01312*np.cos(wd01312/180*math.pi),2) #[ 0.91, 2.71, 0.77, 2.34, 2.75,-3.00  ]

    uu =    np.asarray([5.71,4.83, 3.03 ,3.40,  4.80,  4.31,  4.36,  3.08, 3.95]) #wrf 202306131200 c2z0s2 probe
    vv =    np.asarray([2.77,2.42,-2.89 ,0.90,  1.85,  1.76,  1.78,  0.38,-2.97])
    uu = uu+np.asarray([5.65,4.94, 3.31 ,3.49,  4.83,  4.35,  4.42,  3.36, 4.19]) #wrf 202306131200 c2z0s2 probe + height
    vv = vv+np.asarray([2.91,2.51,-3.12 ,0.95,  1.94,  1.80,  1.86,  0.44,-3.15])
    uu=uu/2
    vv=vv/2
    
    uu =    np.asarray([ 4.41,2.63, 1.22,1.21,  3.64,  2.18,  1.44,  2.29, 5.64]) #wrf 202306131200 c2z0s2 68k
    vv =    np.asarray([ 1.79,1.09, 0.90,0.66,  1.92,  0.43,  1.32, -0.42,-0.38])
    uu = uu+np.asarray([ 4.95,3.22, 4.24,1.84,  4.30,  2.89,  2.23,  2.66, 6.13]) #wrf 202306131200 c2z0s2 68k + height
    vv = vv+np.asarray([ 2.63,1.26,-0.61,0.84,  2.57,  0.26,  1.64, -0.60,-0.26])    
    uu=uu/2
    vv=vv/2

    
    # uu0 = np.asarray([ 0.65,  1.92, 1.78,  1.73,  2.99,  2.37])
    # vv0 = np.asarray([-1.34, -0.55, 0.25, -0.77, -2.33, -1.66])
    
    print('i,U,   V,    WS,   WD   << measurements')
    for i in range(6):
        print(i,round(uu0[i],2),round(vv0[i],2),ws0[i],wd0[i])
        print(i,round(uu[i],2),round(vv[i],2),ws[i],wd[i])
    
    print('i,U,   V,    WS,   WD  << SIMULATION')


    ws=np.zeros_like(uu)
    wd=np.zeros_like(uu)
    for i in range(len(uu)):
        ws[i]=round((uu[i]**2.+vv[i]**2.)**.5,2)
        wd[i]=round(math.atan2(uu[i], vv[i])* 180/math.pi+180,2)
        # wd[i]=round(math.atan(uu[i]/vv[i])*180/math.pi,2)
        # wd[i]=round(math.atan2(uu[i], vv[i]))
        # print(i,round(math.atan2(uu0[i], vv0[i])* 180/math.pi+180,2))
        # if wd[i]<90:
            # wd[i]+=360
        # wd[i] =math.atan2(vv[i], uu[i])*180/math.pi
        # wd[i] =math.atan(vv[i]/uu[i])*180/math.pi
        # if wd[i]<90:
        #     wd[i]+=180
        print(i,round(uu[i],2),round(vv[i],2),round(ws[i],2), round(wd[i],2))

        prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/aerofoil7/postProcessing/probes/0/", nstation=9, ws0=ws01218, wd0=wd01218, uu0=uu01218, vv0=vv01218, line=-1)
        prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2rhosimple/postProcessing/probes/0/", nstation=9, ws0=ws01218, wd0=wd01218, uu0=uu01218, vv0=vv01218, line=-1)


        prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800c2az1s/postProcessing/probes/0/")
        prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2az0s/postProcessing/probes/0/")
        prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simpleb2/postProcessing/probes/0/")
        prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2az0s2/postProcessing/probes/0/")

        prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2az0s2/postProcessing/probes/0/", nstation=9, ws0=ws0, wd0=wd0, uu0=uu0, vv0=vv0)
        prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2az0s2/postProcessing/probes/0/", nstation=6, ws0=ws0, wd0=wd0, uu0=uu0, vv0=vv0)
        prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simple/postProcessing/probes/0/", nstation=60, ws0=ws0, wd0=wd0, uu0=uu0, vv0=vv0)
        prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simpleb2/postProcessing/probes/0/", nstation=60, ws0=ws0, wd0=wd0, uu0=uu0, vv0=vv0)
        prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2az0s3/postProcessing/probes/0/")

        prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200b2az0rhosimple/postProcessing/probes/0/", nstation=9, ws0=ws01312, wd0=wd01312, uu0=uu01312, vv0=vv01312)
        prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200b2az0simple/postProcessing/probes/0/"   , nstation=9, ws0=ws01312, wd0=wd01312, uu0=uu01312, vv0=vv01312)
        
        prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2rhosimple/postProcessing/probes/0/",  nstation=9, ws0=ws01218, wd0=wd01218, uu0=uu01218, vv0=vv01218)
        prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simple/postProcessing/probes/0/", nstation=9, ws0=ws01218, wd0=wd01218, uu0=uu01218, vv0=vv01218)
        prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200c2rhosimple/postProcessing/probes/0/"   , nstation=9, ws0=ws01312, wd0=wd01312, uu0=uu01312, vv0=vv01312)
        prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200c2simple/postProcessing/probes/0/"   , nstation=9, ws0=ws01312, wd0=wd01312, uu0=uu01312, vv0=vv01312)

        
    def prob(file, percent=0.5, verbose=False, ws0=None, wd0=None, line=-1, el=None):
        print(file)
        # test = getProbe(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800c2az1s/postProcessing/probes/0/")
        ws0=np.asarray(ws0)
        wd0=np.asarray(wd0)
        uu0=np.round(-ws0*np.sin(wd0/180*math.pi),2) 
        vv0=np.round(-ws0*np.cos(wd0/180*math.pi),2)

        test = getProbe(file)
        print (datetime.datetime.now(), 'length is ', len(test[1]))
                
        elbkup = el.copy()
        ul=np.zeros(len(elbkup))
        vl=np.zeros(len(elbkup))
        ulh=np.zeros(len(elbkup))
        vlh=np.zeros(len(elbkup))
        # line=-1  # 0, -1
        testvar = [x for x in test[1][line] if x!='']  # remove '' in the list
        for i in range(len(elbkup)):
            while float(testvar[elbkup[i]])<-99999:
                elbkup[i]+=3
            pos = elbkup[i]
            if i==len(el)-1:
               lastpos = len(el) - 1
            else:
               lastpos = elbkup[i+1]
            while float(testvar[pos])==float(testvar[elbkup[i]]) and (pos < lastpos-3):
                pos+=3
                # if pos>len(test[1][line])-2:
                #     pos-=3
                #     break
            ul[i]=float(testvar[elbkup[i]])
            vl[i]=float(testvar[elbkup[i]+1])
            ulh[i]=float(testvar[pos])
            vlh[i]=float(testvar[pos+1])       
        uu=ul*percent+ulh*(1.-percent)
        vv=vl*percent+vlh*(1.-percent)

        ws=np.round((uu**2.+vv**2.)**.5,2)
        np.round(np.arctan2(uu, vv),2)
        wd=np.round(np.arctan2(uu, vv)* 180/math.pi+180,2)
        
        wd[wd>360]=wd[wd>360]-360
        wd[wd-wd0>180] = wd[wd-wd0>180] - 360
        wd[wd0-wd>180] = wd[wd0-wd>180] + 360
        
        print('u:',stat(uu0, uu, kind='r2'),stat(uu0, uu, kind='r'), stat(uu0, uu, kind='rmse'),'/',round(np.mean(uu),4),'+-',round(np.std(uu),4),'>>',round(np.mean(uu0),4),'+-',round(np.std(uu0),4))
        print('v:',stat(vv0, vv, kind='r2'),stat(vv0, vv, kind='r'), stat(vv0, vv, kind='rmse'),'/',round(np.mean(vv),4),'+-',round(np.std(vv),4),'>>',round(np.mean(vv0),4),'+-',round(np.std(vv0),4))
        print('s:',stat(ws0, ws, kind='r2'),stat(ws0, ws, kind='r'), stat(ws0, ws, kind='mae'), stat(ws0, ws, kind='rmse'),'/',round(np.mean(ws),4),'+-',round(np.std(ws),4),'>>',round(np.mean(ws0),4),'+-',round(np.std(ws0),4))
        print('d:',stat(wd0, wd, kind='r2', wind=True),stat(wd0, wd, kind='r', wind=True), stat(wd0, wd, kind='mae', wind=True), stat(wd0, wd, kind='rmse', wind=True),'/',round(np.mean(wd),4),'+-',round(np.std(wd),4),'>>',round(np.mean(wd0),4),'+-',round(np.std(wd0),4))
        if verbose:
            print('u (obs, model)')
            print (uu0)
            print (np.round(uu,2))
            print('v')
            print (vv0)
            print (np.round(vv,2))
            print('ws')
            print (ws0)
            print (np.round(ws,1))
            print('wd')
            print (wd0)
            print (wd.astype(int))
            
            
        distancerx = 0
        distancery = 0
        adistancerx = 0
        adistancery = 0
        maxr = 0
        for i in range(len(wd0)):
            distancerx += uu0[i]-uu[i]
            distancery += vv0[i]-vv[i]
            adistancerx += math.fabs(uu0[i]-uu[i])
            adistancery += math.fabs(vv0[i]-vv[i])
            maxr = max(((uu0[i]-uu[i])**2+(vv0[i]-vv[i])**2)**.5, maxr)
        print('distance (sum, abs, max) = ', (distancerx**2+distancery**2)**.5/len(wd0), (adistancerx**2+adistancery**2)**.5/len(wd0), maxr)

        return 
    
    def probdata(file, percent=0.5, ws0=None, wd0=None, line=-1, el=None):
        # test = getProbe(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800c2az1s/postProcessing/probes/0/")
        ws0=np.asarray(ws0)
        wd0=np.asarray(wd0)
        uu0=np.round(-ws0*np.sin(wd0/180*math.pi),2) 
        vv0=np.round(-ws0*np.cos(wd0/180*math.pi),2)

        test = getProbe(file)
        print (datetime.datetime.now(), 'length is ', len(test[1]))
                
        elbkup = el.copy()
        ul=np.zeros(len(elbkup))
        vl=np.zeros(len(elbkup))
        ulh=np.zeros(len(elbkup))
        vlh=np.zeros(len(elbkup))
        # line=-1  # 0, -1
        testvar = [x for x in test[1][line] if x!='']  # remove '' in the list
        for i in range(len(elbkup)):
            while float(testvar[elbkup[i]])<-99999:
                elbkup[i]+=3
            pos = elbkup[i]
            if i==len(el)-1:
               lastpos = len(el) - 1
            else:
               lastpos = elbkup[i+1]
            while float(testvar[pos])==float(testvar[elbkup[i]]) and (pos < lastpos-3):
                pos+=3
                # if pos>len(test[1][line])-2:
                #     pos-=3
                #     break
            ul[i]=float(testvar[elbkup[i]])
            vl[i]=float(testvar[elbkup[i]+1])
            ulh[i]=float(testvar[pos])
            vlh[i]=float(testvar[pos+1])       
        uu=ul*percent+ulh*(1.-percent)
        vv=vl*percent+vlh*(1.-percent)

        ws=np.round((uu**2.+vv**2.)**.5,2)
        wd=np.round(np.arctan2(uu, vv)* 180/math.pi+180,2)
        
        wd[wd>360]=wd[wd>360]-360
        wd[wd-wd0>180] = wd[wd-wd0>180] - 360
        wd[wd0-wd>180] = wd[wd0-wd>180] + 360
        

        return ws.tolist(),wd.tolist()
    
    print('u (obs, model)')
    print (uu0)
    print (np.round(uu,2))
    print('v')
    print (vv0)
    print (np.round(vv,2))
    print('ws')
    print (ws0)
    print (np.round(ws,1))
    print('wd')
    print (wd0)
    print (wd.astype(int))
    plt.figure()
    plt.plot(uu,uu0,'*b')
    plt.plot(vv,vv0,'*r')
    plt.plot(ws,ws0,'*g')
    xmin=min(uu0.min(),uu.min(),vv0.min(),vv.min())
    xmax=max(uu0.max(),uu.max(),vv0.max(),vv.max())
    xmax=max(ws0.max(),ws.max())
    ln11 = np.linspace(xmin,xmax)
    plt.plot(ln11,ln11,'k')
    plt.xlabel("simulation")
    plt.ylabel("observaion")
    plt.xlim(xmin,xmax)
    plt.ylim(xmin,xmax)
    plt.title('wrf')
    # plt.title('comfort -1.25, -.85, -.73, -1.91')

###start omri
    
    haifa1 = [  
    [202240, 740730, 475], # University
    [202390, 742120, 245], # technion
    [207630, 746612,  14], # Bialik
    [197075, 747605, 113],  #// kiryat shprinzak, ramot school france road 79    
    [199490, 745185, 300],  #//hugim, yair katz 4
    [198988, 743528, 300],  #//ahuze, smolskin st.
    [202315, 743652, 225],  #//hagalil 107, tel hai school
    [204115, 743860,  32],  #// check-post, moshlei yaakov 7
    [205603, 748278,  19]  #//dgania st. 
    ]
    stationss=[]
    stationsd=[]
    ws0=ws1218b.copy()
    wd0=wd1218b.copy()
    for i in range(len(haifa1)):
        stationss.append([haifa1[i][0], haifa1[i][1], ws0[i], haifa1[i][2]])
        stationsd.append([haifa1[i][0], haifa1[i][1], wd0[i], haifa1[i][2]])
    stationss = np.asarray(stationss)
    stationsd = np.asarray(stationsd)
    forecastedws=np.zeros_like(ws0)
    forecastedwd=np.zeros_like(ws0)
    for i in range(len(stationss)):
        loos=[] # leave one out
        lood=[] # leave one out
        for j in range(len(stationss)):
            if i!=j:
                loos.append(stationss[j])
                lood.append(stationsd[j])
        forecastedws[i] = interp(stationss[i,0], stationss[i,1], loos, elev = stationss[i,3])
        forecastedwd[i] = interp(stationsd[i,0], stationsd[i,1], lood, elev = stationsd[i,3])
        ## forecastidw= interp(stations[i,0], stations[i,1], loo)
    # print('d-',stat(wd0, forecastedw, kind='r2'),stat(wd0, forecastedw, kind='r'), stat(wd0, forecastedw, kind='mae'), stat(wd0, forecastedw, kind='rmse'),'/',round(np.mean(forecastedw),4),'+-',round(np.std(forecastedw),4),'>>',round(np.mean(wd0),4),'+-',round(np.std(wd0),4))
    print('s-',stat(ws0, forecastedws, kind='r2'),stat(ws0, forecastedws, kind='r'), stat(ws0, forecastedws, kind='mae'), stat(ws0, forecastedws, kind='rmse'),'/',round(np.mean(forecastedws),4),'+-',round(np.std(forecastedws),4),'>>',round(np.mean(ws0),4),'+-',round(np.std(ws0),4))
    print('d-',stat(wd0, forecastedwd, kind='r2'),stat(wd0, forecastedwd, kind='r'), stat(wd0, forecastedwd, kind='mae'), stat(wd0, forecastedwd, kind='rmse'),'/',round(np.mean(forecastedwd),4),'+-',round(np.std(forecastedwd),4),'>>',round(np.mean(wd0),4),'+-',round(np.std(wd0),4))

    for p in range(3,10):
        allmaxr=0
        allbias=0
        allabs=0
        wsr=0
        wdr=0
        wsr2=0
        wdr2=0
        wsmae=0
        wdmae=0
        wsrmse=0
        wdrmse=0
        
    
        ts0=[]#ws0.copy()
        td0=[]#wd0.copy()
        # ts = forecastedws.copy()
        # td = forecastedwd.copy()
        loops=100000
    
        for j in range(loops):
            distancerx = 0
            distancery = 0
            adistancerx = 0
            adistancery = 0
            maxr = 0
    
    
            
            indices = np.random.choice(ws0.shape[0], p, replace=False)
            mask = np.ones(ws0.shape[0], dtype=bool)
            mask[indices] = False
            stationss=[]
            stationsd=[]
            for i in range(len(haifa1)):
                if ~mask[i]:       
                    stationss.append([haifa1[i][0], haifa1[i][1], ws0[i], haifa1[i][2]])
                    stationsd.append([haifa1[i][0], haifa1[i][1], wd0[i], haifa1[i][2]])
            stationss = np.asarray(stationss)
            stationsd = np.asarray(stationsd)
            forecastedws=np.zeros(len(stationss))
            forecastedwd=np.zeros(len(stationss))
            ts0=[]
            td0=[]
            for i in range(len(stationss)):
                loos=[] # leave one out
                lood=[] # leave one out
                for k in range(len(stationss)):
                    if i!=k:
                        loos.append(stationss[k])
                        lood.append(stationsd[k])
                ts0.append(stationss[i,2])
                td0.append(stationsd[i,2])
                forecastedws[i] = interp(stationss[i,0], stationss[i,1], loos, elev = stationss[i,3])
                forecastedwd[i] = interp(stationsd[i,0], stationsd[i,1], lood, elev = stationsd[i,3])
        
            ts = forecastedws.copy()
            td = forecastedwd.copy()
            ts0=np.asarray(ts0)
            td0=np.asarray(td0)
            
            for i in range(len(ts0)):
                distancerx += (ts0[i]*math.cos(td0[i]/180*math.pi)-ts[i]*math.cos(td[i]/180*math.pi))
                distancery += (ts0[i]*math.sin(td0[i]/180*math.pi)-ts[i]*math.sin(td[i]/180*math.pi))
                adistancerx += math.fabs(ts0[i]*math.cos(td0[i]/180*math.pi)-ts[i]*math.cos(td[i]/180*math.pi))
                adistancery += math.fabs(ts0[i]*math.sin(td0[i]/180*math.pi)-ts[i]*math.sin(td[i]/180*math.pi))
                maxr = max(((ts0[i]*math.cos(td0[i]/180*math.pi)-ts[i]*math.cos(td[i]/180*math.pi))**2 +
                             (ts0[i]*math.sin(td0[i]/180*math.pi)-ts[i]*math.sin(td[i]/180*math.pi))**2)**.5, maxr)
            allmaxr+=maxr
            allbias+=(distancerx**2+distancery**2)**.5
            allabs+=(adistancerx**2+adistancery**2)**.5
            wsr=stat(ts0, ts, kind='r')
            wdr+=stat(td0, td, kind='r')
            wsr2+=stat(ts0, ts, kind='r2')
            wdr2+=stat(td0, td, kind='r2')
            wsmae+=stat(ts0, ts, kind='mae')
            wdmae+=stat(td0, td, kind='mae')
            wsrmse+=stat(ts0, ts, kind='rmse')
            wdrmse+=stat(td0, td, kind='rmse')
            
        print(p,'distance (sum, abs, max) = ', round(allbias/loops/p,2), round(allabs/loops/p,2), round(allmaxr/loops/p,2))
        print(p,'s (r2 r maem rmse)-',round(wsr2/loops,2), round(wsr/loops,2), round(wsmae/loops,2), round(wsrmse/loops,2))
        print(p,'d (r2 r maem rmse)-',round(wdr2/loops,2), round(wdr/loops,2), round(wdmae/loops,2), round(wdrmse/loops,2))
###fin omri

    
    forecastedws[:]=ws0[-1]
    forecastedwd[:]=wd0[-1]
    print('s-',stat(ws0, forecastedws, kind='r2'),stat(ws0, forecastedws, kind='r'), stat(ws0, forecastedws, kind='mae'), stat(ws0, forecastedws, kind='rmse'),'/',round(np.mean(forecastedws),4),'+-',round(np.std(forecastedws),4),'>>',round(np.mean(ws0),4),'+-',round(np.std(ws0),4))
    print('d-',stat(wd0, forecastedwd, kind='r2'),stat(wd0, forecastedwd, kind='r'), stat(wd0, forecastedwd, kind='mae'), stat(wd0, forecastedwd, kind='rmse'),'/',round(np.mean(forecastedwd),4),'+-',round(np.std(forecastedwd),4),'>>',round(np.mean(wd0),4),'+-',round(np.std(wd0),4))
 
    
 
    if 5==7:
        import json
        import requests
        
        url =  "https://api.ims.gov.il/v1/envista/stations/43/data/daily/2023/6/13"
        
        url = "https://api.ims.gov.il/v1/Envista/stations"
        headers = {
            'Authorization': 'ApiToken f058958a-d8bd-47cc-95d7-7ecf98610e47'
                }
        response = requests.request("GET", url, headers=headers)
        data= json.loads(response.text.encode('utf8'))
        print (data)    
        for i in range(len(data)):
            print('>>>',data[i]['stationId'], data[i]['active'],data[i]['location']['latitude'],data[i]['location']['longitude'])
            if (data[i]['active'] is True):
                for j in range(len(data[i]['monitors'])):
                    # print(data[i]['monitors'][j]['name'])
                    if (data[i]['monitors'][j]['name']=='WD'):
                        print ('*WD')
                    if (data[i]['monitors'][j]['name']=='WS'):
                        print ('*WS')
                        url2 = "https://api.ims.gov.il/v1/envista/stations/"+str(data[i]['stationId'])+"/data/latest"
                        response2 = requests.request("GET", url2, headers=headers)
                        data2= json.loads(response2.text.encode('utf8'))
                        dt = data2['data'][0]['datetime']
                        for k in range(len(data2['data'][0]['channels'])):
                            print(data2['data'][0]['channels'][k]['name'],data2['data'][0]['channels'][k]['value'])
                        

ws1218a=np.asarray([1.8,0.8,1.2,1.5,1.0,1.9,0.4,2.5,2.9]) # uni, technion, bialik, shprinzak, galim, ahuza,nosh, checkpost, dgania
ws1218b=np.asarray([1.6,1.1,1.3,1.4,1.2,1.5,0.3,2.9,2.8]) # E-F
ws1218c=np.asarray([1.2,1.0,1.4,1.5,1.6,2.1,1.1,2.9,2.6])
ws1312a=np.asarray([4.9,4.1,2.4,1.9,3.7,2.8,2.8,6.0,3.8])
ws1312b=np.asarray([5.5,3.6,2.0,2.0,3.9,3.1,3.4,5.5,3.8]) # B-C
ws1312c=np.asarray([5.5,4.3,2.4,1.7,4.8,3.2,3.0,5.1,3.8])
ws1316a=np.asarray([2.6,1.5,1.6,1.1,1.0,1.6,1.1,3.2,1.6])
ws1316b=np.asarray([1.7,1.0,1.3,1.1,1.1,1.1,0.3,2.6,1.1]) # A-B
ws1316c=np.asarray([1.7,0.6,0.6,0.8,0.6,0.9,0.2,1.2,0.3])
wd1218a=np.asarray([263,275,277,300,305,273,304,299,287])
wd1218b=np.asarray([264,296,274,306,298,295,319,306,292])
wd1218c=np.asarray([267,309,254,304,296,299,332,321,300])
wd1312a=np.asarray([241,246,312,230,237,248,260,256,300])
wd1312b=np.asarray([237,225,326,235,225,250,237,243,312])
wd1312c=np.asarray([239,249,323,227,225,244,242,245,319])
wd1316a=np.asarray([269,299,245, 30, 14, 62,306,308,296])
wd1316b=np.asarray([242,276,228, 81, 46, 92,326,314,281])
wd1316c=np.asarray([173,122,182, 15,156, 69,349,286,226])
ws1316a=np.asarray([3.5,3.2,0.5,1.1,2.3,1.3,0.3,0.8,0.4])
ws1300b=np.asarray([3.1,3.1,0.6,0.9,2.5,1.3,0.7,0.5,0.3]) # E-F
ws1300c=np.asarray([3.0,2.1,0.5,1.0,2.5,1.2,1.0,0.7,0.3])
wd1300a=np.asarray([228,227, 43,157,206,172,244,226, 41])
wd1300b=np.asarray([238,227, 42,166,207,173,238,158,353])
wd1300c=np.asarray([231,220, 86,144,211,169,232,150,222])


el0 = [0,15,30, 45, 60, 75,90,105,120]  #dynamic
# el4 = [4,19,34, 49, 64, 79,94,109,124]  #dynamic foil pimple
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/aerofoil7/postProcessing/probes/0/", el=el0, ws0=ws1218b, wd0=wd1218b, line= 0)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/aerofoil7/postProcessing/probes/0/", el=el0, ws0=ws1218b, wd0=wd1218b, line=-1)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/aerofoil7/postProcessing/probes/0/", el=el0, ws0=(ws1218a+ws1218b+ws1218c)/3., wd0=(wd1218a+wd1218b+wd1218c)/3., line=1)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/aerofoil7/postProcessing/probes/0/", el=el0, ws0=(ws1218a+ws1218b+ws1218c)/3., wd0=(wd1218a+wd1218b+wd1218c)/3., line=-1)

# prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200b2az0rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1312b, wd0=wd1312b, line= 1, verbose=True)
# prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200b2az0rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1312b, wd0=wd1312b, line=-1, verbose=True)
# prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa//postProcessing/probes/0/", el=el4, ws0=(ws1218a+ws1218b+ws1218c)/3., wd0=(wd1218a+wd1218b+wd1218c)/3., line=1)
# prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa//postProcessing/probes/0/", el=el4, ws0=(ws1218a+ws1218b+ws1218c)/3., wd0=(wd1218a+wd1218b+wd1218c)/3., line=-1)

prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131600a2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1316b, wd0=wd1316b, line= 0)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131600a2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1316b, wd0=wd1316b, line=-1)
# prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131600a2rhosimplea/postProcessing/probes/0/", el=el0, ws0=ws1316b, wd0=wd1316b, line=-1)

# prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131600b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1316b, wd0=wd1316b, line= 1)
# prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131600b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1316b, wd0=wd1316b, line=-1)

prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200c2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1312b, wd0=wd1312b, line= 0)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200c2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1312b, wd0=wd1312b, line=-1)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200a2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1312b, wd0=wd1312b, line=-1)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200a2rhosimplenut/postProcessing/probes/0/", el=el0, ws0=ws1312b, wd0=wd1312b, line=-1)

prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200a2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1312b, wd0=wd1312b, line= 0)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200a2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1312b, wd0=wd1312b, line=-1)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200c2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1312b, wd0=wd1312b, line=-1)

prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306130000a2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1300b, wd0=wd1300b, line= 0)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306130000a2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1300b, wd0=wd1300b, line=-1)

prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1218b, wd0=wd1218b, line=0, verbose=True)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1218b, wd0=wd1218b, line=-1, verbose=True)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1218b, wd0=wd1218b, line=-1, verbose=True)

prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306130000b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1300b, wd0=wd1300b, line=0, verbose=True)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306130000a2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1300b, wd0=wd1300b, line=-1, verbose=True)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306130000b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1300b, wd0=wd1300b, line=-1, verbose=True)

prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1312b, wd0=wd1312b, line=0, verbose=True)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200a2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1312b, wd0=wd1312b, line=-1, verbose=True)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1312b, wd0=wd1312b, line=-1, verbose=True)

prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131600b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1316b, wd0=wd1316b, line=0, verbose=True)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131600a2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1316b, wd0=wd1316b, line=-1, verbose=True)
prob(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131600b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1316b, wd0=wd1316b, line=-1, verbose=True)

wsm = np.asarray(ws1218b.tolist() + ws1300b.tolist() + ws1312b.tolist() + ws1316b.tolist())
wdm = np.asarray(wd1218b.tolist() + wd1300b.tolist() + wd1312b.tolist() + wd1316b.tolist())

ws1 = np.asarray((np.ones(9)*ws1218b[-1]).tolist()+
                 (np.ones(9)*ws1300b[-1]).tolist()+
                 (np.ones(9)*ws1312b[-1]).tolist()+
                 (np.ones(9)*ws1316b[-1]).tolist())
wd1 = np.asarray((np.ones(9)*wd1218b[-1]).tolist()+
                 (np.ones(9)*wd1300b[-1]).tolist()+
                 (np.ones(9)*wd1312b[-1]).tolist()+
                 (np.ones(9)*wd1316b[-1]).tolist())


wss1218r, wds1218r =  probdata(r"/data5/NOBACKUP/nirb/Simulations/Haifa/aerofoil7/postProcessing/probes/0/", el=el0, ws0=ws1218b, wd0=wd1218b, line=-1)
wss1218w, wds1218w =  probdata(r"/data5/NOBACKUP/nirb/Simulations/Haifa/aerofoil7/postProcessing/probes/0/", el=el0, ws0=ws1218b, wd0=wd1218b, line=0)
wss1300r, wds1300r =  probdata(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306130000a2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1300b, wd0=wd1300b, line=-1)
wss1300w, wds1300w =  probdata(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306130000a2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1300b, wd0=wd1300b, line=0)
wss1312r, wds1312r =  probdata(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200a2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1312b, wd0=wd1312b, line=-1)
wss1312w, wds1312w =  probdata(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200a2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1312b, wd0=wd1312b, line=0)
wss1316r, wds1316r =  probdata(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131600a2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1316b, wd0=wd1316b, line=-1)
wss1316w, wds1316w =  probdata(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131600a2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1316b, wd0=wd1316b, line=0)

wss1218r, wds1218r =  probdata(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1218b, wd0=wd1218b, line=-1)
wss1218w, wds1218w =  probdata(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1218b, wd0=wd1218b, line=0)
wss1300r, wds1300r =  probdata(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306130000b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1300b, wd0=wd1300b, line=-1)
wss1300w, wds1300w =  probdata(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306130000b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1300b, wd0=wd1300b, line=0)
wss1312r, wds1312r =  probdata(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1312b, wd0=wd1312b, line=-1)
wss1312w, wds1312w =  probdata(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1312b, wd0=wd1312b, line=0)
wss1316r, wds1316r =  probdata(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131600b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1316b, wd0=wd1316b, line=-1)
wss1316w, wds1316w =  probdata(r"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131600b2rhosimple/postProcessing/probes/0/", el=el0, ws0=ws1316b, wd0=wd1316b, line=0)

wssr=np.asarray(wss1218r+wss1300r+wss1312r+wss1316r)
wssw=np.asarray(wss1218w+wss1300w+wss1312w+wss1316w)
wdsr=np.asarray(wds1218r+wds1300r+wds1312r+wds1316r)
wdsw=np.asarray(wds1218w+wds1300w+wds1312w+wds1316w)

print('s-',stat(wsm, wssr, kind='r2'),stat(wsm, wssr, kind='r'), stat(wsm, wssr, kind='mae'), stat(wsm, wssr, kind='rmse'),'/',round(np.mean(wssr),4),'+-',round(np.std(wssr),4),'>>',round(np.mean(wsm),4),'+-',round(np.std(wsm),4))
print('d-',stat(wdm, wdsr, kind='r2'),stat(wdm, wdsr, kind='r'), stat(wdm, wdsr, kind='mae'), stat(wdm, wdsr, kind='rmse'),'/',round(np.mean(wdsr),4),'+-',round(np.std(wdsr),4),'>>',round(np.mean(wdm),4),'+-',round(np.std(wdm),4))
print('s-',stat(wsm, wssw, kind='r2'),stat(wsm, wssw, kind='r'), stat(wsm, wssw, kind='mae'), stat(wsm, wssw, kind='rmse'),'/',round(np.mean(wssw),4),'+-',round(np.std(wssw),4),'>>',round(np.mean(wsm),4),'+-',round(np.std(wsm),4))
print('d-',stat(wdm, wdsw, kind='r2'),stat(wdm, wdsw, kind='r'), stat(wdm, wdsw, kind='mae'), stat(wdm, wdsw, kind='rmse'),'/',round(np.mean(wdsw),4),'+-',round(np.std(wdsw),4),'>>',round(np.mean(wdm),4),'+-',round(np.std(wdm),4))
print('s-',stat(wsm, ws1, kind='r2'),stat(wsm, ws1, kind='r'), stat(wsm, ws1, kind='mae'), stat(wsm, ws1, kind='rmse'),'/',round(np.mean(ws1),4),'+-',round(np.std(ws1),4),'>>',round(np.mean(wsm),4),'+-',round(np.std(wsm),4))
print('d-',stat(wdm, wd1, kind='r2'),stat(wdm, wd1, kind='r'), stat(wdm, wd1, kind='mae'), stat(wdm, wd1, kind='rmse'),'/',round(np.mean(wd1),4),'+-',round(np.std(wd1),4),'>>',round(np.mean(wdm),4),'+-',round(np.std(wdm),4))

distancerx = 0
distancery = 0
distancewx = 0
distancewy = 0
distance1x = 0
distance1y = 0
adistancerx = 0
adistancery = 0
adistancewx = 0
adistancewy = 0
adistance1x = 0
adistance1y = 0
maxr = 0
maxw = 0
max1 = 0
for i in range(len(wsm)):
    distancerx += (wsm[i]*math.cos(wdm[i]/180*math.pi)-wssr[i]*math.cos(wdsr[i]/180*math.pi))
    distancery += (wsm[i]*math.sin(wdm[i]/180*math.pi)-wssr[i]*math.sin(wdsr[i]/180*math.pi))
    distancewx += (wsm[i]*math.cos(wdm[i]/180*math.pi)-wssw[i]*math.cos(wdsw[i]/180*math.pi))
    distancewy += (wsm[i]*math.sin(wdm[i]/180*math.pi)-wssw[i]*math.sin(wdsw[i]/180*math.pi))
    distance1x += (wsm[i]*math.cos(wdm[i]/180*math.pi)-ws1[i]*math.cos(wd1[i]/180*math.pi))
    distance1y += (wsm[i]*math.sin(wdm[i]/180*math.pi)-ws1[i]*math.sin(wd1[i]/180*math.pi))
    adistancerx += math.fabs(wsm[i]*math.cos(wdm[i]/180*math.pi)-wssr[i]*math.cos(wdsr[i]/180*math.pi))
    adistancery += math.fabs(wsm[i]*math.sin(wdm[i]/180*math.pi)-wssr[i]*math.sin(wdsr[i]/180*math.pi))
    adistancewx += math.fabs(wsm[i]*math.cos(wdm[i]/180*math.pi)-wssw[i]*math.cos(wdsw[i]/180*math.pi))
    adistancewy += math.fabs(wsm[i]*math.sin(wdm[i]/180*math.pi)-wssw[i]*math.sin(wdsw[i]/180*math.pi))
    adistance1x += math.fabs(wsm[i]*math.cos(wdm[i]/180*math.pi)-ws1[i]*math.cos(wd1[i]/180*math.pi))
    adistance1y += math.fabs(wsm[i]*math.sin(wdm[i]/180*math.pi)-ws1[i]*math.sin(wd1[i]/180*math.pi))
    maxr = max(((wsm[i]*math.cos(wdm[i]/180*math.pi)-wssr[i]*math.cos(wdsr[i]/180*math.pi))**2 +
                 (wsm[i]*math.sin(wdm[i]/180*math.pi)-wssr[i]*math.sin(wdsr[i]/180*math.pi))**2)**.5, maxr)
    maxw = max(((wsm[i]*math.cos(wdm[i]/180*math.pi)-wssw[i]*math.cos(wdsw[i]/180*math.pi))**2 +
                 (wsm[i]*math.sin(wdm[i]/180*math.pi)-wssw[i]*math.sin(wdsw[i]/180*math.pi))**2)**.5, maxw)
    max1 = max(((wsm[i]*math.cos(wdm[i]/180*math.pi)-ws1[i]*math.cos(wd1[i]/180*math.pi))**2 +
                 (wsm[i]*math.sin(wdm[i]/180*math.pi)-ws1[i]*math.sin(wd1[i]/180*math.pi))**2)**.5, max1)
print('rho distance (sum, abs, max) = ', (distancerx**2+distancery**2)**.5/len(wsm), (adistancerx**2+adistancery**2)**.5/len(wsm), maxr)
print('wrf distance (sum, abs, max) = ', (distancewx**2+distancewy**2)**.5, (adistancewx**2+adistancewy**2)**.5, maxw)
print('1st distance (sum, abs, max) = ', (distance1x**2+distance1y**2)**.5, (adistance1x**2+adistance1y**2)**.5, max1)


plt.figure() 
plt.plot(wssr,wsm,'*b', label='CFD')
plt.plot(wssw,wsm,'*r', label='WRF')
xmin=min(wsm.min(),wssr.min())
xmax=max(wsm.max(),wssr.max())
ln11 = np.linspace(xmin,xmax)
plt.plot(ln11,ln11,'k')
plt.xlabel("simulation")
plt.ylabel("observaion")
plt.xlim(xmin,xmax)
plt.ylim(xmin,xmax)
plt.legend()
plt.title('Speed')

plt.figure()
plt.plot(wdsr,wdm,'*b', label='CFD')
plt.plot(wdsw,wdm,'*r', label='WRF')
xmin=min(wdm.min(),wdsr.min())
xmax=max(wdm.max(),wdsr.max())
ln11 = np.linspace(xmin,xmax)
plt.plot(ln11,ln11,'k')
plt.xlabel("simulation")
plt.ylabel("observaion")
plt.xlim(xmin,xmax)
plt.ylim(xmin,xmax)
plt.legend()
plt.title('Direction')

distancerx = 0
distancery = 0
adistancerx = 0
adistancery = 0
maxr = 0
ts0=ws1316b.copy()
td0=wd1316b.copy()
ts =np.ones(len(ts0))*ws1316b[-1]
td =np.ones(len(td0))*wd1316b[-1]
# ts = forecastedws.copy()
# td = forecastedwd.copy()
for i in range(len(ws1218b)):
    distancerx += (ts0[i]*math.cos(td0[i]/180*math.pi)-ts[i]*math.cos(td[i]/180*math.pi))
    distancery += (ts0[i]*math.sin(td0[i]/180*math.pi)-ts[i]*math.sin(td[i]/180*math.pi))
    adistancerx += math.fabs(ts0[i]*math.cos(td0[i]/180*math.pi)-ts[i]*math.cos(td[i]/180*math.pi))
    adistancery += math.fabs(ts0[i]*math.sin(td0[i]/180*math.pi)-ts[i]*math.sin(td[i]/180*math.pi))
    maxr = max(((ts0[i]*math.cos(td0[i]/180*math.pi)-ts[i]*math.cos(td[i]/180*math.pi))**2 +
                 (ts0[i]*math.sin(td0[i]/180*math.pi)-ts[i]*math.sin(td[i]/180*math.pi))**2)**.5, maxr)
print('distance (sum, abs, max) = ', (distancerx**2+distancery**2)**.5/len(ws1218b), (adistancerx**2+adistancery**2)**.5/len(ws1218b), maxr)

ua=np.linspace(-5,5,3)
va=np.linspace(-5,5,3)
for i in range(3):
    for j in range(3):
        da = np.arctan2(ua[i],va[j])*180/math.pi + 180
        print(ua[i],va[j],da, np.round(-5*np.sin(da/180*math.pi),2), np.round(-5*np.cos(da/180*math.pi),2) )
