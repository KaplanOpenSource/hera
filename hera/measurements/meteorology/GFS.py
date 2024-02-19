import math
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
    
    uu = np.asarray([1.08  ,  2.07,  1.12,  2.04,  2.09,  2.14]) #wrf2 comforepressure3
    vv = np.asarray([-3.1  , -2.70, -2.41, -2.56, -2.72, -2.25])
    
  
    # uu = np.asarray([]) #p5wrf probe10
    # vv = np.asarray([])
    
    # ws0=np.asarray([1.5,1,1.9,0.2,2.2,1.8]) #day1
    # wd0=np.asarray([58,35,73,258,289,283])
    
    import sklearn.metrics
    import math
    import numpy as np    

    ws0=np.asarray([1.5,2.,1.8,1.9,3.8,2.9]) # 13/06/2023 18GMT
    wd0=np.asarray([334,286,262,294,308,305])
    uu0=np.round(-ws0*np.sin(wd0/180*math.pi),2) #[ 0.66  1.92  1.78  1.74  2.99  2.38]
    vv0=np.round(-ws0*np.cos(wd0/180*math.pi),2) #[-1.35 -0.55  0.25 -0.77 -2.34 -1.66]
    
    # uu0 = np.asarray([ 0.65,  1.92, 1.78,  1.73,  2.99,  2.37])
    # vv0 = np.asarray([-1.34, -0.55, 0.25, -0.77, -2.33, -1.66])
    
    print('i,U,   V,    WS,   WD   << measurements')
    for i in range(6):
        print(i,round(uu0[i],2),round(vv0[i],2),ws0[i],wd0[i])
    
    print('i,U,   V,    WS,   WD  << SIMULATION')

    ws=np.zeros_like(uu)
    wd=np.zeros_like(uu)
    for i in range(len(uu)):
        wd[i]=round(math.atan(uu[i]/vv[i])*180/math.pi,2)
        ws[i]=round((uu[i]**2.+vv[i]**2.)**.5,2)
        if wd[i]<90:
            wd[i]+=360
        # wd[i] =math.atan2(vv[i], uu[i])*180/math.pi
        # wd[i] =math.atan(vv[i]/uu[i])*180/math.pi
        # if wd[i]<90:
        #     wd[i]+=180
        print(i,round(uu[i],2),round(vv[i],2),round(ws[i],2), round(wd[i],2))
        
    print('u:',stat(uu0, uu, kind='r2'),stat(uu0, uu, kind='r'), stat(uu0, uu, kind='rmse'),'/',round(np.mean(uu),4),'+-',round(np.std(uu),4),'>>',round(np.mean(uu0),4),'+-',round(np.std(uu0),4))
    print('v:',stat(vv0, vv, kind='r2'),stat(vv0, vv, kind='r'), stat(vv0, vv, kind='rmse'),'/',round(np.mean(vv),4),'+-',round(np.std(vv),4),'>>',round(np.mean(vv0),4),'+-',round(np.std(vv0),4))
    print('s:',stat(ws0, ws, kind='r2'),stat(ws0, ws, kind='r'), stat(ws0, ws, kind='mae'), stat(ws0, ws, kind='rmse'),'/',round(np.mean(ws),4),'+-',round(np.std(ws),4),'>>',round(np.mean(ws0),4),'+-',round(np.std(ws0),4))
    print('d:',stat(wd0, wd, kind='r2'),stat(wd0, wd, kind='r'), stat(wd0, wd, kind='mae'), stat(wd0, wd, kind='rmse'),'/',round(np.mean(wd),4),'+-',round(np.std(wd),4),'>>',round(np.mean(wd0),4),'+-',round(np.std(wd0),4))
    print('u (obs, model)')
    print (uu0)
    print (uu)
    print('v')
    print (vv0)
    print (vv)
    print('ws')
    print (ws0)
    print (ws)
    print('wd')
    print (wd0)
    print (wd)
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
    # plt.title('comfort -1.25, -.85, -.73, -1.91')
 
    haifa1 = [  
    [197075, 747605, 113],  #// kiryat shprinzak, ramot school france road 79    
    [199490, 745185, 300],  #//hugim, yair katz 4
    [198988, 743528, 300],  #//ahuze, smolskin st.
    [202315, 743652, 225],  #//hagalil 107, tel hai school
    [204115, 743860, 32],  #// check-post, moshlei yaakov 7
    [205603, 748278, 19]  #//dgania st. 
    ]
    stations=[]
    for i in range(len(haifa1)):
        # stations.append([haifa1[i][0], haifa1[i][1], wd0[i], haifa1[i][2]])
        stations.append([haifa1[i][0], haifa1[i][1], ws0[i], haifa1[i][2]])
    stations = np.asarray(stations)
    forecastedw=np.zeros_like(ws0)
    for i in range(len(stations)):
        loo=[] # leave one out
        for j in range(len(stations)):
            if i!=j:
                loo.append(stations[j])
        forecastedw[i] = interp(stations[i,0], stations[i,1], loo, elev = stations[i,3])
        # forecastidw= interp(stations[i,0], stations[i,1], loo)
    # print('d-',stat(wd0, forecastedw, kind='r2'),stat(wd0, forecastedw, kind='r'), stat(wd0, forecastedw, kind='mae'), stat(wd0, forecastedw, kind='rmse'),'/',round(np.mean(forecastedw),4),'+-',round(np.std(forecastedw),4),'>>',round(np.mean(wd0),4),'+-',round(np.std(wd0),4))
    print('s-',stat(ws0, forecastedw, kind='r2'),stat(ws0, forecastedw, kind='r'), stat(ws0, forecastedw, kind='mae'), stat(ws0, forecastedw, kind='rmse'),'/',round(np.mean(forecastedw),4),'+-',round(np.std(forecastedw),4),'>>',round(np.mean(ws0),4),'+-',round(np.std(ws0),4))
   
 
    
 
    if 5==7:
        import json
        import requests
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
                        