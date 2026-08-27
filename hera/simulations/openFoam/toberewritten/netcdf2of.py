#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 14:57:32 2023

@author: nirb
"""
    # This code take netcdf file and put it's data in openFOAM initial condition and boundary condition
    # OF needs coordinates of internal mesh, inlet, outlet, top, topo, buildings
    
    import datetime
    import numpy as np
    import matplotlib.pyplot as plt

    def remove_left(s):
        return s[1:]
    
    
    def remove_right(s):
        return s[:-1]

    directory = r'/data4bk/nirb/Simulations/Haifa/haifa32p2/48000/'  # line 5524
    directory = r'/data4bk/nirb/Simulations/Haifa/haifa32p4/0/'  # line 5524
    directory = r'/data4bk/nirb/Simulations/Haifa/haifa02p5/0/'  # line 5524
    directory = r'/data4bk/nirb/Simulations/Haifa/haifa02p5/0/'  # line 5524
    directory = r'/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200b2/0/'
    
    # postProcess -func writeCellCentres
    
    filex = directory + "Cx"
    filey = directory + "Cy"
    filez = directory + "Cz"
    filec = directory + "C"
    fileu = directory + "U"
    fhx = open(filex, 'r')
    fhy = open(filey, 'r')
    fhz = open(filez, 'r')
    fhu = open(fileu, 'r')
    fhc = open(filec, 'r')
    
    txtx = fhx.read()
    txty = fhy.read()
    txtz = fhz.read()
    txtu = fhu.read()
    txtc = fhc.read()
    fhx.close()
    fhy.close()
    fhz.close()
    fhu.close()
    fhc.close()
    
    txtxsplit = txtx.split()
    txtysplit = txty.split()
    txtzsplit = txtz.split()
    txtusplit = txtu.split()
    txtcsplit = txtc.split()

    def extractvalues(txt, face):
        pos = txt.index(face)
        if face == "internalField":
            leg = 3
        else:
            leg = 7
        if txt[pos+leg].isdigit():
            items = int(txt[pos+leg])
        else:
            items=0
        values = []
        print('start reading U field', items, datetime.datetime.now())
        if (txt[pos+leg+2][0]=='('):
            vector=3
        else:
            vector=1
            
        for i in range(items):
            if vector==3:
                values.append([remove_left(txt[pos+leg+2+i*3]), txt[pos+leg+3+i*3], remove_right(txt[pos+leg+4+i*3])])
            else:
                values.append(txt[pos+leg+2+i])
        # values = np.asarray(values)   
        return values
    
    def insertvalues(txt, face, vals, empty = False):
        pos = txt.index(face)
        if face == "internalField":
            leg = 3
        else:
            leg = 7
        if empty:
            items=len(vals)
        else:
            items = int(txt[pos+leg])
        values = []
        print('start writing U field', items, datetime.datetime.now())
        if (type(vals[0])==list or type(vals[0]==np.ndarray)):
            vector=3
        else:
            vector=1
            
        if empty:
            middle = []
            if vector==3:
                for i in range(items):
                    middle.append('('+str(round(vals[i][0],4)))
                    middle.append(str(round(vals[i][1],4)))
                    middle.append(str(round(vals[i][2],4))+')\n')
            txt = txt[:pos+leg+2]+middle+txt[pos+leg+2:]
        else:
            if vector==3:
                for i in range(items):
                    txt[pos+leg+2+items+0]='('+str(vals[i][0])
                    txt[pos+leg+2+items+1]=str(vals[i][1])
                    txt[pos+leg+2+items+2]=str(vals[i][2])+')\n'
            else:
                txt[pos+leg+2:pos+leg+2+items+1]=vals
        return txt
    

# del inletu
# del txtusplit

    import netCDF4 as nc

    filenc='/data4bk/nirb/Simulations/Haifa/gdas.t18z.atmf000.nc' # 2023061318gmt
    
    import xarray as xr
    import geopandas as gpd
    from shapely.geometry import Point
    dnc = xr.open_dataset(filenc)  

    ugrd=dnc.ugrd[0,:,:,:].data
    vgrd=dnc.vgrd[0,:,:,:].data
    wgrd=dnc.wgrd[0,:,:,:].data
    

    xlatlon = np.load('gfsxlatlon.npy')
    # ylatlon = np.load('gfsylatlon.npy')
    heights = np.load('gfsheights.npy')




    import pyproj
    from pyproj import CRS
    from pyproj import Transformer
    
    crs = CRS.from_epsg(6991) #6991
    crs.to_epsg()
    crs = CRS.from_proj4("+proj=tmerc +lat_0=31.7343936111111 +lon_0=35.2045169444445 +k=1.0000067 +x_0=219529.584 +y_0=626907.39 +ellps=GRS80 +towgs84=-24.002400,-17.103200,-17.844400,-0.33007,-1.852690,1.669690,5.424800 +units=m +no_defs")
    transformer = Transformer.from_crs("EPSG:6991", "EPSG:4326")


    inletx = extractvalues(txtxsplit, "inlet")    
    inlety = extractvalues(txtysplit, "inlet")    
    inletz = extractvalues(txtzsplit, "inlet")    
    inletzmax = np.asarray([float(i) for i in inletz]).max()
    
    inletu = extractvalues(txtusplit, "inlet")    
    for i in range(len(inletu)):
        netxindex = inletx[i]
        netyindex = inlety[i]
        netzindex = float(inletz[i])
        x1=transformer.transform(netxindex,netyindex)
        xindex=np.argmin(np.abs(dnc.grid_yt.data-x1[0]))
        yindex=np.argmin(np.abs(dnc.grid_xt.data-x1[1]))
        zindex=np.argmin(np.abs(heights[:,xindex,yindex]-netzindex))
        if i % 100 ==0:
            print(i, xindex, yindex, zindex)
        inletu[i][0]=ugrd[zindex,xindex,yindex]
        inletu[i][1]=vgrd[zindex,xindex,yindex]
        inletu[i][2]="0."
    txtusplit = insertvalues(txtusplit, "inlet", inletu)

    inletx = extractvalues(txtxsplit, "outlet")    
    inlety = extractvalues(txtysplit, "outlet")    
    inletz = extractvalues(txtzsplit, "outlet")    
    inletu = extractvalues(txtusplit, "outlet")    
    for i in range(len(inletu)):
        netxindex = inletx[i]
        netyindex = inlety[i]
        netzindex = float(inletz[i])
        x1=transformer.transform(netxindex,netyindex)
        xindex=np.argmin(np.abs(dnc.grid_yt.data-x1[0]))
        yindex=np.argmin(np.abs(dnc.grid_xt.data-x1[1]))
        zindex=np.argmin(np.abs(heights[:,xindex,yindex]-netzindex))
        if i % 100 ==0:
            print(i, xindex, yindex, zindex)
        inletu[i][0]=ugrd[zindex,xindex,yindex]
        inletu[i][1]=vgrd[zindex,xindex,yindex]
        inletu[i][2]="0."
    txtusplit = insertvalues(txtusplit, "outlet", inletu)

    inletx = extractvalues(txtxsplit, "internalField")    
    inlety = extractvalues(txtysplit, "internalField")    
    inletz = extractvalues(txtzsplit, "internalField")    
    inletu = extractvalues(txtusplit, "internalField")    
    for i in range(len(inletu)):
        netxindex = inletx[i]
        netyindex = inlety[i]
        netzindex = float(inletz[i])
        x1=transformer.transform(netxindex,netyindex)
        xindex=np.argmin(np.abs(dnc.grid_yt.data-x1[0]))
        yindex=np.argmin(np.abs(dnc.grid_xt.data-x1[1]))
        zindex=np.argmin(np.abs(heights[:,xindex,yindex]-netzindex))
        if i % 100 ==0:
            print(i,'/',len(inletu), xindex, yindex, zindex)
        inletu[i][0]=ugrd[zindex,xindex,yindex]
        inletu[i][1]=vgrd[zindex,xindex,yindex]
        inletu[i][2]="0."
    txtusplit = insertvalues(txtusplit, "internalField", inletu)

    # inletx = extractvalues(txtxsplit, "top")    
    # inlety = extractvalues(txtysplit, "top")    
    # inletz = extractvalues(txtzsplit, "top")
    # if len(inletz)==0:
    #         for i in range(len(inletx)):
    #             inletz.append(inletzmax)
    # inletu = extractvalues(txtusplit, "top")    
    # for i in range(len(inletu)):
    #     netxindex = inletx[i]
    #     netyindex = inlety[i]
    #     netzindex = float(inletz[i])
    #     if i % 100 ==0:
    #         print(i, xindex, yindex, zindex)
    #     x1=transformer.transform(netxindex,netyindex)
    #     xindex=np.argmin(np.abs(dnc.grid_yt.data-x1[0]))
    #     yindex=np.argmin(np.abs(dnc.grid_xt.data-x1[1]))
    #     zindex=np.argmin(np.abs(heights[:,xindex,yindex]-netzindex))
    #     inletu[i][0]=ugrd[zindex,xindex,yindex]
    #     inletu[i][1]=vgrd[zindex,xindex,yindex]
    #     inletu[i][2]="0."
    # txtusplit = insertvalues(txtusplit, "top", inletu)

    # inletx = extractvalues(txtxsplit, "buildings")    
    # inlety = extractvalues(txtysplit, "buildings")    
    # inletz = extractvalues(txtzsplit, "buildings")    
    # inletu = extractvalues(txtusplit, "buildings")    
    # for i in range(len(inletu)):
    #     netxindex = inletx[i]
    #     netyindex = inlety[i]
    #     netzindex = float(inletz[i])
    #     if i % 100 ==0:
    #         print(i, xindex, yindex, zindex)
    #     x1=transformer.transform(netxindex,netyindex)
    #     xindex=np.argmin(np.abs(dnc.grid_yt.data-x1[0]))
    #     yindex=np.argmin(np.abs(dnc.grid_xt.data-x1[1]))
    #     zindex=np.argmin(np.abs(heights[:,xindex,yindex]-netzindex))
    #     inletu[i][0]=ugrd[zindex,xindex,yindex]
    #     inletu[i][1]=vgrd[zindex,xindex,yindex]
    #     inletu[i][2]="0."
    # txtusplit = insertvalues(txtusplit, "buildings", inletu)

    # inletx = extractvalues(txtxsplit, "topo")    
    # inlety = extractvalues(txtysplit, "topo")    
    # inletz = extractvalues(txtzsplit, "topo")    
    # inletu = extractvalues(txtusplit, "topo")    
    # for i in range(len(inletu)):
    #     netxindex = inletx[i]
    #     netyindex = inlety[i]
    #     netzindex = float(inletz[i])
    #     if i % 100 ==0:
    #         print(i, xindex, yindex, zindex)
    #     x1=transformer.transform(netxindex,netyindex)
    #     xindex=np.argmin(np.abs(dnc.grid_yt.data-x1[0]))
    #     yindex=np.argmin(np.abs(dnc.grid_xt.data-x1[1]))
    #     zindex=np.argmin(np.abs(heights[:,xindex,yindex]-netzindex))
    #     inletu[i][0]=ugrd[zindex,xindex,yindex]
    #     inletu[i][1]=vgrd[zindex,xindex,yindex]
    #     inletu[i][2]="0."
    # txtusplit = insertvalues(txtusplit, "topo", inletu)



    utxt2 = ' '.join(txtusplit)
    fhuw = open(r'/data4bk/nirb/Simulations/Haifa/haifa32p2a/0/U', 'w')
    start = 0
    fin = len(txtusplit)
    indextxt = 0
    while indextxt<fin:
        if txtusplit[indextxt][-1]==')':
            utxt3 = ' '.join(txtusplit[start:indextxt+1])        
            fhuw.writelines(utxt3+'\n')
            # fhuw.writelines(utxt2[start:indextxt+1])
            # print('-----',start, indextxt)
            start=indextxt+1
        indextxt+=1
        if indextxt%10000==0:
            print('>>>>',indextxt)

    # fhuw.writelines(utxt2)
    utxt3 = ' '.join(txtusplit[start:indextxt+1])        
    fhuw.writelines(utxt3+'\n')
    fhuw.close()

    
    internalFieldx = extractvalues(txtxsplit, "internalField")    
    outletx = extractvalues(txtxsplit, "outlet")    
    inletx = extractvalues(txtxsplit, "inlet")    
    topx = extractvalues(txtxsplit, "top")    
    buildingsx = extractvalues(txtxsplit, "buildings")    
    topox = extractvalues(txtxsplit, "topo")    
    
    internalFieldy = extractvalues(txtysplit, "internalField")    
    outlety = extractvalues(txtysplit, "outlet")    
    inlety = extractvalues(txtysplit, "inlet")    
    topy = extractvalues(txtysplit, "top")    
    buildingsy = extractvalues(txtysplit, "buildings")    
    topoy = extractvalues(txtysplit, "topo")    

    internalFieldz = extractvalues(txtzsplit, "internalField")    
    outletz = extractvalues(txtzsplit, "outlet")    
    inletz = extractvalues(txtzsplit, "inlet")    
    # topz = extractvalues(txtzsplit, "top")    
    topz=np.zeros_like(topy)
    topz[:]=4000
    buildingsz = extractvalues(txtzsplit, "buildings")    
    topoz = extractvalues(txtzsplit, "topo")    
    
    import netCDF4 as nc

    filenc='/data4bk/nirb/Simulations/Haifa/gdas.t18z.atmf000.nc' # 2023061318gmt
    
    import xarray as xr
    import geopandas as gpd
    from shapely.geometry import Point
    dnc = xr.open_dataset(filenc)  
    # # dnc = dnc.reset_index()
    # geometry = [Point(xy) for xy in zip(dnc.grid_xt, dnc.grid_yt)]
    # geo=gpd.GeoDataFrame(dnc, geometry=geometry)
    # gdf.set_crs(epsg=4326, inplace=True)
    # gdf.to_crs(epsg=3395)
  
    
    # df = dnc.to_dataframe()
    print(dnc['delz'][0,:,66,66].data[4])
    


    crs = CRS.from_epsg(6991) #6991
    crs.to_epsg()
    crs = CRS.from_proj4("+proj=tmerc +lat_0=31.7343936111111 +lon_0=35.2045169444445 +k=1.0000067 +x_0=219529.584 +y_0=626907.39 +ellps=GRS80 +towgs84=-24.002400,-17.103200,-17.844400,-0.33007,-1.852690,1.669690,5.424800 +units=m +no_defs")
    transformer = Transformer.from_crs("EPSG:6991", "EPSG:4326")
    x1=transformer.transform(247814.9186634097,764765.7875774748)


    
    crs = CRS.from_epsg(4326) #6991
    crs.to_epsg()
    crs = CRS.from_proj4("+proj=tmerc +lat_0=31.7343936111111 +lon_0=35.2045169444445 +k=1.0000067 +x_0=219529.584 +y_0=626907.39 +ellps=GRS80 +towgs84=-24.002400,-17.103200,-17.844400,-0.33007,-1.852690,1.669690,5.424800 +units=m +no_defs")
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:6991")
    # transformer.transform(220080, 634451)
    x2=transformer.transform(32.727616, 35.103320)
    xlatlon=np.zeros([len(dnc.grid_yt), len(dnc.grid_xt)])
    ylatlon=np.zeros([len(dnc.grid_yt), len(dnc.grid_xt)])
    for i in range(len(dnc.grid_yt)):
        if i % 100 == 0:
            print('iy',i,len(dnc.grid_yt), datetime.datetime.now())
        for j in range(len(dnc.grid_yt)):
            x2=transformer.transform(dnc.grid_yt.data[i], dnc.grid_xt.data[j])
            xlatlon[i,j]=x2[1]
            ylatlon[i,j]=x2[0]

    np.save('wrfxlatlon',xlatlon)
    np.save('wrfylatlon',ylatlon)
    # xlatlon = np.load('gfsxlatlon.npy')
    # ylatlon = np.load('gfsylatlon.npy')
    
    heights=np.zeros([len(dnc.pfull), len(dnc.grid_yt), len(dnc.grid_xt)])
    delz=dnc.delz[0,:,:,:].data
    ugrd=dnc.ugrd[0,:,:,:].data
    vgrd=dnc.vgrd[0,:,:,:].data
    lpfull=len(dnc.pfull)
    for i in range(len(dnc.grid_yt)):
        if i % 100 == 0:
            print('iy',i,len(dnc.grid_yt), datetime.datetime.now())
        for j in range(len(dnc.grid_yt)):
            for k in range(lpfull):
                heights[k,i,j]=-np.sum(delz[-(lpfull-k):,i,j])
                
    np.save('gfsheights',heights)
    # heights = np.load('gfsheights.npy')
    fhuw = open(r'/data4bk/nirb/Simulations/Haifa/haifa32p1/58000/U2', 'w')
            

    fhuw.writelines(utxt)
    fhuw.close()

    
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
    
    # uu = np.asarray([0.24, 1.62, 1.48, 0.64, 2.78, 2.98])
    # vv = np.asarray([-3.98, -2.36, -1.89, -0.37, -1.43, -0.35])
    # uu = np.asarray([.35,1.21, 1.29, 2.96, 1.65, 1.90])
    # vv = np.asarray([-3.34, -1.36, -1.02, -1.20, -1.22, -0.42])

    # uu = np.asarray([-1.7,1.29,.88, 1.46, 3.22,2.77])
    # vv = np.asarray([-4.74, -2.72, -3.0,-.69,-2.34, -3.25])
    
    # uu = np.asarray([-1.31,-.69, -1.41, .71, 1.20, 1.74])
    # vv = np.asarray([-3.04, -2.47, -3.10, -2.99, -1.43, -1.32])

    # uu = np.asarray([-1.52, -.94,-1.55, 0.48, 0.96, 0.40]) #p2
    # vv = np.asarray([-2.6, -2.24, -2.66, -2.83, -3.48, -3.18])


    # uu = np.asarray([-1.62, -0.73, -1.53, 0.74, 1.05, 0.91]) #p 30525
    # vv = np.asarray([-3.25, -1.99, 2.8, -2.8, -2.27, -2.53])
    
    # ws0=np.asarray([1.5,1,1.9,0.2,2.2,1.8]) #day1
    # wd0=np.asarray([58,35,73,258,289,283])

    # ws0=np.asarray([1.5,2.,1.8,1.9,3.8,2.9]) # 1306202318GMT
    # wd0=np.asarray([334,286,262,294,308,305])
    # uu0=-ws0*np.sin(wd0/180*math.pi)
    # vv0=-ws0*np.cos(wd0/180*math.pi)
    
    # print('i,U,   V,    WS,   WD   << measurements')
    # for i in range(6):
    #     print(i,round(uu0[i],2),round(vv0[i],2),ws0[i],wd0[i])
    
    # ws=np.zeros_like(uu)
    # wd=np.zeros_like(uu)
    # print('i,U,   V,    WS,   WD  << SIMULATION')
    # for i in range(len(uu)):
    #     wd[i] =math.atan2(vv[i], uu[i])*180/math.pi
    #     # wd[i] =math.atan(vv[i]/uu[i])*180/math.pi
    #     wd[i] =math.atan(uu[i]/vv[i])*180/math.pi
    #     ws[i]=(uu[i]**2.+vv[i]**2.)**.5
    #     if wd[i]<90:
    #         wd[i]+=360
    #     # if wd[i]<90:
    #     #     wd[i]+=180
    #     print(i,round(uu[i],2),round(vv[i],2),round(ws[i],2), round(wd[i],2))
        
    # print('u-',stat(uu0, uu, kind='r2'),stat(uu0, uu, kind='r'), stat(uu0, uu, kind='rmse'))
    # print('v-',stat(vv0, vv, kind='r2'),stat(vv0, vv, kind='r'), stat(vv0, vv, kind='rmse'))
    # print('s-',stat(ws0, ws, kind='r2'),stat(ws0, ws, kind='r'), stat(ws0, ws, kind='mae'), stat(ws0, ws, kind='rmse'))
    # print('d-',stat(wd0, wd, kind='r2'),stat(wd0, wd, kind='r'), stat(wd0, wd, kind='mae'), stat(wd0, wd, kind='rmse'))

    filenc='/data4bk/nirb/wrf12062023.npz' # 2023061318gmt
    data = np.load(filenc)
    print(data.files)
    print(data['xlat'], data['xlong'])
    # heights=data['gpt_hgt_M'] # domain
    # heights=data['hgt'] # surface
    for i in range(len(data['Times'])):
        s=''
        for j in range(len(data['Times'][i][:])):
            s+=str(data['Times'][i][j])[2]
        print(i,s)
    si = 36 # chosen time
    heights = data['gpt_hgt_M'][si,:,:,:]
    # for i in range(heights.shape[0]):
        # heights[i,:,:]+=data['hgt']
        
    ugrd=data['Ue'][si,:,:,:]
    vgrd=data['Ve'][si,:,:,:]
    wgrd=data['W'][si,:,:,:]

    crs = CRS.from_epsg(4326) #6991
    crs.to_epsg()
    crs = CRS.from_proj4("+proj=tmerc +lat_0=31.7343936111111 +lon_0=35.2045169444445 +k=1.0000067 +x_0=219529.584 +y_0=626907.39 +ellps=GRS80 +towgs84=-24.002400,-17.103200,-17.844400,-0.33007,-1.852690,1.669690,5.424800 +units=m +no_defs")
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:6991")
    # transformer.transform(220080, 634451)
    x2=transformer.transform(32.727616, 35.103320)
    xlatlon=np.zeros_like(data['xlat'])
    ylatlon=np.zeros_like(data['xlat'])
    for i in range(xlatlon.shape[0]):
        if i % 100 == 0:
            print('iy',i,xlatlon.shape[0], datetime.datetime.now())
        for j in range(xlatlon.shape[1]):
            x2=transformer.transform(data['xlat'][i,j], data['xlong'][i,j])
            xlatlon[i,j]=x2[1]
            ylatlon[i,j]=x2[0]


    crs = CRS.from_epsg(6991) #6991
    crs.to_epsg()
    crs = CRS.from_proj4("+proj=tmerc +lat_0=31.7343936111111 +lon_0=35.2045169444445 +k=1.0000067 +x_0=219529.584 +y_0=626907.39 +ellps=GRS80 +towgs84=-24.002400,-17.103200,-17.844400,-0.33007,-1.852690,1.669690,5.424800 +units=m +no_defs")
    transformer = Transformer.from_crs("EPSG:6991", "EPSG:4326")


    inletx = extractvalues(txtxsplit, "inlet")    
    inlety = extractvalues(txtysplit, "inlet")    
    inletz = extractvalues(txtzsplit, "inlet")    
    inletc = extractvalues(txtcsplit, "inlet")    
    inletzmax = np.asarray([float(i) for i in inletz]).max()
    
    empty = True
    if empty:
        inletu = extractvalues(txtcsplit, "inlet")    #############
        print('empty')
    else:
        inletu = extractvalues(txtusplit, "inlet")    
    for i in range(len(inletu)):
        netxindex = inletx[i]
        netyindex = inlety[i]
        netzindex = float(inletz[i])
        x1=transformer.transform(netxindex,netyindex)
        xindex=np.argmin(np.abs(data['xlat'][:,0]-x1[0]))
        yindex=np.argmin(np.abs(data['xlong'][0,:]-x1[1]))
        zindex=np.argmin(np.abs(heights[:,xindex,yindex]-netzindex))
        if i % 100 ==0:
            print(i,'/', len(inletu),xindex, yindex, zindex)
        inletu[i][0]=ugrd[zindex,xindex,yindex]
        inletu[i][1]=vgrd[zindex,xindex,yindex]
        inletu[i][2]=wgrd[zindex,xindex,yindex]
    txtusplit = insertvalues(txtusplit, "inlet", inletu, empty)

    inletx = extractvalues(txtxsplit, "outlet")    
    inlety = extractvalues(txtysplit, "outlet")    
    inletz = extractvalues(txtzsplit, "outlet")    
    inletc = extractvalues(txtcsplit, "outlet")    
    inletzmax = np.asarray([float(i) for i in inletz]).max()
    
    empty = True
    if empty:
        inletu = extractvalues(txtcsplit, "outlet")    #############
        print('empty')
    else:
        inletu = extractvalues(txtusplit, "outlet")    
    for i in range(len(inletu)):
        netxindex = inletx[i]
        netyindex = inlety[i]
        netzindex = float(inletz[i])
        x1=transformer.transform(netxindex,netyindex)
        xindex=np.argmin(np.abs(data['xlat'][:,0]-x1[0]))
        yindex=np.argmin(np.abs(data['xlong'][0,:]-x1[1]))
        zindex=np.argmin(np.abs(heights[:,xindex,yindex]-netzindex))
        if i % 100 ==0:
            print(i,'/', len(inletu),xindex, yindex, zindex)
        inletu[i][0]=ugrd[zindex,xindex,yindex]
        inletu[i][1]=vgrd[zindex,xindex,yindex]
        inletu[i][2]=wgrd[zindex,xindex,yindex]
    txtusplit = insertvalues(txtusplit, "outlet", inletu, empty)



    inletx = extractvalues(txtxsplit, "top")    
    inlety = extractvalues(txtysplit, "top")    
    inletz = extractvalues(txtzsplit, "top")    
    if empty:
        inletu = extractvalues(txtcsplit, "top")    #############
        print('empty')
    else:
        inletu = extractvalues(txtusplit, "top")    
    for i in range(len(inletu)):
        netxindex = inletx[i]
        netyindex = inlety[i]
        netzindex = float(inletz[i])
        x1=transformer.transform(netxindex,netyindex)
        xindex=np.argmin(np.abs(data['xlat'][:,0]-x1[0]))
        yindex=np.argmin(np.abs(data['xlong'][0,:]-x1[1]))
        zindex=np.argmin(np.abs(heights[:,xindex,yindex]-netzindex))
        if i % 100 ==0:
            print(i, '/',len(inletu),xindex, yindex, zindex)
        # if (zindex>1):
            # if ugrd[zindex,xindex,yindex]<2:
                # print('outlet:',zindex,xindex,yindex, ugrd[zindex,xindex,yindex])
        inletu[i][0]=ugrd[zindex,xindex,yindex]
        inletu[i][1]=vgrd[zindex,xindex,yindex]
        inletu[i][2]=wgrd[zindex,xindex,yindex]
    txtusplit = insertvalues(txtusplit, "top", inletu, empty)

    inletx = extractvalues(txtxsplit, "internalField")    
    inlety = extractvalues(txtysplit, "internalField")    
    inletz = extractvalues(txtzsplit, "internalField")    
    if empty:
        inletu = extractvalues(txtcsplit, "internalField")    #############
        print('empty')
    else:
        inletu = extractvalues(txtusplit, "internalField")    
    
    for i in range(len(inletu)):
        netxindex = inletx[i]
        netyindex = inlety[i]
        netzindex = float(inletz[i])
        x1=transformer.transform(netxindex,netyindex)
        xindex=np.argmin(np.abs(data['xlat'][:,0]-x1[0]))
        yindex=np.argmin(np.abs(data['xlong'][0,:]-x1[1]))
        zindex=np.argmin(np.abs(heights[:,xindex,yindex]-netzindex))
        if i % 100 ==0:
            print(i,'/',len(inletu), xindex, yindex, zindex)
        inletu[i][0]=ugrd[zindex,xindex,yindex]
        inletu[i][1]=vgrd[zindex,xindex,yindex]
        inletu[i][2]=wgrd[zindex,xindex,yindex]
    txtusplit = insertvalues(txtusplit, "internalField", inletu, empty)


    utxt2 = ' '.join(txtusplit)
    
    fhuw = open(r'/data5/nirb/Simulations/Haifa/haifa22p5/0/U2', 'w')
    # start = 0
    # fin = len(txtusplit)
    # indextxt = 0
    # while indextxt<fin:
    #     if txtusplit[indextxt][-1]==')':
    #         utxt3 = ' '.join(txtusplit[start:indextxt+1])        
    #         fhuw.writelines(utxt3+'\n')
    #         # fhuw.writelines(utxt2[start:indextxt+1])
    #         # print('-----',start, indextxt)
    #         start=indextxt+1
    #     indextxt+=1
    #     if indextxt%10000==0:
    #         print('>>>>',indextxt)

    # fhuw.writelines(utxt2)
    # utxt3 = ' '.join(txtusplit[start:indextxt+1])        
    utxt3 = ' '.join(txtusplit)        
    fhuw.writelines(utxt3+'\n')
    fhuw.close()

plt.figure()
plt.imshow(data['Ue'][48,17,:,:]) #19
plt.title(heights[17,100,100])
plt.colorbar()

for i in range(len(inletu)):
    if inletu[i][0]<0.:
        print(i,inletu[i])





    fhuw = open(r'/data4bk/nirb/Simulations/Haifa/haifa32p2a/0/U', 'w')
    # start = 0
    # fin = len(txtusplit)
    # indextxt = 0
    # while indextxt<fin:
    #     if txtusplit[indextxt][-1]==')':
    #         utxt3 = ' '.join(txtusplit[start:indextxt+1])        
    #         fhuw.writelines(utxt3+'\n')
    #         # fhuw.writelines(utxt2[start:indextxt+1])
    #         # print('-----',start, indextxt)
    #         start=indextxt+1
    #     indextxt+=1
    #     if indextxt%10000==0:
    #         print('>>>>',indextxt)

    # fhuw.writelines(utxt2)
    # utxt3 = ' '.join(txtusplit[start:indextxt+1])        
    utxt3 = ' '.join(txtusplit)        
    fhuw.writelines(utxt3+'\n')
    fhuw.close()



    fhu = open(r'/data4bk/nirb/Simulations/Haifa/haifa22p5/0/U2', 'r')
    txtu = fhu.read()
    fhu.close()
    txtusplit = txtu.split()
    for i in range(len(txtusplit)):
        if i % 10000 == 0:
            print (i , '/',  len(txtusplit))
        try:
            txtusplit[i]=str(round(float(txtusplit[i]),4))
        except:
            try:
                txtusplit[i]='('+str(round(float(txtusplit[i][1:]),4))
            except:
                try:
                    txtusplit[i]=str(round(float(txtusplit[i][:-1]),4))+')\n'
                except:
                    pass
            
            

    utxt2 = ' '.join(txtusplit)
    fhuw = open(r'/data4bk/nirb/Simulations/Haifa/haifa22p5/0/U2', 'w')
    # start = 0
    # fin = len(txtusplit)
    # indextxt = 0
    # while indextxt<fin:
    #     if txtusplit[indextxt][-1]==')':
    #         utxt3 = ' '.join(txtusplit[start:indextxt+1])        
    #         fhuw.writelines(utxt3+'\n')
    #         # fhuw.writelines(utxt2[start:indextxt+1])
    #         # print('-----',start, indextxt)
    #         start=indextxt+1
    #     indextxt+=1
    #     if indextxt%10000==0:
    #         print('>>>>',indextxt)

    # fhuw.writelines(utxt2)
    # utxt3 = ' '.join(txtusplit[start:indextxt+1])        
    utxt3 = ' '.join(txtusplit)        
    fhuw.writelines(utxt3+'\n')
    fhuw.close()
    
######### stability start #################    
    filenc='/data4bk/nirb/wrf12062023.npz' # 2023061318gmt
    data = np.load(filenc)
    print(data.files)
    print(data['xlat'], data['xlong'])
    # data['xlat'][52][:]
    # data['xlong'][:][54]
    latindex = 43 # carmel
    longindex = 60 # haifa carmel

    latindex = 52 # haifa sea
    longindex = 54 # haifa sea

    latindex =  50 # haifa sea
    longindex = 56 # haifa sea
    # heights=data['gpt_hgt_M'] # domain
    # heights=data['hgt'] # surface
    print(data['xlat'][latindex,longindex],data['xlong'][latindex,longindex],data['hgt'][latindex,longindex])
    for i in range(len(data['Times'])):
        s=''
        for j in range(len(data['Times'][i][:])):
            s+=str(data['Times'][i][j])[2]
        print(i,s)



    top = 12  #12 # 52
    heights = data['gpt_hgt_M'][si,:top,latindex,longindex]


    plt.figure()
    # for i in range(0,len(data['Times']),12):
    for i in range(44,45):
    
        ttl =''
        for j in range(len(data['Times'][i][:])):
            ttl+=str(data['Times'][i][j])[2]   
        # wd = np.degrees(np.arctan2(data['Ue'][i,10,latindex,longindex], data['Ve'][i,10,latindex,longindex])) % 360 + 180
        print(i, ttl, wd)
        # uv = (data['Ue'][i,:top,latindex,longindex]**2+data['Ve'][i,:top,latindex,longindex]**2)**0.5
        # temp = data['Temp'][i,:top,latindex,longindex]
        temp = data['Theta'][i,:top,latindex,longindex]
        # for i in range(heights.shape[0]):
            # heights[i,:,:]+=data['hgt']
    
        # plt.figure()
        plt.plot(temp, heights, label = ttl)


    plt.legend()
    plt.show()    

    
    
######### stability end #################    
    
