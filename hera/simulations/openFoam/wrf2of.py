#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 11:02:33 2024
Issue 159
@author: nirb
"""
from hera import toolkitHome
import numpy as np
import datetime
import pyproj
from pyproj import CRS
from pyproj import Transformer
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator


def interp2(x,y, x0,x1,y0,y1,d00,d01,d10,d11):
    # d01 is d at x=0, y=1
    d0 = d00*(x1-x)/(x1-x0) + d10*(x-x0)/(x1-x0)
    d1 = d01*(x1-x)/(x1-x0) + d11*(x-x0)/(x1-x0)
    d = d0*(y1-y)/(y1-y0) + d1*(y-y0)/(y1-y0)
    # print(d0,d1,d)
    return d

par = False
# reading of case
tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM)
case = r'/data5/NOBACKUP/nirb/Simulations/Haifa/testcode'
# cc = tk.OFObjectHome.readFieldFromCase(fieldName="cellCenters",flowType=tk.FLOWTYPE_INCOMPRESSIBLE,caseDirectory=case)


U =   tk.OFObjectHome.readFieldFromCase("U",tk.FLOWTYPE_INCOMPRESSIBLE,case, readParallel=par) #"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2/")
Udf = tk.OFObjectHome.readFieldAsDataFrame("U", case, readParallel=par)
T =   tk.OFObjectHome.readFieldFromCase("T",tk.FLOWTYPE_INCOMPRESSIBLE,case, readParallel=par) #"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2/")
Tdf = tk.OFObjectHome.readFieldAsDataFrame("T", case, readParallel=par)
print(Udf)
print(U.componentNames)
A = tk.getMesh(case, readParallel=par)
Adf = A.getDataFrame() 


# read the WRF file
filenc='/data4bk/nirb/wrf12062023.npz' # 2023061318gmt
data = np.load(filenc)
# print(data.files)
# print(data['xlat'], data['xlong'])

# print the times so we can choose
for i in range(len(data['Times'])):
    s=''
    for j in range(len(data['Times'][i][:])):
        s+=str(data['Times'][i][j])[2]
    print(i,s)
si = 36 # chosen time


heights = data['gpt_hgt_M'][si,:,:,:] # altitude
# for i in range(heights.shape[0]):
    # heights[i,:,:]+=data['hgt'] # ground elevation
    
# ugrd=data['Ue'][si,:,:,:]
# vgrd=data['Ve'][si,:,:,:]
# wgrd=data['W'][si,:,:,:]


# transfer coordinates to ITM, this is what OF wants
crs = CRS.from_epsg(4326) #6991
crs.to_epsg()
crs = CRS.from_proj4("+proj=tmerc +lat_0=31.7343936111111 +lon_0=35.2045169444445 +k=1.0000067 +x_0=219529.584 +y_0=626907.39 +ellps=GRS80 +towgs84=-24.002400,-17.103200,-17.844400,-0.33007,-1.852690,1.669690,5.424800 +units=m +no_defs")
transformer = Transformer.from_crs("EPSG:4326", "EPSG:6991")
# transformer.transform(220080, 634451)
x2=transformer.transform(32.727616, 35.103320)


crs = CRS.from_epsg(6991) #6991 #2039
crs.to_epsg()
crs = CRS.from_proj4("+proj=tmerc +lat_0=31.7343936111111 +lon_0=35.2045169444445 +k=1.0000067 +x_0=219529.584 +y_0=626907.39 +ellps=GRS80 +towgs84=-24.002400,-17.103200,-17.844400,-0.33007,-1.852690,1.669690,5.424800 +units=m +no_defs")
transformer = Transformer.from_crs("EPSG:6991", "EPSG:4326")
transformer.transform(209977.17691279756, 737011.5791438483)


# interpolate the WRF coordinates to OF coordinates
    
datax= data['xlat']
datay= data['xlong']
datau = data['Ue'][si,:,:,:]
datav = data['Ve'][si,:,:,:]
Adz = Adf.Cz.values
uxi = []
vxi = []
# f(716373,153153)[0]  #15,15
x0=transformer.transform(Adf.Cx.values, Adf.Cy.values)

for i in range(len(Adz)):
    if i%1000==0: print(i,'/',len(Adf.Cx.values))
    # n0=datetime.datetime.now()
    # it is almost regular data, so we need better method than scatter interpolation, but not regular
    # x0=transformer.transform(Adf.Cx.values[i], Adf.Cy.values[i])
    ilat = np.argmin(np.abs(data['xlat'][:,data['xlat'].shape[1]//2]-x0[0][i]))
    if data['xlat'][ilat,data['xlat'].shape[1]//2]-x0[0][i]>0:
        ilat-=1
        # print('good practice ilat')
    ilong = np.argmin(np.abs(data['xlong'][ilat,:]-x0[1][i]))
    if data['xlong'][ilat,ilong]-x0[1][i]>0:
        ilong-=1
        # print('good practice ilong')
    # print(x0,data['xlat'][ilat,ilong],data['xlat'][ilat+1,ilong+1],data['xlong'][ilat,ilong],data['xlong'][ilat+1,ilong+1])
    level = np.argmin(np.abs(heights[:,ilat,ilong]-Adz[i]))
    if heights[level,ilat,ilong]-Adz[i]>0:
        level-=1
        # print('good practice level')
    # x = np.asarray([data['xlat'][ilat,ilong],data['xlat'][ilat+1,ilong]])
    # y = np.asarray([data['xlong'][ilat,ilong],data['xlong'][ilat,ilong+1]])
    # pts = np.array([x[0], x0[1]])
    # datac =[[data['Ue'][si,level,ilat,ilong],data['Ue'][si,level,ilat+1,ilong]],[data['Ue'][si,level,ilat,ilong+1],data['Ue'][si,level,ilat+1,ilong+1]]]
    ux0 = interp2(x0[0][i], x0[1][i], 
        datax[ilat,ilong], datax[ilat+1,ilong], datay[ilat,ilong], datay[ilat,ilong+1],
        datau[level,ilat,ilong], datau[level,ilat+1,ilong],datau[level,ilat,ilong+1],datau[level,ilat+1,ilong+1])
    ux1 = interp2(x0[0][i], x0[1][i], 
        datax[ilat,ilong], datax[ilat+1,ilong], datay[ilat,ilong], datay[ilat,ilong+1],
        datau[level+1,ilat,ilong], datau[level+1,ilat+1,ilong],datau[level+1,ilat,ilong+1],datau[level+1,ilat+1,ilong+1])
    vx0 = interp2(x0[0][i], x0[1][i], 
        datax[ilat,ilong], datax[ilat+1,ilong], datay[ilat,ilong], datay[ilat,ilong+1],
        datav[level,ilat,ilong], datav[level,ilat+1,ilong],datav[level,ilat,ilong+1],datav[level,ilat+1,ilong+1])
    vx1 = interp2(x0[0][i], x0[1][i], 
        datax[ilat,ilong], datax[ilat+1,ilong], datay[ilat,ilong], datay[ilat,ilong+1],
        datav[level+1,ilat,ilong], datav[level+1,ilat+1,ilong],datav[level+1,ilat,ilong+1],datav[level+1,ilat+1,ilong+1])
    # n1=datetime.datetime.now()
    # print(n1-n0)
    # ux0 =interp(pts)[0]
    # datac =[[data['Ve'][si,level,ilat,ilong],data['Ve'][si,level,ilat+1,ilong]],[data['Ve'][si,level,ilat,ilong+1],data['Ve'][si,level,ilat+1,ilong+1]]]
    # interp = RegularGridInterpolator((x, y), datac)
    # vx0 =interp(pts)[0]
    # datac =[[data['Ue'][si,level+1,ilat,ilong],data['Ue'][si,level+1,ilat+1,ilong]],[data['Ue'][si,level+1,ilat,ilong+1],data['Ue'][si,level+1,ilat+1,ilong+1]]]
    # interp = RegularGridInterpolator((x, y), datac)
    # ux1 =interp(pts)[0]
    # datac =[[data['Ve'][si,level+1,ilat,ilong],data['Ve'][si,level+1,ilat+1,ilong]],[data['Ve'][si,level+1,ilat,ilong+1],data['Ve'][si,level+1,ilat+1,ilong+1]]]
    # interp = RegularGridInterpolator((x, y), datac)
    # vx1 =interp(pts)[0]
    height0 = heights[level,ilat,ilong]
    height1 = heights[level+1,ilat,ilong]
    ux  = ux0*(height1-Adz[i])/(height1-height0)+ux1*(Adz[i]-height0)/(height1-height0)
    vx  = vx0*(height1-Adz[i])/(height1-height0)+vx1*(Adz[i]-height0)/(height1-height0)  
    if (ux>100):
        print('uxbig',ux,ux0,ux1,height0,height1, level, Adf.Cy.values[i],Adf.Cx.values[i])
    uxi.append(ux)
    vxi.append(vx)


# writing on the parallel OF case
# toSet = Adf.assign(Ux=Adf.Cx+1,Uy=0,Uz=0)
toSet = Adf.assign(Ux=uxi,Uy=vxi,Uz=0)#.dropna()
toSet.head()
U.setFieldFromDataFrame(toSet)

U.writeToCase(case,"0")



############# setfieldsdict ###############
import numpy as np
from pyproj import CRS
from pyproj import Transformer
import datetime
from scipy.interpolate import RegularGridInterpolator


def regrid2(data, factor):
    out_x = int(data.shape[0]*factor)
    out_y = int(data.shape[1]*factor)
    m = max(data.shape[0], data.shape[1])
    y = np.linspace(0, 1.0/m, data.shape[0])
    x = np.linspace(0, 1.0/m, data.shape[1])
    interpolating_function = RegularGridInterpolator((y, x), data)

    yv, xv = np.meshgrid(np.linspace(0, 1.0/m, out_y), np.linspace(0, 1.0/m, out_x))

    return interpolating_function((xv, yv))

def regrid3(data, factor):
    # make the data grid more focus (factor more than 1) or blur (factor less than 1)
    out_x = int(data.shape[0]*factor)
    out_y = int(data.shape[1]*factor)
    out_z = int(data.shape[2]*factor)
    m = max(data.shape[0], data.shape[1], data.shape[2])
    x = np.linspace(0, 1.0/m, data.shape[0])
    y = np.linspace(0, 1.0/m, data.shape[1])
    z = np.linspace(0, 1.0/m, data.shape[2])
    interpolating_function = RegularGridInterpolator((x, y, z), data)

    zv, yv, xv = np.meshgrid(np.linspace(0, 1.0/m, out_y), np.linspace(0, 1.0/m, out_x), np.linspace(0, 1.0/m, out_z))
    return interpolating_function((yv, zv, xv))

miny=195100
maxy=209800
minx=737000
maxx=752000
minz=-10
maxz=2490

miny=193100 # minx - blockmesh_dx
maxy=222800 # maxx + blockmesh_dx
minx=735000
maxx=754000
minz=-60
maxz=2690

# miny=195100
# maxy=196900
# minx=737000
# maxx=738800
# minz=-10
# maxz=40

filenc='/data4bk/nirb/wrf12062023.npz' # 2023061318gmt
data = np.load(filenc)

# print(data.files)
# print(data['xlat'], data['xlong'])

# print the times so we can choose
print ('GMT please add 3 hour')
for i in range(len(data['Times'])):
    s=''
    for j in range(len(data['Times'][i][:])):
        s+=str(data['Times'][i][j])[2]
    print(i,s)
si = 36 # chosen time

heights = data['gpt_hgt_M'][si,:,:,:] # altitude
# for i in range(heights.shape[0]):
    # heights[i,:,:]+=data['hgt'] # ground elevation

datau = data['Ue'][si,:,:,:]
datav = data['Ve'][si,:,:,:]
dataw = data['W'][si,:,:,:]
datat = data['Temp'][si,:,:,:]

factor = 2
rdatau = regrid3(datau, factor)
rdatav = regrid3(datav, factor)
rdataw = regrid3(dataw, factor)
rdatat = regrid3(datat, factor)

crs = CRS.from_epsg(4326) #6991
crs.to_epsg()
crs = CRS.from_proj4("+proj=tmerc +lat_0=31.7343936111111 +lon_0=35.2045169444445 +k=1.0000067 +x_0=219529.584 +y_0=626907.39 +ellps=GRS80 +towgs84=-24.002400,-17.103200,-17.844400,-0.33007,-1.852690,1.669690,5.424800 +units=m +no_defs")
transformer = Transformer.from_crs("EPSG:4326", "EPSG:6991")
# transformer.transform(220080, 634451)
# x2=transformer.transform(32.727616, 35.103320)
xlatlon=np.zeros_like(data['xlat'])
ylatlon=np.zeros_like(data['xlat'])
for i in range(xlatlon.shape[0]):
    if i % 100 == 0:
        print('iy',i,xlatlon.shape[0], datetime.datetime.now())
    for j in range(xlatlon.shape[1]):
        x2=transformer.transform(data['xlat'][i,j], data['xlong'][i,j])
        xlatlon[i,j]=x2[1]
        ylatlon[i,j]=x2[0]

rdatax = regrid2(xlatlon, factor)
rdatay = regrid2(ylatlon, factor)

rheights = regrid3(heights, factor)

    
outputdir = r'/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200a2z0'

lines = []
lines.append('// built with wrf2of.py for '+outputdir)
lines.append("/*--------------------------------*- C++ -*----------------------------------*\'")
lines.append('| =========                |                                                  |')
lines.append('| \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox            |')
lines.append('|  \    /   O peration     | Version:  10.0                                   |')
lines.append('|   \  /    A nd           | Web:      www.OpenFOAM.org                       |')
lines.append('|    \/     M anipulation  |                                                  |')
lines.append('\*---------------------------------------------------------------------------*/')
lines.append('FoamFile')
lines.append('{')
lines.append('    version     10.0;')
lines.append('    format      ascii;')
lines.append('    class       dictionary;')
lines.append('    location    "system";')
lines.append('    object      setFieldsDict;')
lines.append('}')
lines.append('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //')
lines.append('')
lines.append('defaultFieldValues')
lines.append('(')
lines.append('    volScalarFieldValue T 345')
lines.append(');')
lines.append('')
lines.append('regions')
lines.append('(')
shp1=rdatau.shape[0]
shp2=rdatau.shape[1]
shp3=rdatau.shape[2]
for i in range(shp2): # x
    print('i:',i,shp2)
    for j in range(shp3): # y
        celly=rdatay[i,j]
        cellx=rdatax[i,j]
        for k in range(shp1): # z
            cellz=rheights[k,i,j]
            if celly<=maxy and celly>=miny and cellx>=minx and cellx<=maxx and cellz>=minz and cellz<=maxz:
                T = rdatat[k,i,j]
                Ux = rdatau[k,i,j]
                Uy = rdatav[k,i,j]
                Uz = rdataw[k,i,j]
                if i==0:
                    dxp = rdatax[i+1,j]-rdatax[i,j]
                    dxm = rdatax[i,j]-minx
                elif i==shp2:
                    dxm = rdatax[i,j]-rdatax[i-1,j]
                    dxp = maxx-rdatax[i,j]
                else:
                    dxm = rdatax[i,j]-rdatax[i-1,j]
                    dxp = rdatax[i+1,j]-rdatax[i,j]
    
                if j==0:
                    dyp = rdatay[i,j+1]-rdatay[i,j]
                    dym = rdatay[i,j] - miny
                elif j==shp3:
                    dym = rdatay[i,j]-rdatay[i,j-1]
                    dyp = maxy-rdatay[i,j]
                else:
                    dym = rdatay[i,j]-rdatay[i,j-1]
                    dyp = rdatay[i,j+1]-rdatay[i,j]
                                
                if k==0:
                    dzp = rheights[k+1,i,j]-rheights[k,i,j]
                    dzm = rheights[k,i,j] - minz
                elif k==shp1:
                    dzm = rheights[k,i,j]-rheights[k-1,i,j]
                    dzp = maxz - rheights[k,i,j] 
                else:
                    dzm = rheights[k,i,j]-rheights[k-1,i,j]
                    dzp = rheights[k+1,i,j]-rheights[k,i,j]
    
                lines.append('boxToCell {box ('+str(celly-dym/2)+' '+str(cellx-dxm/2)+' '+str(cellz-dzm/2)+') ('+str(celly+dyp/2)+' '+str(cellx+dxp/2)+' '+str(cellz+dzp/2)+'); fieldValues (volScalarFieldValue T '+str(T)+' volVectorFieldValue U ('+str(Ux)+' '+str(Uy)+' '+str(Uz)+'));}')
                lines.append('boxToFace {box ('+str(celly-dym/2)+' '+str(cellx-dxm/2)+' '+str(cellz-dzm/2)+') ('+str(celly+dyp/2)+' '+str(cellx+dxp/2)+' '+str(cellz+dzp/2)+'); fieldValues (volScalarFieldValue T '+str(T)+' volVectorFieldValue U ('+str(Ux)+' '+str(Uy)+' '+str(Uz)+'));}')
    
lines.append(');')

with open(outputdir+'/system/setFieldsDict', 'w') as f:
    f.write('\n'.join(lines))
