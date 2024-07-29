#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 11:02:33 2024

@author: nirb
"""
from hera import toolkitHome
import numpy as np
import datetime
import pyproj
from pyproj import CRS
from pyproj import Transformer


# reading of case
tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM)
case = r'/data5/NOBACKUP/nirb/Simulations/Haifa/testcode'
U =   tk.OFObjectHome.readFieldFromCase("U",tk.FLOWTYPE_INCOMPRESSIBLE,case) #"/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2/")
Udf = tk.OFObjectHome.readFieldAsDataFrame("U", case)
print(Udf)
print(U.componentNames)
A = tk.getMesh(case)
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


heights = data['gpt_hgt_M'][si,:,:,:]
# for i in range(heights.shape[0]):
    # heights[i,:,:]+=data['hgt']
    
# ugrd=data['Ue'][si,:,:,:]
# vgrd=data['Ve'][si,:,:,:]
# wgrd=data['W'][si,:,:,:]



# crs = CRS.from_epsg(6991) #6991
# crs.to_epsg()
# crs = CRS.from_proj4("+proj=tmerc +lat_0=31.7343936111111 +lon_0=35.2045169444445 +k=1.0000067 +x_0=219529.584 +y_0=626907.39 +ellps=GRS80 +towgs84=-24.002400,-17.103200,-17.844400,-0.33007,-1.852690,1.669690,5.424800 +units=m +no_defs")
# transformer = Transformer.from_crs("EPSG:6991", "EPSG:4326")


# transfer coordinates to ITM, this is what OF wants
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


# interpolate the WRF coordinates to OF coordinates
# https://stackoverflow.com/questions/33259896/python-interpolation-2d-array-for-huge-arrays
from scipy import interpolate
# making interpolation for each level, the height can be different for each coordinate
fu=[]
fy=[]
for i in range(heights.shape[0]):
    fu.append(interpolate.interp2d(xlatlon,ylatlon,data['Ue'][si,13,:,:],kind='cubic'))
    fv.append(interpolate.interp2d(xlatlon,ylatlon,data['Ve'][si,13,:,:],kind='cubic'))
    
uxi = []
vxi = []
# f(716373,153153)[0]  #15,15
xi=[]
yi=[]
for i in range(len(Adf.Cx.values)):
    # using the middle coordinate, not so accurate, but good enough
    xi.append(np.argmin(np.abs(xlatlon[:,61]-Adf.Cy.values[i])))
    yi.append(np.argmin(np.abs(ylatlon[61,:]-Adf.Cx.values[i])))
for i in range (len(Adf.Cx.values)):
    level = np.argmin(np.abs(heights[:,xi[i],yi[i]]-Adf.Cz.values[i]))
    ??? calculate level and level+1
    uxi.append(fu[level](Adf.Cy.values[i],Adf.Cx.values[i])[0])
    vxi.append(fv[level](Adf.Cy.values[i],Adf.Cx.values[i])[0])
    
    
# I wanted to interpolated 3D, but I failed.    
# from scipy.interpolate import RegularGridInterpolator
# # not the best interpolation because the grid is not regular !!!!
# # xp = xlatlon[:,xlatlon.shape[1]//2] 
# # yp = ylatlon[ylatlon.shape[0]//2,:]
# # zp = 
# meshx = np.zeros_like(heights)
# meshy = np.zeros_like(heights)
# for i in range(heights.shape[0]):
#     meshx[i,:,:]=xlatlon
#     meshy[i,:,:]=ylatlon

# fu = RegularGridInterpolator((meshx,meshy,heights), data['Ue'][si,:,:,:], solver="slinear")
# fu = RegularGridInterpolator((meshx,meshy,heights), data['Ue'][si,:,:,:])
# pts = array([[13,15,15],[13,15,17]])
# print(fu(pts))










# writing on the parallel OF case
toSet = Adf.assign(Ux=Adf.Cx+1,Uy=0,Uz=-1)
toSet.head()
U.setFieldFromDataFrame(toSet)
U.writeToCase("test","0")

U.writeToCase(case,"0")
