import geopandas 
from shapely.geometry import  MultiLineString,LineString
import scipy
from scipy.interpolate import griddata
from numpy import array,cross,sqrt
import numpy 
import numpy as np
import matplotlib.pyplot as plt 
from itertools import product
import argparse
import sys

import glob
import os

#def filterShapefilesWithOgr2Ogr(inputFolder, outputFolder, attributesToSave):
def filterShapefilesWithOgr2Ogr(inputFolder, outputFolder,**kwargs):

    LowerLeft_x=kwargs["LowerLeft_x"]
    LowerLeft_y=kwargs["LowerLeft_y"]
    UpperRight_x=kwargs["UpperRight_x"]
    UpperRight_y=kwargs["UpperRight_y"]

    # inputFolder: the folder where the original shapefiles are located
    # outputFolder: the folder where the filtered shapefiles are to be stored
    # attributesToSave: string of attributes that shall be copied into the output, e.g. 'ID, location', length'
    # Example: ogr2ogr -clipsrc  $LOWERLEFT1 $LOWERLEFT2 $UPPERRIGHT1 $UPPERRIGHT2 $DESTDIR/$newfile $file
    # ,LowerLeft_x=216325,LowerLeft_y=629804,UpperRight_x=221802,UpperRight_y=633448 # 17500
    import subprocess
    # traverse through the input folder
    os.chdir(inputFolder)
    for filename in glob.glob('*.shp'):
        # filter each shapefile from the input folder and save it at output folder
        print('filtering ' + filename)
        #subprocess.call(["ogr2ogr", "-f", "ESRI Shapefile", "-select", attributesToSave, outputFolder + '/' + filename, filename])
        subprocess.call(["ogr2ogr", "-f", "ESRI Shapefile", "-clipsrc", LowerLeft_x, LowerLeft_y, UpperRight_x, UpperRight_y, outputFolder + '/' + filename, filename])
    return
    #ogr2ogr -clipsrc  $LOWERLEFT1 $LOWERLEFT2 $UPPERRIGHT1 $UPPERRIGHT2 $DESTDIR/$newfile $file	


def ReadTopo(filename):
                """
                        Reads the requested file.
                        with the format

                        and returns an array of z(x,y)
                        Example: TOPO
                """
                parser = lambda l: [x for x in l.split(" ") if x != '']
                f = open(filename)

                l = f.readline()
                xmin, xmax, ymin, ymax = map(float,parser(l) )
                l = f.readline()
                nx,ny                  = map(int,parser(l))

#                X = [ xmin + i*(xmax-xmin)/(nx-1) for i in range(0,nx)]   #last grid point: xmax
                X = [ xmin + i*(xmax-xmin)/(nx) for i in range(0,nx)]      #last grid point: xmax-dx
                Y = [ ymin + i*(ymax-ymin)/(ny) for i in range(0,ny)]
                i = 0
                j = 0

                Lines = f.read()

                Zv = map(float,parser(Lines))
                Z = []
                k = 0;
                for i in range(0,len(X)):
                        t = []
                        for j in range(0,len(Y)):
                                if Zv[k]<=0:
                                        Zv[k] =-1.0
                                t.append(Zv[k])
                                k += 1
                        Z.append(t)

                return X,Y,Z,(xmin, xmax, ymin, ymax)

def extrapolate_nans(x, y, v):
    '''  
    Extrapolate the NaNs or masked values in a grid INPLACE using nearest
    value.

    .. warning:: Replaces the NaN or masked values of the original array!

    Parameters:

    * x, y : 1D arrays
        Arrays with the x and y coordinates of the data points.
    * v : 1D array
        Array with the scalar value assigned to the data points.

    Returns:

    * v : 1D array
        The array with NaNs or masked values extrapolated.
    '''

    if np.ma.is_masked(v):
        nans = v.mask
    else:
        nans = np.isnan(v)
    notnans = np.logical_not(nans)
    v[nans] = scipy.interpolate.griddata((x[notnans], y[notnans]), v[notnans],
        (x[nans], y[nans]), method='nearest').ravel()
    return v

if __name__ == "__main__":
###########################################   INIT stdin TRIAL NO. AND RELEASE NO.
    print(len(sys.argv))
    if (len(sys.argv) != 2):
        print "Error: you should provide a *.shp filename"
        exit()
    else:
        fname = sys.argv[1]
    
    '''    
    inputFolder="/home/eyal/tmp/BNTL_MALE_ARZI/BNTL_MALE_ARZI"
    outputFolder=os.getcwd()+"/new"
    kwargs = {'LowerLeft_x': '150000.0', 'LowerLeft_y': '558000.0', 'UpperRight_x': '160000.0', 'UpperRight_y': '568000.0'}
    filterShapefilesWithOgr2Ogr(inputFolder, outputFolder,**kwargs )
    '''
    
    dxdy=80.0
    	#gpandas = geopandas.read_file('JERU-CONTOUR.shp') # fnas is a shape file. 
    	#gpandas = geopandas.read_file('JERU-10x10-CONTOUR.shp') # fnas is a shape file. 
    	#gpandas = geopandas.read_file('ZEELIM-CONTOUR.shp') # fnas is a shape file. 
    	#gpandas = geopandas.read_file('MALA-CONTOUR.shp') # fnas is a shape file. 
    gpandas = geopandas.read_file(fname) # fnas is a shape file. 
    	
     # 2. Convert contour map to regular height map. 
    	# 2.1 get boundaries 
    
    xmin    = gpandas['geometry'].bounds['minx'].min()
    xmax    = gpandas['geometry'].bounds['maxx'].max()
    
    ymin    = gpandas['geometry'].bounds['miny'].min()
    ymax    = gpandas['geometry'].bounds['maxy'].max()
     
    Nx=int(((xmax-xmin)/dxdy) )
    Ny=int(((ymax-ymin)/dxdy ))
    
    
    print("Mesh boundaries x=(%s,%s) ; y=(%s,%s); N=(%s,%s)" % (xmin,xmax,ymin,ymax,Nx,Ny))
    dx=(xmax-xmin)/(Nx)
    dy=(ymax-ymin)/(Ny)
    print("Mesh increments: D=(%s,%s); N=(%s,%s)" % (dx,dy,Nx,Ny))
    		# 2.2 build the mesh. 
    grid_x, grid_y = numpy.mgrid[xmin:xmax:dxdy,ymin:ymax:dxdy]
    print grid_x
    	# 3. Get the points from the geom
    Height = []
    XY     = []
    
    for i,line in enumerate(gpandas.iterrows()):
    		if isinstance(line[1]['geometry'	],LineString): 
    			linecoords = [x for x in line[1]['geometry'].coords]
    			lineheight = [line[1]['HEIGHT']]*len(linecoords)
    			XY	  += linecoords
    			Height    += lineheight
    		else: 
    			for ll in line[1]['geometry']:
    				linecoords = [x for x in ll.coords]
    				lineheight = [line[1]['HEIGHT']]*len(linecoords)
    				XY	  += linecoords
    				Height    += lineheight
    
    
    grid_z2 = griddata(XY, Height, (grid_x, grid_y),method='cubic')
    
    extrapolate_nans(grid_x,grid_y,grid_z2)
    
    	#print grid_z2
    	#plt.contour(grid_x,grid_y,grid_z2)
    	#plt.show()
    
    Gx = grid_x.reshape(grid_x.shape[0]*grid_x.shape[1],1)
    Gy = grid_y.reshape(grid_y.shape[0]*grid_y.shape[1],1)
    Gz = grid_z2.reshape(grid_z2.shape[0]*grid_z2.shape[1],1)
    
    GG = numpy.concatenate([Gx,Gy,Gz],axis=1)
    
    	#I = numpy.isnan(GG[:,2])
    	#plt.tricontour(GG[~I,0],GG[~I,1],GG[~I,2])
    
    with open('TOPO','w') as f:
    		f.write("%s %s %s %s\n" %(xmin,xmax,ymin,ymax))
    		#f.write("%s %s %s %s\n" %(0.0,float(xmax-xmin),0.0,float(ymax-ymin)))
    		f.write("%s %s\n" %(int(Nx),int(Ny)))
    		#numpy.savetxt(f,grid_z2.T) # needs transposed .T? (think not)
    	
    		#print grid_z2.shape[0],grid_z2.shape[1]
    		#import pdb
    		#pdb.set_trace()
    
    		for i in range(0,Nx):
    			for j in range(0,Ny):
    				#print i*int(Nx)+j
    				#f.write("0.0 ")
    				f.write("%s " %(float(grid_z2[j,i])))  # for my LSM so that it will be consistent with haifa and tlv
    			f.write("\n")
    		#print grid_z2[3,10],grid_z2[3,50], grid_z2[3,70]
    	#numpy.savetxt("Topo.txt",GG)
    print GG[:,2]
