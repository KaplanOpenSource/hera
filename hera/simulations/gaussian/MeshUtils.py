import pandas
try:
	import xarray
except ModuleNotFoundError:
	pass

import numpy
from unum.units import *
from scipy.special import erf

class GaussianToMesh(object):
    """
        Gets a list of gaussians as DataFrame
        with the columns:

        x,sigmaX,sigmaY

        and evaluates the gaussians on a grid.

        A regular grid is constructed  according to the parameters of the gaussians.
        The columns names can be different.

    """

    _sigmaXName = None
    _sigmaYName = None
    _dxdy = None

    @property
    def sigmaXName(self):
        return self._sigmaXName

    @property
    def sigmaYName(self):
        return self._sigmaYName

    @sigmaXName.setter
    def sigmaXName(self,value):
        self._sigmaXName = value

    @sigmaYName.setter
    def sigmaYName(self, value):
        self._sigmaYName = value


    @property
    def sigmashifts(self):
        return 6.

    @property 
    def dxdy(self): 
        return self._dxdy 
    @dxdy.setter 
    def dxdy(self,value):
        self._dxdy = value 

    def __init__(self,dxdy = None):
        self.sigmaXName = "sigmaXCorrected"
        self.sigmaYName = "sigmaXCorrected"
        self.QName      = "QName"	
        self._dxdy      = dxdy
    def _defineCoordinates(self,gaussians):
        """
            The mesh coordinates will be defined by the [min X - 6*gaussian std, max X + 6*gaussian std],[min Y - 6*gaussian std, max Y + 6*gaussian std]

            The number of points will be determined by domain width / [min std/5.].


        :return:
            the x,y coordinates.
        """
        # fine the bounds.
        if isinstance(gaussians,pandas.core.series.Series):
            maxX = gaussians['x'] + self.sigmashifts * gaussians[self.sigmaXName]
            minX = gaussians['x'] - self.sigmashifts * gaussians[self.sigmaXName]
            maxY = gaussians['y'] + self.sigmashifts * gaussians[self.sigmaYName]
            minY = gaussians['y'] - self.sigmashifts * gaussians[self.sigmaYName]
        else:
            maxXID = gaussians[gaussians['x']==gaussians['x'].max()].index[0]
            minXID = gaussians[gaussians['x'] == gaussians['x'].min()].index[0]
            maxYID = gaussians[gaussians['y'] == gaussians['y'].max()].index[0]
            minYID = gaussians[gaussians['y'] == gaussians['y'].min()].index[0]

            maxX = gaussians.iloc[maxXID]['x'] + self.sigmashifts * gaussians.iloc[maxXID][self.sigmaXName]
            minX = gaussians.iloc[minXID]['x'] - self.sigmashifts * gaussians.iloc[minXID][self.sigmaXName]
            maxY = gaussians.iloc[maxYID]['y'] + self.sigmashifts * gaussians.iloc[maxYID][self.sigmaYName]
            minY = gaussians.iloc[minYID]['y'] - self.sigmashifts * gaussians.iloc[minYID][self.sigmaYName]

        if self._dxdy is None: 
           dx   = gaussians[self.sigmaXName].min()/5.
           dy   = gaussians[self.sigmaYName].min()/5.
        else: 
           dx   = self._dxdy
           dy   = self._dxdy

        xCoord = numpy.arange(minX, maxX, dx)
        yCoord = numpy.arange(minY, maxY, dy)

        return xCoord,yCoord

    def gaussianToMesh(self,gaussians,groupby=None,QName="Q",QUnits=kg):
        """
		Spreads the list of gaussians on the regular mesh. 
		if groupby is None then sum all the gaussians, 
		if it has the name a column, return the sum of the gaussians 
		grouped by this fields. 
        """
        xCoord, yCoord = self._defineCoordinates(gaussians)
        if groupby is None: 
              ret = None

              if isinstance(gaussians, pandas.core.series.Series):
                 ret = self._gaussianToMesh(gaussians, xCoord, yCoord,QName=QName,QUnits=QUnits)
              else:
                 for lid,igauss in gaussians.iterrows():
                     ret = self._gaussianToMesh(igauss,xCoord,yCoord,addTo=ret,QName=QName,QUnits=QUnits)
      
        else: 
              groupedDim = gaussians[groupby].unique()
              ret = xarray.DataArray(numpy.zeros([xCoord.size,yCoord.size,groupedDim.size]),dims=["x","y",groupby],coords={"x":xCoord[:-1],"y":yCoord[:-1],groupby:groupedDim})
          
              for curDim,group in gaussians.groupby(groupby):
                     for lid,igauss in gaussians.iterrows():
                           curdata = self._gaussianToMesh(igauss,xCoord,yCoord,addTo=None,QName=QName,QUnits=QUnits)                        
                           ret.loc[{groupby:curDim}] += curdata


        return ret


    def _gaussianToMesh(self,gaussian,xCoords,yCoords,addTo=None,QName="Q",QUnits=kg):
        """
            Return a single gaussian that was spread on the mesh.

        :param gaussian:
                A pandas with the columns of a gaussian: x,y and the sigma columns
        :param addTo:
                If not None, add the value of the current gaussian to the addTo xarray.
        :param xCoords:
                A numpy array of the x coordinates.
        :param yCoords:
                A numpy array of the y coordinates.
        :return:
            xarray dataset.
        """
        x = gaussian['x']
        sigx = gaussian[self.sigmaXName]

        y = gaussian['y']
        sigy = gaussian[self.sigmaYName]

        twoPiFactor = (2*numpy.pi)**0.5

        # 1. evaluate the value of the gaussian on the x coord.

        gaussX = xarray.DataArray(1/(twoPiFactor*sigx)*numpy.exp(-0.5*((xCoords-x)/sigx)**2), dims='x',coords={'x':xCoords})
        gaussY = xarray.DataArray(1/(twoPiFactor*sigy)*numpy.exp(-0.5*((yCoords-y)/sigy)**2), dims='y',coords={'y':yCoords})

        fullX, fullY = xarray.broadcast(gaussX, gaussY)

        ret = fullX*fullY*gaussian[QName].asNumber(QUnits)
        if addTo is not None:
            ret = addTo+ ret

        return ret


class GaussianIntegrationToMesh(GaussianToMesh):
    """
       This class calculates the concentration of meterial there is in the mesh by subtracting erf.
       This is more accurate than estimating the gaussian.  	
    """
    def _gaussianToMesh(self,gaussian,xCoords,yCoords,addTo=None,QName="Q",QUnits=kg):
        """
            Return the concentration a single gaussian that was spread on the mesh.

        :param gaussian:
                A pandas with the columns of a gaussian: x,y and the sigma columns
        :param addTo:
                If not None, add the value of the current gaussian to the addTo xarray.
        :param xCoords:
                A numpy array of the x coordinates.
        :param yCoords:
                A numpy array of the y coordinates.
        :return:
            xarray dataset.
        """
        x = gaussian['x']
        sigx = gaussian[self.sigmaXName]

        y = gaussian['y']
        sigy = gaussian[self.sigmaYName]

        TwoSqrt =(2)**(0.5)

        # 1. evaluate the value of the gaussian on the x coord.
        dx = numpy.diff(xCoords)
        dy = numpy.diff(yCoords)  

        gaussX = xarray.DataArray(erf( (xCoords-x)/(TwoSqrt*sigx)), dims='x',coords={'x':xCoords}).diff(dim="x")/dx/2.
        gaussY = xarray.DataArray(erf( (yCoords-y)/(TwoSqrt*sigy)), dims='y',coords={'y':yCoords}).diff(dim="y")/dy/2.

        fullX, fullY = xarray.broadcast(gaussX, gaussY)

        ret = fullX*fullY*gaussian[QName].asNumber(QUnits)
        if addTo is not None:
            ret = addTo+ ret
        return ret




gaussianToMesh = GaussianToMesh()

