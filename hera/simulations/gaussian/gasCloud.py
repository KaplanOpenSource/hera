import pandas
import numpy
from numpy import matlib
import math
from ...utils import *
import xarray
from scipy import special


class abstractGasCloud:

    def __init__(self, sourceQ, sourceHeight, initialCloudSize, sigmaType):
        """

        Parameters
        ----------
        sourceQ : unum, method
            If unum:
                The unit determine the release time.
                [mass] - Instantaneous
                [mass/time] - Continuous
            else
                Continuous (not implementaed yet.)

        sourceHeight : unum
        initialCloudSize : 3-touple unum, the sigmas in each axis.
        sigmaType : The sigma type, for example from Briggs, rural/urban.
        """
        self.sourceHeight = sourceHeight
        self.initialCloudSize = initialCloudSize
        self.sigmaType = sigmaType
        self.sourceQ = sourceQ



    @staticmethod
    def createGasCloud(sourceQ,sourceHeight,initialCloudSize,sigmaType):
        """
            Return the type of the release based on the units of Q
        Parameters
        ----------
        sourceQ : unum, method
            If unum:
                The unit determine the release time.
                [mass] - Instantaneous
                [mass/time] - Continuous
            else
                Continuous (not implementaed yet.)

        sourceHeight : unum
        initialCloudSize : 3-touple unum, the sigmas in each axis.

        Returns
        -------

        """
        try:
            sourceQ.asUnit(mg)
            instantaneous = True
        except:
            try:
                sourceQ.asUnit(mg/min)
                instantaneous = False
            except:
                raise ValueError("Must be mass or mass per time!")

        returnCls = instantaneousReleaseGasCloud if instantaneous else continuousReleaseGasCloud

        return returnCls(sourceQ=sourceQ,sourceHeight=sourceHeight,initialCloudSize=initialCloudSize,sigmaType=sigmaType)


    def _getTXterm(self, stability, u, xcoordRange, tcoordRange):
        """
        Parameters
        ----------
        initialCloudSize : 3-tuple of float/unum (default m)
            The initial cloud size (standard deviation) in the x,y and z dimensions.
        xcoordRange : Tuple in numpy.arange format, unitless.
        tcoordRange : Tuple in numpy.arange format, unitless.

        Returns
        -------
        The X component of the Gaussian concentration formula.
        """
        T, X = numpy.meshgrid(tcoordRange, xcoordRange, indexing='ij')
        sigmaX = self.sigmaType.getSigma(x=X, stability=stability, sigma0=self.initialCloudSize, units=False)['sigmaX']
        downwind = (1 / (numpy.sqrt(2 * numpy.pi) * sigmaX)) * numpy.exp((-(X - u * T) ** 2) / (2 * sigmaX ** 2))
        XR_downwind = xarray.DataArray(downwind, dims=("time", "x"), coords={"time": tcoordRange, "x": xcoordRange})
        return XR_downwind

    def _getXYterm(self, stability, xcoordRange, ycoordRange):
        """

        Parameters
        ----------
        initialCloudSize : 3-tuple of float/unum (default m)
            The initial cloud size (standard deviation) in the x,y and z dimensions.
        xcoordRange : Tuple in numpy.arange format, unitless.
        ycoordRange : Tuple in numpy.arange format, unitless.

        Returns
        -------
        The Y component of the Gaussian concentration formula.
        """
        X, Y = numpy.meshgrid(xcoordRange, ycoordRange, indexing='ij')
        sigmaY = self.sigmaType.getSigma(x=X, stability=stability, sigma0=self.initialCloudSize, units=False)['sigmaY']
        crosswind = (1 / (numpy.sqrt(2 * numpy.pi) * sigmaY)) * numpy.exp((-(Y - 0) ** 2) / (2 * sigmaY ** 2))
        XR_crosswind = xarray.DataArray(crosswind, dims=("x", "y"), coords={"x": xcoordRange, "y": ycoordRange})
        return XR_crosswind

    def _getXZterm(self, stability, inversion, xcoordRange, zcoordRange, numOfReflections):
        """

        Parameters
        ----------
        initialCloudSize : 3-tuple of float/unum (default m)
            The initial cloud size (standard deviation) in the x,y and z dimensions.
        xcoordRange : Tuple in numpy.arange format, unitless.
        zcoordRange : Tuple in numpy.arange format, unitless.
        numOfReflections : The number of reflections of the summation of the Z component.

        Returns
        -------
        The Z component of the Gaussian concentration formula.
        """
        sourceHeight = tonumber(self.sourceHeight, m)
        X, Z = numpy.meshgrid(xcoordRange, zcoordRange, indexing='ij')
        sigmaZ = self.sigmaType.getSigma(x=X, stability=stability, sigma0=self.initialCloudSize, units=False)['sigmaZ']

        nSum = numpy.arange(-numOfReflections, numOfReflections + 1, 1)
        vertical = 0
        for n in nSum:
            vertical += (1 / (numpy.sqrt(2 * numpy.pi) * sigmaZ)) * (
                        numpy.exp((-(Z -(sourceHeight + 2 * n * inversion)) ** 2) / (2 * sigmaZ ** 2)) +
                        numpy.exp((-(Z +(sourceHeight + 2 * n * inversion)) ** 2) / (2 * sigmaZ ** 2)))
        XR_vertical = xarray.DataArray(vertical, dims=("x", "z"), coords={"x": xcoordRange, "z": zcoordRange})
        return XR_vertical



    def _getTXDosage(self, stability, u, xcoordRange, tcoordRange):
        """
        Parameters
        ----------
        initialCloudSize : 3-tuple of float/unum (default m)
            The initial cloud size (standard deviation) in the x,y and z dimensions.
        xcoordRange : Tuple in numpy.arange format, unitless.
        tcoordRange : Tuple in numpy.arange format, unitless.

        Returns
        -------
        The X component of the Gaussian dosage formula.
        """
        T, X = numpy.meshgrid(tcoordRange, xcoordRange, indexing='ij')
        sigmaX = self.sigmaType.getSigma(x=X, stability=stability, sigma0=self.initialCloudSize, units=False)['sigmaX']
        downwind_erf = (1/(2*u))*(special.erf(X/(numpy.sqrt(2)*sigmaX))-special.erf((X-u*T)/(numpy.sqrt(2)*sigmaX)))
        XR_downwind_erf = xarray.DataArray(downwind_erf, dims=( "time", "x"), coords={"time": tcoordRange, "x": xcoordRange})
        return XR_downwind_erf


    def trapezoidal_integration(self, data, dim='time'):
        """
        Performs trapezoidal integration along the specified dimension of an xarray DataArray.

        Args:
            data: Input DataArray.
            dim: Dimension along which to perform integration (default: 'time').

        Returns:
            DataArray containing the integrated values.
        """

        # Calculate time differences
        dt = data[dim].shift({dim: -1}) - data[dim]
        # Handle potential edge cases (e.g., first or last time step)
        dt = dt.fillna(dt.mean(dim=dim))

        # Calculate trapezoidal areas
        areas = 0.5 * (data + data.shift({dim: -1})) * dt

        # Perform cumulative sum along the specified dimension
        integrated_data = areas.cumsum(dim=dim)

        return integrated_data


    def getConcentrationFromMinMaxRange_inst_noQ(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                                        dxdy=10*m, dz=1*m, dt=1*min, numOfReflections=3):

        xcoordRange = numpy.arange(tonumber(minx, m), tonumber(maxx, m), tonumber(dxdy,m))
        ycoordRange = numpy.arange(tonumber(miny,m),tonumber(maxy,m),tonumber(dxdy,m))
        zcoordRange = numpy.arange(tonumber(minz, m), tonumber(maxz, m), tonumber(dz,m))
        tcoordRange = numpy.arange(0,tonumber(timeSpan,s),tonumber(dt,s))

        stability = meteorology.stability
        u = tonumber(meteorology.u10, m/s)
        inversion = tonumber(meteorology.inversion, m)

        TX = self._getTXterm(stability=stability, u=u, xcoordRange=xcoordRange, tcoordRange=tcoordRange)
        XY = self._getXYterm(stability=stability, xcoordRange=xcoordRange, ycoordRange=ycoordRange)
        XZ = self._getXZterm(stability=stability, inversion=inversion, xcoordRange=xcoordRange, zcoordRange=zcoordRange, numOfReflections=numOfReflections)

        return TX*XY*XZ


    def getDosageFromMinMaxRange_inst_noQ(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                           dxdy=10*m, dz=1*m, dt=1*min, numOfReflections=3):
        stability = meteorology.stability
        u = tonumber(meteorology.u10, m/s)
        inversion = tonumber(meteorology.inversion, m)

        xcoordRange = numpy.arange(tonumber(minx, m), tonumber(maxx, m), tonumber(dxdy, m))
        ycoordRange = numpy.arange(tonumber(miny, m), tonumber(maxy, m), tonumber(dxdy, m))
        zcoordRange = numpy.arange(tonumber(minz, m), tonumber(maxz, m), tonumber(dz, m))
        tcoordRange = numpy.arange(0, tonumber(timeSpan, s), tonumber(dt, s))

        TX = self._getTXDosage(stability=stability, u=u, xcoordRange=xcoordRange, tcoordRange=tcoordRange)
        XY = self._getXYterm(stability=meteorology.stability, xcoordRange=xcoordRange, ycoordRange=ycoordRange)
        XZ = self._getXZterm(stability=stability, inversion=inversion, xcoordRange=xcoordRange, zcoordRange=zcoordRange, numOfReflections=numOfReflections)

        return TX*XY*XZ

    def getDosageFromMinMaxRange_inst_NoERF_noQ(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                           dxdy=10*m, dz=1*m, dt=1*min, numOfReflections=3):
        C_without_Q = self.getConcentrationFromMinMaxRange_inst_noQ(meteorology=meteorology, minx=minx, miny=miny, minz=minz,
                                                                    maxx=maxx, maxy=maxy, maxz=maxz, timeSpan=timeSpan,
                                                                    dxdy=dxdy, dz=dz, dt=dt,numOfReflections=numOfReflections)
        D_without_Q = self.trapezoidal_integration(data=C_without_Q)

        return D_without_Q



class instantaneousReleaseGasCloud(abstractGasCloud):


    def getConcentrationFromMinMaxRange_inst(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                                        dxdy=10*m, dz=1*m, dt=1*min, numOfReflections=3):
        C_without_Q = self.getConcentrationFromMinMaxRange_inst_noQ(meteorology=meteorology, minx=minx, miny=miny, minz=minz,
                                                                    maxx=maxx, maxy=maxy, maxz=maxz, timeSpan=timeSpan,
                                                                    dxdy=dxdy, dz=dz, dt=dt, numOfReflections=numOfReflections)

        return tonumber(self.sourceQ, mg)*C_without_Q

    def getConcentrationFromDomain(self, meteorology, xcoord="x", ycoord="y"):
        xcoord = domain.coords[xcoord]
        ycoord = domain.coords[ycoord]

        TX = self._getTXterm()
        XY = self._getXYterm()
        XZ = self._getXZterm()

        pass

    def getDosageFromMinMaxRange_inst(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                           dxdy=10*m, dz=1*m, dt=1*min, numOfReflections=3):
        D_without_Q = self.getDosageFromMinMaxRange_inst_noQ(meteorology=meteorology, minx=minx, miny=miny, minz=minz, maxx=maxx,
                                                             maxy=maxy, maxz=maxz, timeSpan=timeSpan, dxdy=dxdy, dz=dz,
                                                             dt=dt, numOfReflections=numOfReflections)

        return tonumber(self.sourceQ, mg)*D_without_Q/60 #devide by 60 go get units of mg*min/m^3 (rather than mg*sec/m^3)



    def getDosageFromMinMaxRange_inst_NoERF(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                           dxdy=10*m, dz=1*m, dt=1*min, numOfReflections=3):

        D_without_Q = self.getDosageFromMinMaxRange_inst_NoERF_noQ(meteorology=meteorology, minx=minx, miny=miny, minz=minz,
                                                                   maxx=maxx, maxy=maxy,maxz=maxz, timeSpan=timeSpan,
                                                                   dxdy=dxdy, dz=dz, dt=dt, numOfReflections=numOfReflections)

        return tonumber(self.sourceQ, mg)*D_without_Q/60 #devide by 60 go get units of mg*min/m^3 (rather than mg*sec/m^3)



class continuousReleaseGasCloud(abstractGasCloud):

    def getConcentrationFromMinMaxRange_cont(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                                        dxdy=10*m, dz=1*m, dt=1*min, numOfReflections=3):
        """
        Returns
        -------
        An xarray of concentrations at every grid-poit, which is the dosage of the instantaneous release,
        since we assume the release rate is constant.
        Here we take the concentration xarray that was claculated using the error function (erf).
        """
        C_without_Q = self.getDosageFromMinMaxRange_inst_noQ(meteorology=meteorology, minx=minx, miny=miny, minz=minz, maxx=maxx,
                                                        maxy=maxy, maxz=maxz, timeSpan=timeSpan, dxdy=dxdy, dz=dz, dt=dt,
                                                        numOfReflections=numOfReflections)
        return tonumber(self.sourceQ, mg/s)*C_without_Q


    def getConcentrationFromMinMaxRange_cont_NoERF(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                                        dxdy=10*m, dz=1*m, dt=1*min, numOfReflections=3):
        """
        Returns
        -------
        An xarray of concentrations at every grid-poit, which is the dosage of the instantaneous release,
        since we assume the release rate is constant.
        Here we take the concentration xarray that was claculated without the error function (erf).
        """
        C_without_Q = self.getDosageFromMinMaxRange_inst_NoERF_noQ(meteorology=meteorology, minx=minx, miny=miny, minz=minz,
                                                          maxx=maxx, maxy=maxy, maxz=maxz, timeSpan=timeSpan, dxdy=dxdy,
                                                          dz=dz, dt=dt, numOfReflections=numOfReflections)
        return tonumber(self.sourceQ, mg/s)*C_without_Q


    def getDosageFromMinMaxRange_cont_NoERF(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                                        dxdy=10*m, dz=1*m, dt=1*min, numOfReflections=3):

        C_without_Q = self.getConcentrationFromMinMaxRange_cont(meteorology=meteorology, minx=minx, miny=miny, minz=minz,
                                                           maxx=maxx, maxy=maxy, maxz=maxz, timeSpan=timeSpan, dxdy=dxdy,
                                                           dz=dz, dt=dt, numOfReflections=numOfReflections)
        D_without_Q = self.trapezoidal_integration(data=C_without_Q)

        return tonumber(self.sourceQ, mg/s) * D_without_Q / 60  # devide by 60 go get units of mg*min/m^3 (rather than mg*sec/m^3)

    def getDosageFromMinMaxRange_cont_doubleNoERF(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                                        dxdy=10*m, dz=1*m, dt=1*min, numOfReflections=3):

        C_without_Q = self.getConcentrationFromMinMaxRange_cont_NoERF(meteorology=meteorology, minx=minx, miny=miny, minz=minz,
                                                           maxx=maxx, maxy=maxy, maxz=maxz, timeSpan=timeSpan, dxdy=dxdy,
                                                           dz=dz, dt=dt, numOfReflections=numOfReflections)
        D_without_Q = self.trapezoidal_integration(data=C_without_Q)

        return tonumber(self.sourceQ, mg/s) * D_without_Q / 60  # devide by 60 go get units of mg*min/m^3 (rather than mg*sec/m^3)



#-------------------------- Yehuda's Code For Convolution --------------------------

class Continuous(object):
    dt = None
    M  = None
    Timekernel = None
    _FullKernel = None

    def __init__(self,dt,kernelsize,timetofinish=10*min):
        """
        Time to finish.
        the time (min) it take to reach 0.1.

        Now define:
        $\dot{Q} = Aexp^{-\alpha t}$
        Therefore
        $\int_0^{\infty} \dot{Q} = Q \leftarrow A = Q\cdot \alpha

        Therefore in one time step the amount that is released is
        $Q(t,t+dt) = \int_t^{t+\delta t} Aexp^{-\alpha t} = \frac{A}{\alpha}\left[exp^{-alpha t} - exp^{-alpha (t+\delta t)}\right]$
        so the amount is $Q\left[exp^{-alpha t} - exp^{-alpha (t+\delta t)}$
        the rate will be $\frac{Q(t,t+dt)}{\delta t}$

        and that is what we should put in the kernel.
        """
        dt = tounum(dt,min)
        self.dt = tonumber(dt,min)
        self.kernelsize  = kernelsize


        timetofinish = tounum(timetofinish,min)
        alpha = numpy.log(0.1)/(-timetofinish.asNumber(min))

        # build the kernel.
        ts = numpy.arange(kernelsize,-1,-1)*dt.asNumber(min)

        # the kernel
        self.Timekernel = (numpy.exp(-alpha*ts[1:]) - numpy.exp(-alpha*ts[:-1]))/dt.asNumber(min)
        self.Timekernel = self.Timekernel.reshape([kernelsize,1,1,1])


    def _convolve(self,data,axis,FullKernel):
        # convolve.
        return (data*FullKernel[FullKernel.shape[0]-data.shape[0]:,:,:]).sum(axis=axis)*self.dt

    def calc(self,data):
        # build the kernel.
        FullKernel = numpy.tile(self.Timekernel,[1,data.x.size,data.y.size,data.z.size])
        return data.rolling(datetime=self.kernelsize,min_periods=1).reduce(self._convolve,FullKernel=FullKernel)


#-------------------------- End Of Yehuda's Code For Convolution --------------------------









    # def instantaneousReleaseGasCloud(self,minx, miny, minz, maxx, maxy, maxz, timeSpan,
    #                                  numOfIteration=None, xcoord="x", ycoord="y", zcoord="z", tcoord="time"):
    #
    #     if numOfIteration == None:
    #         numOfIteration = 3
    #     else:
    #         numOfIteration = numOfIteration
    #
    #     xcoord = domain.coords[xcoord]
    #     ycoord = domain.coords[ycoord]
    #     zcoord = domain.coords[zcoord]
    #     tcoord = domain.coords[tcoord]
    #
    #     u = self.meteorolgy.getWindVeclocity(self.sourceHeight)
    #
    #     xcoord = np.arange(tonumber(minx[0],m),tonumber(maxx[1],m),1)
    #     ycoord = np.arange(tonumber(miny[0],m),tonumber(maxy[1],m),1)
    #     zcoord = np.arange(tonumber(minz[0],m),tonumber(maxz[1],m),1)
    #     tcoord = np.arange(0,tonumber(timeSpan,s),1)
    #
    #     X, T = np.meshgrid(xcoord, tcoord, indexing='ij')
    #     sigmaX = 1
    #     downwind = (1 / (np.sqrt(2 * np.pi) * sigmaX)) * np.exp(-(X - u * T) ** 2 / (2 * sigmaX ** 2))
    #     XR_downwind = xr.DataArray(downwind, dims=("x", "time"), coords={"x": xcoord, "time": tcoord})
    #
    #     X, Y = np.meshgrid(xcoord, ycoord, indexing='ij')
    #     sigmaY = 1
    #     crosswind = (1 / (np.sqrt(2 * np.pi) * sigmaY)) * np.exp(-(Y - 0) ** 2 / (2 * sigmaX ** 2))
    #     XR_crosswind = xr.DataArray(crosswind, dims=("x", "y"), coords={"x": xcoord, "y": ycoord})
    #
    #     X, Z = np.meshgrid(xcoord, zcoord, indexing='ij')
    #     sigmaZ = 1
    #
    #     nSum = np.arange(-numOfIteration, numOfIteration + 1, 1)
    #     vertical = 0
    #     for n in numOfIteration:
    #         vertical += (1 / (np.sqrt(2 * np.pi) * sigmaZ)) * (
    #                     np.exp(-(Z -(self.sourceHeight + 2 * n * self.inversion)) ** 2 / (2 * sigmaX ** 2)) +
    #                     np.exp(-(Z +(self.sourceHeight + 2 * n * self.inversion)) ** 2 / (2 * sigmaX ** 2)))
    #     XR_vertical = xr.DataArray(vertical, dims=("x", "z"), coords={"x": xcoord, "z": zcoord})
    #
    #     down, cross, vert = xr.broadcast(XR_downwind, XR_crosswind, XR_vertical)
    #
    #     concentrations = down * cross * vert
    #
    #     return concentrations


#
# import pandas as pd
# import numpy as np
# import math
# from ...utils import *
# import xarray as xr
#
# class gasCloud:
#
#     @property
#     def sourceHeight(self):
#         return self.source.height
#
#     @property
#     def Q(self):
#         return self.source.Q
#
#     @property
#     def initialCloudSize(self):
#         return self.source.initialCloudSize
#
#     def __init__(self, sigmaType, meteorology, source):
#
#         self.source = source
#         self.meteorology = meteorology
#         self.sigmaType = sigmaType
#
#
#         # self.sourceHeight = source.getHeight #needs correction
#         # self.Q = source.getQ #needs correction
#         # self.u = meteorology.getWindVeclocity(height) #needs correction
#         # self.inversion = meteorology.inversionHeight  # needs correction
#         # self.stability = meteorology.stability
#
#         pass
#
#
#
#
#
#     def getConcentrationDomainRange(self):
#         xRange = [0 * m, 500 * m]
#         yRange = [-50 * m, 50 * m]
#         zRange = [-20 * m, 20 * m]
#         timeSpan = 5 * min
#
#         xcoord = np.arange(tonumber(xRange[0],m),tonumber(xRange[1],m),1)
#         ycoord = np.arange(tonumber(yRange[0],m),tonumber(yRange[1],m),1)
#         zcoord = np.arange(tonumber(zRange[0],m),tonumber(zRange[1],m),1)
#         tcoord = np.arange(0,tonumber(timeSpan,s),1)
#
#     def getConcentrationXarray(self, domain, time, numOfIteration=None, xcoord="x", ycoord="y", zcoord="z", tcoord="time"):
#         """
#
#         Parameters
#         ----------
#         domain = Three arrays of [min, max] for each coordinate, with units.
#         time = length of time, with units
#
#         Returns
#         -------
#         The concentration at all grid points of the domain, at time "time".
#         """
#
#         if numOfIteration == None:
#             numOfIteration = 3
#         else:
#             numOfIteration = numOfIteration
#
#
#         #create the mesh grid from the coordinats.
#
#         xcoord = domain.coords[xcoord]
#         ycoord = domain.coords[ycoord]
#         zcoord = domain.coords[zcoord]
#         tcoord = domain.coords[tcoord]
#
#         u = self.meteorolgy.getWindVeclocity(self.sourceHeight)
#
#         sigma = getSigma(self, x, stability, sigma0=None)
#
#
#         X,T = np.meshgrid(xcoord, tcoord, indexing='ij')
#         sigmaX = 1
#         downwind = (1/(np.sqrt(2*np.pi)*sigmaX))*np.exp(-(X-u*T)**2 / (2*sigmaX**2))
#         XR_downwind = xr.DataArray(downwind, dims = ("x","time"), coords={"x": xcoord, "time": tcoord})
#
#         X,Y = np.meshgrid(xcoord, ycoord, indexing='ij')
#         sigmaY = 1
#         crosswind = (1/(np.sqrt(2*np.pi)*sigmaY))*np.exp(-(Y-0)**2 / (2*sigmaX**2))
#         XR_crosswind = xr.DataArray(crosswind, dims = ("x","y"), coords={"x": xcoord, "y": ycoord})
#
#         X,Z = np.meshgrid(xcoord, zcoord, indexing='ij')
#         sigmaZ = 1
#
#         nSum = np.arange(-numOfIteration,numOfIteration+1,1)
#         vertical = 0
#         for n in numOfIteration:
#             vertical += (1/(np.sqrt(2*np.pi)*sigmaZ))*(np.exp(-(Z-2*n*self.inversion)**2 / (2*sigmaX**2)) +
#                                                   np.exp(-(Z+2*n*self.inversion)**2 / (2*sigmaX**2)))
#         XR_vertical = xr.DataArray(vertical, dims = ("x","z"), coords={"x": xcoord, "z": zcoord})
#
#         down, cross, vert = xr.broadcast(XR_downwind, XR_crosswind, XR_vertical)
#
#         concentrations = down*cross*vert
#
#         return concentrations
#
#
#
#
#
#
#
#
#
#
#



























