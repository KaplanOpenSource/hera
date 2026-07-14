import pandas
import numpy
from numpy import matlib
import math
from hera.utils import *
from hera.utils.unitHandler import ureg, unumToPint
import xarray
from scipy import special



class abstractGasCloud:

    def __init__(self, sourceQ, sourceHeight, initialCloudSize, sigmaType):
        """

        Parameters
        ----------
        sourceQ : pint Quantity or unum
            If Quantity:
                The unit determine the release time.
                [mass] - Instantaneous
                [mass/time] - Continuous
            else
                Continuous (not implementaed yet.)

        sourceHeight : pint Quantity
        initialCloudSize : 3-touple pint Quantity, the sigmas in each axis.
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
        sourceQ : pint Quantity or unum
            If Quantity:
                The unit determine the release time.
                [mass] - Instantaneous
                [mass/time] - Continuous
            else
                Continuous (not implementaed yet.)

        sourceHeight : pint Quantity
        initialCloudSize : 3-touple pint Quantity, the sigmas in each axis.

        Returns
        -------

        """
        sourceQ = unumToPint(sourceQ)
        try:
            sourceQ.to(ureg.mg)
            instantaneous = True
        except Exception:
            try:
                sourceQ.to(ureg.mg/ureg.min)
                instantaneous = False
            except Exception:
                raise ValueError("Must be mass or mass per time!")

        returnCls = instantaneousReleaseGasCloud if instantaneous else continuousReleaseGasCloud

        return returnCls(sourceQ=sourceQ,sourceHeight=sourceHeight,initialCloudSize=initialCloudSize,sigmaType=sigmaType)


    def _getTXterm(self, stability, u, xcoordRange, tcoordRange):
        """
        Parameters
        ----------
        initialCloudSize : 3-tuple of float/pint Quantity (default m)
            The initial cloud size (standard deviation) in the x,y and z dimensions.
        stability : the stability class
        u: the wind speed / pint Quantity (default m/s)
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
        initialCloudSize : 3-tuple of float/pint Quantity (default m)
            The initial cloud size (standard deviation) in the x,y and z dimensions.
        stability : the stability class
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
        initialCloudSize : 3-tuple of float/pint Quantity (default m)
            The initial cloud size (standard deviation) in the x,y and z dimensions.
        stability : the stability class
        inversion: the inversion height / pint Quantity (default m)
        xcoordRange : Tuple in numpy.arange format, unitless.
        zcoordRange : Tuple in numpy.arange format, unitless.
        numOfReflections : The number of reflections of the summation of the Z component.

        Returns
        -------
        The Z component of the Gaussian concentration formula.
        """
        sourceHeight = tonumber(self.sourceHeight, ureg.m)
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




    def fractions(self, fracVector, minx, miny, minz, maxx, maxy, maxz, timeSpan, dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min):

        """
        This function generates an xarray of fractions of the mass of each isotope at every time-step.
        After generating this xarray, we multiply it by the concentration xarray in order to get the relative
        concentration of a given isotope.
        :param fracVector: Tuple in numpy.arange format, unitless.
                            A vector of the isotope's fractions at every time-step (calculated by Eshed's code).
                            Note - This vector should have the same length as tcoordRange (described within the function)
        :return:
        """

        xcoordRange = numpy.arange(tonumber(minx, ureg.m), tonumber(maxx, ureg.m), tonumber(dxdy,ureg.m))
        ycoordRange = numpy.arange(tonumber(miny,ureg.m),tonumber(maxy,ureg.m),tonumber(dxdy,ureg.m))
        zcoordRange = numpy.arange(tonumber(minz, ureg.m), tonumber(maxz, ureg.m), tonumber(dz,ureg.m))
        tcoordRange = numpy.arange(0,tonumber(timeSpan,ureg.min),tonumber(dt,ureg.min))

        frac, X = numpy.meshgrid(fracVector, xcoordRange, indexing='ij')
        XR_downwind = xarray.DataArray(frac, dims=("time", "x"), coords={"time": tcoordRange, "x": xcoordRange})
        XR_crosswind = xarray.DataArray(1, dims=("x", "y"), coords={"x": xcoordRange, "y": ycoordRange})
        XR_vertical = xarray.DataArray(1, dims=("x", "z"), coords={"x": xcoordRange, "z": zcoordRange})
        FULL = XR_downwind * XR_crosswind * XR_vertical

        return FULL



    def _getTXDosage(self, stability, u, xcoordRange, tcoordRange):
        """
        Parameters
        ----------
        initialCloudSize : 3-tuple of float/pint Quantity (default m)
            The initial cloud size (standard deviation) in the x,y and z dimensions.
        stability : the stability class
        u: the wind speed / pint Quantity (default m/s)
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



#---------------------------------------------DF (depletion factor)------------------------------------------------

    def _getTXterm_ones(self, xcoordRange, tcoordRange):
        """
        Parameters
        ----------
        xcoordRange : Tuple in numpy.arange format, unitless.
        tcoordRange : Tuple in numpy.arange format, unitless.

        Returns
        -------
        Array of 1s.
        """
        # T, X = numpy.meshgrid(tcoordRange, xcoordRange, indexing='ij')
        XR_downwind = xarray.DataArray(1, dims=("time", "x"), coords={"time": tcoordRange, "x": xcoordRange})
        return XR_downwind

    def _getXYterm_ones(self, xcoordRange, ycoordRange):
        """

        Parameters
        ----------
        xcoordRange : Tuple in numpy.arange format, unitless.
        ycoordRange : Tuple in numpy.arange format, unitless.

        Returns
        -------
        Array of 1s.
        """
        # X, Y = numpy.meshgrid(xcoordRange, ycoordRange, indexing='ij')
        XR_crosswind = xarray.DataArray(1, dims=("x", "y"), coords={"x": xcoordRange, "y": ycoordRange})
        return XR_crosswind

    def _getXZterm_ones(self, xcoordRange, zcoordRange):
        """

        Parameters
        ----------
        xcoordRange : Tuple in numpy.arange format, unitless.
        zcoordRange : Tuple in numpy.arange format, unitless.

        Returns
        -------
        Array of 1s.
        """

        # X, Z = numpy.meshgrid(xcoordRange, zcoordRange, indexing='ij')
        XR_vertical = xarray.DataArray(1, dims=("x", "z"), coords={"x": xcoordRange, "z": zcoordRange})

        return XR_vertical


    def _getDF(self, stability, u, xcoordRange, zcoordRange):

        """
        :param stability: The stability state
        :param u: Wind speed
        :param xcoordRange: x-coordinates of the grid
        :param zcoordRange: z-coordinates of the grid
        :return: The Depletion Factor
        """

        v = tonumber(0.003*ureg.m/ureg.s, ureg.m/ureg.min) #Deposition velocity. By default we take this value to be 0.003 [m/s]

        X, Z = numpy.meshgrid(xcoordRange, zcoordRange, indexing='ij')
        sigmaZ = self.sigmaType.getSigma(x=X, stability=stability, sigma0=self.initialCloudSize, units=False)['sigmaZ']
        H = tonumber(self.sourceHeight, ureg.m)
        dx = xcoordRange[1] - xcoordRange[0]
        DF = (numpy.e**(numpy.cumsum(1/(sigmaZ*numpy.e**(0.5*(H/sigmaZ)**2))*dx, axis=0)))**((-v/u)*numpy.sqrt(2/numpy.pi))

        return DF



    def getDF_noQ_xarray(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan, dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min):
        stability = meteorology.stability
        u = tonumber(meteorology.u10, ureg.m/ureg.min)

        xcoordRange = numpy.arange(tonumber(minx, ureg.m), tonumber(maxx, ureg.m), tonumber(dxdy, ureg.m))
        ycoordRange = numpy.arange(tonumber(miny, ureg.m), tonumber(maxy, ureg.m), tonumber(dxdy, ureg.m))
        zcoordRange = numpy.arange(tonumber(minz, ureg.m), tonumber(maxz, ureg.m), tonumber(dz, ureg.m))
        tcoordRange = numpy.arange(0, tonumber(timeSpan, ureg.min), tonumber(dt, ureg.min))

        TX = self._getTXterm_ones(xcoordRange=xcoordRange, tcoordRange=tcoordRange)
        XY = self._getXYterm_ones(xcoordRange=xcoordRange, ycoordRange=ycoordRange)
        XZ = self._getXZterm_ones(xcoordRange=xcoordRange, zcoordRange=zcoordRange)
        DF = self._getDF(stability=stability, u=u, xcoordRange=xcoordRange, zcoordRange=zcoordRange)

        return TX*XY*(XZ*DF)

# ---------------------------------------------END OF DF------------------------------------------------


class instantaneousReleaseGasCloud(abstractGasCloud):


    def getConcentrationFromMinMaxRange_inst_noQ(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                                        dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min, numOfReflections=3):

        xcoordRange = numpy.arange(tonumber(minx, ureg.m), tonumber(maxx, ureg.m), tonumber(dxdy,ureg.m))
        ycoordRange = numpy.arange(tonumber(miny,ureg.m),tonumber(maxy,ureg.m),tonumber(dxdy,ureg.m))
        zcoordRange = numpy.arange(tonumber(minz, ureg.m), tonumber(maxz, ureg.m), tonumber(dz,ureg.m))
        tcoordRange = numpy.arange(0,tonumber(timeSpan,ureg.min),tonumber(dt,ureg.min))

        stability = meteorology.stability
        u = tonumber(meteorology.u10, ureg.m/ureg.min)
        inversion = tonumber(meteorology.inversion, ureg.m)

        TX = self._getTXterm(stability=stability, u=u, xcoordRange=xcoordRange, tcoordRange=tcoordRange)
        XY = self._getXYterm(stability=stability, xcoordRange=xcoordRange, ycoordRange=ycoordRange)
        XZ = self._getXZterm(stability=stability, inversion=inversion, xcoordRange=xcoordRange, zcoordRange=zcoordRange,
                             numOfReflections=numOfReflections)

        ret = TX*XY*XZ
        ret.attrs['Q'] = 1/ureg.m**3

        return ret


    def getDosageFromMinMaxRange_inst_noQ(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                           dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min, numOfReflections=3, DF=False):
        stability = meteorology.stability
        u = tonumber(meteorology.u10, ureg.m/ureg.min)
        inversion = tonumber(meteorology.inversion, ureg.m)

        xcoordRange = numpy.arange(tonumber(minx, ureg.m), tonumber(maxx, ureg.m), tonumber(dxdy, ureg.m))
        ycoordRange = numpy.arange(tonumber(miny, ureg.m), tonumber(maxy, ureg.m), tonumber(dxdy, ureg.m))
        zcoordRange = numpy.arange(tonumber(minz, ureg.m), tonumber(maxz, ureg.m), tonumber(dz, ureg.m))
        tcoordRange = numpy.arange(0, tonumber(timeSpan, ureg.min), tonumber(dt, ureg.min))

        TX = self._getTXDosage(stability=stability, u=u, xcoordRange=xcoordRange, tcoordRange=tcoordRange)
        XY = self._getXYterm(stability=meteorology.stability, xcoordRange=xcoordRange, ycoordRange=ycoordRange)
        XZ = self._getXZterm(stability=stability, inversion=inversion, xcoordRange=xcoordRange, zcoordRange=zcoordRange,
                             numOfReflections=numOfReflections)
        D_F = 1
        if DF:
            D_F = self._getDF(stability=stability, u=u, xcoordRange=xcoordRange, zcoordRange=zcoordRange)

        ret = TX*XY*(XZ*D_F)
        ret.attrs['Q'] = 1*ureg.min/ureg.m**3

        return ret

    def getDosageFromMinMaxRange_inst_NoERF_noQ(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                           dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min, numOfReflections=3, DF=False):
        C_without_Q = self.getConcentrationFromMinMaxRange_inst_noQ(meteorology=meteorology, minx=minx, miny=miny, minz=minz,
                                                                    maxx=maxx, maxy=maxy, maxz=maxz, timeSpan=timeSpan,dxdy=dxdy,
                                                                    dz=dz, dt=dt,numOfReflections=numOfReflections)

        D_without_Q = self.trapezoidal_integration(data=C_without_Q)
        D_F = 1
        if DF:
            D_F = self.getDF_noQ_xarray(meteorology=meteorology, minx=minx, miny=miny, minz=minz, maxx=maxx, maxy=maxy,
                                        maxz=maxz, timeSpan=timeSpan, dxdy=dxdy, dz=dz, dt=dt)

        ret = D_without_Q*D_F
        ret.attrs['Q'] = 1*ureg.min/ureg.m**3

        return ret


class instantaneousReleaseGasCloud(abstractGasCloud):


    def getConcentrationFromMinMaxRange_inst_noQ(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                                        dxdy=10*m, dz=1*m, dt=1*min, numOfReflections=3):

        xcoordRange = numpy.arange(tonumber(minx, m), tonumber(maxx, m), tonumber(dxdy,m))
        ycoordRange = numpy.arange(tonumber(miny,m),tonumber(maxy,m),tonumber(dxdy,m))
        zcoordRange = numpy.arange(tonumber(minz, m), tonumber(maxz, m), tonumber(dz,m))
        tcoordRange = numpy.arange(0,tonumber(timeSpan,min),tonumber(dt,min))

        stability = meteorology.stability
        u = tonumber(meteorology.u10, m/min)
        inversion = tonumber(meteorology.inversion, m)

        TX = self._getTXterm(stability=stability, u=u, xcoordRange=xcoordRange, tcoordRange=tcoordRange)
        XY = self._getXYterm(stability=stability, xcoordRange=xcoordRange, ycoordRange=ycoordRange)
        XZ = self._getXZterm(stability=stability, inversion=inversion, xcoordRange=xcoordRange, zcoordRange=zcoordRange,
                             numOfReflections=numOfReflections)

        ret = TX*XY*XZ
        ret.attrs['Q'] = 1/m**3

        return ret


    def getConcentrationFromMinMaxRange_inst(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                                        dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min, numOfReflections=3, DF=False):
        C_without_Q = self.getConcentrationFromMinMaxRange_inst_noQ(meteorology=meteorology, minx=minx, miny=miny, minz=minz,
                                                                    maxx=maxx, maxy=maxy, maxz=maxz, timeSpan=timeSpan,dxdy=dxdy,
                                                                    dz=dz, dt=dt,numOfReflections=numOfReflections)

        ret = tonumber(self.sourceQ, ureg.mg)*C_without_Q
        ret.attrs['Q'] = 1*ureg.mg/ureg.m**3
        return ret


    def getDosageFromMinMaxRange_inst(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                           dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min, numOfReflections=3, DF=False):
        D_without_Q = self.getDosageFromMinMaxRange_inst_noQ(meteorology=meteorology, minx=minx, miny=miny, minz=minz, maxx=maxx,
                                                             maxy=maxy, maxz=maxz, timeSpan=timeSpan, dxdy=dxdy, dz=dz,
                                                             dt=dt, numOfReflections=numOfReflections, DF=DF)

        ret = tonumber(self.sourceQ, ureg.mg)*D_without_Q
        ret.attrs['Q'] = 1*ureg.mg*ureg.min/ureg.m**3
        return ret



    def getDosageFromMinMaxRange_inst_NoERF(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                           dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min, numOfReflections=3, DF=False):

        D_without_Q = self.getDosageFromMinMaxRange_inst_NoERF_noQ(meteorology=meteorology, minx=minx, miny=miny, minz=minz,
                                                                   maxx=maxx, maxy=maxy,maxz=maxz, timeSpan=timeSpan,dxdy=dxdy,
                                                                   dz=dz, dt=dt, numOfReflections=numOfReflections, DF=DF)

        ret = tonumber(self.sourceQ, ureg.mg)*D_without_Q
        ret.attrs['Q'] = 1*ureg.mg*ureg.min/ureg.m**3
        return ret



    #---------------------------------------------Radiology---------------------------------------------

    def concentrationConversion_mass_to_Bq(self, C, outputUnits, specificActivity):
        units = unumToPint(C.attrs['Q'])
        out_units = unumToPint(outputUnits)
        factor = (units * specificActivity).m_as(out_units)
        C_Bq = C * factor
        C_Bq.attrs['Q'] = out_units
        return C_Bq


    def getTIACFromMinMaxRange_inst(self,specifitActivity, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                           dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min, numOfReflections=3, outputUnits=ureg.Bq*ureg.s/ureg.m**3, DF=False):

        D_without_Q = self.getDosageFromMinMaxRange_inst_noQ(meteorology=meteorology, minx=minx, miny=miny, minz=minz, maxx=maxx,
                                                             maxy=maxy, maxz=maxz, timeSpan=timeSpan, dxdy=dxdy, dz=dz,
                                                             dt=dt, numOfReflections=numOfReflections, DF=DF)

        out_units = unumToPint(outputUnits)
        ret = tonumber(self.sourceQ, ureg.mg)*D_without_Q
        factor = (1*ureg.mg*ureg.min/ureg.m**3 * specifitActivity).m_as(out_units)
        ret *= factor
        ret.attrs['Q'] = out_units
        return ret


    def getTIACFromMinMaxRange_inst_noQ(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                           dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min, numOfReflections=3, outputUnits=ureg.s/ureg.m**3, DF=False):

        D_without_Q = self.getDosageFromMinMaxRange_inst_noQ(meteorology=meteorology, minx=minx, miny=miny, minz=minz, maxx=maxx,
                                                             maxy=maxy, maxz=maxz, timeSpan=timeSpan, dxdy=dxdy, dz=dz,
                                                             dt=dt, numOfReflections=numOfReflections, DF=DF)
        out_units = unumToPint(outputUnits)
        currentUnites = unumToPint(D_without_Q.attrs['Q'])
        factor = currentUnites.m_as(out_units)
        D_without_Q *= factor
        D_without_Q.attrs['Q'] = out_units
        return D_without_Q



    def getTIACFromMinMaxRange_inst_NoERF(self, specifitActivity, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                                    dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min, numOfReflections=3, outputUnits=ureg.Bq*ureg.s/ureg.m**3, DF=False):

        D_without_Q = self.getDosageFromMinMaxRange_inst_NoERF_noQ(meteorology=meteorology, minx=minx, miny=miny, minz=minz,
                                                             maxx=maxx, maxy=maxy, maxz=maxz, timeSpan=timeSpan, dxdy=dxdy,
                                                             dz=dz, dt=dt, numOfReflections=numOfReflections, DF=DF)

        out_units = unumToPint(outputUnits)
        ret = tonumber(self.sourceQ, ureg.mg) * D_without_Q
        factor = (1 * ureg.mg * ureg.min / ureg.m ** 3 * specifitActivity).m_as(out_units)
        ret *= factor
        ret.attrs['Q'] = out_units
        return ret


    def getTIACFromMinMaxRange_inst_NoERF_noQ(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                                    dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min, numOfReflections=3, outputUnits=ureg.s/ureg.m**3, DF=False):

        D_without_Q = self.getDosageFromMinMaxRange_inst_NoERF_noQ(meteorology=meteorology, minx=minx, miny=miny, minz=minz,
                                                             maxx=maxx, maxy=maxy, maxz=maxz, timeSpan=timeSpan, dxdy=dxdy,
                                                             dz=dz, dt=dt, numOfReflections=numOfReflections, DF=DF)

        out_units = unumToPint(outputUnits)
        currentUnites = unumToPint(D_without_Q.attrs['Q'])
        factor = currentUnites.m_as(out_units)
        D_without_Q *= factor
        D_without_Q.attrs['Q'] = out_units
        return D_without_Q




    def getTIACFromConcentration_inst_NoERF(self, C, specifitActivity, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                                    dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min, outputUnits=ureg.Bq*ureg.s/ureg.m**3, DF=False):
        """

        :param C: xarray of concentrations in units of [mass/volume].
        :param specifitActivity: The specific activity of the isotope.
        :param outputUnits: The desired output units [Bq*time/volume].
        :return: TIAC (Time Integrated Air Concentration) in unites of [Bq*time/volume]
        """

        out_units = unumToPint(outputUnits)
        factor = unumToPint(C.attrs['Q']).m_as(ureg.mg/ureg.m**3)
        C_mg_m3 = C*factor
        C_mg_m3.attrs['Q'] = ureg.mg/ureg.m**3

        C_without_Q = C_mg_m3 / tonumber(self.sourceQ, ureg.mg)
        D_without_Q = self.trapezoidal_integration(data = C_without_Q)

        D_F = 1
        if DF:
            D_F = self.getDF_noQ_xarray(meteorology=meteorology, minx=minx, miny=miny, minz=minz, maxx=maxx, maxy=maxy,
                                        maxz=maxz, timeSpan=timeSpan, dxdy=dxdy, dz=dz, dt=dt)

        D_without_Q = D_without_Q * D_F

        ret = tonumber(self.sourceQ, ureg.mg) * D_without_Q
        factor = (1 * ureg.mg * ureg.min / ureg.m ** 3 * specifitActivity).m_as(out_units)
        ret *= factor
        ret.attrs['Q'] = out_units
        return ret


    def getTIACFromConcentration_inst_NoERF_noQ(self, C_noQ, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                                    dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min, outputUnits=ureg.s/ureg.m**3, DF=False):
        """

        :param C_noQ: xarray of concentrations in units of [1/volume].
        :return: TIAC (Time Integrated Air Concentration) in unites of [time/volume]
        """

        D_without_Q = self.trapezoidal_integration(data = C_noQ)

        D_F = 1
        if DF:
            D_F = self.getDF_noQ_xarray(meteorology=meteorology, minx=minx, miny=miny, minz=minz, maxx=maxx, maxy=maxy,
                                        maxz=maxz, timeSpan=timeSpan, dxdy=dxdy, dz=dz, dt=dt)

        D_without_Q = D_without_Q*D_F
        D_without_Q.attrs['Q'] = 1*ureg.min/ureg.m**3 #need to verify that C_noQ was generated with time steps in [min], not [s]

        out_units = unumToPint(outputUnits)
        currentUnites = unumToPint(D_without_Q.attrs['Q'])
        factor = currentUnites.m_as(out_units)
        D_without_Q *= factor
        D_without_Q.attrs['Q'] = out_units

        return D_without_Q



    def get_TIAC_from_dist(self, data, y, z, dist_list):
        """

        :param data: TIAC xarray with coordinates x,y,z,time
        :param y: The y value we wish to get the TIAC values for ([m])
        :param z: The z value we wish to get the TIAC values for ([m])
        :param dist_list: A list of x values we wish to get the TIAC values for ([m])
        :return: A list of tuples (x, TIAC)

        """
        distances = numpy.array(data.squeeze().x)
        TIAC_lastTime_DF = numpy.array(data.sel(y=y, z=z, time=data.time[-1], method='nearest'))
        tuples = list(tuple(zip(distances, TIAC_lastTime_DF)))

        for x in dist_list:
            if x not in distances:
                print(f"x={x} not in Xarray. Try multiples of {distances[1] - distances[0]}")

        return [(x, tiac) for x, tiac in tuples if x in dist_list]




class continuousReleaseGasCloud(instantaneousReleaseGasCloud):

    def getConcentrationFromMinMaxRange_cont(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                                        dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min, numOfReflections=3, DF=False):
        """
        Returns
        -------
        An xarray of concentrations at every grid-poit, which is the dosage of the instantaneous release,
        since we assume the release rate is constant.
        Here we take the concentration xarray that was claculated using the error function (erf).
        """
        C_without_Q = self.getDosageFromMinMaxRange_inst_noQ(meteorology=meteorology, minx=minx, miny=miny, minz=minz, maxx=maxx,
                                                        maxy=maxy, maxz=maxz, timeSpan=timeSpan, dxdy=dxdy, dz=dz, dt=dt,
                                                        numOfReflections=numOfReflections, DF=DF)

        return tonumber(self.sourceQ, ureg.mg/ureg.s)*C_without_Q


    def getConcentrationFromMinMaxRange_cont_NoERF(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                                        dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min, numOfReflections=3, DF=False):
        """
        Returns
        -------
        An xarray of concentrations at every grid-poit, which is the dosage of the instantaneous release,
        since we assume the release rate is constant.
        Here we take the concentration xarray that was claculated without the error function (erf).
        """
        C_without_Q = self.getDosageFromMinMaxRange_inst_NoERF_noQ(meteorology=meteorology, minx=minx, miny=miny, minz=minz,
                                                          maxx=maxx, maxy=maxy, maxz=maxz, timeSpan=timeSpan, dxdy=dxdy,
                                                          dz=dz, dt=dt, numOfReflections=numOfReflections, DF=DF)
        return tonumber(self.sourceQ, ureg.mg/ureg.s)*C_without_Q


    def getDosageFromMinMaxRange_cont_NoERF(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                                        dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min, numOfReflections=3, DF=False):

        C_without_Q = self.getConcentrationFromMinMaxRange_cont(meteorology=meteorology, minx=minx, miny=miny, minz=minz,
                                                           maxx=maxx, maxy=maxy, maxz=maxz, timeSpan=timeSpan, dxdy=dxdy,
                                                           dz=dz, dt=dt, numOfReflections=numOfReflections, DF=DF)
        D_without_Q = self.trapezoidal_integration(data=C_without_Q)

        return tonumber(self.sourceQ, ureg.mg/ureg.min) * D_without_Q

    def getDosageFromMinMaxRange_cont_doubleNoERF(self, meteorology, minx, miny, minz, maxx, maxy, maxz, timeSpan,
                                        dxdy=10*ureg.m, dz=1*ureg.m, dt=1*ureg.min, numOfReflections=3, DF=False):

        C_without_Q = self.getConcentrationFromMinMaxRange_cont_NoERF(meteorology=meteorology, minx=minx, miny=miny, minz=minz,
                                                           maxx=maxx, maxy=maxy, maxz=maxz, timeSpan=timeSpan, dxdy=dxdy,
                                                           dz=dz, dt=dt, numOfReflections=numOfReflections, DF=DF)
        D_without_Q = self.trapezoidal_integration(data=C_without_Q)

        return tonumber(self.sourceQ, ureg.mg/ureg.min) * D_without_Q



#-------------------------- Yehuda's Code For Convolution --------------------------

class Continuous(object):
    dt = None
    M  = None
    Timekernel = None
    _FullKernel = None

    def __init__(self,dt,kernelsize,timetofinish=10*ureg.min):
        r"""
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
        dt = unumToPint(dt).to(ureg.min)
        self.dt = dt.m_as(ureg.min)
        self.kernelsize  = kernelsize


        timetofinish = unumToPint(timetofinish).to(ureg.min)
        alpha = numpy.log(0.1)/(-timetofinish.m_as(ureg.min))

        # build the kernel.
        ts = numpy.arange(kernelsize,-1,-1)*dt.m_as(ureg.min)

        # the kernel
        self.Timekernel = (numpy.exp(-alpha*ts[1:]) - numpy.exp(-alpha*ts[:-1]))/dt.m_as(ureg.min)
        self.Timekernel = self.Timekernel.reshape([kernelsize,1,1,1])


    def _convolve(self,data,axis,FullKernel):
        # convolve.
        return (data*FullKernel[FullKernel.shape[0]-data.shape[0]:,:,:]).sum(axis=axis)*self.dt

    def calc(self,data):
        # build the kernel.
        FullKernel = numpy.tile(self.Timekernel,[1,data.x.size,data.y.size,data.z.size])
        return data.rolling(datetime=self.kernelsize,min_periods=1).reduce(self._convolve,FullKernel=FullKernel)


#-------------------------- End Of Yehuda's Code For Convolution --------------------------



