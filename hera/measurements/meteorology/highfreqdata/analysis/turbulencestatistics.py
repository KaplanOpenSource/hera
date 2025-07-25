import os
import json
import numpy
import pandas
import dask.dataframe
from scipy.constants import g
from .abstractcalculator import AbstractCalculator


class singlePointTurbulenceStatistics(AbstractCalculator):

    def __init__(self, rawData, metadata):
        super(singlePointTurbulenceStatistics, self).__init__(rawData=rawData, metadata=metadata)
        self.data = None

    @property
    def isMissingData(self):
        return self.metaData["isMissingData"]

    def getData(self):
        return self.data

    def fluctuations(self, inMemory=None):
        """
        Calculates the mean of u,v,w,T,wind_dir and the fluctuations u',v',w',T',wind_dir'.

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'up' not in self._RawData.columns:
            avg = self._RawData[['u','v','w','T']]
            if self.SamplingWindow is None:
                avg = avg.mean()
                if self._DataType == 'pandas':
                    avg = pandas.DataFrame(avg).T
                    avg.index = [self._RawData.index[0]]
                else:
                    self._RawData = self._RawData.repartition(npartitions=1)
                    avg = pandas.DataFrame(avg.compute()).T
                    avg.index = self._RawData.head(1).index
                    avg = dask.dataframe.from_pandas(avg, npartitions=1)
            else:
                avg = avg.resample(self.SamplingWindow).mean()

            avg = avg.rename(columns={'u': 'u_bar', 'v': 'v_bar', 'w': 'w_bar', 'T': 'T_bar'})

            avg['wind_dir_bar'] = numpy.arctan2(avg['v_bar'], avg['u_bar'])
            avg['wind_dir_bar'] = numpy.rad2deg(avg['wind_dir_bar'])
            avg['wind_dir_bar'] = (270 - avg['wind_dir_bar']) % 360

            self._TemporaryData = avg
            self._CalculatedParams += [['u_bar',{}], ['v_bar',{}], ['w_bar',{}], ['T_bar',{}], ['wind_dir_bar', {}]]
            if self.isMissingData:
                self._RawData = self._RawData.merge(avg, how='outer', left_index=True, right_index=True)
                self._RawData = self._RawData.dropna(how='all')
                self._RawData[['u_bar', 'v_bar', 'w_bar', 'T_bar', 'wind_dir_bar']] = self._RawData[['u_bar', 'v_bar', 'w_bar', 'T_bar', 'wind_dir_bar']].ffill()
                self._RawData = self._RawData.dropna(how='any')
            else:
                self._RawData = self._RawData.merge(avg, how='left', left_index=True, right_index=True)
                self._RawData = self._RawData.ffill()

            self._RawData['wind_dir'] = numpy.arctan2(self._RawData['v'], self._RawData['u'])
            self._RawData['wind_dir'] = numpy.rad2deg(self._RawData['wind_dir'])
            self._RawData['wind_dir'] = (270 - self._RawData['wind_dir']) % 360


            self._RawData['up'] = self._RawData['u'] - self._RawData['u_bar']
            self._RawData['vp'] = self._RawData['v'] - self._RawData['v_bar']
            self._RawData['wp'] = self._RawData['w'] - self._RawData['w_bar']
            self._RawData['Tp'] = self._RawData['T'] - self._RawData['T_bar']
            self._RawData['wind_dir_p'] = (180-(180-(self._RawData['wind_dir'] - self._RawData['wind_dir_bar']).abs()).abs()).abs()

            self.data = self._RawData
        return self

    def sigma(self, inMemory=None):
        """


        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'sigmaU' not in self._TemporaryData.columns:
            self.fluctuations()

            if self.SamplingWindow is None:
                sigmaU = self._RawData['u'].std()
                sigmaV = self._RawData['v'].std()
                sigmaW = self._RawData['w'].std()
            else:
                sigmaU = self._RawData['u'].resample(self.SamplingWindow).std()
                sigmaV = self._RawData['v'].resample(self.SamplingWindow).std()
                sigmaW = self._RawData['w'].resample(self.SamplingWindow).std()

            self._TemporaryData['sigmaU'] = sigmaU
            self._CalculatedParams.append(['sigmaU',{}])

            self._TemporaryData['sigmaV'] = sigmaV
            self._CalculatedParams.append(['sigmaV',{}])

            self._TemporaryData['sigmaW'] = sigmaW
            self._CalculatedParams.append(['sigmaW',{}])

        return self

    def sigmaH(self, inMemory=None):
        """


        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'sigmaH' not in self._TemporaryData.columns:
            self.sigma()
            sigmaH = numpy.hypot(self._TemporaryData['sigmaU'], self._TemporaryData['sigmaV']) / numpy.sqrt(2)
            self._TemporaryData['sigmaH'] = sigmaH
            self._CalculatedParams.append(['sigmaH',{}])

        return self

    def sigmaHOverUstar(self, inMemory=None):
        """


        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'sigmaHOverUstar' not in self._TemporaryData.columns:
            self.sigmaH()
            self.Ustar()
            sigmaHOverUstar = self._TemporaryData['sigmaH']/self._TemporaryData['Ustar']
            self._TemporaryData['sigmaHOverUstar'] = sigmaHOverUstar
            self._CalculatedParams.append(['sigmaHOverUstar',{}])

        return self

    def sigmaWOverUstar(self, inMemory=None):
        """



        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'sigmaWOverUstar' not in self._TemporaryData.columns:
            self.sigma()
            self.Ustar()
            sigmaWOverUstar = self._TemporaryData['sigmaW']/self._TemporaryData['Ustar']
            self._TemporaryData['sigmaWOverUstar'] = sigmaWOverUstar
            self._CalculatedParams.append(['sigmaWOverUstar',{}])

        return self

    def horizontalSpeed(self, inMemory=None):
        """
        Calculates the mean and the std of the horizontal speed.

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'horizontal_speed_bar' not in self._TemporaryData.columns:
            self.fluctuations()
            horizontal_speed_bar = numpy.hypot(self._TemporaryData['u_bar'], self._TemporaryData['v_bar'])
            self._TemporaryData['horizontal_speed_bar'] = horizontal_speed_bar
            self._CalculatedParams.append(['horizontal_speed_bar',{}])

        return self



    def wind_dir_std(self, inMemory=None):
        """
        Calculates the std of the wind direction.

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'wind_dir_std' not in self._TemporaryData.columns:
            self.fluctuations()

            std = numpy.square(self._RawData['wind_dir_p'])
            std = std if self.SamplingWindow is None else std.resample(self.SamplingWindow)
            std = numpy.sqrt(std.mean())

            self._TemporaryData['wind_dir_std'] = std
            self._CalculatedParams.append(['wind_dir_std',{}])

        return self

    def sigmaHOverWindSpeed(self, inMemory=None):
        """


        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'sigmaHOverWindSpeed' not in self._TemporaryData.columns:
            self.sigmaH()
            self.horizontalSpeed()
            sigmaHOverWindSpeed = self._TemporaryData['sigmaH']/self._TemporaryData['horizontal_speed_bar']
            self._TemporaryData['sigmaHOverWindSpeed'] = sigmaHOverWindSpeed
            self._CalculatedParams.append(['sigmaHOverWindSpeed',{}])

        return self

    def sigmaWOverWindSpeed(self, inMemory=None):
        """


        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'sigmaWOverWindSpeed' not in self._TemporaryData.columns:
            self.sigma()
            self.horizontalSpeed()
            sigmaWOverWindSpeed = self._TemporaryData['sigmaW']/self._TemporaryData['horizontal_speed_bar']
            self._TemporaryData['sigmaWOverWindSpeed'] = sigmaWOverWindSpeed
            self._CalculatedParams.append(['sigmaWOverWindSpeed',{}])

        return self

    def w3OverSigmaW3(self, inMemory=None):
        """


        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'w3OverSigmaW3' not in self._TemporaryData.columns:
            self.w3()
            self.sigma()
            w3OverSigmaW3 = self._TemporaryData['w3']/self._TemporaryData['sigmaW']**3
            self._TemporaryData['w3OverSigmaW3'] = w3OverSigmaW3
            self._CalculatedParams.append(['w3OverSigmaW3',{}])

        return self

    def uStarOverWindSpeed(self, inMemory=None):
        """

uStarOverWindSpeed
        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'uStarOverWindSpeed' not in self._TemporaryData.columns:
            self.Ustar()
            self.horizontalSpeed()
            uStarOverWindPeed = self._TemporaryData['Ustar']/self._TemporaryData['horizontal_speed_bar']
            self._TemporaryData['uStarOverWindSpeed'] = uStarOverWindPeed
            self._CalculatedParams.append(['uStarOverWindSpeed',{}])

        return self

    def secondMoments(self, inMemory=None):
        """


        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """

        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        moments_list = ["uu","uv","uw","uT","vv","vw","vT","ww","wT","TT"]

        for moment in moments_list:
            getattr(self, moment)(inMemory= inMemory)

        return self

    def uu(self, inMemory=None):
        """
        Calculates the mean of u'*u'.

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'uu' not in self._TemporaryData.columns:
            self.fluctuations()
            if self.SamplingWindow is None:
                uu = (self._RawData['up'] * self._RawData['up']).mean()
            else:
                uu = (self._RawData['up'] * self._RawData['up']).resample(self.SamplingWindow).mean()
            self._TemporaryData['uu'] = uu
            self._CalculatedParams.append(['uu',{}])

        return self

    def vv(self, inMemory=None):
        """
        Calculates the mean of v'*v'.

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'vv' not in self._TemporaryData.columns:
            self.fluctuations()
            if self.SamplingWindow is None:
                vv = (self._RawData['vp'] * self._RawData['vp']).mean()
            else:
                vv = (self._RawData['vp'] * self._RawData['vp']).resample(self.SamplingWindow).mean()
            self._TemporaryData['vv'] = vv
            self._CalculatedParams.append(['vv',{}])

        return self

    def ww(self, inMemory=None):
        """
        Calculates the mean of w'*w'.

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'ww' not in self._TemporaryData.columns:
            self.fluctuations()
            if self.SamplingWindow is None:
                ww = (self._RawData['wp'] * self._RawData['wp']).mean()
            else:
                ww = (self._RawData['wp'] * self._RawData['wp']).resample(self.SamplingWindow).mean()
            self._TemporaryData['ww'] = ww
            self._CalculatedParams.append(['ww',{}])

        return self

    def TT(self, inMemory=None):
        """
        Calculates the mean of T'*T'.

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'TT' not in self._TemporaryData.columns:
            self.fluctuations()
            if self.SamplingWindow is None:
                TT = (self._RawData['Tp'] * self._RawData['Tp']).mean()
            else:
                TT = (self._RawData['Tp'] * self._RawData['Tp']).resample(self.SamplingWindow).mean()
            self._TemporaryData['TT'] = TT
            self._CalculatedParams.append(['TT',{}])

        return self

    def uT(self, inMemory=None):
        """
        Calculates the mean of u'*T'.

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'uT' not in self._TemporaryData.columns:
            self.fluctuations()
            if self. SamplingWindow is None:
                uT = (self._RawData['up'] * self._RawData['Tp']).mean()
            else:
                uT = (self._RawData['up'] * self._RawData['Tp']).resample(self.SamplingWindow).mean()
            self._TemporaryData['uT'] = uT
            self._CalculatedParams.append(['uT',{}])

        return self

    def vT(self, inMemory=None):
        """
        Calculates the mean of v'*T'.

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'vT' not in self._TemporaryData.columns:
            self.fluctuations()
            if self. SamplingWindow is None:
                vT = (self._RawData['vp'] * self._RawData['Tp']).mean()
            else:
                vT = (self._RawData['vp'] * self._RawData['Tp']).resample(self.SamplingWindow).mean()
            self._TemporaryData['vT'] = vT
            self._CalculatedParams.append(['vT',{}])

        return self

    def wT(self, inMemory=None):
        """
        Calculates the mean of w'*T'.

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'wT' not in self._TemporaryData.columns:
            self.fluctuations()
            if self. SamplingWindow is None:
                wT = (self._RawData['wp'] * self._RawData['Tp']).mean()
            else:
                wT = (self._RawData['wp'] * self._RawData['Tp']).resample(self.SamplingWindow).mean()
            self._TemporaryData['wT'] = wT
            self._CalculatedParams.append(['wT',{}])

        return self

    def uv(self, inMemory=None):
        """
        Calculates the mean of u'*v'.

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'uv' not in self._TemporaryData.columns:
            self.fluctuations()
            if self.SamplingWindow is None:
                uv = (self._RawData['up'] * self._RawData['vp']).mean()
            else:
                uv = (self._RawData['up'] * self._RawData['vp']).resample(self.SamplingWindow).mean()
            self._TemporaryData['uv'] = uv
            self._CalculatedParams.append(['uv',{}])

        return self

    def uw(self, inMemory=None):
        """
        Calculates the mean of u'*w'.

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'uw' not in self._TemporaryData.columns:
            self.fluctuations()
            if self.SamplingWindow is None:
                uw = (self._RawData['up'] * self._RawData['wp']).mean()
            else:
                uw = (self._RawData['up'] * self._RawData['wp']).resample(self.SamplingWindow).mean()
            self._TemporaryData['uw'] = uw
            self._CalculatedParams.append(['uw',{}])

        return self

    def vw(self, inMemory=None):
        """
        Calculates the mean of v'*w'.

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'vw' not in self._TemporaryData.columns:
            self.fluctuations()
            if self.SamplingWindow is None:
                vw = (self._RawData['vp'] * self._RawData['wp']).mean()
            else:
                vw = (self._RawData['vp'] * self._RawData['wp']).resample(self.SamplingWindow).mean()
            self._TemporaryData['vw'] = vw
            self._CalculatedParams.append(['vw',{}])

        return self

    def w3(self, inMemory=None):
        """
        Calculates the mean of w'^3.

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'w3' not in self._TemporaryData.columns:
            self.fluctuations()
            if self.SamplingWindow is None:
                www = (self._RawData['wp'] ** 3).mean()
            else:
                www = (self._RawData['wp'] ** 3).resample(self.SamplingWindow).mean()
            self._TemporaryData['w3'] = www
            self._CalculatedParams.append(['w3',{}])

        return self

    def w4(self, inMemory=None):
        """
        Calculates the mean of w'^4.

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'w4' not in self._TemporaryData.columns:
            self.fluctuations()
            if self.SamplingWindow is None:
                wwww = (self._RawData['wp'] ** 4).mean()
            else:
                wwww = (self._RawData['wp'] ** 4).resample(self.SamplingWindow).mean()
            self._TemporaryData['w4'] = wwww
            self._CalculatedParams.append(['w4',{}])

        return self

    def TKE(self, inMemory=None):
        """
        Calculates the turbulence kinetic energy.

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'TKE' not in self._TemporaryData.columns:
            self.uu().vv().ww()
            TKE = 0.5 * (self._TemporaryData['uu'] + self._TemporaryData['vv'] + self._TemporaryData['ww'])
            self._TemporaryData['TKE'] = TKE
            self._CalculatedParams.append(['TKE',{}])

        return self

    def wTKE(self, inMemory=None):
        """

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'wTKE' not in self._TemporaryData.columns:
            self.fluctuations()
            uu = self._RawData['up'] ** 2
            vv = self._RawData['vp'] ** 2
            ww = self._RawData['wp'] ** 2
            wp = self._RawData['wp']
            if self.SamplingWindow is None:
                wTKE = (0.5 * (uu + vv + ww) * wp).mean()
            else:
                wTKE = (0.5 * (uu + vv + ww) * wp).resample(self.SamplingWindow).mean()
            self._TemporaryData['wTKE'] = wTKE
            self._CalculatedParams.append(['wTKE',{}])

        return self

    def Ustar(self, inMemory=None):
        """

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'Ustar' not in self._TemporaryData.columns:
            self.uw().vw()
            Ustar = (self._TemporaryData['uw'] ** 2 + self._TemporaryData['vw'] ** 2) ** 0.25
            self._TemporaryData['Ustar'] = Ustar
            self._CalculatedParams.append(['Ustar',{}])

        return self

    def Rvw(self, inMemory=None):
        """

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'Rvw' not in self._TemporaryData.columns:
            self.vw().vv().ww()
            Rvw = self._TemporaryData['vw'] / numpy.sqrt(self._TemporaryData['vv'] * self._TemporaryData['ww'])
            self._TemporaryData['Rvw'] = Rvw
            self._CalculatedParams.append(['Rvw',{}])

        return self

    def Ruw(self, inMemory=None):
        """

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'Ruw' not in self._TemporaryData.columns:
            self.uw().uu().ww()
            Ruw = self._TemporaryData['uw'] / numpy.sqrt(self._TemporaryData['uu'] * self._TemporaryData['ww'])
            self._TemporaryData['Ruw'] = Ruw
            self._CalculatedParams.append(['Ruw',{}])

        return self

    def MOLength_Sonic(self, inMemory=None):
        """
        Calculates the Monin-Obukhov length. The mean temperature used is the one calculated from the Sonic raw data
        (less accurate than using TRH).

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'L_Sonic' not in self._TemporaryData.columns:
            self.wT().Ustar()
            L = -(self._TemporaryData['T_bar']+273.15) * self._TemporaryData['Ustar'] ** 3 / (
                        self.Karman * g * self._TemporaryData['wT'])
            self._TemporaryData['L_Sonic'] = L
            self._CalculatedParams.append(['L_Sonic',{}])

        return self

    def zoL_Sonic(self, zmd, inMemory=None):
        """

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        zmd: float
            Height.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        i = 1

        while 'zoL_Sonic%s' % i in self._TemporaryData.columns:
            if ['zoL_Sonic%s' % i, {'zmd': zmd}] in self._AllCalculatedParams:
                return self
            i += 1

        self.MOLength_Sonic()
        zoL = zmd / self._TemporaryData['L_Sonic']
        self._TemporaryData['zoL_Sonic%s' % i] = zoL
        self._CalculatedParams.append(['zoL_Sonic%s' % i, {'zmd': zmd}])

        return self

    def zOverL_Sonic(self, inMemory=None):
        """

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'zOverL_Sonic' not in self._TemporaryData.columns:
            self.MOLength_Sonic()

            H = int(self.metaData['buildingHeight'])
            instrumentHeight = int(self.metaData['height'])
            averagedHeight = int(self.metaData['averagedHeight'])
            effectivez = instrumentHeight + H - 0.7 * averagedHeight
            zOverL = effectivez / self._TemporaryData['L_Sonic']
            self._TemporaryData['zOverL_Sonic'] = zOverL
            self._CalculatedParams.append(['zOverL_Sonic',{}])

        return self

    def Lminus1_masked_Sonic(self, inMemory=None):
        """

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'Lminus1_masked_Sonic' not in self._TemporaryData.columns:
            self.MOLength_Sonic()
            mask = ((numpy.abs(self._TemporaryData['wT']) > 0.05) & (numpy.abs(self._TemporaryData['Ustar']) > 0.15))
            maskedData = self._TemporaryData[mask]
            Lminus1_masked = -self.Karman * (g / maskedData['T_bar']) * maskedData['wT'] / maskedData['Ustar'] ** 3
            self._TemporaryData['Lminus1_masked_Sonic'] = Lminus1_masked
            self._CalculatedParams.append('Lminus1_masked_Sonic')

        return self

    def StabilityMOLength_Sonic(self, inMemory=None):
        """
        Calculates the MOlength stability.

        Parameters
        ----------
        inMemory : boolean
            Default value is None.

        Returns
        -------
        singlePointTurbulenceStatistics
            The object himself.
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'StabilityMOLength_Sonic' not in self._TemporaryData.columns:
            self.MOLength_Sonic()
            stability = self._TemporaryData['L_Sonic'].apply(self._ClassifyStability) if self._DataType == 'pandas' \
                else self._TemporaryData['L_Sonic'].apply(self._ClassifyStability, meta='str')
            self._TemporaryData['StabilityMOLength_Sonic'] = stability
            self._CalculatedParams.append(['StabilityMOLength_Sonic',{}])

        return self

    def _ClassifyStability(self, L):
        """
            According to 1/L categories:
            0 - Very Unstable
            1 - Unstable
            2 - Near Neutral
            3 - Stable
            4 - Very Stable
        """

        # For Z_0=1 (Irwin1979 Table 1)
        ret = 0
        if L is None:
            return "No Stability"
        if numpy.isnan(L):
            return "No Stability"
        if 1. / L < -.0875:
            ret = "very unstable"  # very un stable (A)
        elif 1. / L < -0.0081:
            ret = "unstable"  # un stable (C,B)
        elif 1. / L < 0.0081:
            ret = "neutral/near neutral"  # Neutral/Near Neutral (D)
        elif 1. / L < 0.25:  # (Mahrt1999: z/L>O(1)) #(z-d)/L<0.1667 from Delft Conference
            ret = "stable"  # stable (E,F)
        else:
            ret = "very stable"  # very stable (G)

        return ret

    def StrucFunDir(self, tau_range = None, dir1_data = None, u_dir1 = "u_dir1", v_dir1 = "v_dir1", w_dir1 = "w_dir1",
                    dir2_data = None, u_dir2 = "u_dir2", v_dir2 = "v_dir2", w_dir2 = "w_dir2", title = "", inMemory = None):

        """
        Calculates the 2nd order structure function <(ui(t+tau)-ui(t)) * (uj(t+tau)-uj(t))> for a given list of tau values,
        on given directions for ui, uj (which might be different for each temporal window).

        :param tau_range: list of tau values (floats)

        :param dir1_data: pandas dataframe containing the direction vectors on which ui is calculated (index not necessarily
                            with the same frequency as the sampling window, but the sampling window is assumed to divide
                            the index frequency of dir1_data), vectors not necesarily normalized

        :param u_dir1: str name of column in dir1_data containing the x direction of the direction vectors

        :param v_dir1: str name of column in dir1_data containing the y direction of the direction vectors

        :param w_dir1: str name of column in dir1_data containing the z direction of the direction vectors

        :param dir2_data: same for uj - if not given, same directions as in dir1_data will be used

        :param u_dir2: same for uj

        :param v_dir2: same for uj

        :param w_dir2: same for uj

        :param title: title to add for the structure function column in the returned data (see col_names below)

        :param inMemory:
        :return:
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if dir2_data is None:
            dir2_data = dir1_data
            u_dir2 = u_dir1
            v_dir2 = v_dir1
            w_dir2 = w_dir1

        col_names = {tau:"D" + title + "_" + str(tau) + "s" for tau in tau_range}
        if set(col_names.values()).issubset(set(self._TemporaryData.columns)):
            return self

        self.fluctuations()

        # Extracting the supplied direction data of ui,uj:
        dir1_data_new = dir1_data[[u_dir1,v_dir1,w_dir1]]\
                       .rename(columns = {u_dir1: "u_dir1", v_dir1: "v_dir1", w_dir1: "w_dir1"})\
                       .loc[(dir1_data.index >= self.metaData["start"]) & (dir1_data.index < self.metaData["end"])]

        dir2_data_new = dir2_data[[u_dir2, v_dir2, w_dir2]] \
            .rename(columns={u_dir2: "u_dir2", v_dir2: "v_dir2", w_dir2: "w_dir2"}) \
            .loc[(dir2_data.index >= self.metaData["start"]) & (dir2_data.index < self.metaData["end"])]

        # Computing the magnitude of the direction vectors and normalizing them:
        dir1_data_new["dir1_mag"] = (dir1_data_new["u_dir1"] ** 2 + dir1_data_new["v_dir1"] ** 2 + dir1_data_new["w_dir1"] ** 2) ** 0.5

        dir1_data_new.loc[:,['u_dir1','v_dir1','w_dir1']] = dir1_data_new.loc[:,['u_dir1','v_dir1','w_dir1']].div(dir1_data_new["dir1_mag"], axis=0)

        dir1_data_new = dir1_data_new.drop(columns = 'dir1_mag')

        dir2_data_new["dir2_mag"] = (dir2_data_new["u_dir2"] ** 2 + dir2_data_new["v_dir2"] ** 2 + dir2_data_new["w_dir2"] ** 2) ** 0.5

        dir2_data_new.loc[:,['u_dir2','v_dir2','w_dir2']] = dir2_data_new.loc[:,['u_dir2','v_dir2','w_dir2']].div(dir2_data_new["dir2_mag"], axis=0)

        dir2_data_new = dir2_data_new.drop(columns = 'dir2_mag')

        # Creating the temporary raw data on which calculations are performed (u,v,w + direction data):
        united_data = self._RawData.merge(dir1_data_new, how = "outer", left_index = True, right_index = True)\
                           .merge(dir2_data_new, how = "outer", left_index = True, right_index = True)\
                           .dropna(how='all')

        united_data[["u_dir1", "v_dir1", "w_dir1", "u_dir2", "v_dir2", "w_dir2"]] = \
            united_data[["u_dir1", "v_dir1", "w_dir1", "u_dir2", "v_dir2", "w_dir2"]].ffill()

        united_data = united_data.dropna(how='any')

        # Computing ui,uj as projections of u,v,w on the direction vectors:
        united_data["ui"] = 0
        united_data["uj"] = 0
        for component in ["u","v","w"]:
            united_data["ui"] += united_data[component] * united_data["%s_dir1" % component]
            united_data["uj"] += united_data[component] * united_data["%s_dir2" % component]

        for tau in tau_range:
            # Computation of <(ui(t+tau)-ui(t)) * (uj(t+tau)-uj(t))>:
            if col_names[tau] not in self._TemporaryData.columns:
                # Computing the data of ui,uj(t+tau):
                data_tau = united_data[["ui","uj"]].reset_index()
                data_tau["Time"] -= pandas.Timedelta(tau, unit="s")
                data_tau = data_tau.set_index("Time").rename(columns = {"ui":"ui_shifted","uj":"uj_shifted"})

                # Drop potentially used columns from last tau:
                to_drop = set(united_data.columns) & {"ui_shifted", "uj_shifted"}
                if bool(to_drop):
                    united_data = united_data.drop(columns = to_drop)

                # Merge data(t) with data(t+tau), repartition (otherwise gets crazy):
                united_data = united_data.merge(data_tau, how = "left", left_index = True, right_index = True).repartition(freq = "1W")
                # Final calculation:
                self._TemporaryData[col_names[tau]] = ((united_data["ui_shifted"] - united_data["ui"]) *\
                                                 (united_data["uj_shifted"] - united_data["uj"])).resample(self.SamplingWindow).mean()
                self._CalculatedParams.append([col_names[tau],{}])

        return self

    def StrucFun(self, tau_range = None, ubar_data = None, u_bar = "u_bar", v_bar = "v_bar", w_bar = "w_bar",
                     mode = "MeanDir", title_additions = "", inMemory = None):

        """
        Calculates (at least 1 component of) the 2nd order structure function <(ui(t+tau)-ui(t)) * (uj(t+tau)-uj(t))>
        for a given list of tau values, with axes aligned with the 3d mean velocity direction.

        :param tau_range: range of tau values

        :param ubar_data: pandas dataframe of externally supplied mean velocity data. Default is the same mean velocity
                            computed using fluctuations. The frequency of ubar_data.index is assumed to be divisible by
                            the sampling window.

        :param u_bar: str name of column in ubar_data that contains the x component of the mean velocity

        :param v_bar: str name of column in ubar_data that contains the y component of the mean velocity

        :param w_bar: str name of col   umn in ubar_data that contains the z component of the mean velocity

        :param mode: either default "MeanDir" or "3dMeanDir"
                        MeanDir - calculate only D11, where the 1 component is aligned with the mean velocity 3d direction
                        3dMeanDir - calculate all 6 components Dij, i,j=1,2,3, where the axes are defined as:
                                    1 - mean velocity direction
                                    2 - a horizontal direction perpendicular to the mean velocity (z cross mean vel dir)
                                    3 - a direction perpendicular to 1,2 (dir 1 cross dir 2)

        :param title_additions: str additions to the column titles of the resulting computations

        :param inMemory:
        :return:
        """
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if ubar_data is None:
            # Default ubar_data
            self.fluctuations()
            ubar_data = self._TemporaryData[["u_bar","v_bar","w_bar"]].compute()

        if mode == "MeanDir":
            self.StrucFunDir(tau_range=tau_range, dir1_data=ubar_data, u_dir1=u_bar, v_dir1=v_bar, w_dir1=w_bar, title = "11" + title_additions)

        elif mode == "3dMeanDir":
            # Create direction data for all 3 directions (directions are not normalized since StrucFunDir takes care of that)
            ubar_new = ubar_data[[u_bar,v_bar,w_bar]].rename(columns = {u_bar:"x_hat1",v_bar:"x_hat2",w_bar:"x_hat3"})
            ubar_new["y_hat1"] = - ubar_new["x_hat2"]
            ubar_new["y_hat2"] = ubar_new["x_hat1"]
            ubar_new["y_hat3"] = 0
            ubar_new["z_hat1"] = - ubar_new["x_hat1"] * ubar_new["x_hat3"]
            ubar_new["z_hat2"] = - ubar_new["x_hat2"] * ubar_new["x_hat3"]
            ubar_new["z_hat3"] = ubar_new["x_hat1"] ** 2 +  ubar_new["x_hat2"] ** 2

            # Compute for each ij
            self.StrucFunDir(tau_range=tau_range, dir1_data=ubar_new, u_dir1="x_hat1", v_dir1="x_hat2", w_dir1="x_hat3", title="11" + title_additions)
            self.StrucFunDir(tau_range=tau_range, dir1_data=ubar_new, u_dir1="y_hat1", v_dir1="y_hat2", w_dir1="y_hat3", title="22" + title_additions)
            self.StrucFunDir(tau_range=tau_range, dir1_data=ubar_new, u_dir1="z_hat1", v_dir1="z_hat2", w_dir1="z_hat3", title="33" + title_additions)
            self.StrucFunDir(tau_range=tau_range, dir1_data=ubar_new, u_dir1="x_hat1", v_dir1="x_hat2", w_dir1="x_hat3",
                             dir2_data = ubar_new, u_dir2="y_hat1", v_dir2="y_hat2", w_dir2="y_hat3",title="12" + title_additions)
            self.StrucFunDir(tau_range=tau_range, dir1_data=ubar_new, u_dir1="x_hat1", v_dir1="x_hat2", w_dir1="x_hat3",
                             dir2_data = ubar_new, u_dir2="z_hat1", v_dir2="z_hat2", w_dir2="z_hat3",title="13" + title_additions)
            self.StrucFunDir(tau_range=tau_range, dir1_data=ubar_new, u_dir1="y_hat1", v_dir1="y_hat2", w_dir1="y_hat3",
                             dir2_data = ubar_new, u_dir2="z_hat1", v_dir2="z_hat2", w_dir2="z_hat3",title="23" + title_additions)
        else:
            raise("mode must be either MeanDir or 3dMeanDir")

        if "u_mag" not in self.TemporaryData.columns:
            self._TemporaryData["u_mag" + title_additions] = ((ubar_data[u_bar] ** 2 + ubar_data[v_bar] ** 2 + ubar_data[w_bar] ** 2) ** 0.5).loc[(ubar_data.index >=
                                    self.metaData["start"]) & (ubar_data.index < self.metaData["end"])]
            self._TemporaryData["u_mag" + title_additions] = self._TemporaryData["u_mag" + title_additions].ffill()
            self._CalculatedParams.append(["u_mag" + title_additions,{}])
        return self

    def ThirdStrucFun(self, tau_range = None, ubar_data = None, u_bar = "u_bar", v_bar = "v_bar", w_bar = "w_bar",
                     title_additions = "", inMemory = None):

        """
        Calculates the 3nd order structure function <(u1(t+tau)-u1(t))^3> for a given list of tau values, where u1 is the
        velocity component along the mean velocity direction (which might be different for each temporal window).

        :param tau_range: list of tau values (floats)

        :param ubar_data: pandas dataframe of externally supplied mean velocity data. Default is the same mean velocity
                            computed using fluctuations. The frequency of ubar_data.index is assumed to be divisible by
                            the sampling window.

        :param u_bar: str name of column in ubar_data that contains the x component of the mean velocity

        :param v_bar: str name of column in ubar_data that contains the y component of the mean velocity

        :param w_bar: str name of column in ubar_data that contains the z component of the mean velocity

        :param title_additions: str additions to the column titles of the resulting computations

        :param inMemory:

        :return:
        """

        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        self.fluctuations()
        if ubar_data is None:
            # Default ubar_data
            ubar_data = self._TemporaryData[["u_bar","v_bar","w_bar"]].compute()

        col_names = {tau:"D111" + title_additions + "_" + str(tau) + "s" for tau in tau_range}
        if set(col_names.values()).issubset(set(self._TemporaryData.columns)):
            return self

        # Extracting the supplied mean velocity data:
        dir_data = ubar_data[[u_bar,v_bar,w_bar]].rename(columns={u_bar: "u_dir", v_bar: "v_dir", w_bar: "w_dir"}).loc[
            (ubar_data.index >= self.metaData["start"]) & (ubar_data.index < self.metaData["end"])]

        # Computing the velocity magnitude and normalizing the velocity vectors:
        dir_data["u_mag"] = (dir_data["u_dir"] ** 2 + dir_data["v_dir"] ** 2 + dir_data["w_dir"] ** 2) ** 0.5
        dir_data.loc[:, ['u_dir', 'v_dir', 'w_dir']] = dir_data.loc[:, ['u_dir', 'v_dir', 'w_dir']].div(dir_data["u_mag"], axis=0)

        # Creating the temporary raw data on which calculations are performed (u,v,w + direction data):
        united_data = self._RawData.merge(dir_data[['u_dir', 'v_dir', 'w_dir']], how="outer", left_index=True, right_index=True)\
                .dropna(how='all')
        united_data[["u_dir", "v_dir", "w_dir"]] = united_data[["u_dir", "v_dir", "w_dir"]].ffill()
        united_data = united_data.dropna(how='any')

        # Computing u1 as a projection of u,v,w on the direction vectors:
        united_data["u1"] = 0
        for component in ["u", "v", "w"]:
            united_data["u1"] += united_data[component] * united_data["%s_dir" % component]

        for tau in tau_range:
            # Computation of <(u1(t+tau)-u1(t))^3>:
            if col_names[tau] not in self._TemporaryData.columns:
                # Computing the data of u1(t+tau):
                data_tau = united_data[["u1"]].reset_index()
                data_tau["Time"] -= pandas.Timedelta(tau, unit="s")
                data_tau = data_tau.set_index("Time").rename(columns = {"u1":"u1_shifted"})

                # Drop potentially used columns from last tau:
                if "u1_shifted" in united_data.columns:
                    united_data = united_data.drop(columns = "u1_shifted")

                # Merge data(t) with data(t+tau), repartition (otherwise gets crazy):
                united_data = united_data.merge(data_tau, how = "left", left_index = True, right_index = True).repartition(freq = "1W")
                # Final calculation:
                self._TemporaryData[col_names[tau]] = ((united_data["u1_shifted"] - united_data["u1"]) ** 3)\
                                                        .resample(self.SamplingWindow).mean()
                self._CalculatedParams.append([col_names[tau],{}])

        if "u_mag" not in self.TemporaryData.columns:
            self._TemporaryData["u_mag" + title_additions] = dir_data["u_mag"]
            self._TemporaryData["u_mag" + title_additions] = self._TemporaryData["u_mag" + title_additions].ffill()
            self._CalculatedParams.append(["u_mag" + title_additions,{}])

        return self


class SinglePointStatisticsSpark(singlePointTurbulenceStatistics):

    def fluctuations(self, inMemory=None):
        if self._InMemoryAvgRef is None:
            self._InMemoryAvgRef = inMemory

        if 'up' not in self._RawData.columns:
            avg = self._RawData
            if self.SamplingWindow is None:
                avg = avg.mean()
                if self._DataType == 'pandas':
                    avg = pandas.DataFrame(avg).T
                    avg.index = [self._RawData.index[0]]
                else:
                    avg = pandas.DataFrame(avg.compute()).T
                    avg.index = self._RawData.head(1).index
                    npartitions = self._RawData.npartitions
                    avg = dask.dataframe.from_pandas(avg, npartitions=npartitions)
            else:
                avg = avg.resample(self.SamplingWindow).mean()

            avg = avg.rename(columns={'u': 'u_bar', 'v': 'v_bar', 'w': 'w_bar', 'T': 'T_bar'})

            avg['wind_dir_bar'] = numpy.arctan2(avg['v_bar'], avg['u_bar'])
            avg['wind_dir_bar'] = (2 * numpy.pi + avg['wind_dir_bar']) % (2 * numpy.pi)
            avg['wind_dir_bar'] = numpy.rad2deg(avg['wind_dir_bar'])

            avg['wind_dir_bar'] = avg['wind_dir_bar'].apply(lambda x: 270 - x if 270 - x >= 0 else 630 - x)

            self._TemporaryData = avg
            self._CalculatedParams += [['u_bar',{}], ['v_bar',{}], ['w_bar',{}], ['T_bar',{}]]

            # correcting the first index to be the same as the avg.
            self._RawData = self._RawData.reset_index()
            self._RawData.at[0,'Time'] = avg.index[0]
            self._RawData = self._RawData.set_index("Time")

            self._RawData = self._RawData.merge(avg, how='left', left_index=True, right_index=True)
            self._RawData = self._RawData.ffill()

            self._RawData['wind_dir'] = numpy.arctan2(self._RawData['v'], self._RawData['u'])
            self._RawData['wind_dir'] = (2 * numpy.pi + self._RawData['wind_dir']) % (2 * numpy.pi)
            self._RawData['wind_dir'] = numpy.rad2deg(self._RawData['wind_dir'])
            self._RawData['wind_dir'] = self._RawData['wind_dir'].apply(lambda x: 270 - x if 270 - x >= 0 else 630 - x)

            self._RawData['up'] = self._RawData['u'] - self._RawData['u_bar']
            self._RawData['vp'] = self._RawData['v'] - self._RawData['v_bar']
            self._RawData['wp'] = self._RawData['w'] - self._RawData['w_bar']
            self._RawData['Tp'] = self._RawData['T'] - self._RawData['T_bar']
            self._RawData['wind_dir_p'] = (180 - (180 - (self._RawData['wind_dir'] - self._RawData['wind_dir_bar']).abs()).abs()).abs()

        return self


class InMemoryRawData(pandas.DataFrame):
    _Attrs = None

    def __init__(self, data=None, index=None, columns=None, dtype=None, copy=False):
        super(InMemoryRawData, self).__init__(data=data, index=index, columns=columns, dtype=dtype, copy=copy)
        self._Attrs = {}

    def append(self, other, ignore_index=False, verify_integrity=False):
        ret = super(InMemoryRawData, self).append(other, ignore_index=ignore_index, verify_integrity=verify_integrity)
        ret = InMemoryRawData(ret)
        ret._Attrs = other._Attrs
        ret._Attrs.update(self._Attrs)

        return ret

    @classmethod
    def read_hdf(cls, path_or_buf, key=None, **kwargs):
        ret = InMemoryRawData(pandas.read_hdf(path_or_buf, key, **kwargs))
        path_or_buf = '%s%s' % (path_or_buf.rpartition('.')[0], '.json')

        if os.path.isfile(path_or_buf):
            with open(path_or_buf, 'r') as jsonFile:
                ret._Attrs = json.load(jsonFile)

        return ret

    def to_hdf(self, path_or_buf, key, **kwargs):
        pandasCopy = self.copy()
        path_or_buf = '%s%s' % (path_or_buf.rpartition('.')[0], '.hdf')
        pandasCopy.to_hdf(path_or_buf, key, **kwargs)
        path_or_buf = '%s%s' % (path_or_buf.rpartition('.')[0], '.json')
        attrsToSave = self._Attrs

        if len(self._Attrs) > 0:
            if os.path.isfile(path_or_buf):
                with open(path_or_buf, 'r') as jsonFile:
                    attrsFile = json.load(jsonFile)
                    attrsFile.update(attrsToSave)
                    attrsToSave = attrsFile

            with open(path_or_buf, 'w') as jsonFile:
                json.dump(attrsToSave, jsonFile, indent=4, sort_keys=True)


class InMemoryAvgData(InMemoryRawData):
    _TurbulenceCalculator = None

    def __init__(self, data = None, index = None, columns = None, dtype = None, copy = False, turbulenceCalculator = None):
        super(InMemoryAvgData, self).__init__(data = data, index = index, columns = columns, dtype = dtype, copy = copy)
        self._TurbulenceCalculator = turbulenceCalculator
        self._Attrs['samplingWindow'] = turbulenceCalculator.SamplingWindow

    def __getattr__(self, item):
        if self._TurbulenceCalculator is None:
            raise AttributeError("The attribute '_TurbulenceCalculator' is None.")
        elif not item in dir(self._TurbulenceCalculator):
            raise NotImplementedError("The attribute '%s' is not implemented." % item)
        elif item == 'compute':
            ret = getattr(self._TurbulenceCalculator, item)
        else:
            ret = lambda *args, **kwargs: getattr(self._TurbulenceCalculator, item)(inMemory = self, *args, **kwargs)

        return ret

