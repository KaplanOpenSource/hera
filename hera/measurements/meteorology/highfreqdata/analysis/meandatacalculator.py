import numpy
import pandas
import dask
from copy import deepcopy
from scipy.constants import g
from .abstractcalculator import AbstractCalculator
from .turbulencestatistics import singlePointTurbulenceStatistics
from hera.utils.filter_immediate import Filter

class AveragingCalculator(AbstractCalculator):
    """Calculator that resamples raw data to compute time-averaged means.

    Resamples the input data at the sampling window interval and computes
    the mean of each column, appending ``_bar`` to column names.
    """

    def __init__(self, rawData, metadata):
        """Initialise and immediately compute windowed averages.

        Parameters
        ----------
        rawData : pandas.DataFrame or dask.dataframe.DataFrame
            High-frequency time-series data.
        metadata : dict
            Must contain ``'samplingWindow'`` (e.g. ``'30min'``).
        """
        super(AveragingCalculator, self).__init__(rawData=rawData, metadata=metadata)

        self._TemporaryData = self._RawData.resample(self.SamplingWindow).mean().rename(
            columns={col: col + "_bar" for col in self._RawData.columns})
        for col in self._TemporaryData.columns:
            self._CalculatedParams.append([col, {}])

        self.data = self._TemporaryData

    def getData(self):
        """Return the time-averaged DataFrame.

        Returns
        -------
        pandas.DataFrame
            Columns named ``<original>_bar`` containing windowed means.
        """
        return self.data


class MeanDataCalculator:
    """Computes mean-field turbulence statistics from second moments.

    Combines turbulence calculator output (second moments) with optionally
    time-averaged mean fields. Provides a fluent API for computing derived
    quantities (friction velocity, stability parameters, anisotropy, etc.)
    via method chaining.
    """

    def __init__(self, TurbCalcOrData = None, compute_mode_turb = 'not_from_db_and_not_save', AverageCalcOrData = None,
                 compute_mode_AverageCalc = None, **metadata):
        """Initialise with turbulence data and optional mean fields.

        Parameters
        ----------
        TurbCalcOrData : singlePointTurbulenceStatistics or pandas.DataFrame or dask.DataFrame
            Source of second-moment data.
        compute_mode_turb : str
            Compute mode for the turbulence calculator.
        AverageCalcOrData : AveragingCalculator or pandas.DataFrame or dask.DataFrame, optional
            Source of time-averaged mean fields (``u_bar``, ``v_bar``, etc.).
        compute_mode_AverageCalc : str, optional
            Compute mode for the averaging calculator. Defaults to *compute_mode_turb*.
        **metadata
            Experiment metadata (``start``, ``end``, ``height``, etc.).
        """

        if compute_mode_AverageCalc is None:
            compute_mode_AverageCalc = compute_mode_turb

        self._Karman = 0.4

        if type(TurbCalcOrData) == pandas:
            self.metaData = deepcopy(metadata)
            self.MeanData = TurbCalcOrData.copy()
        elif type(TurbCalcOrData) == dask:
            self.metaData = deepcopy(metadata)
            self.MeanData = TurbCalcOrData[self.metaData["start"]:self.metaData["end"]].compute()
        elif type(TurbCalcOrData) == singlePointTurbulenceStatistics:
            self.TurbCalc = TurbCalcOrData
            self.metaData = deepcopy(self.TurbCalc.metaData)
            self.metaData.update(deepcopy(metadata))
            self.MeanData = self.TurbCalc.secondMoments().compute(mode=compute_mode_turb)
        else:
            raise ValueError("TurbCalcOrData must be either a singlePointTurbulenceStatistics instance or a pandas/dask dataframe")

        if type(AverageCalcOrData) == pandas:
            AverageData = AverageCalcOrData
        elif type(AverageCalcOrData) == dask:
            AverageData = AverageCalcOrData[self.metaData["start"]:self.metaData["end"]].compute()
        elif type(AverageCalcOrData) == AveragingCalculator:
            AverageData = AverageCalcOrData.compute(mode=compute_mode_AverageCalc)
        elif AverageCalcOrData is None:
            AverageData = None
        else:
            raise ValueError(
                "AverageCalcOrData must be either a singlePointTurbulenceStatistics instance or a pandas/dask dataframe")

        if AverageData is not None:
            self.MeanData = self.MeanData.join(AverageData)

        self.MeanData = self.MeanData.loc[(self.MeanData.index >= self.metaData["start"]) &
                                          (self.MeanData.index < self.metaData["end"])]

    def thresholds(self, threshold_list, inplace = False):
        """

        :param threshold_list: Format - [("u_bar","lt",20), ... , (field,preposition,bound)]
        :param inplace:
        :return:
        """

        filter_obj = Filter(data = self.MeanData, inplace = True)

        for cut in threshold_list:
            filter_obj.threshold(column = cut[0], preposition= cut[1], bound= cut[2])

        if inplace:
            self.MeanData = filter_obj.data
            return self
        else:
            return MeanDataCalculator(filter_obj.data, **self.metaData)

    def filterDates(self, start, end, inplace = False):
        """Filter data to a date range.

        Parameters
        ----------
        start : str or pandas.Timestamp
            Start of the date range.
        end : str or pandas.Timestamp
            End of the date range.
        inplace : bool
            If ``True``, modify this instance. Otherwise return a new one.

        Returns
        -------
        MeanDataCalculator
        """
        filter_obj = Filter(data = self.MeanData, inplace = True)

        filter_obj.outsideInterval(lower_bound = start, upper_bound = end)

        if inplace:
            self.MeanData = filter_obj.data
            return self
        else:
            return MeanDataCalculator(filter_obj.data, **self.metaData)


    def hour(self):
        """Add an ``hour`` column extracted from the time index.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if "hour" not in self.MeanData.columns:
            self.MeanData["hour"] = self.MeanData.index.hour

        return self

    def timeWithinDay(self):
        """Add a ``timeWithinDay`` column as fractional hours (hour + min/60 + sec/3600).

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if "timeWithinDay" not in self.MeanData.columns:
            self.MeanData["timeWithinDay"] = self.MeanData.index.hour + self.MeanData.index.minute / 60 + \
                                             self.MeanData.index.second / 3600

        return self

    def horizontalSpeed(self):
        """Compute mean horizontal wind speed from u_bar and v_bar using ``hypot``.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if 'horizontal_speed_bar' not in self.MeanData.columns:
            self.MeanData['horizontal_speed_bar'] = numpy.hypot(self.MeanData['u_bar'], self.MeanData['v_bar'])

        return self

    def _UV_to_SpdDir(self,U, V):
        return (U ** 2 + V ** 2) ** 0.5, (-numpy.degrees(numpy.arctan2(V, U)) + 90) % 360

    def alignedStress(self):
        """Rotate the Reynolds stress tensor to align with the mean wind direction.

        Adds columns: ``uu_aligned``, ``uv_aligned``, ``vv_aligned``,
        ``uw_aligned``, ``vw_aligned``.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if "uu_aligned" not in self.MeanData.columns:
            self.horizontalSpeed()
            cos_theta = self.MeanData["u_bar"] / self.MeanData["horizontal_speed_bar"]
            sin_theta = self.MeanData["v_bar"] / self.MeanData["horizontal_speed_bar"]
            self.MeanData["uu_aligned"] = (cos_theta ** 2) * self.MeanData["uu"] + 2 * sin_theta * cos_theta * self.MeanData["uv"] + (
                        sin_theta ** 2) * self.MeanData["vv"]
            self.MeanData["uv_aligned"] = - sin_theta * cos_theta * self.MeanData["uu"] + (cos_theta ** 2 - sin_theta ** 2) * \
                                      self.MeanData["uv"] + sin_theta * cos_theta * self.MeanData["vv"]
            self.MeanData["vv_aligned"] = (sin_theta ** 2) * self.MeanData["uu"] - 2 * sin_theta * cos_theta * self.MeanData["uv"] + (
                        cos_theta ** 2) * self.MeanData["vv"]
            self.MeanData["uw_aligned"] = cos_theta * self.MeanData["uw"] + sin_theta * self.MeanData["vw"]
            self.MeanData["vw_aligned"] = - sin_theta * self.MeanData["uw"] + cos_theta * self.MeanData["vw"]

        return self

    def sigma(self):
        """Compute velocity standard deviations: ``sigmaU``, ``sigmaV``, ``sigmaW``.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if 'sigmaU' not in self.MeanData.columns:
            self.MeanData['sigmaU'] = numpy.sqrt(self.MeanData["uu"])
            self.MeanData['sigmaV'] = numpy.sqrt(self.MeanData["vv"])
            self.MeanData['sigmaW'] = numpy.sqrt(self.MeanData["ww"])

        return self

    def sigmaAligned(self):
        """Compute wind-aligned standard deviations: ``sigmaU_aligned``, ``sigmaV_aligned``.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if "sigmaU_aligned" not in self.MeanData.columns:
            self.sigma()
            self.alignedStress()
            self.MeanData['sigmaU_aligned'] = numpy.sqrt(self.MeanData["uu_aligned"])
            self.MeanData['sigmaV_aligned'] = numpy.sqrt(self.MeanData["vv_aligned"])

        return self

    def sigmaH(self):
        """Compute horizontal velocity standard deviation ``sigmaH = hypot(sigmaU, sigmaV) / sqrt(2)``.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if 'sigmaH' not in self.MeanData.columns:
            self.sigma()
            self.MeanData['sigmaH'] = numpy.hypot(self.MeanData['sigmaU'], self.MeanData['sigmaV']) / numpy.sqrt(2)

        return self

    def Ustar(self):
        """Compute friction velocity ``u* = (uw² + vw²)^0.25``.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if 'Ustar' not in self.MeanData.columns:
            self.MeanData['Ustar'] = (self.MeanData['uw'] ** 2 + self.MeanData['vw'] ** 2) ** 0.25

        return self

    def sigmaHOverUstar(self):
        """Compute dimensionless ratio ``sigmaH / u*``.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if 'sigmaHOverUstar' not in self.MeanData.columns:
            self.sigmaH()
            self.Ustar()
            self.MeanData['sigmaHOverUstar'] = self.MeanData['sigmaH']/self.MeanData['Ustar']

        return self

    def sigmaUOverUstar(self):
        """Compute dimensionless ratio ``sigmaU_aligned / u*``.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if 'sigmaUOverUstar' not in self.MeanData.columns:
            self.sigmaAligned()
            self.Ustar()
            self.MeanData['sigmaUOverUstar'] = self.MeanData['sigmaU_aligned']/self.MeanData['Ustar']

        return self

    def sigmaVOverUstar(self):
        """Compute dimensionless ratio ``sigmaV_aligned / u*``.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if 'sigmaVOverUstar' not in self.MeanData.columns:
            self.sigmaAligned()
            self.Ustar()
            self.MeanData['sigmaVOverUstar'] = self.MeanData['sigmaV_aligned']/self.MeanData['Ustar']

        return self

    def sigmaWOverUstar(self):
        """Compute dimensionless ratio ``sigmaW / u*``.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if 'sigmaWOverUstar' not in self.MeanData.columns:
            self.sigma()
            self.Ustar()
            self.MeanData['sigmaWOverUstar'] = self.MeanData['sigmaW']/self.MeanData['Ustar']

        return self

    def sigmaHOverWindSpeed(self):
        """Compute turbulence intensity ``sigmaH / horizontal_speed_bar``.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if 'sigmaHOverWindSpeed' not in self.MeanData.columns:
            self.sigmaH()
            self.horizontalSpeed()
            self.MeanData['sigmaHOverWindSpeed'] = self.MeanData['sigmaH']/self.MeanData['horizontal_speed_bar']

        return self

    def sigmaWOverWindSpeed(self):
        """Compute vertical turbulence intensity ``sigmaW / horizontal_speed_bar``.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if 'sigmaWOverWindSpeed' not in self.MeanData.columns:
            self.sigma()
            self.horizontalSpeed()
            self.MeanData['sigmaWOverWindSpeed'] = self.MeanData['sigmaW']/self.MeanData['horizontal_speed_bar']

        return self

    def absWOverSigmaW(self):
        """Compute ratio ``|w_bar| / sigmaW``.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if "absWOverSigmaW" not in self.MeanData.columns:
            self.sigma()
            self.MeanData["absWOverSigmaW"] = abs(self.MeanData["w_bar"]) / self.MeanData["sigmaW"]

        return self

    def uStarOverWindSpeed(self):
        """Compute ratio ``u* / horizontal_speed_bar``.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if 'uStarOverWindSpeed' not in self.MeanData.columns:
            self.Ustar()
            self.horizontalSpeed()
            self.MeanData['uStarOverWindSpeed'] = self.MeanData['Ustar']/self.MeanData['horizontal_speed_bar']

        return self

    def TKE(self):
        """
        Calculates the turbulence kinetic energy.
        :return:
        """
        if 'TKE' not in self.MeanData.columns:
            self.MeanData['TKE'] = 0.5 * (self.MeanData['uu'] + self.MeanData['vv'] + self.MeanData['ww'])

        return self

    def Rvw(self):
        """Compute correlation coefficient between v' and w' fluctuations.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if 'Rvw' not in self.MeanData.columns:
            self.MeanData['Rvw'] = self.MeanData['vw'] / numpy.sqrt(self.MeanData['vv'] * self.MeanData['ww'])

        return self

    def Ruw(self):
        """Compute correlation coefficient between u' and w' fluctuations.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if 'Ruw' not in self.MeanData.columns:
            self.MeanData['Ruw'] = self.MeanData['uw'] / numpy.sqrt(self.MeanData['uu'] * self.MeanData['ww'])

        return self

    def MOLength(self):
        """
        Calculates the Monin-Obukhov length.
        """

        if 'L' not in self.MeanData.columns:
            self.Ustar()
            L = -(self.MeanData['TC_T_bar']+273.15) * self.MeanData['Ustar'] ** 3 / (
                        self._Karman * g * self.MeanData['wT'])
            self.MeanData['L'] = L

        return self

    # def zoL(self, zmd):
    #     """
    #     Parameters
    #     ----------
    #     zmd: float
    #         Height.
    #     """
    #
    #     i = 1
    #
    #     while 'zoL%s' % i in self.MeanData.columns:
    #         if ['zoL%s' % i, {'zmd': zmd}] in self._AllCalculatedParams:
    #             return self
    #         i += 1
    #
    #     self.MOLength()
    #     zoL = zmd / self.MeanData['L']
    #     self.MeanData['zoL%s' % i] = zoL
    #     self._CalculatedParams.append(['zoL%s' % i, {'zmd': zmd}])
    #
    #     return self

    def effectivez(self):
        """Compute effective measurement height accounting for buildings.

        ``z_eff = instrumentHeight + buildingHeight - 0.7 * averagedHeight``

        Stores the result in ``self.metaData['effectivez']``.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if "effectivez" not in self.metaData.keys():
            H = int(self.metaData['buildingHeight'])
            instrumentHeight = int(self.metaData['height'])
            averagedHeight = int(self.metaData['averagedHeight'])
            self.metaData["effectivez"] = instrumentHeight + H - 0.7 * averagedHeight

        return self

    def zOverL(self):
        """Compute dimensionless stability parameter ``z/L`` (effective height / Monin-Obukhov length).

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if 'zOverL' not in self.MeanData.columns:
            self.MOLength().effectivez()
            self.MeanData['zOverL'] = self.metaData["effectivez"] / self.MeanData['L']

        return self

    def StabilityMOLength(self):
        """
        Calculates the MOlength stability.
        """

        if 'StabilityMOLength' not in self.MeanData.columns:
            self.MOLength()
            self.MeanData['StabilityMOLength'] = self.MeanData['L'].apply(self._ClassifyStability)

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

    def _eig(series):
        if series["uu"] + series["vv"] + series["ww"] == 0:
            eig_ser = pandas.Series(data=numpy.nan, index=["lambda_1", "lambda_2", "lambda_3"])
        else:
            A = numpy.array([[series["uu"], series["uv"], series["uw"]], [series["uv"], series["vv"], series["vw"]],
                          [series["uw"], series["vw"], series["ww"]]]) / (
                            series["uu"] + series["vv"] + series["ww"]) - numpy.identity(3) / 3
            if numpy.isnan(A).any():
                eig_ser = pandas.Series(data=numpy.nan, index=["lambda_1", "lambda_2", "lambda_3"])
            else:
                eig_ser = pandas.Series(data=list(numpy.linalg.eigvalsh(A))[::-1], index=["lambda_1", "lambda_2", "lambda_3"])
        return eig_ser

    def anisotropyEigs(self):
        """Compute eigenvalues of the Reynolds stress anisotropy tensor.

        Adds columns: ``lambda_1``, ``lambda_2``, ``lambda_3``, ``x_B``, ``y_B``
        (Lumley triangle coordinates).

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if "lambda_1" not in self.MeanData.columns:
            eig_data = self.MeanData.apply(self._eig, axis=1)
            for col in eig_data.columns:
                self.MeanData[col] = eig_data[col]

            self.MeanData["x_B"] = self.MeanData["lambda_1"] - self.MeanData["lambda_2"] + 3 / 2 * self.MeanData[
                "lambda_3"] + 1 / 2
            self.MeanData["y_B"] = numpy.sqrt(3) / 2 * (3 * self.MeanData["lambda_3"] + 1)

        return self

    def anisotropyCats(self):
        """Classify turbulence anisotropy into categories on the Lumley triangle.

        Categories: ``'2-component axisymmetric'``, ``'isotropic'``,
        ``'1-component'``, or ``'non-pure'``.

        Returns
        -------
        MeanDataCalculator
            Self (for method chaining).
        """
        if "isotropy_cat" not in self.MeanData.columns:
            self.anisotropyEigs()

            self.MeanData["isotropy_cat"] = "non-pure"
            self.MeanData.loc[(self.MeanData["x_B"] <= 0.35) & (
                        self.MeanData["y_B"] <= -1 / numpy.sqrt(3) * self.MeanData["x_B"] + 7 / 30 * numpy.sqrt(
                    3)), "isotropy_cat"] = "2-component axisymmetric"
            self.MeanData.loc[(self.MeanData["y_B"] >= -1 / numpy.sqrt(3) * self.MeanData["x_B"] + 13 / 30 * numpy.sqrt(3)) & (
                        self.MeanData["y_B"] >= 1 / numpy.sqrt(3) * self.MeanData["x_B"] + 1 / 10 * numpy.sqrt(
                    3)), "isotropy_cat"] = "isotropic"
            self.MeanData.loc[(self.MeanData["x_B"] >= 0.65) & (
                        self.MeanData["y_B"] <= 1 / numpy.sqrt(3) * self.MeanData["x_B"] - 1 / 10 * numpy.sqrt(
                    3)), "isotropy_cat"] = "1-component"

        return self

    def StrucFun_eps(self, tau_range = None, title_additions = "", title_additions_eps = "",
                     rmin = 0, rmax = 10, max = False, plus_minus = 2):
        """

        :param tau_range:
        :param title_additions:
        :param title_additions_eps:
        :param rmin:
        :param rmax:
        :param max:
        :return:
        """

        a = 0.52
        col_names = {tau:"D11" + title_additions + "_" + str(tau) + "s" for tau in tau_range}
        col_names_reversed = {("D11" + title_additions + "_" + str(tau) + "s"):tau for tau in tau_range}
        data = self.MeanData[list(col_names.values()) + ["u_mag" + title_additions]]
        estimations = pandas.DataFrame(index=data.index, columns=col_names.values())

        if max:
            max_tau = data[list(col_names.values())].idxmax(axis=1).map(col_names_reversed)

        for tau in tau_range:
            data_temp = ((a * data[col_names[tau]]) ** (3 / 2)) / (tau * data["u_mag" + title_additions])
            if max:
                mask = (tau * data["u_mag" + title_additions] < max_tau * data["u_mag" + title_additions] + plus_minus) \
                       & (tau * data["u_mag" + title_additions] > max_tau * data["u_mag" + title_additions] - plus_minus)
            else:
                mask = (tau * data["u_mag" + title_additions] < rmax) & (tau * data["u_mag" + title_additions] > rmin)
            # estimations[col_names[tau]] = ((((a * self.MeanData[col_names[tau]]) ** (3 / 2)) / (tau * self.MeanData["u_mag"])).compute())\
            #     .loc[((tau * self.MeanData["u_mag" + title_additions] < rmax) & (tau * self.MeanData["u_mag" + title_additions] > rmin)).compute()]
            estimations[col_names[tau]] = data_temp.loc[mask]

        self.MeanData["eps_D11" + title_additions_eps] = estimations.mean(axis=1)

        return self

    def ThirdStrucFun_eps(self, tau_range = None, title_additions = "", title_additions_eps = "",
                          rmin = 0, rmax = 10, max = False, plus_minus = 2):
        """

        :param tau_range:
        :param title_additions:
        :param title_additions_eps:
        :param rmin:
        :param rmax:
        :return:
        """

        col_names = {tau:"D111" + title_additions + "_" + str(tau) + "s" for tau in tau_range}
        col_names_reversed = {("D111" + title_additions + "_" + str(tau) + "s"):tau for tau in tau_range}
        data = self.MeanData[list(col_names.values()) + ["u_mag" + title_additions]]
        estimations = pandas.DataFrame(index=data.index, columns=col_names.values())

        if max:
            max_tau = data[list(col_names.values())].idxmax(axis=1).map(col_names_reversed)

        for tau in tau_range:
            data_temp = 1.25 * data[col_names[tau]] / (tau * data["u_mag" + title_additions])

            if max:
                mask = (tau * data["u_mag" + title_additions] < max_tau * data["u_mag" + title_additions] + plus_minus) \
                       & (tau * data["u_mag" + title_additions] > max_tau * data["u_mag" + title_additions] - plus_minus)
            else:
                mask = (tau * data["u_mag" + title_additions] < rmax) & (tau * data["u_mag" + title_additions] > rmin)

            estimations[col_names[tau]] = data_temp.loc[mask]

        self.MeanData["eps_D111" + title_additions_eps] = estimations.mean(axis=1)

        return self

    def compute(self):
        """Return the accumulated mean-data DataFrame.

        Returns
        -------
        pandas.DataFrame
            All computed columns joined into one DataFrame.
        """
        return self.MeanData