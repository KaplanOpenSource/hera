import numpy
import pandas
from scipy.constants import g
from .abstractcalculator import AbstractCalculator

class AveragingCalculator(AbstractCalculator):
    def __init__(self, rawData, metadata):
        super(AveragingCalculator, self).__init__(rawData=rawData, metadata=metadata)

        self._TemporaryData = self._RawData.resample(self.SamplingWindow).mean().rename(
            columns={col:col+"_bar" for col in self._RawData.columns})
        for col in self._TemporaryData.columns:
            self._CalculatedParams.append([col, {}])


class MeanDataCalculator:
    def __init__(self, TurbCalc = None, compute_mode_turb = 'not_from_db_and_not_save', AverageCalc = None,
                 compute_mode_AverageCalc = None, data = None):
        """

        :param query_fields:
        :param functions:
        :param compute_mode:
        """

        if compute_mode_AverageCalc is None:
            compute_mode_AverageCalc = compute_mode_turb

        self._Karman = 0.4

        self.TurbCalc = TurbCalc
        self.AverageCalc = AverageCalc

        self.metaData = self.TurbCalc.metaData

        self.MeanData = self.TurbCalc.secondMoments().compute(mode = compute_mode_turb)
        AverageData = self.AverageCalc.compute(mode = compute_mode_AverageCalc)
        self.MeanData = self.MeanData.join(AverageData)

        self.MeanData = self.MeanData[self.metaData["start"]:self.metaData["end"]]

    def thresholds(self, threshold_list, inplace = False):
        """

        :param threshold_list: Format - [("u_bar","lt",20), ... , (field,preposition,bound)]
        :param inplace:
        :return:
        """

        mask = pandas.Series(data = True, index = self.MeanData.index)
        for cut in threshold_list:
            if cut[1] == "lt":
                mask = mask & (self.MeanData[cut[0]] < cut[2])
            elif cut[1] == "lte":
                mask = mask & (self.MeanData[cut[0]] <= cut[2])
            elif cut[1] == "gt":
                mask = mask & (self.MeanData[cut[0]] > cut[2])
            elif cut[1] == "gte":
                mask = mask & (self.MeanData[cut[0]] >= cut[2])
            elif cut[1] == "abs_lt":
                mask = mask & (abs(self.MeanData[cut[0]]) < cut[2])
            elif cut[1] == "abs_lte":
                mask = mask & (abs(self.MeanData[cut[0]]) <= cut[2])
            elif cut[1] == "abs_gt":
                mask = mask & (abs(self.MeanData[cut[0]]) > cut[2])
            elif cut[1] == "abs_gte":
                mask = mask & (abs(self.MeanData[cut[0]]) >= cut[2])
            elif cut[1] == "eq":
                mask = mask & (self.MeanData[cut[0]] == cut[2])
            elif cut[1] == "neq":
                mask = mask & (self.MeanData[cut[0]] != cut[2])

        data_cut = self.MeanData.loc[mask]

        if inplace:
            self.MeanData = data_cut
            return self
        else:
            return data_cut

    def hour(self):

        if "hour" not in self.MeanData.columns:
            self.MeanData["hour"] = self.MeanData.index.hour

        return self

    def timeWithinDay(self):

        if "timeWithinDay" not in self.MeanData.columns:
            self.MeanData["timeWithinDay"] = self.MeanData.index.hour + self.MeanData.index.minute / 60 + \
                                             self.MeanData.index.second / 3600

        return self

    def horizontalSpeed(self):

        if 'horizontal_speed_bar' not in self.MeanData.columns:
            self.MeanData['horizontal_speed_bar'] = numpy.hypot(self.MeanData['u_bar'], self.MeanData['v_bar'])

        return self

    def _UV_to_SpdDir(self,U, V):
        return (U ** 2 + V ** 2) ** 0.5, (-np.degrees(np.arctan2(V, U)) + 90) % 360

    def alignedStress(self):

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

        if 'sigmaU' not in self.MeanData.columns:
            self.MeanData['sigmaU'] = numpy.sqrt(self.MeanData["uu"])
            self.MeanData['sigmaV'] = numpy.sqrt(self.MeanData["vv"])
            self.MeanData['sigmaW'] = numpy.sqrt(self.MeanData["ww"])

        return self

    def sigmaAligned(self):

        if "sigmaU_aligned" not in self.MeanData.columns:
            self.sigma()
            self.alignedStress()
            self.MeanData['sigmaU_aligned'] = numpy.sqrt(self.MeanData["uu_aligned"])
            self.MeanData['sigmaV_aligned'] = numpy.sqrt(self.MeanData["vv_aligned"])

        return self

    def sigmaH(self):

        if 'sigmaH' not in self.MeanData.columns:
            self.sigma()
            self.MeanData['sigmaH'] = numpy.hypot(self.MeanData['sigmaU'], self.MeanData['sigmaV']) / numpy.sqrt(2)

        return self

    def Ustar(self):

        if 'Ustar' not in self.MeanData.columns:
            self.MeanData['Ustar'] = (self.MeanData['uw'] ** 2 + self.MeanData['vw'] ** 2) ** 0.25

        return self

    def sigmaHOverUstar(self):

        if 'sigmaHOverUstar' not in self.MeanData.columns:
            self.sigmaH()
            self.Ustar()
            self.MeanData['sigmaHOverUstar'] = self.MeanData['sigmaH']/self.MeanData['Ustar']

        return self

    def sigmaUOverUstar(self):

        if 'sigmaUOverUstar' not in self.MeanData.columns:
            self.sigmaAligned()
            self.Ustar()
            self.MeanData['sigmaUOverUstar'] = self.MeanData['sigmaU_aligned']/self.MeanData['Ustar']

        return self

    def sigmaVOverUstar(self):

        if 'sigmaVOverUstar' not in self.MeanData.columns:
            self.sigmaAligned()
            self.Ustar()
            self.MeanData['sigmaVOverUstar'] = self.MeanData['sigmaV_aligned']/self.MeanData['Ustar']

        return self

    def sigmaWOverUstar(self):

        if 'sigmaWOverUstar' not in self.MeanData.columns:
            self.sigma()
            self.Ustar()
            self.MeanData['sigmaWOverUstar'] = self.MeanData['sigmaW']/self.MeanData['Ustar']

        return self

    def sigmaHOverWindSpeed(self):

        if 'sigmaHOverWindSpeed' not in self.MeanData.columns:
            self.sigmaH()
            self.horizontalSpeed()
            self.MeanData['sigmaHOverWindSpeed'] = self.MeanData['sigmaH']/self.MeanData['horizontal_speed_bar']

        return self

    def sigmaWOverWindSpeed(self):

        if 'sigmaWOverWindSpeed' not in self.MeanData.columns:
            self.sigma()
            self.horizontalSpeed()
            self.MeanData['sigmaWOverWindSpeed'] = self.MeanData['sigmaW']/self.MeanData['horizontal_speed_bar']

        return self

    def absWOverSigmaW(self):

        if "absWOverSigmaW" not in self.MeanData.columns:
            self.sigma()
            self.MeanData["absWOverSigmaW"] = abs(self.MeanData["w_bar"]) / self.MeanData["sigmaW"]

        return self

    def uStarOverWindSpeed(self):

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

        if 'Rvw' not in self.MeanData.columns:
            self.MeanData['Rvw'] = self.MeanData['vw'] / numpy.sqrt(self.MeanData['vv'] * self.MeanData['ww'])

        return self

    def Ruw(self):

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

    def zOverL(self):

        if 'zOverL' not in self.MeanData.columns:
            self.MOLength()

            H = int(self.metaData['buildingHeight'])
            instrumentHeight = int(self.metaData['height'])
            averagedHeight = int(self.metaData['averagedHeight'])
            effectivez = instrumentHeight + H - 0.7 * averagedHeight
            self.MeanData['zOverL'] = effectivez / self.MeanData['L']

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
            eig_ser = pd.Series(data=np.nan, index=["lambda_1", "lambda_2", "lambda_3"])
        else:
            A = np.array([[series["uu"], series["uv"], series["uw"]], [series["uv"], series["vv"], series["vw"]],
                          [series["uw"], series["vw"], series["ww"]]]) / (
                            series["uu"] + series["vv"] + series["ww"]) - np.identity(3) / 3
            if np.isnan(A).any():
                eig_ser = pd.Series(data=np.nan, index=["lambda_1", "lambda_2", "lambda_3"])
            else:
                eig_ser = pd.Series(data=list(np.linalg.eigvalsh(A))[::-1], index=["lambda_1", "lambda_2", "lambda_3"])
        return eig_ser

    def anisotropyEigs(self):

        if "lambda_1" not in self.MeanData.columns:
            eig_data = self.MeanData.apply(eig, axis=1)
            for col in eig_data.columns:
                self.MeanData[col] = eig_data[col]

            self.MeanData["x_B"] = self.MeanData["lambda_1"] - self.MeanData["lambda_2"] + 3 / 2 * self.MeanData[
                "lambda_3"] + 1 / 2
            self.MeanData["y_B"] = np.sqrt(3) / 2 * (3 * self.MeanData["lambda_3"] + 1)

        return self

    def anisotropyCats(self):

        if "isotropy_cat" not in self.MeanData.columns:
            self.anisotropyEigs()

            self.MeanData["isotropy_cat"] = "non-pure"
            self.MeanData.loc[(self.MeanData["x_B"] <= 0.35) & (
                        self.MeanData["y_B"] <= -1 / np.sqrt(3) * self.MeanData["x_B"] + 7 / 30 * np.sqrt(
                    3)), "isotropy_cat"] = "2-component axisymmetric"
            self.MeanData.loc[(self.MeanData["y_B"] >= -1 / np.sqrt(3) * self.MeanData["x_B"] + 13 / 30 * np.sqrt(3)) & (
                        self.MeanData["y_B"] >= 1 / np.sqrt(3) * self.MeanData["x_B"] + 1 / 10 * np.sqrt(
                    3)), "isotropy_cat"] = "isotropic"
            self.MeanData.loc[(self.MeanData["x_B"] >= 0.65) & (
                        self.MeanData["y_B"] <= 1 / np.sqrt(3) * self.MeanData["x_B"] - 1 / 10 * np.sqrt(
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
        return self.MeanData