import numpy
import pandas
import dask.diagnostics as diag
from .analysislayer import RawdataAnalysis

class MeanDataCalculator:
    def __init__(self, query_fields, functions, compute_mode):
        """

        :param query_fields:
        :param functions:
        :param compute_mode:
        """
        self._TurbCalc = RawdataAnalysis.singlePointStatistics(**query_fields)

        with diag.ProgressBar():
            for func in functions:
                getattr(self._TurbCalc, func)(**functions[func])

            self._MeanData = self._TurbCalc.compute(mode = compute_mode)

        if "start" in query_fields.keys():
            self._MeanData = self._MeanData.loc[self._MeanData.index >= query_fields["start"]]
        if "end" in query_fields.keys():
            self._MeanData = self._MeanData.loc[self._MeanData.index < query_fields["end"]]

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
        data = self._MeanData[list(col_names.values()) + ["u_mag" + title_additions]]
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
            # estimations[col_names[tau]] = ((((a * self._TemporaryData[col_names[tau]]) ** (3 / 2)) / (tau * self._TemporaryData["u_mag"])).compute())\
            #     .loc[((tau * self._TemporaryData["u_mag" + title_additions] < rmax) & (tau * self._TemporaryData["u_mag" + title_additions] > rmin)).compute()]
            estimations[col_names[tau]] = data_temp.loc[mask]

        self._MeanData["eps_D11" + title_additions_eps] = estimations.mean(axis=1)

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
        data = self._MeanData[list(col_names.values()) + ["u_mag" + title_additions]]
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

        self._MeanData["eps_D111" + title_additions_eps] = estimations.mean(axis=1)

        return self

    def compute(self):
        return self._MeanData