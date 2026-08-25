"""MeanDataCalculator: the fluent second-moment statistics chain, and
AveragingCalculator's windowed means.

Two defects surfaced while probing construction:

* B70: ``__init__`` dispatches on ``type(TurbCalcOrData) == pandas`` and
  ``type(TurbCalcOrData) == dask`` -- comparing the argument's type to the
  *module objects* ``pandas``/``dask`` themselves, which no object's type
  is ever equal to. Both branches (and the matching ones for
  ``AverageCalcOrData``) can never match, so passing a plain
  ``pandas.DataFrame`` or ``dask.DataFrame`` directly -- both explicitly
  documented and type-hinted as accepted -- always raises ``ValueError``.
  Only a ``singlePointTurbulenceStatistics``/``AveragingCalculator``
  instance (or ``None``, for ``AverageCalcOrData``) actually works.
* B71: ``anisotropyEigs`` calls ``self.MeanData.apply(self._eig, axis=1)``,
  but ``_eig`` is defined as ``def _eig(series):`` -- missing ``self`` --
  so the bound-method call always passes two arguments into a
  one-parameter function and raises ``TypeError``. ``anisotropyCats``
  calls ``anisotropyEigs`` first, so it fails the same way.

Because of B70, the tests below build ``MeanDataCalculator.MeanData``
directly via ``__new__`` rather than going through the real constructor --
the same technique batch 10 used for evaporationModels (B50).
"""
import numpy
import pandas
import pytest

from hera.measurements.meteorology.highfreqdata.analysis.meandatacalculator import (
    AveragingCalculator,
    MeanDataCalculator,
)


def _mean_data_calculator(**columns):
    mc = MeanDataCalculator.__new__(MeanDataCalculator)
    mc._Karman = 0.4
    mc.metaData = {}
    mc.MeanData = pandas.DataFrame(columns)
    return mc


def _second_moments_row():
    return dict(
        u_bar=[2.0], v_bar=[1.0], w_bar=[0.1],
        uu=[0.5], vv=[0.3], ww=[0.1],
        uv=[0.05], uw=[0.02], vw=[0.01],
    )


@pytest.mark.unit
class TestConstructionRejectsPlainDataFrames:
    @pytest.mark.xfail(
        strict=True,
        reason="B70: type(x) == pandas compares against the pandas MODULE, "
               "which no DataFrame's type ever equals, so this documented "
               "input path can never succeed. See the consolidated "
               "findings issue.",
    )
    def test_a_plain_pandas_dataframe_should_be_accepted(self):
        MeanDataCalculator(pandas.DataFrame(_second_moments_row()), start="2020", end="2021")

    def test_a_plain_pandas_dataframe_currently_raises(self):
        """Characterisation of B70."""
        with pytest.raises(ValueError, match="singlePointTurbulenceStatistics"):
            MeanDataCalculator(pandas.DataFrame(_second_moments_row()), start="2020", end="2021")

    def test_an_unrecognised_average_calc_or_data_raises(self):
        """The AverageCalcOrData branch has the same bug (type(x) == pandas
        never matches), reachable independently of TurbCalcOrData: passing
        {} matches none of pandas/dask/AveragingCalculator/None."""
        turb_calc = _real_turb_calc()
        with pytest.raises(ValueError, match="pandas/dask dataframe"):
            MeanDataCalculator(turb_calc, AverageCalcOrData={}, start="2020-01-01", end="2020-01-02")


def _real_turb_calc():
    from hera.measurements.meteorology.highfreqdata.analysis.turbulencestatistics import (
        singlePointTurbulenceStatistics,
    )

    times = pandas.date_range("2020-01-01", periods=5, freq="30s")
    raw = pandas.DataFrame(
        {"u": [2.0] * 5, "v": [1.0] * 5, "w": [0.1] * 5, "T": [20.0] * 5}, index=times
    )
    metadata = dict(isMissingData=False, samplingWindow=None, start=times[0], end=times[-1])
    return singlePointTurbulenceStatistics(raw, metadata)


@pytest.mark.unit
class TestFluentChainOnSyntheticData:
    def test_horizontal_speed_is_the_hypotenuse_of_u_bar_and_v_bar(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.horizontalSpeed()
        assert mc.MeanData["horizontal_speed_bar"].iloc[0] == pytest.approx(numpy.hypot(2.0, 1.0))

    def test_sigma_takes_the_square_root_of_the_diagonal_moments(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.sigma()
        assert mc.MeanData["sigmaU"].iloc[0] == pytest.approx(numpy.sqrt(0.5))
        assert mc.MeanData["sigmaV"].iloc[0] == pytest.approx(numpy.sqrt(0.3))
        assert mc.MeanData["sigmaW"].iloc[0] == pytest.approx(numpy.sqrt(0.1))

    def test_ustar_is_the_fourth_root_of_the_summed_squared_fluxes(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.Ustar()
        assert mc.MeanData["Ustar"].iloc[0] == pytest.approx((0.02**2 + 0.01**2) ** 0.25)

    def test_sigma_h_over_ustar_chains_sigma_h_and_ustar(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.sigmaHOverUstar()
        assert mc.MeanData["sigmaHOverUstar"].iloc[0] == pytest.approx(
            mc.MeanData["sigmaH"].iloc[0] / mc.MeanData["Ustar"].iloc[0]
        )

    def test_sigma_u_over_ustar_uses_the_aligned_sigma_u(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.sigmaUOverUstar()
        assert mc.MeanData["sigmaUOverUstar"].iloc[0] == pytest.approx(
            mc.MeanData["sigmaU_aligned"].iloc[0] / mc.MeanData["Ustar"].iloc[0]
        )

    def test_sigma_v_over_ustar_uses_the_aligned_sigma_v(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.sigmaVOverUstar()
        assert mc.MeanData["sigmaVOverUstar"].iloc[0] == pytest.approx(
            mc.MeanData["sigmaV_aligned"].iloc[0] / mc.MeanData["Ustar"].iloc[0]
        )

    def test_sigma_w_over_ustar_divides_the_two(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.sigmaWOverUstar()
        assert mc.MeanData["sigmaWOverUstar"].iloc[0] == pytest.approx(
            numpy.sqrt(0.1) / mc.MeanData["Ustar"].iloc[0]
        )

    def test_aligned_stress_rotates_into_the_mean_wind_direction(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.alignedStress()
        # trace of the stress tensor is rotation-invariant
        trace_before = 0.5 + 0.3 + 0.1
        trace_after = (
            mc.MeanData["uu_aligned"].iloc[0] + mc.MeanData["vv_aligned"].iloc[0]
        )
        # uu_aligned + vv_aligned == uu + vv exactly (2D rotation of the horizontal block)
        assert trace_after == pytest.approx(0.5 + 0.3)
        assert trace_before > 0  # sanity: fixture isn't degenerate

    def test_absw_over_sigmaw_is_the_absolute_vertical_mean_over_sigma_w(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.absWOverSigmaW()
        assert mc.MeanData["absWOverSigmaW"].iloc[0] == pytest.approx(abs(0.1) / numpy.sqrt(0.1))

    def test_sigma_h_over_wind_speed_is_turbulence_intensity(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.sigmaHOverWindSpeed()
        assert mc.MeanData["sigmaHOverWindSpeed"].iloc[0] == pytest.approx(
            mc.MeanData["sigmaH"].iloc[0] / mc.MeanData["horizontal_speed_bar"].iloc[0]
        )

    def test_ustar_over_wind_speed_divides_the_two(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.uStarOverWindSpeed()
        assert mc.MeanData["uStarOverWindSpeed"].iloc[0] == pytest.approx(
            mc.MeanData["Ustar"].iloc[0] / mc.MeanData["horizontal_speed_bar"].iloc[0]
        )

    def test_rvw_is_the_normalised_vw_covariance(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.Rvw()
        assert mc.MeanData["Rvw"].iloc[0] == pytest.approx(0.01 / numpy.sqrt(0.3 * 0.1))

    def test_ruw_is_the_normalised_uw_covariance(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.Ruw()
        assert mc.MeanData["Ruw"].iloc[0] == pytest.approx(0.02 / numpy.sqrt(0.5 * 0.1))


@pytest.mark.unit
class TestStabilityAndAnisotropy:
    def test_effectivez_combines_instrument_and_building_height(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.metaData = {"buildingHeight": "5", "height": "10", "averagedHeight": "3"}
        mc.effectivez()
        assert mc.metaData["effectivez"] == pytest.approx(10 + 5 - 0.7 * 3)

    def test_zoverl_combines_mo_length_and_effective_height(self):
        mc = _mean_data_calculator(TC_T_bar=[20.0], wT=[0.05], **_second_moments_row())
        mc.metaData = {"buildingHeight": "5", "height": "10", "averagedHeight": "3"}
        mc.zOverL()
        assert mc.MeanData["zOverL"].iloc[0] == pytest.approx(
            mc.metaData["effectivez"] / mc.MeanData["L"].iloc[0]
        )

    def test_stability_mo_length_classifies_every_row(self):
        mc = _mean_data_calculator(TC_T_bar=[20.0], wT=[0.05], **_second_moments_row())
        mc.StabilityMOLength()
        assert mc.MeanData["StabilityMOLength"].iloc[0] in {
            "very unstable", "unstable", "neutral/near neutral", "stable", "very stable", "No Stability",
        }

    @pytest.mark.xfail(
        strict=True,
        reason="B71: _eig is defined without `self` "
               "(`def _eig(series):`), but is invoked as the bound method "
               "self._eig via DataFrame.apply, which always passes both "
               "self and the row -- one argument too many for a "
               "one-parameter function. See the consolidated findings issue.",
    )
    def test_anisotropy_eigs_should_add_lambda_columns(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.anisotropyEigs()

    def test_anisotropy_eigs_currently_raises(self):
        """Characterisation of B71."""
        mc = _mean_data_calculator(**_second_moments_row())
        with pytest.raises(TypeError, match="_eig"):
            mc.anisotropyEigs()

    def test_anisotropy_cats_fails_the_same_way_since_it_calls_eigs_first(self):
        mc = _mean_data_calculator(**_second_moments_row())
        with pytest.raises(TypeError, match="_eig"):
            mc.anisotropyCats()


@pytest.mark.unit
class TestAveragingCalculator:
    def test_get_data_returns_the_resampled_bar_columns(self):
        times = pandas.date_range("2020-01-01", periods=10, freq="1s")
        raw = pandas.DataFrame({"u": numpy.arange(10.0)}, index=times)
        avg = AveragingCalculator(raw, {"samplingWindow": "5s"})
        data = avg.getData()
        assert list(data.columns) == ["u_bar"]
        assert data["u_bar"].iloc[0] == pytest.approx(numpy.arange(5).mean())
