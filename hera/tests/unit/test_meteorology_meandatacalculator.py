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
class TestSimpleDerivedColumns:
    def test_hour_extracts_the_hour_from_the_index(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.MeanData.index = pandas.date_range("2020-01-01 03:00", periods=1, freq="30min")
        mc.hour()
        assert mc.MeanData["hour"].iloc[0] == 3

    def test_timewithinday_is_fractional_hours(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.MeanData.index = pandas.date_range("2020-01-01 03:30", periods=1, freq="30min")
        mc.timeWithinDay()
        assert mc.MeanData["timeWithinDay"].iloc[0] == pytest.approx(3.5)

    def test_tke_is_half_the_trace_of_the_stress_tensor(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.TKE()
        assert mc.MeanData["TKE"].iloc[0] == pytest.approx(0.5 * (0.5 + 0.3 + 0.1))

    def test_sigmaw_over_windspeed_divides_sigmaw_by_horizontal_speed(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.sigmaWOverWindSpeed()
        assert mc.MeanData["sigmaWOverWindSpeed"].iloc[0] == pytest.approx(
            mc.MeanData["sigmaW"].iloc[0] / mc.MeanData["horizontal_speed_bar"].iloc[0]
        )

    def test_uv_to_spddir_converts_components_to_speed_and_bearing(self):
        mc = _mean_data_calculator(**_second_moments_row())
        speed, direction = mc._UV_to_SpdDir(1.0, 1.0)
        assert speed == pytest.approx(numpy.sqrt(2))
        assert direction == pytest.approx(45.0)

    def test_compute_returns_the_accumulated_mean_data(self):
        mc = _mean_data_calculator(**_second_moments_row())
        mc.TKE()
        assert mc.compute() is mc.MeanData


@pytest.mark.unit
class TestEigDirect:
    def test_eig_called_directly_on_the_class_computes_eigenvalues(self):
        """B71 is about the bound-method call through DataFrame.apply;
        calling the underlying function directly (as it must be called,
        since it has no `self` parameter) works fine and lets it be
        exercised on its own."""
        from hera.measurements.meteorology.highfreqdata.analysis.meandatacalculator import (
            MeanDataCalculator,
        )

        row = pandas.Series({"uu": 0.5, "vv": 0.3, "ww": 0.1, "uv": 0.05, "uw": 0.02, "vw": 0.01})
        result = MeanDataCalculator._eig(row)
        assert set(result.index) == {"lambda_1", "lambda_2", "lambda_3"}
        assert not result.isna().any()

    def test_eig_returns_nan_when_the_diagonal_sums_to_zero(self):
        from hera.measurements.meteorology.highfreqdata.analysis.meandatacalculator import (
            MeanDataCalculator,
        )

        row = pandas.Series({"uu": 0.0, "vv": 0.0, "ww": 0.0, "uv": 0.0, "uw": 0.0, "vw": 0.0})
        result = MeanDataCalculator._eig(row)
        assert result.isna().all()


@pytest.mark.unit
class TestThresholdsAndFilterDatesInplace:
    """Non-inplace mode hits B70 (same as the constructor bug documented
    above): it builds a fresh MeanDataCalculator from the filtered plain
    DataFrame, and that path can never succeed. inplace=True mutates the
    existing instance instead and sidesteps the constructor entirely."""

    def test_thresholds_inplace_filters_and_returns_self(self):
        mc = _mean_data_calculator(u_bar=[2.0, 3.0], v_bar=[1.0, 0.5])
        result = mc.thresholds([("u_bar", "gt", 2.5)], inplace=True)
        assert result is mc
        assert mc.MeanData["u_bar"].tolist() == [3.0]

    def test_thresholds_not_inplace_currently_raises(self):
        """Characterisation of B70 resurfacing in thresholds()."""
        mc = _mean_data_calculator(u_bar=[2.0, 3.0], v_bar=[1.0, 0.5])
        with pytest.raises(ValueError, match="singlePointTurbulenceStatistics"):
            mc.thresholds([("u_bar", "gt", 2.5)])

    def test_filterdates_inplace_narrows_to_the_interval(self):
        mc = _mean_data_calculator(u_bar=[2.0, 3.0], v_bar=[1.0, 0.5])
        mc.MeanData.index = pandas.date_range("2020-01-01", periods=2, freq="30min")
        result = mc.filterDates(start=mc.MeanData.index[0], end=mc.MeanData.index[1], inplace=True)
        assert result is mc
        assert len(mc.MeanData) == 1

    def test_filterdates_not_inplace_currently_raises(self):
        """Characterisation of B70 resurfacing in filterDates()."""
        mc = _mean_data_calculator(u_bar=[2.0, 3.0], v_bar=[1.0, 0.5])
        mc.MeanData.index = pandas.date_range("2020-01-01", periods=2, freq="30min")
        with pytest.raises(ValueError, match="singlePointTurbulenceStatistics"):
            mc.filterDates(start=mc.MeanData.index[0], end=mc.MeanData.index[1])


@pytest.mark.unit
class TestStructureFunctionDissipation:
    def test_strucfun_eps_averages_the_masked_estimations_per_tau(self):
        mc = _mean_data_calculator(D11_1s=[0.1, 0.2, 0.3], D11_2s=[0.05, 0.1, 0.15], u_mag=[1.0, 2.0, 3.0])
        mc.StrucFun_eps(tau_range=[1, 2], rmin=0, rmax=10)
        assert "eps_D11" in mc.MeanData.columns
        assert not mc.MeanData["eps_D11"].isna().any()

    def test_thirdstrucfun_eps_averages_the_masked_estimations_per_tau(self):
        mc = _mean_data_calculator(D111_1s=[0.1, 0.2, 0.3], D111_2s=[0.05, 0.1, 0.15], u_mag=[1.0, 2.0, 3.0])
        mc.ThirdStrucFun_eps(tau_range=[1, 2], rmin=0, rmax=10)
        assert "eps_D111" in mc.MeanData.columns
        assert not mc.MeanData["eps_D111"].isna().any()


@pytest.mark.unit
class TestAveragingCalculator:
    def test_get_data_returns_the_resampled_bar_columns(self):
        times = pandas.date_range("2020-01-01", periods=10, freq="1s")
        raw = pandas.DataFrame({"u": numpy.arange(10.0)}, index=times)
        avg = AveragingCalculator(raw, {"samplingWindow": "5s"})
        data = avg.getData()
        assert list(data.columns) == ["u_bar"]
        assert data["u_bar"].iloc[0] == pytest.approx(numpy.arange(5).mean())
