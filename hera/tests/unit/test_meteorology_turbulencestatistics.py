"""singlePointTurbulenceStatistics: the fluent turbulence-statistics chain.

Batch 6 and the integration suite already exercise ``fluctuations``,
``sigma``, ``uu``/``vv``/``ww``/``uw``/``vw``, ``Ustar`` and
``_ClassifyStability`` with real high-frequency data. This file targets the
methods layered on top of those -- ratios, higher moments, and
Monin-Obukhov stability -- plus the small ``InMemoryRawData`` container.

``StrucFunDir``/``StrucFun``/``ThirdStrucFun`` are left untested here: they
need a dask DataFrame (``.repartition(freq=...)`` is dask-only, real
pandas has no such method) with an index literally named ``"Time"``, plus
separate direction-vector DataFrames for a tau-lagged self-join. That is an
integration-shaped test, not a hermetic one.

One defect surfaced in ``InMemoryRawData``:

* B72: ``append`` calls ``super().append(...)``, but
  ``pandas.DataFrame.append`` was removed in pandas 2.0 (this environment
  runs 3.0.2). Every call raises ``AttributeError``.
* B73: ``zoL_Sonic``'s dedup guard checks whether ``[column, {'zmd': zmd}]``
  is already in ``self._AllCalculatedParams`` before adding a new column --
  but nothing anywhere in this class ever appends to
  ``_AllCalculatedParams`` (only to ``_CalculatedParams``, a same-named but
  different list). The guard can never fire, so calling ``zoL_Sonic`` twice
  with the same ``zmd`` adds a second, identical column instead of being a
  no-op.
"""
import numpy
import pandas
import pytest

from hera.measurements.meteorology.highfreqdata.analysis.turbulencestatistics import (
    InMemoryRawData,
    singlePointTurbulenceStatistics,
)


def _raw_data(n=60, freq="1s", seed=0):
    rng = numpy.random.default_rng(seed)
    times = pandas.date_range("2020-01-01", periods=n, freq=freq)
    return pandas.DataFrame({
        "u": 3.0 + rng.normal(0, 0.5, n),
        "v": 1.0 + rng.normal(0, 0.5, n),
        "w": 0.0 + rng.normal(0, 0.2, n),
        "T": 20.0 + rng.normal(0, 0.1, n),
    }, index=times)


@pytest.fixture()
def ts():
    metadata = dict(
        isMissingData=False,
        samplingWindow="30s",
        start=pandas.Timestamp("2020-01-01"),
        end=pandas.Timestamp("2020-01-01 00:02:00"),
        buildingHeight=5,
        height=10,
        averagedHeight=3,
    )
    return singlePointTurbulenceStatistics(_raw_data(), metadata)


@pytest.mark.unit
class TestGetData:
    def test_get_data_is_none_before_any_computation(self, ts):
        assert ts.getData() is None

    def test_get_data_returns_the_augmented_raw_data_after_fluctuations(self, ts):
        ts.fluctuations()
        assert ts.getData() is ts._RawData
        assert "up" in ts.getData().columns


@pytest.mark.unit
class TestDerivedRatios:
    def test_sigma_h_is_the_root_mean_square_of_sigma_u_and_v_over_sqrt_2(self, ts):
        ts.sigmaH()
        expected = numpy.hypot(ts._TemporaryData["sigmaU"], ts._TemporaryData["sigmaV"]) / numpy.sqrt(2)
        assert numpy.allclose(ts._TemporaryData["sigmaH"], expected)

    def test_sigma_h_over_ustar_divides_the_two(self, ts):
        ts.sigmaHOverUstar()
        assert numpy.allclose(
            ts._TemporaryData["sigmaHOverUstar"],
            ts._TemporaryData["sigmaH"] / ts._TemporaryData["Ustar"],
        )

    def test_sigma_w_over_ustar_divides_the_two(self, ts):
        ts.sigmaWOverUstar()
        assert numpy.allclose(
            ts._TemporaryData["sigmaWOverUstar"],
            ts._TemporaryData["sigmaW"] / ts._TemporaryData["Ustar"],
        )

    def test_wind_dir_std_is_non_negative(self, ts):
        ts.wind_dir_std()
        assert (ts._TemporaryData["wind_dir_std"].dropna() >= 0).all()

    def test_sigma_h_over_wind_speed_is_turbulence_intensity(self, ts):
        ts.sigmaHOverWindSpeed()
        assert numpy.allclose(
            ts._TemporaryData["sigmaHOverWindSpeed"],
            ts._TemporaryData["sigmaH"] / ts._TemporaryData["horizontal_speed_bar"],
        )

    def test_sigma_w_over_wind_speed_is_vertical_turbulence_intensity(self, ts):
        ts.sigmaWOverWindSpeed()
        assert numpy.allclose(
            ts._TemporaryData["sigmaWOverWindSpeed"],
            ts._TemporaryData["sigmaW"] / ts._TemporaryData["horizontal_speed_bar"],
        )

    def test_ustar_over_wind_speed_divides_the_two(self, ts):
        ts.uStarOverWindSpeed()
        assert numpy.allclose(
            ts._TemporaryData["uStarOverWindSpeed"],
            ts._TemporaryData["Ustar"] / ts._TemporaryData["horizontal_speed_bar"],
        )


@pytest.mark.unit
class TestHigherMoments:
    def test_w3_is_the_mean_cubed_vertical_fluctuation(self, ts):
        ts.fluctuations()
        ts.w3()
        expected = (ts._RawData["wp"] ** 3).resample("30s").mean()
        assert numpy.allclose(ts._TemporaryData["w3"].dropna(), expected.dropna())

    def test_w4_is_the_mean_fourth_power_vertical_fluctuation(self, ts):
        ts.fluctuations()
        ts.w4()
        expected = (ts._RawData["wp"] ** 4).resample("30s").mean()
        assert numpy.allclose(ts._TemporaryData["w4"].dropna(), expected.dropna())

    def test_w3_over_sigma_w3_normalises_the_third_moment(self, ts):
        ts.w3OverSigmaW3()
        assert numpy.allclose(
            ts._TemporaryData["w3OverSigmaW3"],
            ts._TemporaryData["w3"] / ts._TemporaryData["sigmaW"] ** 3,
        )

    def test_wtke_is_the_vertical_flux_of_turbulent_kinetic_energy(self, ts):
        ts.fluctuations()
        ts.wTKE()
        uu, vv, ww = ts._RawData["up"] ** 2, ts._RawData["vp"] ** 2, ts._RawData["wp"] ** 2
        expected = (0.5 * (uu + vv + ww) * ts._RawData["wp"]).resample("30s").mean()
        assert numpy.allclose(ts._TemporaryData["wTKE"].dropna(), expected.dropna())


@pytest.mark.unit
class TestMomentumCorrelations:
    def test_rvw_is_bounded_between_minus_one_and_one(self, ts):
        ts.Rvw()
        assert ((ts._TemporaryData["Rvw"].dropna() >= -1.0001) & (ts._TemporaryData["Rvw"].dropna() <= 1.0001)).all()

    def test_ruw_is_bounded_between_minus_one_and_one(self, ts):
        ts.Ruw()
        assert ((ts._TemporaryData["Ruw"].dropna() >= -1.0001) & (ts._TemporaryData["Ruw"].dropna() <= 1.0001)).all()

    def test_rvw_is_the_normalised_vw_covariance(self, ts):
        ts.Rvw()
        expected = ts._TemporaryData["vw"] / numpy.sqrt(ts._TemporaryData["vv"] * ts._TemporaryData["ww"])
        assert numpy.allclose(ts._TemporaryData["Rvw"].dropna(), expected.dropna())


@pytest.mark.unit
class TestMoninObukhovStability:
    def test_zol_sonic_divides_the_probe_height_by_the_mo_length(self, ts):
        ts.zoL_Sonic(zmd=2.0)
        assert numpy.allclose(
            ts._TemporaryData["zoL_Sonic1"].dropna(), (2.0 / ts._TemporaryData["L_Sonic"]).dropna()
        )

    @pytest.mark.xfail(
        strict=True,
        reason="B73: the dedup guard reads self._AllCalculatedParams, which "
               "nothing ever appends to (only _CalculatedParams is "
               "written), so it never recognises a repeat call and always "
               "adds another column. See the consolidated findings issue.",
    )
    def test_calling_zol_sonic_twice_with_the_same_zmd_does_not_duplicate_the_column(self, ts):
        ts.zoL_Sonic(zmd=2.0)
        ts.zoL_Sonic(zmd=2.0)
        assert "zoL_Sonic2" not in ts._TemporaryData.columns

    def test_calling_zol_sonic_twice_with_the_same_zmd_currently_duplicates_it(self, ts):
        """Characterisation of B73."""
        ts.zoL_Sonic(zmd=2.0)
        ts.zoL_Sonic(zmd=2.0)
        assert numpy.allclose(
            ts._TemporaryData["zoL_Sonic1"].dropna(), ts._TemporaryData["zoL_Sonic2"].dropna()
        )

    def test_calling_zol_sonic_with_a_different_zmd_adds_a_new_column(self, ts):
        ts.zoL_Sonic(zmd=2.0)
        ts.zoL_Sonic(zmd=5.0)
        assert "zoL_Sonic2" in ts._TemporaryData.columns
        assert numpy.allclose(
            ts._TemporaryData["zoL_Sonic2"].dropna(), (5.0 / ts._TemporaryData["L_Sonic"]).dropna()
        )

    def test_z_over_l_sonic_uses_the_effective_height_from_metadata(self, ts):
        ts.zOverL_Sonic()
        effectivez = 10 + 5 - 0.7 * 3
        assert numpy.allclose(
            ts._TemporaryData["zOverL_Sonic"].dropna(), (effectivez / ts._TemporaryData["L_Sonic"]).dropna()
        )

    def test_lminus1_masked_sonic_only_keeps_rows_past_the_wt_and_ustar_thresholds(self, ts):
        ts.Lminus1_masked_Sonic()
        kept = ts._TemporaryData["Lminus1_masked_Sonic"].dropna()
        for idx in kept.index:
            assert abs(ts._TemporaryData.loc[idx, "wT"]) > 0.05
            assert abs(ts._TemporaryData.loc[idx, "Ustar"]) > 0.15

    def test_stability_mo_length_sonic_classifies_every_row(self, ts):
        ts.StabilityMOLength_Sonic()
        categories = set(ts._TemporaryData["StabilityMOLength_Sonic"].dropna())
        assert categories <= {
            "very unstable", "unstable", "neutral/near neutral", "stable", "very stable", "No Stability",
        }


@pytest.mark.unit
class TestInMemoryRawData:
    def test_it_starts_with_an_empty_attrs_dict(self):
        data = InMemoryRawData(pandas.DataFrame({"a": [1, 2]}))
        assert data._Attrs == {}

    @pytest.mark.xfail(
        strict=True,
        reason="B72: append() calls super().append(...), but "
               "pandas.DataFrame.append was removed in pandas 2.0 (this "
               "environment runs 3.0.2). Every call raises AttributeError. "
               "See the consolidated findings issue.",
    )
    def test_append_should_combine_two_frames_and_merge_attrs(self):
        first = InMemoryRawData(pandas.DataFrame({"a": [1]}))
        first._Attrs = {"x": 1}
        second = InMemoryRawData(pandas.DataFrame({"a": [2]}))
        second._Attrs = {"y": 2}
        combined = first.append(second)
        assert list(combined["a"]) == [1, 2]
        assert combined._Attrs == {"x": 1, "y": 2}

    def test_append_currently_raises(self):
        """Characterisation of B72."""
        first = InMemoryRawData(pandas.DataFrame({"a": [1]}))
        second = InMemoryRawData(pandas.DataFrame({"a": [2]}))
        with pytest.raises(AttributeError, match="append"):
            first.append(second)

    def test_to_hdf_then_read_hdf_round_trips_the_data_and_attrs(self, tmp_path):
        pytest.importorskip("tables")
        data = InMemoryRawData(pandas.DataFrame({"a": [1, 2, 3]}))
        data._Attrs = {"height": 10}
        target = str(tmp_path / "out.hdf")
        data.to_hdf(target, key="raw")

        loaded = InMemoryRawData.read_hdf(target, key="raw")
        assert list(loaded["a"]) == [1, 2, 3]
        assert loaded._Attrs == {"height": 10}

    def test_read_hdf_without_a_json_sidecar_leaves_attrs_empty(self, tmp_path):
        pytest.importorskip("tables")
        plain = pandas.DataFrame({"a": [1, 2]})
        target = str(tmp_path / "plain.hdf")
        plain.to_hdf(target, key="raw")

        loaded = InMemoryRawData.read_hdf(target, key="raw")
        assert loaded._Attrs == {}
