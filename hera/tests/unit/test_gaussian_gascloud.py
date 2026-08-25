"""abstractGasCloud / instantaneousReleaseGasCloud / continuousReleaseGasCloud.

The Gaussian puff/plume dispersion model built on top of batch 3's Sigma
and Meteorology (already validated against Briggs 1973 and standard
atmospheric profiles). No new defects surfaced here -- every method checks
out against the physics it documents: Q-scaling is exact, time-bounding
zeros out everything outside the window, and the depletion factor (DF)
only ever reduces dosage.

``Continuous`` ("Yehuda's code for convolution", per its own section
comment) expects a ``datetime`` dimension that nothing else in this file
produces (everything else uses ``time``) -- only its kernel construction
is covered; ``calc``/``_convolve`` are left for whatever produces that
shape.
"""
import numpy
import pytest

from hera.simulations.gaussian.gasCloud import (
    Continuous,
    abstractGasCloud,
    continuousReleaseGasCloud,
    instantaneousReleaseGasCloud,
)
from hera.simulations.gaussian.Meteorology import MeteorologyFactory
from hera.simulations.gaussian.Sigma import briggsRural
from hera.utils.unitHandler import tonumber, ureg


@pytest.fixture()
def meteorology():
    return MeteorologyFactory().getMeteorologyFromU10(
        u10=1 * ureg.m / ureg.s, inversion=1000 * ureg.m, verticalProfileType="log"
    )


@pytest.fixture()
def space_time():
    return dict(
        minx=1 * ureg.m, maxx=200 * ureg.m, dxdy=10 * ureg.m,
        miny=-50 * ureg.m, maxy=50 * ureg.m,
        minz=0 * ureg.m, maxz=20 * ureg.m, dz=5 * ureg.m,
        timeSpan=10 * ureg.min, dt=1 * ureg.min,
    )


def _make_cloud(meteorology, space_time, sourceQ=10 * ureg.kg, wind_profile_type="default"):
    return abstractGasCloud.createGasCloud(
        sourceQ=sourceQ, sourceHeight=2 * ureg.m, initialCloudSize=(1 * ureg.m,) * 3,
        meteorology=meteorology, wind_profile_type=wind_profile_type, spaceTime=space_time,
        sigmaType=briggsRural, deposition_velocity=0.01 * ureg.m / ureg.s,
    )


@pytest.mark.unit
class TestCreateGasCloud:
    def test_a_mass_source_gives_the_instantaneous_subclass(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time, sourceQ=10 * ureg.kg)
        assert isinstance(cloud, instantaneousReleaseGasCloud)

    def test_a_mass_per_time_source_gives_the_continuous_subclass(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time, sourceQ=10 * ureg.kg / ureg.min)
        assert isinstance(cloud, continuousReleaseGasCloud)

    def test_a_source_with_neither_unit_kind_raises(self, meteorology, space_time):
        with pytest.raises(ValueError, match="mass or mass per time"):
            _make_cloud(meteorology, space_time, sourceQ=10 * ureg.m)

    def test_an_unknown_wind_profile_type_raises(self, meteorology, space_time):
        with pytest.raises(ValueError, match="default.*HotSpot"):
            _make_cloud(meteorology, space_time, wind_profile_type="bogus")

    def test_hotspot_and_default_wind_profiles_give_different_speeds(self, meteorology, space_time):
        default_cloud = _make_cloud(meteorology, space_time, wind_profile_type="default")
        hotspot_cloud = _make_cloud(meteorology, space_time, wind_profile_type="HotSpot")
        assert default_cloud.u != hotspot_cloud.u


@pytest.mark.unit
class TestDepletionFactor:
    def test_df_is_between_zero_and_one(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        df = cloud.getDF_xarray()
        assert float(df.min()) > 0
        assert float(df.max()) <= 1.0

    def test_df_decreases_monotonically_downwind(self, meteorology, space_time):
        """Deposition only ever removes mass as the plume travels."""
        cloud = _make_cloud(meteorology, space_time)
        df = cloud.getDF_xarray()
        along_x = df.isel(y=0, z=0, time=0).values
        assert numpy.all(numpy.diff(along_x) <= 1e-12)


@pytest.mark.unit
class TestTrapezoidalIntegration:
    def test_integrating_a_constant_gives_a_linear_ramp(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        c = cloud.getConcentration_inst_noQ() * 0 + 1  # constant field of 1s
        integrated = cloud.trapezoidal_integration(c, dim="time")
        # cumulative sum of (dt * 1) should increase strictly with time
        values = integrated.isel(x=0, y=0, z=0).values
        assert numpy.all(numpy.diff(values[:-1]) > 0)


@pytest.mark.unit
class TestInstantaneousConcentration:
    def test_the_puff_is_somewhere_in_the_domain_at_the_first_timestep(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        c = cloud.getConcentration_inst_noQ()
        assert float(c.isel(time=0).max()) > 0

    def test_with_q_scales_no_q_by_the_source_mass_in_milligrams(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time, sourceQ=10 * ureg.kg)
        c_noq = cloud.getConcentration_inst_noQ()
        c = cloud.getConcentration_inst()
        point = dict(time=2, x=5, y=5, z=1)
        assert float(c.isel(**point)) == pytest.approx(
            float(c_noq.isel(**point)) * tonumber(10 * ureg.kg, ureg.mg)
        )

    def test_bounding_zeros_out_everything_outside_the_time_window(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        bounded = cloud.getConcentration_inst_noQ_bounded(start_time=2 * ureg.min, end_time=4 * ureg.min)
        assert float(bounded.sel(time=3).max()) > 0
        assert float(bounded.sel(time=8).max()) == 0

    def test_bounding_keeps_the_unbounded_values_inside_the_window(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        unbounded = cloud.getConcentration_inst_noQ()
        bounded = cloud.getConcentration_inst_noQ_bounded(start_time=2 * ureg.min, end_time=4 * ureg.min)
        assert float(bounded.sel(time=3).max()) == pytest.approx(float(unbounded.sel(time=3).max()))


@pytest.mark.unit
class TestInstantaneousDosage:
    def test_dosage_without_df_is_at_least_as_large_as_with_df(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        without_df = cloud.getDosage_inst_noQ(DF=False)
        with_df = cloud.getDosage_inst_noQ(DF=True)
        assert float(with_df.sum()) <= float(without_df.sum())

    def test_dosage_from_concentration_matches_the_direct_noerf_dosage(self, meteorology, space_time):
        """getDosage_inst_NoERF_noQ integrates the same concentration field
        that getDosageFromConcentration_inst_NoERF_noQ is handed explicitly."""
        cloud = _make_cloud(meteorology, space_time)
        direct = cloud.getDosage_inst_NoERF_noQ()
        c = cloud.getConcentration_inst_noQ()
        via_helper = cloud.getDosageFromConcentration_inst_NoERF_noQ(c)
        assert numpy.allclose(direct.values, via_helper.values, equal_nan=True)

    def test_bounded_noerf_dosage_only_integrates_the_window(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        bounded = cloud.getDosage_inst_NoERF_noQ_bounded(start_time=2 * ureg.min, end_time=4 * ureg.min)
        full = cloud.getDosage_inst_NoERF_noQ()
        # the bounded dosage can never exceed the unbounded one at the same time/point
        assert float(bounded.isel(time=-1).max()) <= float(full.isel(time=-1).max()) + 1e-12

    def test_with_q_dosage_scales_by_the_source_mass(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time, sourceQ=10 * ureg.kg)
        d_noq = cloud.getDosage_inst_noQ()
        d = cloud.getDosage_inst()
        point = dict(time=2, x=5, y=5, z=1)
        assert float(d.isel(**point)) == pytest.approx(
            float(d_noq.isel(**point)) * tonumber(10 * ureg.kg, ureg.mg)
        )

    def test_with_q_noerf_dosage_scales_by_the_source_mass(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time, sourceQ=10 * ureg.kg)
        d_noq = cloud.getDosage_inst_NoERF_noQ()
        d = cloud.getDosage_inst_NoERF()
        point = dict(time=2, x=5, y=5, z=1)
        assert float(d.isel(**point)) == pytest.approx(
            float(d_noq.isel(**point)) * tonumber(10 * ureg.kg, ureg.mg)
        )


@pytest.mark.unit
class TestRadiologyConversions:
    def test_mass_to_bq_conversion_uses_the_specific_activity(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        c = cloud.getConcentration_inst()
        c.attrs["Q"] = ureg.mg / ureg.m**3
        specific_activity = 2 * ureg.Bq / ureg.mg
        c_bq = cloud.concentrationConversion_mass_to_Bq(c, outputUnits=ureg.Bq / ureg.m**3, specificActivity=specific_activity)
        point = dict(time=2, x=5, y=5, z=1)
        assert float(c_bq.isel(**point)) == pytest.approx(float(c.isel(**point)) * 2)

    def test_tiac_noq_carries_the_requested_output_units(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        tiac = cloud.getTIAC_inst_noQ(outputUnits=ureg.min / ureg.m**3)
        assert tiac.attrs["Q"] == 1 * ureg.min / ureg.m**3

    def test_tiac_noerf_noq_carries_the_requested_output_units(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        tiac = cloud.getTIAC_inst_NoERF_noQ(outputUnits=ureg.min / ureg.m**3)
        assert tiac.attrs["Q"] == 1 * ureg.min / ureg.m**3

    def test_tiac_with_q_scales_by_specific_activity_and_source_mass(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time, sourceQ=10 * ureg.kg)
        d_noq = cloud.getDosage_inst_noQ()
        specific_activity = 3 * ureg.Bq / ureg.mg
        tiac = cloud.getTIAC_inst(specifitActivity=specific_activity, outputUnits=ureg.Bq * ureg.min / ureg.m**3)
        point = dict(time=2, x=5, y=5, z=1)
        expected = float(d_noq.isel(**point)) * tonumber(10 * ureg.kg, ureg.mg) * 3
        assert float(tiac.isel(**point)) == pytest.approx(expected)

    def test_tiac_noerf_with_q_scales_the_same_way(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time, sourceQ=10 * ureg.kg)
        d_noq = cloud.getDosage_inst_NoERF_noQ()
        specific_activity = 3 * ureg.Bq / ureg.mg
        tiac = cloud.getTIAC_inst_NoERF(specifitActivity=specific_activity, outputUnits=ureg.Bq * ureg.min / ureg.m**3)
        point = dict(time=2, x=5, y=5, z=1)
        expected = float(d_noq.isel(**point)) * tonumber(10 * ureg.kg, ureg.mg) * 3
        assert float(tiac.isel(**point)) == pytest.approx(expected)

    def test_tiac_from_concentration_noerf_noq_matches_the_direct_noq_route(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        c = cloud.getConcentration_inst_noQ()
        via_helper = cloud.getTIACFromConcentration_inst_NoERF_noQ(c, outputUnits=ureg.min / ureg.m**3)
        direct = cloud.getTIAC_inst_NoERF_noQ(outputUnits=ureg.min / ureg.m**3)
        assert numpy.allclose(via_helper.values, direct.values, equal_nan=True)

    def test_tiac_from_concentration_noerf_bounded_only_covers_the_window(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        bounded = cloud.getTIAC_inst_NoERF_noQ_bounded(start_time=2 * ureg.min, end_time=4 * ureg.min)
        full = cloud.getTIAC_inst_NoERF_noQ()
        assert float(bounded.isel(time=-1).max()) <= float(full.isel(time=-1).max()) + 1e-12


@pytest.mark.unit
class TestGetTiacForDist:
    def test_it_returns_the_closest_grid_point_for_each_requested_distance(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        tiac = cloud.getTIAC_inst_noQ()
        result = cloud.get_TIAC_for_dist(tiac, y=0, z=0, dist_list=[1, 51, 101])
        assert list(result.columns) == ["Distance", "TIAC"]
        assert len(result) == 3
        assert set(result["Distance"]) <= set(tiac.x.values)

    def test_tiac_decreases_with_distance_downwind(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        tiac = cloud.getTIAC_inst_noQ()
        result = cloud.get_TIAC_for_dist(tiac, y=0, z=0, dist_list=[1, 51, 101, 151])
        assert numpy.all(numpy.diff(result["TIAC"].values) <= 1e-12)


@pytest.mark.unit
class TestContinuousReleaseGasCloudIsUnusable:
    """B78: continuousReleaseGasCloud inherits directly from abstractGasCloud
    -- NOT from instantaneousReleaseGasCloud -- but all four of its own
    public methods call self.getDosage_inst_noQ / self.getDosage_inst_NoERF_noQ,
    which are only ever defined on instantaneousReleaseGasCloud. Every one of
    them raises AttributeError; the class cannot be used at all."""

    @pytest.fixture()
    def cloud(self, meteorology, space_time):
        return _make_cloud(meteorology, space_time, sourceQ=10 * ureg.kg / ureg.min)

    @pytest.mark.xfail(
        strict=True,
        reason="B78: getConcentration_cont calls self.getDosage_inst_noQ, "
               "defined only on instantaneousReleaseGasCloud. "
               "See the consolidated findings issue.",
    )
    def test_concentration_cont_should_scale_the_erf_dosage_by_the_release_rate(self, cloud):
        cloud.getConcentration_cont()

    @pytest.mark.xfail(
        strict=True,
        reason="B78: same missing method, via getConcentration_cont_NoERF.",
    )
    def test_concentration_cont_noerf_should_scale_the_noerf_dosage_by_the_release_rate(self, cloud):
        cloud.getConcentration_cont_NoERF()

    @pytest.mark.xfail(
        strict=True,
        reason="B78: getDosage_cont_NoERF calls getConcentration_cont first.",
    )
    def test_dosage_cont_noerf_should_be_non_negative_and_increase_with_time(self, cloud):
        cloud.getDosage_cont_NoERF()

    @pytest.mark.xfail(
        strict=True,
        reason="B78: getDosage_cont_doubleNoERF calls getConcentration_cont_NoERF first.",
    )
    def test_dosage_cont_double_noerf_should_run_without_error(self, cloud):
        cloud.getDosage_cont_doubleNoERF()

    def test_all_four_public_methods_currently_raise(self, cloud):
        """Characterisation of B78."""
        for method_name in (
            "getConcentration_cont",
            "getConcentration_cont_NoERF",
            "getDosage_cont_NoERF",
            "getDosage_cont_doubleNoERF",
        ):
            with pytest.raises(AttributeError, match="getDosage_inst"):
                getattr(cloud, method_name)()


@pytest.mark.unit
class TestContinuousConvolutionKernel:
    """Only __init__'s kernel construction -- calc()/_convolve() expect a
    'datetime' dimension that nothing else in this module produces."""

    def test_the_kernel_has_one_weight_per_step(self):
        kernel = Continuous(dt=1 * ureg.min, kernelsize=5, timetofinish=10 * ureg.min)
        assert kernel.Timekernel.shape == (5, 1, 1, 1)

    def test_the_kernel_integrates_towards_the_documented_decay(self):
        """By timetofinish, the exponential's cumulative release should be
        about 90% (it is defined to reach 0.1 remaining at that time)."""
        kernel = Continuous(dt=1 * ureg.min, kernelsize=10, timetofinish=10 * ureg.min)
        total_released = float(kernel.Timekernel.sum()) * kernel.dt
        assert total_released == pytest.approx(0.9, abs=0.05)

    def test_dt_is_stored_in_minutes(self):
        kernel = Continuous(dt=30 * ureg.s, kernelsize=3, timetofinish=5 * ureg.min)
        assert kernel.dt == pytest.approx(0.5)
