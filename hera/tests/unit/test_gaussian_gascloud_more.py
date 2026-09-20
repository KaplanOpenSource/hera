"""gaussian/gasCloud.py: the three members left over from
test_gaussian_gascloud.py --
``instantaneousReleaseGasCloud.getTIACFromConcentration_inst_NoERF`` and
``Continuous._convolve`` / ``Continuous.calc``.

``getTIACFromConcentration_inst_NoERF`` is the with-Q counterpart of the
already-covered ``..._noQ`` variant: it takes a finished concentration
field in [mass/volume], strips the source strength back out, integrates in
time and re-applies both the source strength and the isotope's specific
activity.  Because the source strength is divided out and multiplied back
in, the answer must not depend on it at all -- that identity, and
agreement with the direct ``getTIAC_inst_NoERF`` route, are what the tests
below check.

``Continuous`` is the convolution helper in the "Yehuda's code" section at
the bottom of the file.  Its kernel construction is already covered;
``_convolve`` is exercised here directly (window on axis 0, which is what
it assumes), and ``calc`` is pinned as broken.

Bug pinned here:
  * B271: ``Continuous.calc`` and ``Continuous._convolve`` disagree
    about the shape xarray hands a rolling reducer, so ``calc`` always
    raises ValueError.
"""
import numpy
import pytest
import xarray

from hera.simulations.gaussian.Meteorology import MeteorologyFactory
from hera.simulations.gaussian.Sigma import briggsRural
from hera.simulations.gaussian.gasCloud import Continuous, abstractGasCloud
from hera.utils.unitHandler import ureg

SPECIFIC_ACTIVITY = 1e6 * ureg.Bq / ureg.mg


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


def _make_cloud(meteorology, space_time, sourceQ=10 * ureg.kg):
    return abstractGasCloud.createGasCloud(
        sourceQ=sourceQ,
        sourceHeight=10 * ureg.m,
        initialCloudSize=(1 * ureg.m,) * 3,
        meteorology=meteorology,
        wind_profile_type="default",
        spaceTime=space_time,
        sigmaType=briggsRural,
        deposition_velocity=0.01 * ureg.m / ureg.s,
    )


# ---------------------------------------------------------------------------
# getTIACFromConcentration_inst_NoERF
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestTIACFromConcentrationWithQ:
    def test_it_matches_the_direct_route_from_the_same_cloud(
        self, meteorology, space_time
    ):
        """Feeding back the cloud's own concentration must reproduce the
        dosage it would have computed internally."""
        cloud = _make_cloud(meteorology, space_time)
        concentration = cloud.getConcentration_inst()

        viaConcentration = cloud.getTIACFromConcentration_inst_NoERF(
            concentration, specifitActivity=SPECIFIC_ACTIVITY
        )
        direct = cloud.getTIAC_inst_NoERF(specifitActivity=SPECIFIC_ACTIVITY)

        assert float(viaConcentration.max()) == pytest.approx(float(direct.max()))
        assert numpy.asarray(viaConcentration) == pytest.approx(
            numpy.asarray(direct), rel=1e-9
        )

    def test_the_result_carries_the_requested_output_units(
        self, meteorology, space_time
    ):
        cloud = _make_cloud(meteorology, space_time)
        result = cloud.getTIACFromConcentration_inst_NoERF(
            cloud.getConcentration_inst(),
            specifitActivity=SPECIFIC_ACTIVITY,
            outputUnits=ureg.Bq * ureg.min / ureg.m**3,
        )
        assert result.attrs["Q"].units == ureg.Bq * ureg.min / ureg.m**3

    def test_the_default_output_units_are_becquerel_seconds_per_cubic_metre(
        self, meteorology, space_time
    ):
        cloud = _make_cloud(meteorology, space_time)
        result = cloud.getTIACFromConcentration_inst_NoERF(
            cloud.getConcentration_inst(), specifitActivity=SPECIFIC_ACTIVITY
        )
        assert result.attrs["Q"].units == ureg.Bq * ureg.s / ureg.m**3

    def test_asking_for_minutes_instead_of_seconds_divides_by_sixty(
        self, meteorology, space_time
    ):
        cloud = _make_cloud(meteorology, space_time)
        concentration = cloud.getConcentration_inst()
        inSeconds = cloud.getTIACFromConcentration_inst_NoERF(
            concentration, specifitActivity=SPECIFIC_ACTIVITY,
            outputUnits=ureg.Bq * ureg.s / ureg.m**3,
        )
        inMinutes = cloud.getTIACFromConcentration_inst_NoERF(
            concentration, specifitActivity=SPECIFIC_ACTIVITY,
            outputUnits=ureg.Bq * ureg.min / ureg.m**3,
        )
        assert float(inSeconds.max()) == pytest.approx(60 * float(inMinutes.max()))

    def test_it_is_linear_in_the_specific_activity(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        concentration = cloud.getConcentration_inst()
        single = cloud.getTIACFromConcentration_inst_NoERF(
            concentration, specifitActivity=SPECIFIC_ACTIVITY
        )
        tenfold = cloud.getTIACFromConcentration_inst_NoERF(
            concentration, specifitActivity=10 * SPECIFIC_ACTIVITY
        )
        assert float(tenfold.max()) == pytest.approx(10 * float(single.max()))

    def test_it_does_not_depend_on_the_clouds_own_source_strength(
        self, meteorology, space_time
    ):
        """The source strength is divided out of C and multiplied back in, so
        two clouds of different Q must agree on the same concentration field."""
        light = _make_cloud(meteorology, space_time, sourceQ=1 * ureg.kg)
        heavy = _make_cloud(meteorology, space_time, sourceQ=1000 * ureg.kg)
        concentration = light.getConcentration_inst()

        fromLight = light.getTIACFromConcentration_inst_NoERF(
            concentration, specifitActivity=SPECIFIC_ACTIVITY
        )
        fromHeavy = heavy.getTIACFromConcentration_inst_NoERF(
            concentration, specifitActivity=SPECIFIC_ACTIVITY
        )
        assert numpy.asarray(fromLight) == pytest.approx(
            numpy.asarray(fromHeavy), rel=1e-9
        )

    def test_it_is_linear_in_the_concentration_field(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        concentration = cloud.getConcentration_inst()
        doubled = concentration * 2
        doubled.attrs["Q"] = concentration.attrs["Q"]

        single = cloud.getTIACFromConcentration_inst_NoERF(
            concentration, specifitActivity=SPECIFIC_ACTIVITY
        )
        twice = cloud.getTIACFromConcentration_inst_NoERF(
            doubled, specifitActivity=SPECIFIC_ACTIVITY
        )
        assert float(twice.max()) == pytest.approx(2 * float(single.max()))

    def test_depletion_never_increases_the_integrated_exposure(
        self, meteorology, space_time
    ):
        cloud = _make_cloud(meteorology, space_time)
        concentration = cloud.getConcentration_inst()
        without = cloud.getTIACFromConcentration_inst_NoERF(
            concentration, specifitActivity=SPECIFIC_ACTIVITY, DF=False
        )
        with_ = cloud.getTIACFromConcentration_inst_NoERF(
            concentration, specifitActivity=SPECIFIC_ACTIVITY, DF=True
        )
        assert float(with_.sum()) <= float(without.sum())

    def test_the_result_is_non_negative_everywhere(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        result = cloud.getTIACFromConcentration_inst_NoERF(
            cloud.getConcentration_inst(), specifitActivity=SPECIFIC_ACTIVITY
        )
        assert float(result.min()) >= 0

    def test_the_exposure_accumulates_over_time(self, meteorology, space_time):
        """A time-integrated quantity computed by cumulative sum can only grow."""
        cloud = _make_cloud(meteorology, space_time)
        result = cloud.getTIACFromConcentration_inst_NoERF(
            cloud.getConcentration_inst(), specifitActivity=SPECIFIC_ACTIVITY
        )
        totals = [float(result.isel(time=index).sum()) for index in range(result.time.size)]
        assert totals == sorted(totals)

    def test_the_grid_of_the_input_field_is_preserved(self, meteorology, space_time):
        cloud = _make_cloud(meteorology, space_time)
        concentration = cloud.getConcentration_inst()
        result = cloud.getTIACFromConcentration_inst_NoERF(
            concentration, specifitActivity=SPECIFIC_ACTIVITY
        )
        assert set(result.dims) == set(concentration.dims)
        assert result.sizes == concentration.sizes

    def test_a_concentration_field_in_the_wrong_units_is_refused(
        self, meteorology, space_time
    ):
        """The input is documented as [mass/volume]; a dosage is not."""
        cloud = _make_cloud(meteorology, space_time)
        wrong = cloud.getConcentration_inst()
        wrong.attrs["Q"] = 1 * ureg.mg * ureg.min / ureg.m**3

        with pytest.raises(Exception):
            cloud.getTIACFromConcentration_inst_NoERF(
                wrong, specifitActivity=SPECIFIC_ACTIVITY
            )

    def test_a_concentration_field_without_a_q_attribute_is_refused(
        self, meteorology, space_time
    ):
        cloud = _make_cloud(meteorology, space_time)
        bare = cloud.getConcentration_inst()
        bare.attrs.pop("Q")

        with pytest.raises(KeyError, match="Q"):
            cloud.getTIACFromConcentration_inst_NoERF(
                bare, specifitActivity=SPECIFIC_ACTIVITY
            )


# ---------------------------------------------------------------------------
# Continuous._convolve / Continuous.calc
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestConvolveReducer:
    """_convolve is written for a window on axis 0, matching the kernel it
    is handed; called that way it is a plain kernel-weighted sum times dt."""

    @pytest.fixture()
    def kernel(self):
        return Continuous(dt=1 * ureg.min, kernelsize=3, timetofinish=10 * ureg.min)

    def test_a_full_window_is_the_kernel_weighted_sum_times_dt(self, kernel):
        fullKernel = numpy.tile(kernel.Timekernel, [1, 1, 1, 1])
        data = numpy.array([2.0, 3.0, 4.0]).reshape(3, 1, 1, 1)

        result = kernel._convolve(data, axis=0, FullKernel=fullKernel)

        weights = kernel.Timekernel.ravel()
        expected = sum(w * d for w, d in zip(weights, [2.0, 3.0, 4.0])) * kernel.dt
        assert float(result.ravel()[0]) == pytest.approx(expected)

    def test_a_short_window_is_aligned_with_the_end_of_the_kernel(self, kernel):
        """A window shorter than the kernel (the leading edge of a rolling
        reduction) drops the *oldest* weights, not the newest."""
        fullKernel = numpy.tile(kernel.Timekernel, [1, 1, 1, 1])
        data = numpy.array([2.0, 3.0]).reshape(2, 1, 1, 1)

        result = kernel._convolve(data, axis=0, FullKernel=fullKernel)

        weights = kernel.Timekernel.ravel()[-2:]
        expected = sum(w * d for w, d in zip(weights, [2.0, 3.0])) * kernel.dt
        assert float(result.ravel()[0]) == pytest.approx(expected)

    def test_it_reduces_the_axis_it_is_told_to(self, kernel):
        fullKernel = numpy.tile(kernel.Timekernel, [1, 2, 2, 2])
        data = numpy.ones((3, 2, 2, 2))
        assert kernel._convolve(data, axis=0, FullKernel=fullKernel).shape == (2, 2, 2)

    def test_it_scales_linearly_with_the_data(self, kernel):
        fullKernel = numpy.tile(kernel.Timekernel, [1, 1, 1, 1])
        data = numpy.ones((3, 1, 1, 1))
        single = float(kernel._convolve(data, axis=0, FullKernel=fullKernel).ravel()[0])
        double = float(
            kernel._convolve(2 * data, axis=0, FullKernel=fullKernel).ravel()[0]
        )
        assert double == pytest.approx(2 * single)

    def test_it_multiplies_by_the_timestep_in_minutes(self):
        oneMinute = Continuous(dt=1 * ureg.min, kernelsize=1, timetofinish=10 * ureg.min)
        twoMinutes = Continuous(dt=2 * ureg.min, kernelsize=1, timetofinish=10 * ureg.min)
        assert twoMinutes.dt == pytest.approx(2 * oneMinute.dt)


@pytest.mark.unit
class TestContinuousCalcIsBroken:
    @pytest.fixture()
    def kernel(self):
        return Continuous(dt=1 * ureg.min, kernelsize=3, timetofinish=10 * ureg.min)

    @pytest.fixture()
    def data(self):
        return xarray.DataArray(
            numpy.ones((5, 2, 2, 2)),
            dims=("datetime", "x", "y", "z"),
            coords={
                "datetime": numpy.arange(5),
                "x": [0, 1], "y": [0, 1], "z": [0, 1],
            },
        )

    @pytest.mark.xfail(
        strict=True,
        reason="B271: Continuous.calc builds a kernel shaped "
               "(kernelsize, nx, ny, nz) and hands it to _convolve through "
               "xarray's rolling(...).reduce(). xarray gives a reducer an "
               "array whose *last* axis is the rolling window -- shape "
               "(nt, nx, ny, nz, kernelsize) -- but _convolve indexes the "
               "kernel on axis 0 using data.shape[0] (the number of time "
               "steps, not the window length) and slices only three of the "
               "four kernel dimensions. The multiplication is therefore "
               "always shape-incompatible and calc() can never return. "
               "See the consolidated findings issue.",
    )
    def test_it_should_convolve_the_release_kernel_over_time(self, kernel, data):
        result = kernel.calc(data)
        assert set(result.dims) == {"datetime", "x", "y", "z"}
        assert numpy.isfinite(numpy.asarray(result)).all()

    def test_it_currently_raises_a_broadcasting_error(self, kernel, data):
        """Characterisation of B271."""
        with pytest.raises(ValueError, match="could not be broadcast"):
            kernel.calc(data)

    def test_the_mismatch_is_between_the_window_axis_and_the_kernel_axis(
        self, kernel, data
    ):
        """Characterisation of B271: the two shapes named in the error are
        the windowed data (window last) and the kernel slice (window first)."""
        with pytest.raises(ValueError) as raised:
            kernel.calc(data)
        message = str(raised.value)
        assert "(5,2,2,2,3)" in message
        assert "(2,2,2,2)" in message

    def test_a_kernel_as_long_as_the_series_does_not_rescue_it(self, kernel, data):
        """Characterisation of B271: it is not an off-by-one in the slice."""
        sized = Continuous(dt=1 * ureg.min, kernelsize=5, timetofinish=10 * ureg.min)
        with pytest.raises(ValueError, match="could not be broadcast"):
            sized.calc(data)
