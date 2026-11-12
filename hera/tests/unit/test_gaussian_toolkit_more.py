"""gaussianToolkit.getSpaceTime and gaussianToolkit.getGasCloud.

These are the two remaining members of hera/simulations/gaussian/toolkit.py
(the sigma registry, the meteorology pass-throughs and the presentation
wiring are covered by test_gaussian_toolkit.py).

`getSpaceTime` turns a handful of physical requirements into the grid
dictionary every gasCloud method consumes.  Its two nested helpers,
`find_x_for_timeSpan` and `get_function`, are only reachable through the
root-finding it performs, so they are exercised by asserting the defining
property the code comments state: the time span must be long enough for
the requested downwind extent to sit three alongwind sigmas behind the
cloud centre.  Every expected value below is derived independently -- from
the meteorology profile and the Briggs coefficients (both validated in
test_gaussian_meteorology.py / test_gaussian_sigma.py) plus scipy's own
root finder -- never read back out of getSpaceTime's result.

`getGasCloud` is pure dispatch: pick the sigma model by name, then let
abstractGasCloud.createGasCloud choose instantaneous vs continuous from
the units of Q.

No new defects surfaced in either function.  One hard limit is worth
recording rather than filing: `find_x_for_timeSpan` brackets its root
search at a fixed [0, 1e6] m, so a domain deeper than roughly 976 km
raises "f(a) and f(b) must have different signs" from scipy instead of a
domain error of its own.  That is asserted below as a limitation, not
pinned as a bug.
"""
import numpy
import pytest
import scipy.optimize

from hera import toolkitHome
from hera.simulations.gaussian.Meteorology import MeteorologyFactory
from hera.simulations.gaussian.Sigma import BriggsRural
from hera.simulations.gaussian.gasCloud import (
    continuousReleaseGasCloud,
    instantaneousReleaseGasCloud,
)
from hera.utils.unitHandler import ureg

SOURCE_HEIGHT = 10 * ureg.m
INITIAL_CLOUD_SIZE = (1 * ureg.m,) * 3


@pytest.fixture()
def toolkit(unit_toolkit_factory):
    return unit_toolkit_factory(toolkitHome.GAUSSIANDISPERSION)


@pytest.fixture()
def meteorology():
    return MeteorologyFactory().getMeteorologyFromU10(
        u10=3 * ureg.m / ureg.s,
        inversion=1000 * ureg.m,
        verticalProfileType="powerLaw",
    )


def _spaceTime(toolkit, meteorology, **overrides):
    arguments = dict(
        meteorology=meteorology,
        sourceHeight=SOURCE_HEIGHT,
        wind_profile_type="default",
        maxx=500 * ureg.m,
        maxz=20 * ureg.m,
        dt=1 * ureg.min,
        dz=5 * ureg.m,
        dxdy_multiplier=1,
        minimal_maxy=50 * ureg.m,
        initialCloudSize=INITIAL_CLOUD_SIZE,
    )
    arguments.update(overrides)
    return toolkit.getSpaceTime(**arguments)


def _alongwindSigma(x, stability, sigma0=INITIAL_CLOUD_SIZE):
    """sigma_x from the Briggs table, evaluated independently."""
    return BriggsRural().getSigma(x=x, stability=stability, sigma0=sigma0, units=False)[
        "sigmaX"
    ][0]


def _threeSigmaRoot(maxx_m, stability):
    """The x at which maxx sits three alongwind sigmas behind the centre."""
    return scipy.optimize.brentq(
        lambda x: x - 3 * _alongwindSigma(x, stability) - maxx_m, 1.0, 1e6
    )


# ---------------------------------------------------------------------------
# getSpaceTime
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestSpaceTimeWindProfileSelection:
    def test_an_unknown_wind_profile_type_is_refused(self, toolkit, meteorology):
        with pytest.raises(ValueError, match="'default' or 'HotSpot'"):
            _spaceTime(toolkit, meteorology, wind_profile_type="powerLaw")

    def test_the_default_profile_uses_the_meteorologys_own_wind_velocity(
        self, toolkit, meteorology
    ):
        windSpeed = meteorology.getWindVelocity(height=SOURCE_HEIGHT)
        spaceTime = _spaceTime(toolkit, meteorology, wind_profile_type="default")
        expected = (1 * ureg.min * windSpeed * 1).to(ureg.m)
        assert spaceTime["dxdy"].m_as(ureg.m) == pytest.approx(expected.m_as(ureg.m))

    def test_the_hotspot_profile_uses_the_hotspot_wind_velocity(
        self, toolkit, meteorology
    ):
        windSpeed = meteorology.getWindVelocity_hotSpot(height=SOURCE_HEIGHT)
        spaceTime = _spaceTime(toolkit, meteorology, wind_profile_type="HotSpot")
        expected = (1 * ureg.min * windSpeed * 1).to(ureg.m)
        assert spaceTime["dxdy"].m_as(ureg.m) == pytest.approx(expected.m_as(ureg.m))

    def test_the_two_profiles_disagree_at_a_height_away_from_the_reference(
        self, toolkit, meteorology
    ):
        default = _spaceTime(
            toolkit, meteorology, wind_profile_type="default", sourceHeight=50 * ureg.m
        )
        hotspot = _spaceTime(
            toolkit, meteorology, wind_profile_type="HotSpot", sourceHeight=50 * ureg.m
        )
        assert default["dxdy"].m_as(ureg.m) != pytest.approx(
            hotspot["dxdy"].m_as(ureg.m)
        )


@pytest.mark.unit
class TestSpaceTimeHorizontalGrid:
    def test_the_horizontal_step_is_one_timestep_of_travel(self, toolkit, meteorology):
        """dxdy = dt * u * multiplier, so the puff centre lands on a node."""
        windSpeed = meteorology.getWindVelocity(height=SOURCE_HEIGHT)
        spaceTime = _spaceTime(toolkit, meteorology, dt=2 * ureg.min)
        expected = (2 * ureg.min * windSpeed).m_as(ureg.m)
        assert spaceTime["dxdy"].m_as(ureg.m) == pytest.approx(expected)

    @pytest.mark.parametrize("multiplier", [1, 2, 5])
    def test_the_multiplier_scales_the_horizontal_step(
        self, toolkit, meteorology, multiplier
    ):
        windSpeed = meteorology.getWindVelocity(height=SOURCE_HEIGHT)
        spaceTime = _spaceTime(toolkit, meteorology, dxdy_multiplier=multiplier)
        expected = (1 * ureg.min * windSpeed * multiplier).m_as(ureg.m)
        assert spaceTime["dxdy"].m_as(ureg.m) == pytest.approx(expected)

    def test_the_horizontal_step_is_reported_in_metres(self, toolkit, meteorology):
        assert _spaceTime(toolkit, meteorology)["dxdy"].units == ureg.m

    def test_the_crosswind_extent_is_rounded_up_to_whole_steps(
        self, toolkit, meteorology
    ):
        spaceTime = _spaceTime(toolkit, meteorology, minimal_maxy=50 * ureg.m)
        dxdy = spaceTime["dxdy"].m_as(ureg.m)
        expected = numpy.ceil(50.0 / dxdy) * dxdy
        assert -spaceTime["miny"].m_as(ureg.m) == pytest.approx(expected)

    def test_the_crosswind_extent_is_never_smaller_than_requested(
        self, toolkit, meteorology
    ):
        for minimal in (10, 50, 200, 1000):
            spaceTime = _spaceTime(toolkit, meteorology, minimal_maxy=minimal * ureg.m)
            assert -spaceTime["miny"].m_as(ureg.m) >= minimal

    def test_the_crosswind_extent_is_a_whole_number_of_steps(
        self, toolkit, meteorology
    ):
        spaceTime = _spaceTime(toolkit, meteorology, minimal_maxy=250 * ureg.m)
        steps = -spaceTime["miny"].m_as(ureg.m) / spaceTime["dxdy"].m_as(ureg.m)
        assert steps == pytest.approx(round(steps))

    def test_the_upper_crosswind_bound_carries_one_extra_metre(
        self, toolkit, meteorology
    ):
        """numpy.arange excludes its stop, so the +1 m keeps +maxy on the grid."""
        spaceTime = _spaceTime(toolkit, meteorology)
        assert spaceTime["maxy"].m_as(ureg.m) == pytest.approx(
            -spaceTime["miny"].m_as(ureg.m) + 1
        )

    def test_the_downwind_span_starts_at_the_source(self, toolkit, meteorology):
        assert _spaceTime(toolkit, meteorology)["minx"].m_as(ureg.m) == 0


@pytest.mark.unit
class TestSpaceTimePassThroughs:
    def test_the_requested_downwind_extent_is_returned_untouched(
        self, toolkit, meteorology
    ):
        maxx = 500 * ureg.m
        assert _spaceTime(toolkit, meteorology, maxx=maxx)["maxx"] is maxx

    def test_the_requested_height_and_vertical_step_are_returned_untouched(
        self, toolkit, meteorology
    ):
        maxz = 20 * ureg.m
        dz = 5 * ureg.m
        spaceTime = _spaceTime(toolkit, meteorology, maxz=maxz, dz=dz)
        assert spaceTime["maxz"] is maxz
        assert spaceTime["dz"] is dz
        assert spaceTime["minz"].m_as(ureg.m) == 0

    def test_the_requested_timestep_is_returned_untouched(self, toolkit, meteorology):
        dt = 1 * ureg.min
        assert _spaceTime(toolkit, meteorology, dt=dt)["dt"] is dt

    def test_every_documented_key_is_present(self, toolkit, meteorology):
        assert set(_spaceTime(toolkit, meteorology)) == {
            "minx", "maxx", "miny", "maxy", "minz", "maxz",
            "dxdy", "dz", "timeSpan", "dt",
        }


@pytest.mark.unit
class TestSpaceTimeTimeSpan:
    """The time span is set by find_x_for_timeSpan/get_function, which solve
    x - 3*sigma_x(x) = maxx and convert the answer to a travel time."""

    def test_the_time_span_is_reported_in_whole_minutes(self, toolkit, meteorology):
        timeSpan = _spaceTime(toolkit, meteorology)["timeSpan"]
        assert timeSpan.units == ureg.min
        assert timeSpan.magnitude == pytest.approx(round(timeSpan.magnitude))

    def test_the_cloud_centre_reaches_three_sigmas_beyond_the_domain(
        self, toolkit, meteorology
    ):
        """The property the nested root finder exists to guarantee."""
        maxx_m = 500.0
        spaceTime = _spaceTime(toolkit, meteorology, maxx=maxx_m * ureg.m)
        windSpeed = meteorology.getWindVelocity(height=SOURCE_HEIGHT)
        travelled = (windSpeed * spaceTime["timeSpan"]).m_as(ureg.m)
        sigma = _alongwindSigma(travelled, meteorology.stability)
        assert travelled - 3 * sigma >= maxx_m

    def test_one_minute_less_would_not_have_been_enough(self, toolkit, meteorology):
        """Rounding up is by a single step, not an arbitrary safety margin."""
        maxx_m = 500.0
        spaceTime = _spaceTime(toolkit, meteorology, maxx=maxx_m * ureg.m)
        windSpeed = meteorology.getWindVelocity(height=SOURCE_HEIGHT)
        shorter = spaceTime["timeSpan"] - 1 * ureg.min
        travelled = (windSpeed * shorter).m_as(ureg.m)
        sigma = _alongwindSigma(travelled, meteorology.stability)
        assert travelled - 3 * sigma < maxx_m

    def test_it_matches_an_independently_solved_root(self, toolkit, meteorology):
        maxx_m = 500.0
        root = _threeSigmaRoot(maxx_m, meteorology.stability)
        windSpeed = meteorology.getWindVelocity(height=SOURCE_HEIGHT).m_as(
            ureg.m / ureg.s
        )
        expected = numpy.ceil(root / (60.0 * windSpeed))
        spaceTime = _spaceTime(toolkit, meteorology, maxx=maxx_m * ureg.m)
        assert spaceTime["timeSpan"].m_as(ureg.min) == pytest.approx(expected)

    @pytest.mark.parametrize("maxx_m", [200.0, 1000.0, 5000.0])
    def test_it_matches_the_independent_root_at_several_domain_depths(
        self, toolkit, meteorology, maxx_m
    ):
        root = _threeSigmaRoot(maxx_m, meteorology.stability)
        windSpeed = meteorology.getWindVelocity(height=SOURCE_HEIGHT).m_as(
            ureg.m / ureg.s
        )
        expected = numpy.ceil(root / (60.0 * windSpeed))
        spaceTime = _spaceTime(toolkit, meteorology, maxx=maxx_m * ureg.m)
        assert spaceTime["timeSpan"].m_as(ureg.min) == pytest.approx(expected)

    def test_a_deeper_domain_needs_a_longer_run(self, toolkit, meteorology):
        near = _spaceTime(toolkit, meteorology, maxx=500 * ureg.m)["timeSpan"]
        far = _spaceTime(toolkit, meteorology, maxx=5000 * ureg.m)["timeSpan"]
        assert far.m_as(ureg.min) > near.m_as(ureg.min)

    def test_a_faster_wind_needs_a_shorter_run(self, toolkit):
        slow = MeteorologyFactory().getMeteorologyFromU10(
            u10=1 * ureg.m / ureg.s, inversion=1000 * ureg.m,
            verticalProfileType="powerLaw",
        )
        fast = MeteorologyFactory().getMeteorologyFromU10(
            u10=10 * ureg.m / ureg.s, inversion=1000 * ureg.m,
            verticalProfileType="powerLaw",
        )
        assert _spaceTime(toolkit, fast, maxx=2000 * ureg.m)["timeSpan"].m_as(
            ureg.min
        ) < _spaceTime(toolkit, slow, maxx=2000 * ureg.m)["timeSpan"].m_as(ureg.min)

    def test_the_domain_depth_may_be_given_in_kilometres(self, toolkit, meteorology):
        """maxx only ever reaches the solver through .m_as(m), so units convert."""
        inMetres = _spaceTime(toolkit, meteorology, maxx=2000 * ureg.m)["timeSpan"]
        inKilometres = _spaceTime(toolkit, meteorology, maxx=2 * ureg.km)["timeSpan"]
        assert inKilometres.m_as(ureg.min) == pytest.approx(inMetres.m_as(ureg.min))

    def test_a_domain_beyond_the_hardcoded_bracket_cannot_be_solved(
        self, toolkit, meteorology
    ):
        """Limitation, not a defect: the root search is bracketed at 1e6 m."""
        with pytest.raises(ValueError, match="different signs"):
            _spaceTime(toolkit, meteorology, maxx=2000 * ureg.km)

    def test_a_larger_initial_cloud_lengthens_the_run(self, toolkit, meteorology):
        """sigma0 shifts the virtual source upwind, so sigma_x is larger at
        every x and the three-sigma condition is only met further downwind."""
        small = _spaceTime(
            toolkit, meteorology, initialCloudSize=(1 * ureg.m,) * 3,
            maxx=2000 * ureg.m,
        )["timeSpan"]
        large = _spaceTime(
            toolkit, meteorology, initialCloudSize=(50 * ureg.m,) * 3,
            maxx=2000 * ureg.m,
        )["timeSpan"]
        assert large.m_as(ureg.min) >= small.m_as(ureg.min)
        # rounding up to whole minutes can hide the difference, so check the
        # underlying condition too, solved independently.
        stability = meteorology.stability
        pointRoot = scipy.optimize.brentq(
            lambda x: x - 3 * _alongwindSigma(x, stability, (1 * ureg.m,) * 3) - 2000.0,
            1.0, 1e6,
        )
        cloudRoot = scipy.optimize.brentq(
            lambda x: x - 3 * _alongwindSigma(x, stability, (50 * ureg.m,) * 3) - 2000.0,
            1.0, 1e6,
        )
        assert cloudRoot > pointRoot


# ---------------------------------------------------------------------------
# getGasCloud
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetGasCloud:
    @pytest.fixture()
    def spaceTime(self, toolkit, meteorology):
        return _spaceTime(toolkit, meteorology)

    def _build(self, toolkit, meteorology, spaceTime, **overrides):
        arguments = dict(
            sourceQ=10 * ureg.kg,
            sourceHeight=SOURCE_HEIGHT,
            initialCloudSize=INITIAL_CLOUD_SIZE,
            meteorology=meteorology,
            wind_profile_type="default",
            spaceTime=spaceTime,
            deposition_velocity=0.01 * ureg.m / ureg.s,
        )
        arguments.update(overrides)
        return toolkit.getGasCloud(**arguments)

    def test_a_mass_source_gives_an_instantaneous_cloud(
        self, toolkit, meteorology, spaceTime
    ):
        cloud = self._build(toolkit, meteorology, spaceTime, sourceQ=10 * ureg.kg)
        assert isinstance(cloud, instantaneousReleaseGasCloud)

    def test_a_mass_rate_source_gives_a_continuous_cloud(
        self, toolkit, meteorology, spaceTime
    ):
        cloud = self._build(
            toolkit, meteorology, spaceTime, sourceQ=10 * ureg.kg / ureg.min
        )
        assert isinstance(cloud, continuousReleaseGasCloud)

    def test_a_source_with_neither_unit_kind_is_refused(
        self, toolkit, meteorology, spaceTime
    ):
        with pytest.raises(ValueError, match="mass or mass per time"):
            self._build(toolkit, meteorology, spaceTime, sourceQ=10 * ureg.m)

    def test_the_default_sigma_model_is_briggs_rural(
        self, toolkit, meteorology, spaceTime
    ):
        cloud = self._build(toolkit, meteorology, spaceTime)
        assert isinstance(cloud.sigmaType, BriggsRural)

    def test_the_sigma_model_is_an_instance_not_the_class(
        self, toolkit, meteorology, spaceTime
    ):
        cloud = self._build(toolkit, meteorology, spaceTime)
        assert cloud.sigmaType is not BriggsRural
        assert cloud.sigmaType.getSigma(100.0, "D", units=False)["sigmaY"][0] > 0

    def test_an_unknown_sigma_model_name_lists_the_valid_ones(
        self, toolkit, meteorology, spaceTime
    ):
        with pytest.raises(ValueError, match="briggsRural"):
            self._build(toolkit, meteorology, spaceTime, sigmaTypeName="briggsUrban")

    def test_the_grid_dictionary_is_handed_over_unchanged(
        self, toolkit, meteorology, spaceTime
    ):
        cloud = self._build(toolkit, meteorology, spaceTime)
        assert cloud.spaceTime is spaceTime

    def test_the_meteorology_is_handed_over_unchanged(
        self, toolkit, meteorology, spaceTime
    ):
        cloud = self._build(toolkit, meteorology, spaceTime)
        assert cloud.meteorology is meteorology

    def test_the_source_terms_are_handed_over_unchanged(
        self, toolkit, meteorology, spaceTime
    ):
        deposition = 0.02 * ureg.m / ureg.s
        cloud = self._build(
            toolkit, meteorology, spaceTime, deposition_velocity=deposition
        )
        assert cloud.deposition_velocity is deposition
        assert cloud.sourceHeight is SOURCE_HEIGHT
        assert cloud.initialCloudSize is INITIAL_CLOUD_SIZE

    def test_the_source_strength_is_normalised_to_a_pint_quantity(
        self, toolkit, meteorology, spaceTime
    ):
        cloud = self._build(toolkit, meteorology, spaceTime, sourceQ=10 * ureg.kg)
        assert cloud.sourceQ.m_as(ureg.mg) == pytest.approx(1e7)

    def test_the_cloud_speed_follows_the_requested_wind_profile(
        self, toolkit, meteorology, spaceTime
    ):
        default = self._build(
            toolkit, meteorology, spaceTime, wind_profile_type="default",
            sourceHeight=50 * ureg.m,
        )
        hotspot = self._build(
            toolkit, meteorology, spaceTime, wind_profile_type="HotSpot",
            sourceHeight=50 * ureg.m,
        )
        assert default.u.m_as(ureg.m / ureg.s) == pytest.approx(
            meteorology.getWindVelocity(height=50 * ureg.m).m_as(ureg.m / ureg.s)
        )
        assert hotspot.u.m_as(ureg.m / ureg.s) == pytest.approx(
            meteorology.getWindVelocity_hotSpot(height=50 * ureg.m).m_as(
                ureg.m / ureg.s
            )
        )

    def test_an_unknown_wind_profile_type_is_refused(
        self, toolkit, meteorology, spaceTime
    ):
        with pytest.raises(ValueError, match="'default' or 'HotSpot'"):
            self._build(toolkit, meteorology, spaceTime, wind_profile_type="log")


@pytest.mark.unit
class TestSpaceTimeAndGasCloudTogether:
    """The two functions are designed to be used back to back."""

    def test_a_grid_from_getSpaceTime_produces_a_usable_concentration_field(
        self, toolkit, meteorology
    ):
        spaceTime = _spaceTime(toolkit, meteorology)
        cloud = toolkit.getGasCloud(
            sourceQ=1 * ureg.kg,
            sourceHeight=SOURCE_HEIGHT,
            initialCloudSize=INITIAL_CLOUD_SIZE,
            meteorology=meteorology,
            wind_profile_type="default",
            spaceTime=spaceTime,
            deposition_velocity=0.01 * ureg.m / ureg.s,
        )
        concentration = cloud.getConcentration_inst()

        assert set(concentration.dims) == {"time", "x", "y", "z"}
        assert numpy.isfinite(numpy.asarray(concentration)).all()
        assert float(concentration.max()) > 0

    def test_the_grid_covers_the_requested_domain(self, toolkit, meteorology):
        spaceTime = _spaceTime(toolkit, meteorology, maxx=500 * ureg.m)
        cloud = toolkit.getGasCloud(
            sourceQ=1 * ureg.kg,
            sourceHeight=SOURCE_HEIGHT,
            initialCloudSize=INITIAL_CLOUD_SIZE,
            meteorology=meteorology,
            wind_profile_type="default",
            spaceTime=spaceTime,
            deposition_velocity=0.01 * ureg.m / ureg.s,
        )
        concentration = cloud.getConcentration_inst()

        assert float(concentration.x.min()) >= 0
        assert float(concentration.x.max()) < 500
        assert float(concentration.y.min()) == pytest.approx(
            spaceTime["miny"].m_as(ureg.m)
        )
