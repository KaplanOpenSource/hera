"""gaussian/DropletCloud.py: the ``dropletList`` property and ``getGround``,
the two members test_gaussian_dropletcloud.py could not reach.

Every class in this file builds ``FallingNonEvaporatingDroplets`` instances
internally, and that constructor is unconditionally broken (B96, already
documented) -- which is why the existing test file can only assert that
construction raises.  Here ``FallingNonEvaporatingDroplets`` is replaced,
in the DropletCloud module's own namespace, with a recording stand-in.
That makes the discretisation logic and ``getGround`` reachable without
touching the defective constructor, and it is the right seam anyway: what
these classes own is *how the lognormal droplet spectrum is split into
sub-clouds and where each sub-cloud is placed*, not the settling physics.

Every expected value is derived independently from scipy's lognormal
distribution (the geometric bin edges, the mass in each bin, the
normalisation for the clipped variants), never read back out of the
implementation.

``FixedPointClippedDropletCloud`` stays unreachable even with the stand-in
in place: B100 (already pinned in test_gaussian_dropletcloud.py) makes its
constructor dispatch to its own ``_initDropletPosition`` override before
assigning ``self.clippedDiameter``.  Its clipped discretisation is
therefore exercised through ``CirclePositionClippedDropletsCloud``, which
sets ``clippedDiameter`` first and inherits the same override.

Bugs pinned here:
  * B272: ``CirclePositionClippedDropletsCloud`` builds the crosswind
    coordinate of every sub-cloud from ``position[0]`` instead of
    ``position[1]``, so the ring is centred on (x0, x0) rather than
    (x0, y0).
  * B273: the same constructor's ``if radius is None`` branch is dead --
    ``tounit(radius, ureg.m)`` has already run on the line above and
    raises TypeError for None.
  * B274: its ring angles come from ``numpy.linspace(0, 2*pi, n)``,
    which includes both endpoints, so the 0-degree position is duplicated:
    n positions occupy n-1 locations and one of them gets a double share.
"""
import numpy
import pandas
import pytest
from scipy.stats import lognorm

import hera.simulations.gaussian.DropletCloud as dropletCloud_mod
from hera.simulations.gaussian.DropletCloud import (
    CirclePositionClippedDropletsCloud,
    FixedPositionDropletsCloud,
    LinePositionDropletsCloud,
)
from hera.utils import ureg
from hera.utils.unitHandler import tonumber

MMD = 200 * ureg.um
GEOMETRIC_STD = 1.5
POSITION = (100 * ureg.m, 250 * ureg.m, 10 * ureg.m)
TOTAL_Q = 4 * ureg.kg


class FakeFallingDroplets:
    """Records what the cloud handed it; stands in for the B96-broken class."""

    def __init__(self, particleDiameter, Q, position, meteorologyName=None,
                 met_kwargs=None, **kwargs):
        self.particleDiameter = particleDiameter
        self.Q = Q
        self.position = position
        self.meteorologyName = meteorologyName
        self.met_kwargs = met_kwargs
        self.extraKwargs = kwargs
        self.solveToTimeCalls = []
        # Two values a caller can recognise, derived from the diameter so
        # that each stand-in is distinguishable.
        self.N = tonumber(particleDiameter, ureg.m) * 1e6
        self.AreaOnSurface = tonumber(particleDiameter, ureg.m) * 1e12 * ureg.um**2

    def solveToTime(self, T):
        self.solveToTimeCalls.append(T)
        landing = tonumber(self.particleDiameter, ureg.m) * 1e3
        return pandas.DataFrame(
            {"time": [0.0, 1.0, 2.0], "z": [10.0, 5.0, 0.0], "x": [0.0, 1.0, landing]}
        )


@pytest.fixture(autouse=True)
def _fake_droplets(monkeypatch):
    monkeypatch.setattr(
        dropletCloud_mod, "FallingNonEvaporatingDroplets", FakeFallingDroplets
    )


def _distribution(mmd=MMD, geometricstd=GEOMETRIC_STD):
    return lognorm(numpy.log(geometricstd), scale=tonumber(mmd, ureg.m))


def _expectedBinEdges(clouds, mmd=MMD, geometricstd=GEOMETRIC_STD, upper=None):
    """The geometric bin edges the implementation is specified to use."""
    rv = _distribution(mmd, geometricstd)
    lower = rv.ppf(1e-4)
    upper = rv.ppf(1 - 1e-4) if upper is None else upper
    return numpy.logspace(numpy.log(lower), numpy.log(upper), clouds, base=numpy.e)


def _expectedDiameters(clouds, **kwargs):
    edges = _expectedBinEdges(clouds, **kwargs)
    step = numpy.diff(numpy.log(edges))[0]
    return numpy.exp(numpy.log(edges[:-1]) + step / 2.0)


# ---------------------------------------------------------------------------
# dropletList
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestDropletListProperty:
    @pytest.fixture()
    def cloud(self):
        return FixedPositionDropletsCloud(
            mmd=MMD, geometricstd=GEOMETRIC_STD, position=POSITION,
            Q=TOTAL_Q, clouds=6,
        )

    def test_it_exposes_the_internal_list_itself(self, cloud):
        assert cloud.dropletList is cloud._dropletList

    def test_there_is_one_sub_cloud_per_bin(self, cloud):
        """Six bin edges make five bins, so five sub-clouds."""
        assert len(cloud.dropletList) == 5

    @pytest.mark.parametrize("clouds", [3, 6, 12])
    def test_the_count_is_one_less_than_the_requested_cloud_number(self, clouds):
        cloud = FixedPositionDropletsCloud(
            mmd=MMD, geometricstd=GEOMETRIC_STD, position=POSITION,
            Q=TOTAL_Q, clouds=clouds,
        )
        assert len(cloud.dropletList) == clouds - 1

    def test_it_starts_empty_and_is_filled_by_the_constructor(self, cloud):
        assert FixedPositionDropletsCloud._dropletList is None
        assert cloud.dropletList != []

    def test_the_diameters_are_the_geometric_bin_midpoints(self, cloud):
        expected = _expectedDiameters(6)
        actual = [
            tonumber(droplet.particleDiameter, ureg.m) for droplet in cloud.dropletList
        ]
        assert actual == pytest.approx(expected)

    def test_the_diameters_increase_across_the_spectrum(self, cloud):
        diameters = [
            tonumber(droplet.particleDiameter, ureg.m) for droplet in cloud.dropletList
        ]
        assert diameters == sorted(diameters)

    def test_the_mass_of_each_sub_cloud_is_its_share_of_the_distribution(self, cloud):
        rv = _distribution()
        expected = numpy.diff(rv.cdf(_expectedBinEdges(6))) * tonumber(TOTAL_Q, ureg.kg)
        actual = [tonumber(droplet.Q, ureg.kg) for droplet in cloud.dropletList]
        assert actual == pytest.approx(expected)

    def test_the_masses_add_up_to_the_release_less_the_two_tails(self, cloud):
        """Each tail beyond the 1e-4 quantile is discarded, so 0.02% is lost."""
        total = sum(tonumber(droplet.Q, ureg.kg) for droplet in cloud.dropletList)
        assert total == pytest.approx(tonumber(TOTAL_Q, ureg.kg) * (1 - 2e-4), rel=1e-9)

    def test_every_sub_cloud_starts_at_the_requested_position(self, cloud):
        for droplet in cloud.dropletList:
            assert droplet.position == POSITION

    def test_the_meteorology_name_is_forwarded(self, cloud):
        assert {droplet.meteorologyName for droplet in cloud.dropletList} == {
            "logNormal"
        }

    def test_the_default_roughness_length_is_ten_centimetres(self, cloud):
        for droplet in cloud.dropletList:
            assert droplet.met_kwargs["z0"] == 10 * ureg.cm

    def test_caller_supplied_met_kwargs_override_the_default(self):
        cloud = FixedPositionDropletsCloud(
            mmd=MMD, geometricstd=GEOMETRIC_STD, position=POSITION, Q=TOTAL_Q,
            clouds=3, met_kwargs=dict(z0=5 * ureg.cm, u10=3 * ureg.m / ureg.s),
        )
        for droplet in cloud.dropletList:
            assert droplet.met_kwargs["z0"] == 5 * ureg.cm
            assert droplet.met_kwargs["u10"] == 3 * ureg.m / ureg.s

    def test_extra_keyword_arguments_reach_every_droplet(self):
        cloud = FixedPositionDropletsCloud(
            mmd=MMD, geometricstd=GEOMETRIC_STD, position=POSITION, Q=TOTAL_Q,
            clouds=3, dragCoefficientName="Fan",
        )
        for droplet in cloud.dropletList:
            assert droplet.extraKwargs["dragCoefficientName"] == "Fan"

    def test_a_wider_spectrum_spans_a_wider_range_of_diameters(self):
        narrow = FixedPositionDropletsCloud(
            mmd=MMD, geometricstd=1.2, position=POSITION, Q=TOTAL_Q, clouds=6
        )
        wide = FixedPositionDropletsCloud(
            mmd=MMD, geometricstd=2.5, position=POSITION, Q=TOTAL_Q, clouds=6
        )

        def _span(cloud):
            values = [
                tonumber(d.particleDiameter, ureg.m) for d in cloud.dropletList
            ]
            return values[-1] / values[0]

        assert _span(wide) > _span(narrow)


@pytest.mark.unit
class TestDropletListForTheLineSource:
    @pytest.fixture()
    def cloud(self):
        return LinePositionDropletsCloud(
            mmd=MMD, geometricstd=GEOMETRIC_STD, position=POSITION, Q=TOTAL_Q,
            linelength=30 * ureg.m, clouds=4, linepositions=3,
        )

    def test_the_line_holds_one_spectrum_per_position(self, cloud):
        assert len(cloud.dropletList) == 3 * (4 - 1)

    def test_the_positions_are_evenly_spread_along_the_line(self, cloud):
        crosswind = sorted(
            {tonumber(droplet.position[1], ureg.m) for droplet in cloud.dropletList}
        )
        assert crosswind == pytest.approx([0.0, 15.0, 30.0])

    def test_the_downwind_and_vertical_coordinates_are_untouched(self, cloud):
        for droplet in cloud.dropletList:
            assert droplet.position[0] == POSITION[0]
            assert droplet.position[2] == POSITION[2]

    def test_the_release_is_split_evenly_between_the_positions(self, cloud):
        total = sum(tonumber(droplet.Q, ureg.kg) for droplet in cloud.dropletList)
        assert total == pytest.approx(tonumber(TOTAL_Q, ureg.kg) * (1 - 2e-4), rel=1e-9)

    def test_each_position_carries_an_equal_share(self, cloud):
        byPosition = {}
        for droplet in cloud.dropletList:
            key = tonumber(droplet.position[1], ureg.m)
            byPosition[key] = byPosition.get(key, 0.0) + tonumber(droplet.Q, ureg.kg)
        shares = list(byPosition.values())
        assert shares == pytest.approx([shares[0]] * len(shares))


@pytest.mark.unit
class TestDropletListForTheClippedSpectrum:
    """The clipped discretisation is only reachable through the circle
    subclass: FixedPointClippedDropletCloud cannot be constructed at all
    (B100, already pinned in test_gaussian_dropletcloud.py -- it dispatches
    to its own _initDropletPosition override before assigning
    self.clippedDiameter), while CirclePositionClippedDropletsCloud sets
    clippedDiameter first and inherits the same override."""

    CLIP = 300 * ureg.um
    CLOUDS = 4
    CIRCLE_POSITIONS = 4

    @pytest.fixture()
    def cloud(self):
        return CirclePositionClippedDropletsCloud(
            mmd=MMD, geometricstd=GEOMETRIC_STD, position=POSITION, Q=TOTAL_Q,
            clippedDiameter=self.CLIP, clouds=self.CLOUDS,
            radius=10 * ureg.m, circlepositions=self.CIRCLE_POSITIONS,
        )

    def test_no_droplet_exceeds_the_clip(self, cloud):
        for droplet in cloud.dropletList:
            assert tonumber(droplet.particleDiameter, ureg.m) <= 300e-6

    def test_the_diameters_are_the_midpoints_of_the_clipped_bins(self, cloud):
        expected = list(_expectedDiameters(self.CLOUDS, upper=300e-6))
        actual = [
            tonumber(droplet.particleDiameter, ureg.m) for droplet in cloud.dropletList
        ]
        # the same clipped spectrum is repeated once per ring position
        assert actual == pytest.approx(expected * self.CIRCLE_POSITIONS)

    def test_the_clipped_mass_is_renormalised_to_the_whole_release(self, cloud):
        """Dividing by cdf(clip) puts each position's whole share back into
        the retained bins, instead of losing the clipped tail."""
        rv = _distribution()
        edges = _expectedBinEdges(self.CLOUDS, upper=300e-6)
        share = tonumber(TOTAL_Q, ureg.kg) / self.CIRCLE_POSITIONS
        expected = list(numpy.diff(rv.cdf(edges)) * share / rv.cdf(edges[-1]))
        actual = [tonumber(droplet.Q, ureg.kg) for droplet in cloud.dropletList]
        assert actual == pytest.approx(expected * self.CIRCLE_POSITIONS)

    def test_the_retained_mass_is_almost_the_whole_share(self, cloud):
        """Only the sub-1e-4-quantile tail is dropped after renormalisation."""
        rv = _distribution()
        edges = _expectedBinEdges(self.CLOUDS, upper=300e-6)
        expected = tonumber(TOTAL_Q, ureg.kg) * (1 - rv.cdf(edges[0]) / rv.cdf(edges[-1]))
        total = sum(tonumber(droplet.Q, ureg.kg) for droplet in cloud.dropletList)
        assert total == pytest.approx(expected, rel=1e-9)

    def test_a_tighter_clip_keeps_a_narrower_spectrum(self):
        def _largest(clip):
            cloud = CirclePositionClippedDropletsCloud(
                mmd=MMD, geometricstd=GEOMETRIC_STD, position=POSITION, Q=TOTAL_Q,
                clippedDiameter=clip, clouds=4, radius=10 * ureg.m,
                circlepositions=2,
            )
            return max(
                tonumber(d.particleDiameter, ureg.m) for d in cloud.dropletList
            )

        assert _largest(150 * ureg.um) < _largest(400 * ureg.um)

    def test_the_clipped_variant_forwards_no_met_kwargs(self, cloud):
        """Its own _initDropletPosition drops the z0 default its parent sets."""
        for droplet in cloud.dropletList:
            assert droplet.met_kwargs is None

    def test_the_meteorology_name_is_forwarded(self, cloud):
        assert {droplet.meteorologyName for droplet in cloud.dropletList} == {
            "StandardMeteorolgyConstant"
        }


@pytest.mark.unit
class TestDropletListForTheCircleSource:
    def _build(self, **overrides):
        arguments = dict(
            mmd=MMD, geometricstd=GEOMETRIC_STD, position=POSITION, Q=TOTAL_Q,
            clippedDiameter=300 * ureg.um, clouds=4, radius=10 * ureg.m,
            circlepositions=4,
        )
        arguments.update(overrides)
        return CirclePositionClippedDropletsCloud(**arguments)

    @staticmethod
    def _locations(cloud):
        return {
            (
                round(tonumber(droplet.position[0], ureg.m), 6),
                round(tonumber(droplet.position[1], ureg.m), 6),
            )
            for droplet in cloud.dropletList
        }

    def test_the_ring_holds_one_spectrum_per_position(self):
        assert len(self._build().dropletList) == 4 * (4 - 1)

    def test_the_vertical_coordinate_is_untouched(self):
        for droplet in self._build().dropletList:
            assert droplet.position[2] == POSITION[2]

    @pytest.mark.xfail(
        strict=True,
        reason="B272: CirclePositionClippedDropletsCloud builds each "
               "sub-cloud position as (position[0] + r*cos, position[0] + "
               "r*sin, position[2]) -- the crosswind coordinate is taken "
               "from position[0] instead of position[1]. The ring is "
               "therefore centred on (x0, x0), so any source whose x and y "
               "differ is placed somewhere else entirely (here 150 m away "
               "from the requested centre). Its sibling "
               "LinePositionDropletsCloud gets this right. See the "
               "consolidated findings issue.",
    )
    def test_every_position_should_sit_on_a_circle_around_the_source(self):
        x0 = tonumber(POSITION[0], ureg.m)
        y0 = tonumber(POSITION[1], ureg.m)
        for x, y in self._locations(self._build()):
            assert numpy.hypot(x - x0, y - y0) == pytest.approx(10.0)

    def test_the_ring_is_currently_centred_on_the_downwind_coordinate_twice(self):
        """Characterisation of B272."""
        x0 = tonumber(POSITION[0], ureg.m)
        for x, y in self._locations(self._build()):
            assert numpy.hypot(x - x0, y - x0) == pytest.approx(10.0)

    def test_the_ring_currently_sits_far_from_the_requested_centre(self):
        """Characterisation of B272: 150 m off for this source."""
        x0 = tonumber(POSITION[0], ureg.m)
        y0 = tonumber(POSITION[1], ureg.m)
        offsets = [
            numpy.hypot(x - x0, y - y0) for x, y in self._locations(self._build())
        ]
        assert min(offsets) > abs(y0 - x0) - 10.0

    @pytest.mark.xfail(
        strict=True,
        reason="B274: the ring angles are numpy.linspace(0, 2*pi, "
               "circlepositions), which includes both endpoints -- 0 and "
               "2*pi are the same point. A ring of n positions therefore "
               "occupies only n-1 distinct locations, and the one at 0 "
               "receives a double share of the release while the angular "
               "spacing is 2*pi/(n-1) rather than 2*pi/n. endpoint=False is "
               "what a closed ring needs. See the consolidated findings "
               "issue.",
    )
    def test_a_ring_of_four_should_occupy_four_distinct_locations(self):
        assert len(self._locations(self._build(circlepositions=4))) == 4

    @pytest.mark.parametrize("circlepositions", [3, 4, 5, 8])
    def test_a_ring_currently_occupies_one_location_fewer(self, circlepositions):
        """Characterisation of B274."""
        cloud = self._build(circlepositions=circlepositions)
        assert len(self._locations(cloud)) == circlepositions - 1

    def test_the_duplicated_location_currently_carries_a_double_share(self):
        """Characterisation of B274."""
        cloud = self._build(circlepositions=4)
        byLocation = {}
        for droplet in cloud.dropletList:
            key = (
                round(tonumber(droplet.position[0], ureg.m), 6),
                round(tonumber(droplet.position[1], ureg.m), 6),
            )
            byLocation[key] = byLocation.get(key, 0.0) + tonumber(droplet.Q, ureg.kg)
        shares = sorted(byLocation.values())
        assert len(shares) == 3
        assert shares[-1] == pytest.approx(2 * shares[0])

    @pytest.mark.xfail(
        strict=True,
        reason="B273: the `if radius is None: curpos = position` branch is "
               "unreachable -- the line above it already evaluates "
               "tounit(radius, ureg.m), which raises TypeError('Invalid "
               "magnitude for Quantity: None') for None. The documented "
               "'no radius' case cannot be requested. See the consolidated "
               "findings issue.",
    )
    def test_omitting_the_radius_should_place_every_cloud_at_the_source(self):
        cloud = self._build(radius=None)
        for droplet in cloud.dropletList:
            assert droplet.position == POSITION

    def test_omitting_the_radius_currently_raises_a_typeerror(self):
        """Characterisation of B273."""
        with pytest.raises(TypeError, match="Invalid magnitude"):
            self._build(radius=None)


# ---------------------------------------------------------------------------
# getGround
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetGround:
    @pytest.fixture()
    def cloud(self):
        return FixedPositionDropletsCloud(
            mmd=MMD, geometricstd=GEOMETRIC_STD, position=POSITION,
            Q=TOTAL_Q, clouds=4,
        )

    def test_it_returns_one_row_per_sub_cloud(self, cloud):
        ground = cloud.getGround(60 * ureg.s)
        assert len(ground) == len(cloud.dropletList)

    def test_the_result_is_a_dataframe_with_a_fresh_index(self, cloud):
        ground = cloud.getGround(60 * ureg.s)
        assert isinstance(ground, pandas.DataFrame)
        assert list(ground.index) == list(range(len(cloud.dropletList)))

    def test_the_trajectory_columns_survive(self, cloud):
        ground = cloud.getGround(60 * ureg.s)
        assert {"time", "z", "x"} <= set(ground.columns)

    def test_the_particle_count_and_footprint_are_appended(self, cloud):
        ground = cloud.getGround(60 * ureg.s)
        assert "N" in ground.columns
        assert "DropletArea" in ground.columns

    def test_each_row_carries_its_own_droplets_particle_count(self, cloud):
        ground = cloud.getGround(60 * ureg.s)
        expected = [droplet.N for droplet in cloud.dropletList]
        assert list(ground["N"]) == pytest.approx(expected)

    def test_the_footprint_is_reported_in_square_micrometres(self, cloud):
        ground = cloud.getGround(60 * ureg.s)
        expected = [
            droplet.AreaOnSurface.m_as(ureg.um**2) for droplet in cloud.dropletList
        ]
        assert list(ground["DropletArea"]) == pytest.approx(expected)

    def test_only_the_last_moment_of_each_trajectory_is_kept(self, cloud):
        """getGround is the ground deposition, i.e. the end of the fall."""
        ground = cloud.getGround(60 * ureg.s)
        assert list(ground["time"]) == pytest.approx([2.0] * len(cloud.dropletList))
        assert list(ground["z"]) == pytest.approx([0.0] * len(cloud.dropletList))

    def test_the_landing_distance_is_taken_from_the_droplets_own_solution(self, cloud):
        ground = cloud.getGround(60 * ureg.s)
        expected = [
            droplet.solveToTime(None)["x"].iloc[-1] for droplet in cloud.dropletList
        ]
        assert list(ground["x"]) == pytest.approx(expected)

    def test_every_droplet_is_solved_exactly_once_for_the_requested_time(self, cloud):
        target = 90 * ureg.s
        cloud.getGround(target)
        for droplet in cloud.dropletList:
            assert droplet.solveToTimeCalls == [target]

    def test_progress_is_reported_while_solving(self, cloud, capsys):
        cloud.getGround(60 * ureg.s)
        out = capsys.readouterr().out
        assert "solving droplets" in out
        assert f"/{len(cloud.dropletList)}" in out

    def test_a_line_source_reports_every_position(self):
        cloud = LinePositionDropletsCloud(
            mmd=MMD, geometricstd=GEOMETRIC_STD, position=POSITION, Q=TOTAL_Q,
            linelength=30 * ureg.m, clouds=3, linepositions=2,
        )
        assert len(cloud.getGround(60 * ureg.s)) == 2 * (3 - 1)
