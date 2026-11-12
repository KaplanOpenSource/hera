"""simulations/utils/canopyWindProfile.py: ``urbanLogExponentProfile``
and its nested ``calcU``.

The module builds a mass-consistency initial guess for an urban wind
field: from a morphometric grid (frontal and plan area densities
``lambdaF``/``lambdaP`` and canopy height ``hc``) it derives the
Coceal & Belcher canopy length scales

    Lc  = hc (1 - lambdaP) / lambdaF
    l   = 2 beta^3 Lc                      (beta = 0.2)
    z0  = (l / kappa) exp(-kappa / beta)   (kappa = 0.41)
    d   =  l / kappa

then puts an exponential profile inside the canopy and a logarithmic one
above it,

    z <  hc : U(z) = Uh exp(beta (z - hc) / l)
    z >  hc : U(z) = (Uh beta / kappa) ln((z - hc + d) / z0)

with ``Uh`` fixed by matching the log branch to the station's geostrophic
speed at 800 m.  Every assertion below is derived from those published
relations, and two of them are consequences the code cannot dodge:

* the two branches are *constructed* to meet at ``z = hc``, because
  ``d / z0 = exp(kappa / beta)`` makes the logarithm collapse to
  ``kappa / beta`` there and the log branch return exactly ``Uh``;
* the log branch must reproduce the station speed at ``z = 800 m``,
  which is the condition that defined ``Uh`` in the first place.

Two defects surfaced:

* B244: ``calcU`` tests ``if z < hc`` then ``elif (z > hc) and
  (z < surfaceLayerHeight)``.  Both are strict, so ``z == hc`` -- the
  canopy top, the one height where the two profiles are built to agree --
  falls through to the ``else`` branch, which *overwrites* z with the
  2000 m surface-layer height.  A cell exactly at roof level is given the
  wind speed of the top of the boundary layer: 11.10 m/s where the
  profile is 2.44 m/s, a factor of 4.6, and identical to the value
  returned for z = 3000 m.
* B245: the interpolation write-back uses chained assignment,
  ``data['Ux'].iloc[i] = ...``.  pandas already emits
  ``FutureWarning: ChainedAssignmentError`` for it, and under
  copy-on-write -- the pandas 3.0 default -- the write silently does
  nothing: 60 of 80 cells come back NaN, i.e. the whole point of the
  interpolation step is lost without any error.

Not covered: the ``interpolation == 'linear1'`` and ``'Elevation'``
branches.  ``interpolation`` is a local literal set to ``'linear'`` three
lines above, so neither branch is reachable (and the 'Elevation' one is
an admittedly unfinished sketch that only prints).
"""
import numpy
import pandas
import pytest

from hera.simulations.utils.canopyWindProfile import karman, urbanLogExponentProfile

BETA = 0.2
GEOHEIGHT = 800.0
SURFACE_LAYER_HEIGHT = 2000.0

XS = [2.5, 7.5, 12.5, 17.5]      # domain cells, dx = 5, inside 10 m lambda cells
ZS = [1.0, 10.0, 20.0, 800.0, 3000.0]
SPEEDGEO = 10.0


def _lambdaGrid(lambdaF=0.2, lambdaP=0.3, hc=10.0):
    """A 2x2 grid of 10 m lambda cells; only iloc[0].geometry.bounds is read,
    and the merge onto the domain cells is on the (i0, j0) pair."""
    import geopandas
    from shapely.geometry import box

    return geopandas.GeoDataFrame(
        dict(
            i0=[0.0, 0.0, 1.0, 1.0],
            j0=[0.0, 1.0, 0.0, 1.0],
            lambdaF=[lambdaF] * 4,
            lambdaP=[lambdaP] * 4,
            hc=[hc] * 4,
            geometry=[
                box(0, 0, 10, 10),
                box(0, 10, 10, 20),
                box(10, 0, 20, 10),
                box(10, 10, 20, 20),
            ],
        )
    )


def _cellCenters(zs=None):
    zs = ZS if zs is None else zs
    return pandas.DataFrame(
        [dict(x=x, y=y, z=z) for x in XS for y in XS for z in zs]
    )


def _stations(speedgeo=SPEEDGEO, meteorologicalDirection=270.0):
    return pandas.DataFrame(
        dict(speedgeo=[speedgeo], MeteorologicalDirection=[meteorologicalDirection])
    )


def _run(lambdaF=0.2, lambdaP=0.3, hc=10.0, zs=None, direction=270.0):
    return urbanLogExponentProfile(
        _cellCenters(zs), _lambdaGrid(lambdaF, lambdaP, hc), _stations(
            meteorologicalDirection=direction
        )
    )


def _lengthScales(lambdaF=0.2, lambdaP=0.3, hc=10.0):
    """The documented Coceal & Belcher relations, independently evaluated."""
    Lc = hc * (1 - lambdaP) / lambdaF
    ll = 2 * BETA**3 * Lc
    zz0 = ll / karman * numpy.exp(-karman / BETA)
    dd = ll / karman
    return Lc, ll, zz0, dd


def _uh(lambdaF=0.2, lambdaP=0.3, hc=10.0, speedgeo=SPEEDGEO):
    """Uh is whatever makes the log branch equal the geostrophic speed at 800 m."""
    _, _, zz0, dd = _lengthScales(lambdaF, lambdaP, hc)
    return (speedgeo * karman) / (
        BETA * numpy.log((GEOHEIGHT - hc + dd) / zz0)
    )


def _profile(z, lambdaF=0.2, lambdaP=0.3, hc=10.0):
    """The documented two-branch canopy profile."""
    _, ll, zz0, dd = _lengthScales(lambdaF, lambdaP, hc)
    Uh = _uh(lambdaF, lambdaP, hc)
    if z < hc:
        return Uh * numpy.exp(BETA * (z - hc) / ll)
    capped = min(z, SURFACE_LAYER_HEIGHT)
    return (Uh * BETA / karman) * numpy.log((capped - hc + dd) / zz0)


@pytest.fixture(scope="module")
def profile():
    """The reference run. Module-scoped: the interpolation step is the
    expensive part and nothing here mutates the result."""
    return _run()


def _centre(frame):
    """The cells calcU actually writes a value into: the middle of each
    lambda cell, identified by i01 == i01max//2 and j01 == j01max//2."""
    return frame[
        (frame.i01 == frame.i01max // 2) & (frame.j01 == frame.j01max // 2)
    ]


@pytest.mark.unit
class TestCanopyLengthScales:
    def test_the_length_scales_follow_their_definitions(self, profile):
        """Lc = hc(1-lambdaP)/lambdaF, l = 2 beta^3 Lc, z0 = l/kappa
        exp(-kappa/beta), d = l/kappa."""
        Lc, ll, zz0, dd = _lengthScales()
        assert profile.Lc.unique() == pytest.approx([Lc])
        assert profile.ll.unique() == pytest.approx([ll])
        assert profile.zz0.unique() == pytest.approx([zz0])
        assert profile.dd.unique() == pytest.approx([dd])

    def test_the_displacement_height_exceeds_the_roughness_length(self, profile):
        """d/z0 = exp(kappa/beta) > 1 by construction."""
        assert profile.dd.iloc[0] / profile.zz0.iloc[0] == pytest.approx(
            numpy.exp(karman / BETA)
        )

    def test_a_denser_canopy_gives_a_shorter_mixing_length(self):
        """Lc falls as the frontal area density rises."""
        sparse = _run(lambdaF=0.1).ll.iloc[0]
        dense = _run(lambdaF=0.4).ll.iloc[0]
        assert dense < sparse

    def test_the_frontal_density_is_clipped_at_four_tenths(self):
        """lambdaF > 0.4 is replaced by 0.4 before Lc is formed."""
        clipped = _run(lambdaF=0.9)
        assert clipped.lambdaF.unique() == pytest.approx([0.4])
        assert clipped.Lc.unique() == pytest.approx(
            [_lengthScales(lambdaF=0.4)[0]]
        )

    def test_the_plan_density_is_clipped_at_six_tenths(self):
        clipped = _run(lambdaP=0.95)
        assert clipped.lambdaP.unique() == pytest.approx([0.6])

    def test_a_canopy_shorter_than_two_metres_is_replaced_wholesale(self):
        """hc < 2 substitutes hc = 2 and lambdaF = lambdaP = 0.25, so a bare
        cell still gets a defined roughness."""
        short = _run(hc=1.0)
        assert short.hc.unique() == pytest.approx([2.0])
        assert short.lambdaF.unique() == pytest.approx([0.25])
        assert short.lambdaP.unique() == pytest.approx([0.25])
        assert short.Lc.unique() == pytest.approx(
            [_lengthScales(lambdaF=0.25, lambdaP=0.25, hc=2.0)[0]]
        )


@pytest.mark.unit
class TestProfileFrame:
    def test_it_returns_one_row_per_domain_cell(self, profile):
        assert len(profile) == len(XS) ** 2 * len(ZS)

    def test_it_reports_the_velocity_components_and_the_magnitude(self, profile):
        assert {"Umag", "Ux", "Uy"} <= set(profile.columns)

    def test_it_preserves_the_original_cell_numbering(self, profile):
        """The merge is sorted back by 'index' so OpenFOAM cell order
        survives."""
        assert profile["index"].tolist() == sorted(profile["index"].tolist())

    def test_each_domain_cell_is_placed_in_a_lambda_cell(self, profile):
        """dx = dy = 10 for the lambda grid, so x = 2.5 and 7.5 share i0 = 0
        while 12.5 and 17.5 share i0 = 1."""
        assert sorted(profile.i0.unique().tolist()) == [0.0, 1.0]
        assert sorted(profile.j0.unique().tolist()) == [0.0, 1.0]

    def test_only_the_middle_cell_of_each_lambda_cell_carries_a_profile(
        self, profile
    ):
        """calcU returns NaN unless i01 == i01max//2 and j01 == j01max//2, so
        the profile is seeded on a sparse subset and interpolated outward."""
        centre = _centre(profile)
        assert len(centre) == 4 * len(ZS)
        assert not centre.Umag.isna().any()
        assert profile.drop(centre.index).Umag.isna().all()


@pytest.mark.unit
class TestCanopyProfileValues:
    """calcU against the documented two-branch profile."""

    def test_inside_the_canopy_it_decays_exponentially(self, profile):
        """U(z) = Uh exp(beta (z - hc) / l) for z < hc."""
        inside = _centre(profile)[lambda f: f.z == 1.0]
        assert inside.Umag.values == pytest.approx([_profile(1.0)] * len(inside))

    def test_above_the_canopy_it_follows_the_log_law(self, profile):
        """U(z) = (Uh beta / kappa) ln((z - hc + d) / z0) for z > hc."""
        above = _centre(profile)[lambda f: f.z == 20.0]
        assert above.Umag.values == pytest.approx([_profile(20.0)] * len(above))

    def test_the_profile_increases_with_height(self, profile):
        heights = _centre(profile).groupby("z").Umag.mean()
        assert heights.loc[1.0] < heights.loc[20.0] < heights.loc[800.0]

    def test_it_reproduces_the_station_speed_at_the_geostrophic_height(
        self, profile
    ):
        """Uh was defined by matching the log branch to the station's
        speedgeo at 800 m, so the profile must return it there exactly."""
        atGeoHeight = _centre(profile)[lambda f: f.z == 800.0]
        assert atGeoHeight.Umag.values == pytest.approx([SPEEDGEO] * len(atGeoHeight))

    def test_it_is_capped_at_the_surface_layer_height(self, profile):
        """Above 2000 m the else branch clamps z, so 3000 m returns the 2000 m
        value rather than extrapolating."""
        aloft = _centre(profile)[lambda f: f.z == 3000.0]
        assert aloft.Umag.values == pytest.approx(
            [_profile(SURFACE_LAYER_HEIGHT)] * len(aloft)
        )

    def test_a_faster_geostrophic_wind_scales_the_whole_profile(self):
        """Uh is linear in the station speed, and both branches are linear in
        Uh, so doubling speedgeo doubles the profile everywhere."""
        base = _centre(_run()).groupby("z").Umag.mean()
        faster = _centre(
            urbanLogExponentProfile(
                _cellCenters(), _lambdaGrid(), _stations(speedgeo=2 * SPEEDGEO)
            )
        ).groupby("z").Umag.mean()
        assert (faster / base).values == pytest.approx([2.0] * len(base))


@pytest.mark.unit
class TestContinuityAtCanopyTop:
    """B244: see the module docstring."""

    @pytest.mark.xfail(
        strict=True,
        reason="B244: calcU's two branches are `if z < hc` and `elif (z > hc) "
               "and (z < surfaceLayerHeight)`, both strict, so z == hc falls "
               "through to the else branch, which replaces z with the 2000 m "
               "surface-layer height. The canopy top -- the single height where "
               "the exponential and logarithmic branches are constructed to "
               "agree, both giving Uh -- instead returns the wind speed at the "
               "top of the boundary layer. See the consolidated findings issue.",
    )
    def test_the_profile_at_canopy_top_is_the_canopy_top_speed(self, profile):
        """d/z0 = exp(kappa/beta) makes the log branch collapse to Uh at
        z = hc, matching the exponential branch there. A cell at roof level
        must get that value."""
        atRoof = _centre(profile)[lambda f: f.z == 10.0]
        assert len(atRoof) == 4
        assert atRoof.Umag.values == pytest.approx([_uh()] * len(atRoof))

    def test_both_branches_agree_at_canopy_top_analytically(self):
        """The premise of the xfail above, checked without the code: the
        exponential branch gives Uh exp(0) = Uh and the log branch gives
        (Uh beta/kappa) ln(exp(kappa/beta)) = Uh."""
        _, ll, zz0, dd = _lengthScales()
        Uh = _uh()
        fromBelow = Uh * numpy.exp(BETA * (10.0 - 10.0) / ll)
        fromAbove = (Uh * BETA / karman) * numpy.log((10.0 - 10.0 + dd) / zz0)
        assert fromBelow == pytest.approx(Uh)
        assert fromAbove == pytest.approx(Uh)

    def test_canopy_top_currently_returns_the_surface_layer_value(self, profile):
        """Characterisation of B244: z = hc and z = 3000 m come back
        identical, and both are the 2000 m speed."""
        centre = _centre(profile)
        atRoof = centre[centre.z == 10.0].Umag.values
        aloft = centre[centre.z == 3000.0].Umag.values
        assert atRoof == pytest.approx(aloft)
        assert atRoof == pytest.approx(
            [_profile(SURFACE_LAYER_HEIGHT)] * len(atRoof)
        )

    def test_the_error_at_canopy_top_is_a_factor_of_four_and_a_half(self, profile):
        """Characterisation of B244, sized."""
        centre = _centre(profile)
        atRoof = centre[centre.z == 10.0].Umag.mean()
        assert atRoof / _uh() == pytest.approx(4.554, rel=1e-3)


@pytest.mark.unit
class TestWindComponents:
    def test_the_components_are_consistent_with_the_magnitude(self):
        """Ux = cos(theta) Umag, Uy = sin(theta) Umag, so the hypotenuse is
        the magnitude for any direction."""
        result = _centre(_run(direction=225.0))
        assert numpy.hypot(result.Ux, result.Uy).values == pytest.approx(
            result.Umag.abs().values
        )

    def test_a_westerly_wind_blows_along_the_x_axis(self):
        """Meteorological 270 deg is a wind FROM the west, i.e. towards +x, so
        the whole magnitude lands in Ux and none in Uy."""
        result = _centre(_run(direction=270.0))
        assert result.Ux.values == pytest.approx(result.Umag.values)
        assert result.Uy.values == pytest.approx([0.0] * len(result), abs=1e-12)

    def test_a_southerly_wind_blows_along_the_y_axis(self):
        """Meteorological 180 deg is a wind FROM the south, towards +y."""
        result = _centre(_run(direction=180.0))
        assert result.Uy.values == pytest.approx(result.Umag.values)
        assert result.Ux.values == pytest.approx([0.0] * len(result), abs=1e-12)

    def test_the_direction_is_uniform_over_the_domain(self, profile):
        """One station, so one direction for every cell."""
        angles = numpy.arctan2(profile.Uy, profile.Ux)
        assert angles.values == pytest.approx([angles.iloc[0]] * len(angles))


@pytest.mark.unit
class TestInterpolationFillsTheGaps:
    def test_every_cell_ends_up_with_a_velocity(self, profile):
        """The griddata step is what turns the sparse seeding into a field."""
        assert not profile.Ux.isna().any()
        assert not profile.Uy.isna().any()

    def test_the_interpolated_field_reproduces_the_analytic_profile(self, profile):
        """The seeded profile depends on z alone, so a linear interpolation in
        (x, y, z) over a uniform grid has to return it everywhere -- which is
        also the mass-consistency property the step exists to provide."""
        for z in ZS:
            layer = profile[profile.z == z]
            expected = _centre(profile)[lambda f: f.z == z].Umag.iloc[0]
            assert layer.Ux.values == pytest.approx(
                [expected] * len(layer), rel=1e-9
            )

    def test_the_field_is_horizontally_uniform(self, profile):
        """No horizontal gradient can appear from a horizontally uniform
        morphometry."""
        spread = profile.groupby("z").Ux.std()
        assert spread.values == pytest.approx([0.0] * len(spread), abs=1e-9)


@pytest.mark.unit
class TestChainedAssignmentWriteBack:
    """B245: see the module docstring."""

    @pytest.mark.xfail(
        strict=True,
        reason="B245: the interpolation results are written back with chained "
               "assignment, `data['Ux'].iloc[i] = ...`. pandas already warns "
               "(FutureWarning: ChainedAssignmentError) and under copy-on-write, "
               "the pandas 3.0 default, the assignment targets a temporary and is "
               "silently discarded: 60 of the 80 cells come back NaN and the "
               "interpolation step becomes a no-op with no error raised. "
               "See the consolidated findings issue.",
    )
    def test_the_interpolation_survives_copy_on_write(self):
        with pandas.option_context("mode.copy_on_write", True):
            result = _run()
        assert not result.Ux.isna().any()

    def test_under_copy_on_write_only_the_seeded_cells_survive(self):
        """Characterisation of B245: exactly the 20 cells calcU wrote keep a
        value; the 60 the interpolation was supposed to fill stay NaN."""
        with pandas.option_context("mode.copy_on_write", True):
            result = _run()
        assert result.Ux.isna().sum() == 60
        assert len(_centre(result)) == 20
        assert not _centre(result).Umag.isna().any()

    def test_the_write_back_is_flagged_by_pandas_today(self):
        """Characterisation of B245: it is not a latent style issue, pandas
        already emits a FutureWarning for every one of the four writes."""
        with pytest.warns(FutureWarning, match="ChainedAssignmentError"):
            _run()
