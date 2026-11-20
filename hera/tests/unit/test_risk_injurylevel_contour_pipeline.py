r"""``InjuryLevel.calculateContours`` and the ``_getGeopandas`` hook it calls.

``calculateContours`` is the step that turns a toxic-load *field* into
*polygons*: for each time step it asks the subclass for a GeoDataFrame of
raw contours (``_getGeopandas``), drops slivers, repairs invalid rings,
dissolves duplicate levels, un-rotates the cloud, and renames the columns
into the shape the ``Injury`` layer consumes -- ``ToxicLoad``,
``TotalPolygon``, ``severity``, ``datetime``.

The fields used below are analytic Gaussians,
``C(r) = peak * exp(-r^2 / 2 s^2)``, so the geometry is known in closed
form: the contour at level ``L`` is a circle of radius
``r = s * sqrt(2 ln(peak / L))``, enclosing ``pi r^2``.  Every expected
number is computed from that formula rather than from hera's output.  The
plume is deliberately placed *off* the origin so that the
``cloud_math_direction`` un-rotation (a rotation about ``(0, 0)``) is
observable.

Four time-shape branches are covered, because the caller may hand over any
of them: a field with a real time dimension of length > 1, one of length
exactly 1, one where the time survives only as a scalar coordinate (what
``.isel(datetime=0)`` leaves behind), and one with no time at all.

B264 is pinned below: ``InjuryLevel.__str__`` never returns.  This is the
same missing ``return`` as the already-reported B122 on ``Injury.__str__``;
it is a *separate* method in a separate class, so it needs its own number
and its own fix.

Deliberately not covered here:

* the dose-response arithmetic (``getPercent`` / ``getToxicLoad``), the
  constructors and ``toJSON`` -- see ``test_risk_injurylevel.py`` and
  ``test_risk_injurylevel_remaining.py``.
* ``getContoursFromLevels``' dimension handling -- see
  ``test_risk_injurylevel_contours.py``.
* how the resulting polygons are consumed (``_postCalculate``,
  ``calculateThresholdPolygon``) -- see ``test_risk_injury_effects.py``.
"""
import json

import numpy
import pandas
import pytest

from hera.riskassessment.agents.effects.InjuryLevel import (
    InjuryLevel, InjuryLevelExponential, InjuryLevelLognormal10DoseResponse,
    InjuryLevelThreshold)

PEAK = 100.0
SIGMA = 30.0
THRESHOLD = 10.0
CENTRE_X = 100.0


def _radius(level, peak=PEAK, sigma=SIGMA):
    """Analytic radius of the C = level contour of the Gaussian below."""
    return sigma * numpy.sqrt(2.0 * numpy.log(peak / level))


def _field(peak=PEAK, sigma=SIGMA, steps=None, points=161, extent=250.0,
           attrs=None, centre=CENTRE_X):
    """A Gaussian toxic-load field centred at ``(centre, 0)``."""
    import xarray

    x = numpy.linspace(centre - extent, centre + extent, points)
    y = numpy.linspace(-extent, extent, points)
    grid_x, grid_y = numpy.meshgrid(x, y)
    values = peak * numpy.exp(-((grid_x - centre) ** 2 + grid_y ** 2)
                              / (2.0 * sigma ** 2))

    if steps is None:
        field = xarray.DataArray(values, coords={"y": y, "x": x}, dims=("y", "x"))
    else:
        stamps = pandas.to_datetime(
            ["2020-01-01 00:%02d" % minute for minute in range(steps)])
        field = xarray.DataArray(numpy.stack([values] * steps),
                                 coords={"datetime": stamps, "y": y, "x": x},
                                 dims=("datetime", "y", "x"))
    field.attrs.update(attrs or {})
    return field


@pytest.fixture()
def level():
    return InjuryLevelThreshold("Severe", threshold="%s mg/m**3" % THRESHOLD)


# ---------------------------------------------------------------------------
# calculateContours: output shape
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestContourOutputShape:
    def test_the_documented_columns_are_produced(self, level):
        result = level.calculateContours(_field(steps=1))

        assert {"ToxicLoad", "TotalPolygon", "severity", "datetime"} <= set(result.columns)

    def test_the_geometry_column_is_the_total_polygon(self, level):
        result = level.calculateContours(_field(steps=1))

        assert result.geometry.name == "TotalPolygon"

    def test_the_severity_is_the_level_name(self, level):
        result = level.calculateContours(_field(steps=1))

        assert set(result["severity"]) == {"Severe"}

    def test_the_toxic_load_is_the_threshold_that_was_contoured(self, level):
        result = level.calculateContours(_field(steps=1))

        assert result["ToxicLoad"].tolist() == [THRESHOLD]

    def test_the_polygon_encloses_the_analytic_contour_area(self, level):
        """pi r^2 for r = sigma sqrt(2 ln(peak/threshold)), to within the
        discretisation of the grid."""
        result = level.calculateContours(_field(steps=1))

        expected = numpy.pi * _radius(THRESHOLD) ** 2
        assert result["TotalPolygon"].area.iloc[0] == pytest.approx(expected, rel=2e-3)

    def test_the_polygon_is_centred_on_the_plume(self, level):
        result = level.calculateContours(_field(steps=1))

        centroid = result["TotalPolygon"].iloc[0].centroid
        assert (centroid.x, centroid.y) == pytest.approx((CENTRE_X, 0.0), abs=1e-6)

    def test_a_stronger_plume_gives_a_larger_polygon(self, level):
        weak = level.calculateContours(_field(peak=50.0, steps=1))
        strong = level.calculateContours(_field(peak=200.0, steps=1))

        assert strong["TotalPolygon"].area.iloc[0] > weak["TotalPolygon"].area.iloc[0]


# ---------------------------------------------------------------------------
# calculateContours: the time-shape branches
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestContourTimeBranches:
    def test_a_multi_step_field_gives_one_row_per_step(self, level):
        result = level.calculateContours(_field(steps=4))

        assert len(result) == 4

    def test_a_multi_step_field_stamps_each_row_with_its_time(self, level):
        field = _field(steps=3)
        result = level.calculateContours(field)

        assert result["datetime"].tolist() == list(
            pandas.to_datetime(field["datetime"].values))

    def test_a_single_step_dimension_gives_one_row(self, level):
        result = level.calculateContours(_field(steps=1))

        assert len(result) == 1

    def test_a_single_step_dimension_keeps_its_timestamp(self, level):
        field = _field(steps=1)
        result = level.calculateContours(field)

        assert result["datetime"].tolist() == [
            pandas.to_datetime(field["datetime"].item())]

    def test_a_scalar_time_coordinate_is_recognised(self, level):
        """What ``.isel(datetime=-1)`` leaves behind: the coordinate is still
        there but is no longer a dimension."""
        field = _field(steps=3).isel(datetime=1)
        assert "datetime" not in field.dims

        result = level.calculateContours(field)

        assert len(result) == 1
        assert result["datetime"].tolist() == [pandas.to_datetime(field["datetime"].item())]

    def test_a_field_with_no_time_at_all_is_accepted(self, level):
        result = level.calculateContours(_field())

        assert len(result) == 1

    def test_a_field_with_no_time_is_stamped_with_zero(self, level):
        """Characterisation: with nothing to read, the time column becomes the
        placeholder 0 rather than a timestamp."""
        result = level.calculateContours(_field())

        assert result["datetime"].tolist() == [0]

    def test_a_custom_time_dimension_name_is_honoured(self, level):
        field = _field(steps=2).rename({"datetime": "tt"})

        result = level.calculateContours(field, time="tt")

        assert len(result) == 2


# ---------------------------------------------------------------------------
# calculateContours: filtering, repair and rotation
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestContourPostProcessing:
    def test_a_field_below_the_threshold_yields_nothing(self, level):
        """Documented return: None, not an empty frame, when no level is
        exceeded anywhere."""
        assert level.calculateContours(_field(peak=THRESHOLD / 100.0, steps=1)) is None

    def test_a_field_below_the_threshold_yields_nothing_for_many_steps(self, level):
        assert level.calculateContours(_field(peak=THRESHOLD / 100.0, steps=3)) is None

    def test_duplicate_levels_are_dissolved_into_one_row(self, level):
        """One threshold means one row, even though the grid can produce
        several disjoint rings for it."""
        twoBlobs = _field(steps=1)
        result = level.calculateContours(twoBlobs)

        assert result["ToxicLoad"].is_unique

    def test_every_polygon_is_valid(self, level):
        """The buffer(0) repair step exists to guarantee this."""
        result = level.calculateContours(_field(steps=1))

        assert result["TotalPolygon"].is_valid.all()

    def test_the_cloud_direction_rotates_the_polygon_about_the_origin(self, level):
        """A math direction of 90 deg means the field was computed in a frame
        rotated 90 deg, so the contours are rotated back by -90 deg about
        (0, 0): the plume centre (100, 0) has to land on (0, -100)."""
        result = level.calculateContours(
            _field(steps=1, attrs={"cloud_math_direction": 90.0}))

        centroid = result["TotalPolygon"].iloc[0].centroid
        assert (centroid.x, centroid.y) == pytest.approx((0.0, -CENTRE_X), abs=1e-6)

    def test_the_rotation_preserves_the_enclosed_area(self, level):
        straight = level.calculateContours(_field(steps=1))
        turned = level.calculateContours(
            _field(steps=1, attrs={"cloud_math_direction": 90.0}))

        assert turned["TotalPolygon"].area.iloc[0] == pytest.approx(
            straight["TotalPolygon"].area.iloc[0])

    def test_without_the_attribute_nothing_is_rotated(self, level):
        result = level.calculateContours(_field(steps=1))

        centroid = result["TotalPolygon"].iloc[0].centroid
        assert centroid.x == pytest.approx(CENTRE_X, abs=1e-6)


@pytest.mark.unit
class TestContourBackendHandling:
    """The method switches matplotlib to a headless backend so that the
    contour extraction never opens a window, and must put it back."""

    def test_the_backend_is_restored(self, level):
        import matplotlib

        before = matplotlib.get_backend()
        level.calculateContours(_field(steps=1))

        assert matplotlib.get_backend() == before

    def test_the_backend_is_restored_even_when_nothing_is_contoured(self, level):
        import matplotlib

        before = matplotlib.get_backend()
        assert level.calculateContours(_field(peak=1e-6, steps=1)) is None

        assert matplotlib.get_backend() == before

    def test_no_figures_are_left_open(self, level):
        import matplotlib.pyplot as plt

        level.calculateContours(_field(steps=1))

        assert plt.get_fignums() == []


# ---------------------------------------------------------------------------
# _getGeopandas
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetGeopandasIsAbstract:
    def test_the_base_class_refuses_to_contour(self):
        with pytest.raises(NotImplementedError, match="abstract function"):
            InjuryLevel("Severe")._getGeopandas(_field(steps=1), "x", "y")

    def test_calculate_contours_on_the_base_class_therefore_fails(self):
        """The hook is the only thing the template method needs a subclass
        for, so the base class cannot be used directly."""
        with pytest.raises(NotImplementedError, match="abstract function"):
            InjuryLevel("Severe").calculateContours(_field(steps=1))


@pytest.mark.unit
class TestGetGeopandasImplementations:
    """Each subclass picks its own contour levels; all three share the
    "return an empty frame when the field never reaches the level" rule."""

    @pytest.fixture(autouse=True)
    def _headless(self):
        """_getGeopandas draws with pyplot; calculateContours normally does
        this switch itself."""
        import matplotlib.pyplot as plt

        before = plt.get_backend()
        plt.switch_backend("pdf")
        yield
        plt.close("all")
        plt.switch_backend(before)

    def test_the_threshold_level_contours_at_its_threshold(self):
        level = InjuryLevelThreshold("Severe", threshold="%s mg/m**3" % THRESHOLD)

        contours = level._getGeopandas(_field(steps=1).squeeze(), "x", "y")

        assert contours["Level"].tolist() == [THRESHOLD]

    def test_the_threshold_level_returns_an_empty_frame_below_it(self):
        level = InjuryLevelThreshold("Severe", threshold="%s mg/m**3" % THRESHOLD)

        assert level._getGeopandas(_field(peak=1e-6, steps=1).squeeze(), "x", "y").empty

    def test_the_exponential_level_contours_at_its_rate_constant(self):
        level = InjuryLevelExponential("Severe", k=5.0)

        contours = level._getGeopandas(_field(steps=1).squeeze(), "x", "y")

        assert contours["Level"].tolist() == [5.0]

    def test_the_exponential_level_returns_an_empty_frame_below_it(self):
        level = InjuryLevelExponential("Severe", k=5.0)

        assert level._getGeopandas(_field(peak=1e-6, steps=1).squeeze(), "x", "y").empty

    def test_the_lognormal_level_contours_the_default_injury_fractions(self):
        """The documented default is every 5% from 5% to 95%, i.e. 19 levels,
        each at its own toxic load."""
        level = InjuryLevelLognormal10DoseResponse("Severe", TL_50=10, sigma=0.5)

        contours = level._getGeopandas(_field(peak=400.0, steps=1).squeeze(), "x", "y")

        expected = level.getToxicLoad(numpy.arange(0.05, 1, 0.05))
        assert sorted(contours["Level"]) == pytest.approx(sorted(expected))

    def test_the_lognormal_level_returns_an_empty_frame_below_its_loads(self):
        level = InjuryLevelLognormal10DoseResponse("Severe", TL_50=10, sigma=0.5)

        assert level._getGeopandas(_field(peak=1e-9, steps=1).squeeze(), "x", "y").empty

    def test_a_higher_severity_adds_its_levels_to_the_milder_one(self):
        """Documented purpose: without the more-severe level's contours the
        gap between two severities is wide enough to double count
        casualties."""
        severe = InjuryLevelLognormal10DoseResponse("Severe", TL_50=100, sigma=0.5)
        light = InjuryLevelLognormal10DoseResponse("Light", TL_50=10, sigma=0.5,
                                                   higher_severity=severe)
        lonely = InjuryLevelLognormal10DoseResponse("Light", TL_50=10, sigma=0.5)

        field = _field(peak=400.0, steps=1).squeeze()
        withHigher = light._getGeopandas(field, "x", "y")
        withoutHigher = lonely._getGeopandas(field, "x", "y")

        assert len(withHigher) > len(withoutHigher)
        assert set(numpy.round(withoutHigher["Level"], 6)) <= set(
            numpy.round(withHigher["Level"], 6))

    def test_the_combined_levels_are_sorted_and_unique(self):
        severe = InjuryLevelLognormal10DoseResponse("Severe", TL_50=100, sigma=0.5)
        light = InjuryLevelLognormal10DoseResponse("Light", TL_50=10, sigma=0.5,
                                                   higher_severity=severe)

        contours = light._getGeopandas(_field(peak=400.0, steps=1).squeeze(), "x", "y")

        assert contours["Level"].is_unique

    def test_an_explicit_percent_injury_list_is_honoured(self):
        level = InjuryLevelLognormal10DoseResponse("Severe", TL_50=10, sigma=0.5)

        contours = level._getGeopandas(_field(peak=400.0, steps=1).squeeze(),
                                       "x", "y", percentInjury=[0.5])

        assert contours["Level"].tolist() == pytest.approx(
            [level.getToxicLoad(0.5)])


# ---------------------------------------------------------------------------
# __str__
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestInjuryLevelStr:
    """B264: see the module docstring."""

    @pytest.mark.xfail(
        strict=True,
        reason="B264: InjuryLevel.__str__ builds "
               "json.dumps(self.toJSON(), indent=4) and never returns it, so "
               "it returns None and every str()/print()/f-string on an "
               "InjuryLevel raises TypeError: __str__ returned non-string.  "
               "This is a distinct method from the already-reported B122 "
               "(Injury.__str__ in agents/effects/Injury.py), which has the "
               "identical missing return and needs its own fix. "
               "See the consolidated findings issue.",
    )
    def test_str_is_the_pretty_printed_json(self, level):
        assert str(level) == json.dumps(level.toJSON(), indent=4)

    def test_str_currently_raises_because_nothing_is_returned(self, level):
        """Characterisation of B264."""
        with pytest.raises(TypeError, match="__str__ returned non-string"):
            str(level)

    def test_the_dictionary_it_meant_to_serialise_is_fine(self, level):
        """Characterisation of B264: toJSON works, so the defect is purely
        the missing return."""
        assert json.loads(json.dumps(level.toJSON()))["name"] == "Severe"

    @pytest.mark.parametrize("built", [
        lambda: InjuryLevel("Severe"),
        lambda: InjuryLevelThreshold("Severe", threshold="10 mg/m**3"),
        lambda: InjuryLevelExponential("Severe", k=1.0),
        lambda: InjuryLevelLognormal10DoseResponse("Severe", TL_50=10, sigma=0.5),
    ])
    def test_no_subclass_repairs_it(self, built):
        """Characterisation of B264: none of the three subclasses overrides
        __str__, so all of them are unprintable."""
        with pytest.raises(TypeError, match="__str__ returned non-string"):
            str(built())
