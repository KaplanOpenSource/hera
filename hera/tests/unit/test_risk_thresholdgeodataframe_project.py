"""``thresholdGeoDataFrame._project``: the worker that lays injury isopleths
over a demographic layer.

``project`` (already covered in ``test_risk_thresholdgeodataframe.py``) only
dispatches over one or many angles; ``_project`` is where the work happens.
For every ``(severity, datetime)`` group it shifts and rotates that row's
threshold polygon to the release point, intersects it with the demography
after reprojecting to ITM, scales each population column by the row's
``percentEffected``, and stamps the severity, timestamp, toxic load and
percentage onto the result.  With nothing to report it returns ``None``
rather than an empty frame.

Only part of that is reachable today.  The intersection step calls
``_calculatePopulationInPolygon``, which evaluates
``demography.loc[not demography["geometry"].intersection(poly).is_empty]``
-- ``not`` on a pandas Series -- and so raises ``ValueError: The truth value
of a Series is ambiguous`` for any input.  That is the already-reported
B108, pinned in ``test_risk_thresholdgeodataframe.py``, and this file does
not duplicate the pin.  What it does instead is cover the parts of
``_project`` that sit *before* and *around* that call and are therefore
fully exercisable:

* the two documented ``None`` returns -- an empty frame, and a frame whose
  polygons are all invalid;
* the CRS conversion, the shift/rotate of the polygons, and the grouping,
  all observed through a recording stand-in for the intersection step;
* the ``population`` argument's string-or-list handling.

Where a stand-in is used it replaces exactly one static method
(``_calculatePopulationInPolygon``) on the class under test via
``monkeypatch``, so ``_project``'s own arithmetic is the real code
throughout.

Deliberately not covered here: ``shiftLocationAndAngle`` / ``_shiftPolygons``
and the angle dispatch in ``project`` (with the already-reported B97 and
B65 notes), and ``_calculatePopulationInPolygon`` itself -- all in
``test_risk_thresholdgeodataframe.py``.
"""
import pandas
import pytest
from shapely.geometry import Polygon

from hera.measurements.GIS import ITM, WSG84
from hera.riskassessment.agents.effects.thresholdGeoDataFrame import thresholdGeoDataFrame

SEVERE_PERCENT = 0.4
LIGHT_PERCENT = 0.1


def _rect(x0, y0, width, height):
    return Polygon([(x0, y0), (x0 + width, y0), (x0 + width, y0 + height),
                    (x0, y0 + height)])


def _bowtie():
    """A self-intersecting ring: shapely reports it as invalid."""
    return Polygon([(0, 0), (10, 10), (10, 0), (0, 10)])


def _isopleths(rows):
    frame = pandas.DataFrame(rows)
    frame["datetime"] = pandas.to_datetime(frame["datetime"])
    return thresholdGeoDataFrame(frame, geometry="ThresholdPolygon")


def _two_severities():
    return _isopleths([
        {"severity": "Severe", "datetime": "2020-01-01 00:00",
         "ThresholdPolygon": _rect(0, -25, 50, 50),
         "percentEffected": SEVERE_PERCENT, "ToxicLoad": 100.0},
        {"severity": "Light", "datetime": "2020-01-01 00:00",
         "ThresholdPolygon": _rect(0, -50, 100, 100),
         "percentEffected": LIGHT_PERCENT, "ToxicLoad": 10.0},
    ])


#: The single-row frame used wherever exactly one polygon must be observed.
#: ``_project`` groups by (severity, datetime), and pandas visits the groups
#: in sorted order, so a two-severity frame would present "Light" first.
SINGLE_POLYGON = _rect(0, -25, 50, 50)


def _one_severity():
    return _isopleths([
        {"severity": "Severe", "datetime": "2020-01-01 00:00",
         "ThresholdPolygon": SINGLE_POLYGON,
         "percentEffected": SEVERE_PERCENT, "ToxicLoad": 100.0},
    ])


def _demography(crs=ITM):
    import geopandas

    frame = geopandas.GeoDataFrame(
        {"total_pop": [1000.0], "children": [200.0],
         "geometry": [_rect(0, 0, 1000, 1000)]},
        geometry="geometry", crs=ITM)
    return frame.to_crs(crs)


@pytest.fixture()
def recorder(monkeypatch):
    """Replaces only the intersection step (blocked by B108) and records what
    ``_project`` hands it, so the surrounding arithmetic can be asserted."""
    import geopandas

    calls = []

    def _fake(demography, poly, populationTypes):
        calls.append({"demography": demography, "poly": poly,
                      "populationTypes": list(populationTypes)})
        result = geopandas.GeoDataFrame(
            {"geometry": [poly], "areaFraction": [1.0]}, geometry="geometry")
        for name in populationTypes:
            if name in demography:
                result[name] = [float(demography[name].sum())]
        return result

    monkeypatch.setattr(type(_two_severities()), "_calculatePopulationInPolygon",
                        staticmethod(_fake))
    return calls


# ---------------------------------------------------------------------------
# The documented None returns
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestProjectReturnsNone:
    def test_an_empty_isopleth_frame_projects_to_nothing(self):
        """Documented return: None, not an empty frame."""
        import geopandas

        empty = thresholdGeoDataFrame(
            {"severity": pandas.Series([], dtype=object),
             "datetime": pandas.Series([], dtype="datetime64[ns]"),
             "ThresholdPolygon": geopandas.GeoSeries([], dtype=object),
             "percentEffected": pandas.Series([], dtype=float),
             "ToxicLoad": pandas.Series([], dtype=float)},
            geometry="ThresholdPolygon")

        assert empty._project(_demography(), loc=(0, 0), mathematical_angle=0.0) is None

    def test_an_invalid_polygon_is_skipped(self):
        """A self-intersecting ring cannot be intersected meaningfully, so the
        guard drops it -- and with nothing left the answer is None."""
        broken = _isopleths([
            {"severity": "Severe", "datetime": "2020-01-01 00:00",
             "ThresholdPolygon": _bowtie(),
             "percentEffected": 1.0, "ToxicLoad": 100.0},
        ])

        assert broken._project(_demography(), loc=(0, 0),
                               mathematical_angle=0.0) is None

    def test_the_guard_is_the_shifted_polygon_not_the_original(self, recorder):
        """The validity test happens after the shift/rotate, so a frame whose
        rows are valid still projects even though a sibling row is not."""
        mixed = _isopleths([
            {"severity": "Severe", "datetime": "2020-01-01 00:00",
             "ThresholdPolygon": _bowtie(),
             "percentEffected": 1.0, "ToxicLoad": 100.0},
            {"severity": "Light", "datetime": "2020-01-01 00:00",
             "ThresholdPolygon": _rect(0, -50, 100, 100),
             "percentEffected": LIGHT_PERCENT, "ToxicLoad": 10.0},
        ])

        result = mixed._project(_demography(), loc=(0, 0), mathematical_angle=0.0)

        assert result is not None
        assert result["severity"].tolist() == ["Light"]

    def test_an_empty_intersection_yields_nothing(self, monkeypatch):
        """When no population is caught, the documented answer is None."""
        import geopandas

        monkeypatch.setattr(
            type(_two_severities()), "_calculatePopulationInPolygon",
            staticmethod(lambda demography, poly, populationTypes:
                         geopandas.GeoDataFrame({"geometry": []}, geometry="geometry")))

        assert _two_severities()._project(_demography(), loc=(0, 0),
                                          mathematical_angle=0.0) is None


# ---------------------------------------------------------------------------
# What reaches the intersection step
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestProjectDrivesTheIntersection:
    def test_one_intersection_per_isopleth_row(self, recorder):
        _two_severities()._project(_demography(), loc=(0, 0), mathematical_angle=0.0)

        assert len(recorder) == 2

    def test_the_demography_is_converted_to_itm(self, recorder):
        """Documented reason: ITM is in metres, so the area-weighted
        population arithmetic is in m^2."""
        _two_severities()._project(_demography(crs=WSG84), loc=(0, 0),
                                   mathematical_angle=0.0)

        assert recorder[0]["demography"].crs.to_epsg() == ITM

    def test_an_already_projected_demography_is_left_in_itm(self, recorder):
        _two_severities()._project(_demography(crs=ITM), loc=(0, 0),
                                   mathematical_angle=0.0)

        assert recorder[0]["demography"].crs.to_epsg() == ITM

    def test_the_polygon_is_translated_to_the_release_point(self, recorder):
        _one_severity()._project(_demography(), loc=(300, 400),
                                 mathematical_angle=0.0)

        shifted = recorder[0]["poly"]
        assert shifted.bounds[0] == pytest.approx(SINGLE_POLYGON.bounds[0] + 300.0)
        assert shifted.bounds[1] == pytest.approx(SINGLE_POLYGON.bounds[1] + 400.0)

    def test_a_zero_angle_only_translates(self, recorder):
        _one_severity()._project(_demography(), loc=(100, 0),
                                 mathematical_angle=0.0)

        moved = recorder[0]["poly"]
        assert moved.bounds == pytest.approx(
            (SINGLE_POLYGON.bounds[0] + 100, SINGLE_POLYGON.bounds[1],
             SINGLE_POLYGON.bounds[2] + 100, SINGLE_POLYGON.bounds[3]))

    def test_the_angle_rotates_the_polygon_about_the_origin(self, recorder):
        """A 90 degree mathematical angle turns the polygon spanning
        x in [0, 50], y in [-25, 25] into one spanning x in [-25, 25],
        y in [0, 50]."""
        _one_severity()._project(_demography(), loc=(0, 0),
                                 mathematical_angle=90.0)

        turned = recorder[0]["poly"]
        assert turned.bounds == pytest.approx((-25.0, 0.0, 25.0, 50.0))

    def test_the_rotation_preserves_the_polygon_area(self, recorder):
        _one_severity()._project(_demography(), loc=(0, 0),
                                 mathematical_angle=37.0)

        assert recorder[0]["poly"].area == pytest.approx(SINGLE_POLYGON.area)

    def test_a_meteorological_angle_is_accepted(self, recorder):
        """meteorological 270 is mathematical 0, i.e. pure translation."""
        _one_severity()._project(_demography(), loc=(0, 0),
                                 meteorological_angle=270.0)

        assert recorder[0]["poly"].bounds == pytest.approx(SINGLE_POLYGON.bounds)

    def test_without_any_angle_it_refuses(self):
        with pytest.raises(ValueError, match="either met_angle or math_angle"):
            _two_severities()._project(_demography(), loc=(0, 0))

    def test_the_default_population_column_is_total_pop(self, recorder):
        _two_severities()._project(_demography(), loc=(0, 0), mathematical_angle=0.0)

        assert recorder[0]["populationTypes"] == ["total_pop"]

    def test_a_string_population_is_wrapped_in_a_list(self, recorder):
        _two_severities()._project(_demography(), loc=(0, 0), mathematical_angle=0.0,
                                   population="children")

        assert recorder[0]["populationTypes"] == ["children"]

    def test_a_list_of_populations_is_passed_through(self, recorder):
        _two_severities()._project(_demography(), loc=(0, 0), mathematical_angle=0.0,
                                   population=["total_pop", "children"])

        assert recorder[0]["populationTypes"] == ["total_pop", "children"]


# ---------------------------------------------------------------------------
# The assembled result
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestProjectResult:
    @pytest.fixture()
    def projected(self, recorder):
        return _two_severities()._project(_demography(), loc=(0, 0),
                                          mathematical_angle=0.0)

    def test_the_documented_metadata_columns_are_stamped_on(self, projected):
        assert {"severity", "datetime", "percentEffected", "ToxicLoad"} <= \
            set(projected.columns)

    def test_one_row_per_severity(self, projected):
        assert sorted(projected["severity"]) == ["Light", "Severe"]

    def test_the_timestamp_of_the_group_is_recorded(self, projected):
        assert projected["datetime"].unique().tolist() == [
            pandas.to_datetime("2020-01-01 00:00")]

    def test_the_rows_toxic_load_is_carried_through(self, projected):
        byLevel = projected.set_index("severity")
        assert byLevel.loc["Severe", "ToxicLoad"] == 100.0
        assert byLevel.loc["Light", "ToxicLoad"] == 10.0

    def test_the_rows_injured_fraction_is_carried_through(self, projected):
        byLevel = projected.set_index("severity")
        assert byLevel.loc["Severe", "percentEffected"] == SEVERE_PERCENT
        assert byLevel.loc["Light", "percentEffected"] == LIGHT_PERCENT

    def test_an_effected_column_is_added_per_population(self, projected):
        """The naming convention the presentation layer relies on:
        'effected' + the population column name."""
        assert "effectedtotal_pop" in projected.columns

    def test_the_effected_population_is_scaled_by_the_injured_fraction(self, projected):
        """1000 people in the polygon, 40% of them severely injured."""
        byLevel = projected.set_index("severity")
        assert byLevel.loc["Severe", "effectedtotal_pop"] == pytest.approx(
            1000.0 * SEVERE_PERCENT)
        assert byLevel.loc["Light", "effectedtotal_pop"] == pytest.approx(
            1000.0 * LIGHT_PERCENT)

    def test_the_raw_population_column_is_kept_alongside(self, projected):
        assert projected["total_pop"].tolist() == pytest.approx([1000.0, 1000.0])

    def test_each_requested_population_gets_its_own_effected_column(self, recorder):
        result = _two_severities()._project(
            _demography(), loc=(0, 0), mathematical_angle=0.0,
            population=["total_pop", "children"])

        assert {"effectedtotal_pop", "effectedchildren"} <= set(result.columns)
        byLevel = result.set_index("severity")
        assert byLevel.loc["Severe", "effectedchildren"] == pytest.approx(
            200.0 * SEVERE_PERCENT)

    def test_several_time_steps_all_appear(self, recorder):
        rows = []
        for minute in (0, 1, 2):
            rows.append({"severity": "Severe",
                         "datetime": "2020-01-01 00:%02d" % minute,
                         "ThresholdPolygon": _rect(0, -25, 50, 50),
                         "percentEffected": SEVERE_PERCENT, "ToxicLoad": 100.0})

        result = _isopleths(rows)._project(_demography(), loc=(0, 0),
                                           mathematical_angle=0.0)

        assert len(result["datetime"].unique()) == 3

    def test_the_original_frame_is_left_unchanged(self, recorder):
        isopleths = _two_severities()
        before = isopleths["ThresholdPolygon"].tolist()

        isopleths._project(_demography(), loc=(500, 500), mathematical_angle=45.0)

        assert isopleths["ThresholdPolygon"].tolist() == before

    def test_the_supplied_demography_is_left_unchanged(self, recorder):
        demography = _demography(crs=WSG84)
        before = demography.crs

        _two_severities()._project(demography, loc=(0, 0), mathematical_angle=0.0)

        assert demography.crs == before
