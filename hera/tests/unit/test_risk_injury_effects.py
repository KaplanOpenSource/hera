r"""``InjuryFactory`` and the ``Injury`` family: building an effect from its
JSON descriptor, and the level algebra layered on top of it.

An ``Injury`` is a *stack* of ``InjuryLevel`` objects ordered most-severe
first.  Two rules follow from that ordering and drive most of this file:

* ``getPercent(level, TL)`` reports the *differential* fraction -- the
  fraction affected at that level minus the fraction already affected at
  the next-more-severe one -- so that summing over the levels never double
  counts a casualty.  For the log-normal levels used here the reference
  value is computable straight from the documented probit,
  ``Phi(log10(TL/TL_50)/sigma)``, which is what the assertions compare
  against rather than against hera's own output.
* ``calculateThresholdPolygon`` performs the same subtraction in space:
  contours of a monotonically decaying field are nested, so each ring is
  the polygon minus the next-higher-load polygon, and the rings tile the
  outermost polygon exactly.

The end-to-end ``calculateRegionOfInjured`` tests run on an analytic
Gaussian plume, C(r) = C0 exp(-r^2/2s^2), integrated by the ten Berge
calculator with n = 1 and dt = 1 min.  The dose after k steps is therefore
(k+1)*C(r), so the threshold contour sits at r^2 = 2 s^2 ln((k+1) C0 / TL),
and the enclosed area pi r^2 is known in closed form -- no reliance on the
implementation for the expected numbers.

Deliberately not covered here:

* ``calculateToxicLoads`` / ``calculatePointWiseFractionInjured`` on pandas
  input, and the deprecated ``calculate``'s pandas path -- all of that,
  plus bugs B62-B64, already lives in ``test_risk_injury_pandas_path.py``.
* ``InjuryLevel.calculateContours`` / ``_getGeopandas`` internals and the
  probit itself -- see ``test_risk_injurylevel*.py``.  The timestamps the
  contour code stamps onto its output are its business, so the tests below
  assert row counts and geometry rather than datetime values.
* ``Injury.getPercent`` on a *threshold* injury is only characterised, not
  asserted: it dies inside ``InjuryLevelThreshold.getPercent`` on the
  already-recorded B39/B12 cross-registry ``tounit`` defect.

Bugs found while writing this file are pinned below as xfail(strict)/
characterisation pairs: B122 (``Injury.__str__`` returns None),
B123 (the "injuries found" list is mangled), B124 (an empty
``calculator`` mapping does not raise the documented ValueError) and
B125 (a descriptor-level ``units`` string becomes a Quantity where a
pint Unit is documented).
"""
import json

import numpy
import pandas
import pytest
from pint import Quantity, Unit
from scipy.stats import norm
from shapely.geometry import Polygon

from hera.riskassessment.agents.effects import injuryfactory
from hera.riskassessment.agents.effects.Calculator import CalculatorTenBerge
from hera.riskassessment.agents.effects.Injury import (
    Injury,
    InjuryExponential,
    InjuryLognormal10,
    InjuryThreshold,
)
from hera.riskassessment.agents.effects.thresholdGeoDataFrame import thresholdGeoDataFrame
from hera.utils.unitHandler import ureg

SIGMA = 0.5
SEVERE_TL50 = 100.0
LIGHT_TL50 = 10.0


def _probit(toxicLoad, TL_50, sigma=SIGMA):
    """The documented base-10 log-normal dose response, derived here from
    the standard normal CDF rather than from hera's implementation."""
    return norm.cdf(numpy.log10(toxicLoad / TL_50) / sigma)


def _lognormal_descriptor(**extra):
    """The effect-descriptor shape ``Agent``/``RiskToolkit`` build inline,
    with two log-normal levels ordered most-severe first."""
    descriptor = {
        "type": "Lognormal10",
        "calculator": {"TenBerge": {"tenbergeCoefficient": 1}},
        "parameters": {
            "type": "Lognormal10DoseResponse",
            "levels": ["Severe", "Light"],
            "parameters": {
                "Severe": {"TL_50": SEVERE_TL50, "sigma": SIGMA},
                "Light": {"TL_50": LIGHT_TL50, "sigma": SIGMA},
            },
        },
    }
    descriptor.update(extra)
    return descriptor


def _threshold_descriptor(levels=None, **extra):
    levels = levels or {"Severe": "10*mg/m**3"}
    descriptor = {
        "type": "Threshold",
        "calculator": {"TenBerge": {"tenbergeCoefficient": 1}},
        "parameters": {
            "type": "Threshold",
            "levels": list(levels),
            "parameters": {name: {"threshold": value} for name, value in levels.items()},
        },
    }
    descriptor.update(extra)
    return descriptor


def _exponential_descriptor(k=0.1, **extra):
    descriptor = {
        "type": "Exponential",
        "calculator": {"TenBerge": {"tenbergeCoefficient": 1}},
        "parameters": {
            "type": "Exponential",
            "levels": ["Severe"],
            "parameters": {"Severe": {"k": k}},
        },
    }
    descriptor.update(extra)
    return descriptor


@pytest.fixture()
def lognormal_injury():
    return injuryfactory.getInjury("inhalation", _lognormal_descriptor())


@pytest.fixture()
def threshold_injury():
    return injuryfactory.getInjury("inhalation", _threshold_descriptor())


# ---------------------------------------------------------------------------
# InjuryFactory.getInjury
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestInjuryFactoryDispatch:
    @pytest.mark.parametrize(
        "descriptor,expected",
        [
            (_threshold_descriptor(), InjuryThreshold),
            (_lognormal_descriptor(), InjuryLognormal10),
            (_exponential_descriptor(), InjuryExponential),
        ],
        ids=["Threshold", "Lognormal10", "Exponential"],
    )
    def test_the_type_key_selects_the_injury_class(self, descriptor, expected):
        assert type(injuryfactory.getInjury("inhalation", descriptor)) is expected

    def test_the_calculator_key_selects_the_calculator_class(self, lognormal_injury):
        assert isinstance(lognormal_injury.calculator, CalculatorTenBerge)

    def test_the_calculator_block_configures_the_calculator(self):
        injury = injuryfactory.getInjury(
            "inhalation", _lognormal_descriptor(calculator={"TenBerge": {"tenbergeCoefficient": 2.5}})
        )
        assert injury.calculator.n == pytest.approx(2.5)

    def test_additional_parameters_are_passed_on_to_the_calculator(self):
        """The documented purpose of ``**additionalparameters`` -- this is
        how ``Agent`` hands its effectParameters (the ten Berge
        coefficient) down to every effect."""
        descriptor = _lognormal_descriptor(calculator={"TenBerge": {}})
        injury = injuryfactory.getInjury("inhalation", descriptor, tenbergeCoefficient=3)
        assert injury.calculator.n == 3

    def test_an_additional_parameter_that_the_block_also_sets_collides(self):
        """Robustness note (not pinned as a bug): the two parameter sources
        are splatted into one call, so a key present in both is a duplicate
        keyword argument rather than one overriding the other."""
        descriptor = _lognormal_descriptor(calculator={"TenBerge": {"tenbergeCoefficient": 1}})
        with pytest.raises(TypeError, match="tenbergeCoefficient"):
            injuryfactory.getInjury("inhalation", descriptor, tenbergeCoefficient=3)

    def test_the_effect_name_is_stored_on_the_injury(self):
        """``Injury`` keeps the effect name in ``_name`` but -- unlike
        ``InjuryLevel`` and ``InjuryFactory`` -- exposes no ``name``
        property for it, so the private attribute is the only way in."""
        injury = injuryfactory.getInjury("inhalation", _threshold_descriptor())
        assert injury._name == "inhalation"
        assert not hasattr(type(injury), "name")


@pytest.mark.unit
class TestInjuryFactoryValidation:
    def test_a_descriptor_without_a_type_is_rejected(self):
        descriptor = _threshold_descriptor()
        del descriptor["type"]
        with pytest.raises(ValueError, match="Injury type is not defined"):
            injuryfactory.getInjury("inhalation", descriptor)

    def test_a_descriptor_without_a_calculator_is_rejected(self):
        descriptor = _threshold_descriptor()
        del descriptor["calculator"]
        with pytest.raises(ValueError, match="Calculator not defined"):
            injuryfactory.getInjury("inhalation", descriptor)

    @pytest.mark.xfail(
        strict=True,
        reason="B124: an empty `calculator` mapping is exactly the "
               "'Calculator not defined' case the code means to report, but "
               "the guard only catches KeyError -- taking [0] of the empty "
               "items() list raises IndexError instead, so callers cannot "
               "distinguish a malformed descriptor from an internal error. "
               "See the consolidated findings issue.",
    )
    def test_an_empty_calculator_mapping_reports_an_undefined_calculator(self):
        with pytest.raises(ValueError, match="Calculator not defined"):
            injuryfactory.getInjury("inhalation", _threshold_descriptor(calculator={}))

    def test_an_empty_calculator_mapping_currently_raises_indexerror(self):
        """Characterisation of B124."""
        with pytest.raises(IndexError, match="list index out of range"):
            injuryfactory.getInjury("inhalation", _threshold_descriptor(calculator={}))

    def test_an_unknown_injury_type_is_reported_as_not_implemented(self):
        with pytest.raises(NotImplementedError, match="The injury Bogus is not defined"):
            injuryfactory.getInjury("inhalation", _threshold_descriptor(type="Bogus"))

    @staticmethod
    def _injuries_found():
        with pytest.raises(NotImplementedError) as raised:
            injuryfactory.getInjury("inhalation", _threshold_descriptor(type="Bogus"))
        listing = str(raised.value).split("Injuries found:")[1]
        return {name.strip() for name in listing.split(",")}

    @pytest.mark.xfail(
        strict=True,
        reason="B123: the 'Injuries found' list is built by stripping six "
               "characters off every dir() entry starting with 'Injury', so "
               "it advertises 'Factory' (from InjuryFactory, not an injury "
               "at all) and an empty name (from the abstract Injury base), "
               "while omitting nothing else -- the hint names things the "
               "caller cannot use. See the consolidated findings issue.",
    )
    def test_the_hint_lists_exactly_the_implemented_injuries(self):
        assert self._injuries_found() == {"Lognormal10", "Threshold", "Exponential"}

    def test_the_hint_currently_advertises_the_factory_and_the_base_class(self):
        """Characterisation of B123."""
        found = self._injuries_found()
        assert "Factory" in found
        assert "" in found

    def test_an_unknown_calculator_name_currently_fails_opaquely(self):
        """Robustness note (not pinned as a bug): the injury type gets a
        helpful NotImplementedError, but an unresolvable calculator name is
        simply called as None."""
        with pytest.raises(TypeError, match="NoneType"):
            injuryfactory.getInjury("inhalation", _threshold_descriptor(calculator={"Bogus": {}}))


# ---------------------------------------------------------------------------
# Injury.__init__ and the level accessors
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestInjuryConstruction:
    def test_the_levels_are_built_in_the_declared_order(self, lognormal_injury):
        assert lognormal_injury.levelNames == ["Severe", "Light"]

    def test_every_level_is_reachable_by_name(self, lognormal_injury):
        assert lognormal_injury["Severe"] is lognormal_injury.levels[0]
        assert lognormal_injury["Light"] is lognormal_injury.levels[1]

    def test_an_unknown_level_name_raises_keyerror(self, lognormal_injury):
        with pytest.raises(KeyError, match="Fatal"):
            lognormal_injury["Fatal"]

    def test_the_level_parameters_reach_the_level_objects(self, lognormal_injury):
        assert lognormal_injury["Severe"].TL_50.m_as(ureg.mg / ureg.m**3) == pytest.approx(SEVERE_TL50)
        assert lognormal_injury["Light"].sigma == pytest.approx(SIGMA)

    def test_the_most_severe_level_has_nothing_above_it(self, lognormal_injury):
        assert lognormal_injury["Severe"]._parameters["higher_severity"] is None

    def test_each_later_level_points_at_the_level_before_it(self, lognormal_injury):
        """The chain is what stops the milder level's contours from double
        counting the casualties of the severe one."""
        assert lognormal_injury["Light"]._parameters["higher_severity"] is lognormal_injury["Severe"]

    def test_the_descriptor_is_not_mutated_by_construction(self):
        """``higher_severity`` is injected into a copy of each level's
        parameter dict, so the caller's JSON stays serialisable."""
        descriptor = _lognormal_descriptor()
        injuryfactory.getInjury("inhalation", descriptor)
        assert descriptor["parameters"]["parameters"]["Light"] == {"TL_50": LIGHT_TL50, "sigma": SIGMA}

    def test_a_missing_level_type_is_rejected(self):
        descriptor = _lognormal_descriptor()
        del descriptor["parameters"]["type"]
        # the message carries a "nor defined" typo for "not defined"
        with pytest.raises(ValueError, match="InjuryLevel type is nor defined"):
            injuryfactory.getInjury("inhalation", descriptor)

    def test_a_level_without_parameters_is_rejected_by_name(self):
        descriptor = _lognormal_descriptor()
        del descriptor["parameters"]["parameters"]["Light"]
        with pytest.raises(ValueError, match="parameters for level Light not found"):
            injuryfactory.getInjury("inhalation", descriptor)

    def test_a_missing_levels_list_currently_fails_opaquely(self):
        """Robustness note (not pinned as a bug): the sibling checks raise
        ValueError, but an absent ``levels`` key defaults to None and is
        iterated."""
        descriptor = _lognormal_descriptor()
        del descriptor["parameters"]["levels"]
        with pytest.raises(TypeError, match="not iterable"):
            injuryfactory.getInjury("inhalation", descriptor)

    def test_an_unknown_level_type_currently_fails_opaquely(self):
        """Robustness note: pydoc.locate returns None and None is called."""
        descriptor = _lognormal_descriptor()
        descriptor["parameters"]["type"] = "Bogus"
        with pytest.raises(TypeError, match="NoneType"):
            injuryfactory.getInjury("inhalation", descriptor)


@pytest.mark.unit
class TestDescriptorUnits:
    """B125: the descriptor's optional ``units`` string."""

    def test_without_the_key_a_threshold_level_gets_its_default_unit(self, threshold_injury):
        assert threshold_injury["Severe"].units == ureg.mg / ureg.m**3
        assert isinstance(threshold_injury["Severe"].units, Unit)

    @pytest.mark.xfail(
        strict=True,
        reason="B125: Injury.__init__ converts the descriptor's units "
               "string with ureg(units), which yields a Quantity of "
               "magnitude 1 rather than the pint.Unit that InjuryLevel.units "
               "documents and that hera.utils.tounit accepts (Unit or str "
               "only) -- so a declared unit reaches every level as the wrong "
               "type and leaks into toJSON as '1.0 milligram / meter ** 3'. "
               "See the consolidated findings issue.",
    )
    def test_a_declared_unit_reaches_the_levels_as_a_pint_unit(self):
        injury = injuryfactory.getInjury("inhalation", _threshold_descriptor(units="mg/m**3"))
        assert isinstance(injury["Severe"].units, Unit)

    def test_a_declared_unit_currently_arrives_as_a_quantity(self):
        """Characterisation of B125."""
        injury = injuryfactory.getInjury("inhalation", _threshold_descriptor(units="mg/m**3"))
        units = injury["Severe"].units
        assert isinstance(units, Quantity)
        assert units.magnitude == pytest.approx(1.0)
        assert injury.toJSON()["levels"]["Severe"]["units"].startswith("1.0 ")

    def test_the_threshold_value_is_still_converted_to_the_declared_unit(self):
        injury = injuryfactory.getInjury(
            "inhalation", _threshold_descriptor(levels={"Severe": "0.01*g/m**3"}, units="mg/m**3")
        )
        assert injury["Severe"].threshold.magnitude == pytest.approx(10.0)


# ---------------------------------------------------------------------------
# Injury.getPercent
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetPercent:
    def test_the_most_severe_level_reports_its_own_dose_response(self, lognormal_injury):
        """Nothing is more severe, so no subtraction happens: at TL_50 half
        the population is severely injured."""
        assert lognormal_injury.getPercent("Severe", SEVERE_TL50) == pytest.approx(0.5)

    def test_a_milder_level_excludes_those_already_injured_more_severely(self, lognormal_injury):
        """At TL = 100 the light response is Phi(log10(10)/0.5) = Phi(2) =
        0.97725 and the severe one is 0.5, so 0.47725 of the population has
        light injuries as its worst outcome."""
        expected = _probit(SEVERE_TL50, LIGHT_TL50) - 0.5
        assert lognormal_injury.getPercent("Light", SEVERE_TL50) == pytest.approx(expected)
        assert lognormal_injury.getPercent("Light", SEVERE_TL50) == pytest.approx(0.47725, abs=1e-4)

    @pytest.mark.parametrize("toxicLoad", [1.0, 10.0, 100.0, 500.0])
    def test_the_levels_partition_the_affected_population(self, lognormal_injury, toxicLoad):
        """The differential fractions telescope: their sum is the mildest
        level's own cumulative response, never more."""
        total = sum(lognormal_injury.getPercent(name, toxicLoad) for name in lognormal_injury.levelNames)
        assert total == pytest.approx(_probit(toxicLoad, LIGHT_TL50))

    def test_the_differential_fraction_of_a_milder_level_vanishes_at_a_huge_dose(self, lognormal_injury):
        """Both responses saturate at 1, so their difference goes to 0 --
        at an overwhelming dose everyone is in the severe bracket."""
        assert lognormal_injury.getPercent("Light", 1e9) == pytest.approx(0.0, abs=1e-6)
        assert lognormal_injury.getPercent("Severe", 1e9) == pytest.approx(1.0)

    def test_every_differential_fraction_is_a_probability(self, lognormal_injury):
        for toxicLoad in (1e-3, 1.0, 10.0, 100.0, 1e4):
            for name in lognormal_injury.levelNames:
                assert 0.0 <= lognormal_injury.getPercent(name, toxicLoad) <= 1.0

    def test_a_level_that_is_not_configured_raises_valueerror(self, lognormal_injury):
        with pytest.raises(ValueError, match="not in list"):
            lognormal_injury.getPercent("Fatal", 100.0)

    def test_a_threshold_injury_cannot_report_a_percent_at_all(self, threshold_injury):
        """Characterisation of the recorded B39/B12 defect, reached here
        through the Injury API: InjuryLevelThreshold.getPercent builds its
        comparison value with the application-registry pint.Quantity while
        the threshold lives in hera's own registry."""
        with pytest.raises(ValueError, match="different registries"):
            threshold_injury.getPercent("Severe", 20.0)


# ---------------------------------------------------------------------------
# calculateThresholdPolygon
# ---------------------------------------------------------------------------

def _square(side, timestamp=0, severity="X", toxicLoad=1.0):
    return {
        "datetime": timestamp,
        "severity": severity,
        "ToxicLoad": toxicLoad,
        "TotalPolygon": Polygon([(0, 0), (side, 0), (side, side), (0, side)]),
    }


@pytest.mark.unit
class TestCalculateThresholdPolygon:
    """Contours of a decaying field are nested, so the ring belonging to a
    load is that polygon minus the next-higher-load one.  Concentric
    squares of side 10, 4 and 2 make the arithmetic exact: 100, 16 and 4."""

    @staticmethod
    def _nested(timestamp=0):
        return pandas.DataFrame([
            _square(10, timestamp, "Light", 1.0),
            _square(4, timestamp, "Severe", 5.0),
            _square(2, timestamp, "Fatal", 9.0),
        ])

    def test_the_highest_load_keeps_its_whole_polygon(self, threshold_injury):
        result = threshold_injury.calculateThresholdPolygon(self._nested(), "datetime")
        rings = result.set_index("severity")["ThresholdPolygon"]
        assert rings["Fatal"].area == pytest.approx(4.0)

    def test_every_other_load_keeps_only_its_own_ring(self, threshold_injury):
        result = threshold_injury.calculateThresholdPolygon(self._nested(), "datetime")
        rings = result.set_index("severity")["ThresholdPolygon"]
        assert rings["Severe"].area == pytest.approx(16.0 - 4.0)
        assert rings["Light"].area == pytest.approx(100.0 - 16.0)

    def test_the_rings_tile_the_outermost_polygon_exactly(self, threshold_injury):
        result = threshold_injury.calculateThresholdPolygon(self._nested(), "datetime")
        assert sum(poly.area for poly in result["ThresholdPolygon"]) == pytest.approx(100.0)

    def test_the_original_rows_and_columns_survive(self, threshold_injury):
        data = self._nested()
        result = threshold_injury.calculateThresholdPolygon(data, "datetime")
        assert len(result) == len(data)
        assert set(data.columns) <= set(result.columns)

    def test_each_time_step_is_differenced_on_its_own(self, threshold_injury):
        """A second time step with a different geometry must not be
        differenced against the first one's polygons."""
        data = pandas.concat(
            [self._nested(timestamp=0), pandas.DataFrame([
                _square(20, 1, "Light", 1.0), _square(10, 1, "Severe", 5.0),
            ])],
            ignore_index=True,
        )
        result = threshold_injury.calculateThresholdPolygon(data, "datetime")
        later = result[result["datetime"] == 1].set_index("severity")["ThresholdPolygon"]
        assert later["Severe"].area == pytest.approx(100.0)
        assert later["Light"].area == pytest.approx(400.0 - 100.0)


# ---------------------------------------------------------------------------
# _postCalculate: the per-subclass fraction bookkeeping
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPostCalculateDispatch:
    def test_the_base_class_refuses_to_post_calculate(self, threshold_injury):
        """The percent bookkeeping depends on the dose-response law, so the
        base class declares it abstract."""
        with pytest.raises(NotImplementedError, match="Abstract class"):
            Injury._postCalculate(threshold_injury, [], "datetime")

    def test_the_point_wise_post_calculation_is_a_no_op(self, threshold_injury):
        assert threshold_injury._postCalculatePointWise([]) is None


@pytest.mark.unit
class TestInjuryThresholdPostCalculate:
    def test_everyone_inside_a_threshold_contour_counts_as_injured(self, threshold_injury):
        """A threshold law has no partial response: the fraction is 1
        wherever the dose exceeds the threshold."""
        frame = pandas.DataFrame([_square(10, 0, "Severe", 10.0)])
        result = threshold_injury._postCalculate([frame], "datetime")
        assert result["percentEffected"].tolist() == [1.0]

    def test_the_levels_are_differenced_into_rings(self):
        injury = injuryfactory.getInjury(
            "inhalation", _threshold_descriptor(levels={"Severe": "50*mg/m**3", "Light": "10*mg/m**3"})
        )
        frames = [
            pandas.DataFrame([_square(4, 0, "Severe", 50.0)]),
            pandas.DataFrame([_square(10, 0, "Light", 10.0)]),
        ]
        result = injury._postCalculate(frames, "datetime")
        rings = result.set_index("severity")["ThresholdPolygon"]
        assert rings["Severe"].area == pytest.approx(16.0)
        assert rings["Light"].area == pytest.approx(100.0 - 16.0)

    def test_no_contours_at_all_yields_no_result(self, threshold_injury):
        assert threshold_injury._postCalculate([], "datetime") is None


@pytest.mark.unit
class TestInjuryLognormal10PostCalculate:
    def test_the_severe_rows_carry_their_own_dose_response(self, lognormal_injury):
        frame = pandas.DataFrame([_square(4, 0, "Severe", SEVERE_TL50)])
        result = lognormal_injury._postCalculate([frame], "datetime")
        assert result["percentEffected"].iloc[0] == pytest.approx(0.5)

    def test_the_light_rows_carry_only_the_fraction_not_already_severe(self, lognormal_injury):
        frame = pandas.DataFrame([_square(10, 0, "Light", LIGHT_TL50)])
        result = lognormal_injury._postCalculate([frame], "datetime")
        expected = 0.5 - _probit(LIGHT_TL50, SEVERE_TL50)
        assert result["percentEffected"].iloc[0] == pytest.approx(expected)

    def test_the_rows_of_every_level_are_concatenated(self, lognormal_injury):
        frames = [
            pandas.DataFrame([_square(4, 0, "Severe", SEVERE_TL50)]),
            pandas.DataFrame([_square(10, 0, "Light", LIGHT_TL50)]),
        ]
        result = lognormal_injury._postCalculate(frames, "datetime")
        assert result["severity"].tolist() == ["Severe", "Light"]


@pytest.mark.unit
class TestInjuryExponentialPostCalculate:
    def test_the_percent_effected_follows_the_exponential_law(self):
        k, toxicLoad = 0.1, 5.0
        injury = injuryfactory.getInjury("inhalation", _exponential_descriptor(k=k))
        frame = pandas.DataFrame([_square(10, 0, "Severe", toxicLoad)])
        result = injury._postCalculate([frame], "datetime")
        assert result["percentEffected"].iloc[0] == pytest.approx(1 - numpy.exp(-k * toxicLoad))


# ---------------------------------------------------------------------------
# calculateRegionOfInjured, end to end on an analytic plume
# ---------------------------------------------------------------------------

PEAK = 100.0
SPREAD = 40.0


def _gaussian_plume(peak=PEAK, spread=SPREAD, extent=200.0, points=81, steps=2):
    """C(r) = peak * exp(-r^2 / 2 spread^2), constant in time."""
    import xarray

    axis = numpy.linspace(-extent, extent, points)
    xx, yy = numpy.meshgrid(axis, axis)
    field = peak * numpy.exp(-(xx**2 + yy**2) / (2 * spread**2))
    dataset = xarray.Dataset(
        {"C": (("datetime", "y", "x"), numpy.stack([field] * steps))},
        coords={
            "datetime": pandas.date_range("2020-01-01", periods=steps, freq="1min"),
            "x": axis,
            "y": axis,
        },
    )
    dataset.attrs["C"] = 1 * ureg.mg / ureg.m**3
    dataset.attrs["dt"] = 1 * ureg.min
    return dataset


def _expected_area(threshold, step=0, peak=PEAK, spread=SPREAD):
    """Area enclosed by the dose contour after ``step`` one-minute steps.

    The ten Berge integral with n = 1 gives a dose of (step+1)*C(r), so the
    contour is at r^2 = 2 spread^2 ln((step+1) peak / threshold).
    """
    return numpy.pi * 2 * spread**2 * numpy.log((step + 1) * peak / threshold)


@pytest.mark.unit
class TestCalculateRegionOfInjured:
    def test_the_contour_encloses_the_area_where_the_dose_exceeds_the_threshold(self, threshold_injury):
        result = threshold_injury.calculateRegionOfInjured(_gaussian_plume(steps=1), field="C")
        assert result["TotalPolygon"].iloc[0].area == pytest.approx(_expected_area(10.0), rel=0.02)

    def test_the_contour_sits_at_the_configured_threshold(self, threshold_injury):
        result = threshold_injury.calculateRegionOfInjured(_gaussian_plume(steps=1), field="C")
        assert result["ToxicLoad"].tolist() == [10.0]
        assert result["severity"].tolist() == ["Severe"]

    def test_the_injured_area_grows_as_the_dose_accumulates(self, threshold_injury):
        result = threshold_injury.calculateRegionOfInjured(_gaussian_plume(steps=2), field="C")
        areas = [poly.area for poly in result["TotalPolygon"]]
        assert len(areas) == 2
        assert areas[0] == pytest.approx(_expected_area(10.0, step=0), rel=0.02)
        assert areas[1] == pytest.approx(_expected_area(10.0, step=1), rel=0.02)

    def test_the_result_is_a_threshold_geo_data_frame(self, threshold_injury):
        result = threshold_injury.calculateRegionOfInjured(_gaussian_plume(steps=1), field="C")
        assert isinstance(result, thresholdGeoDataFrame)

    def test_a_dose_that_never_reaches_the_threshold_yields_an_empty_frame(self, threshold_injury):
        result = threshold_injury.calculateRegionOfInjured(
            _gaussian_plume(peak=1.0, steps=1), field="C"
        )
        assert result.empty

    def test_isel_restricts_the_calculation_to_the_selected_steps(self, threshold_injury):
        result = threshold_injury.calculateRegionOfInjured(
            _gaussian_plume(steps=3), field="C", isel={"datetime": [0]}
        )
        assert len(result) == 1
        assert result["TotalPolygon"].iloc[0].area == pytest.approx(_expected_area(10.0), rel=0.02)

    def test_two_levels_come_back_as_nested_totals_and_disjoint_rings(self):
        injury = injuryfactory.getInjury(
            "inhalation", _threshold_descriptor(levels={"Severe": "50*mg/m**3", "Light": "10*mg/m**3"})
        )
        result = injury.calculateRegionOfInjured(_gaussian_plume(steps=1), field="C").set_index("severity")
        assert result.loc["Severe", "TotalPolygon"].area == pytest.approx(
            _expected_area(50.0), rel=0.02
        )
        assert result.loc["Light", "TotalPolygon"].area == pytest.approx(
            _expected_area(10.0), rel=0.02
        )
        assert result.loc["Light", "ThresholdPolygon"].area == pytest.approx(
            _expected_area(10.0) - _expected_area(50.0), rel=0.02
        )

    def test_a_lognormal_injury_reports_a_dose_response_per_contour(self):
        """Every contour is drawn at the dose giving a chosen fraction, so
        the severe rows must report exactly those fractions back."""
        injury = injuryfactory.getInjury("inhalation", _lognormal_descriptor())
        result = injury.calculateRegionOfInjured(_gaussian_plume(peak=300.0, steps=1), field="C")
        severe = result[result["severity"] == "Severe"]
        assert len(severe) > 1
        expected = _probit(severe["ToxicLoad"].values, SEVERE_TL50)
        assert severe["percentEffected"].values == pytest.approx(expected, abs=1e-6)

    def test_a_lognormal_milder_level_never_double_counts_the_severe_fraction(self):
        injury = injuryfactory.getInjury("inhalation", _lognormal_descriptor())
        result = injury.calculateRegionOfInjured(_gaussian_plume(peak=300.0, steps=1), field="C")
        light = result[result["severity"] == "Light"]
        expected = _probit(light["ToxicLoad"].values, LIGHT_TL50) - _probit(
            light["ToxicLoad"].values, SEVERE_TL50
        )
        assert light["percentEffected"].values == pytest.approx(expected, abs=1e-6)


@pytest.mark.unit
class TestDeprecatedCalculate:
    def test_it_warns_that_it_is_obsolete(self, threshold_injury, monkeypatch):
        monkeypatch.setattr(
            type(threshold_injury), "calculateRegionOfInjured", lambda self, **kwargs: "done"
        )
        with pytest.warns(UserWarning, match="obselete"):
            threshold_injury.calculate(_gaussian_plume(steps=1), field="C")

    def test_it_forwards_every_argument_to_the_replacement(self, threshold_injury, monkeypatch):
        captured = {}

        def _spy(self, **kwargs):
            captured.update(kwargs)
            return "done"

        monkeypatch.setattr(type(threshold_injury), "calculateRegionOfInjured", _spy)
        with pytest.warns(UserWarning):
            returned = threshold_injury.calculate(
                "the-field", field="C", time="tt", x="xx", y="yy",
                breathingRate=20 * ureg.L / ureg.min, sel={"a": 1}, isel={"b": 2},
            )

        assert returned == "done"
        assert captured["concentrationField"] == "the-field"
        assert captured["field"] == "C"
        assert (captured["time"], captured["x"], captured["y"]) == ("tt", "xx", "yy")
        assert captured["breathingRate"] == 20 * ureg.L / ureg.min
        assert (captured["sel"], captured["isel"]) == ({"a": 1}, {"b": 2})

    def test_the_default_breathing_rate_is_a_man_at_rest(self, threshold_injury, monkeypatch):
        captured = {}
        monkeypatch.setattr(
            type(threshold_injury), "calculateRegionOfInjured",
            lambda self, **kwargs: captured.update(kwargs),
        )
        with pytest.warns(UserWarning):
            threshold_injury.calculate("the-field", field="C")
        assert captured["breathingRate"] == 10 * ureg.L / ureg.min


# ---------------------------------------------------------------------------
# Serialisation
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestSerialisation:
    def test_to_json_holds_one_entry_per_level(self, lognormal_injury):
        serialised = lognormal_injury.toJSON()
        assert set(serialised["levels"]) == {"Severe", "Light"}
        assert serialised["levels"]["Severe"]["type"] == "logNormalBase10"

    def test_to_json_holds_the_calculator(self, lognormal_injury):
        serialised = lognormal_injury.toJSON()
        assert serialised["calculator"]["type"] == "tenBerge"
        assert serialised["calculator"]["n"] == "1"

    def test_to_json_is_json_serialisable(self, lognormal_injury):
        """The point of the dict form: it has to survive json.dumps."""
        assert json.loads(json.dumps(lognormal_injury.toJSON()))["levels"]["Light"]["name"] == "Light"

    @pytest.mark.xfail(
        strict=True,
        reason="B122: Injury.__str__ builds json.dumps(self.toJSON(), "
               "indent=4) but never returns it, so it returns None and "
               "every str()/print()/f-string on an Injury raises TypeError: "
               "__str__ returned non-string. (InjuryLevel.__str__ has the "
               "identical missing return.) See the consolidated findings "
               "issue.",
    )
    def test_str_is_the_pretty_printed_json(self, lognormal_injury):
        assert str(lognormal_injury) == json.dumps(lognormal_injury.toJSON(), indent=4)

    def test_str_currently_raises_because_nothing_is_returned(self, lognormal_injury):
        """Characterisation of B122."""
        with pytest.raises(TypeError, match="__str__ returned non-string"):
            str(lognormal_injury)
