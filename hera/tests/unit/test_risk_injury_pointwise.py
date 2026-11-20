r"""``Injury.calculatePointWiseFractionInjured``: the per-device injury
fraction over time.

This is the *time-series* counterpart of ``calculateRegionOfInjured``.  It
takes a ``pandas.DataFrame`` whose index is the time and whose columns are
measurement points, integrates each column into a toxic load, and returns a
long frame of one row per ``(time, device, level)`` carrying the
*differential* injury fraction at that level.

Why the field-name workaround below.  ``calculateToxicLoads`` calls the
calculator without an ``inUnits`` argument and without a ``field``, so the
calculator falls into ``inUnits = concentrationField.attrs[field]`` -- with
``field=None`` -- whenever the input has an ``attrs`` mapping.  Since pandas
1.0 a DataFrame *does*, so the lookup raises ``KeyError: None``.  That is
the already-reported B63, characterised in
``test_risk_injury_pandas_path.py``, and it blocks this method's body
entirely.  Rather than repeat that pin, these tests put the concentration
units under the key the calculator actually asks for (``attrs[None]``) and
so reach the loop underneath, which nothing has exercised before.

The expected numbers are derived from the documented model, not from hera:
the ten Berge calculator with ``n = 1`` reduces to Haber's law, so a
constant concentration ``C`` sampled every minute gives cumulative toxic
loads ``C, 2C, 3C, ...`` (the final sample is consumed by the differencing
and does not produce a row).  A base-10 log-normal level with median
``TL_50`` and ``sigma`` then reports
``Phi(log10(TL / TL_50) / sigma)``, and ``Injury.getPercent`` subtracts the
next-more-severe level's fraction so that summing over levels never double
counts a casualty.

Deliberately not covered here:

* the non-pandas rejection, the B62-B64 unit defects and the deprecated
  ``calculate`` wrapper -- all in ``test_risk_injury_pandas_path.py``.
* ``_postCalculatePointWise``, which is a documented no-op hook, and the
  xarray path, which the method itself refuses ("Still not implemented").
* ``InjuryFactory``/``Injury`` construction and the level algebra -- see
  ``test_risk_injury_effects.py``.
"""
import numpy
import pandas
import pytest
from scipy.stats import norm

from hera.riskassessment.agents.effects import injuryfactory
from hera.utils import ureg

SEVERE_TL50 = 100.0
LIGHT_TL50 = 10.0
SIGMA = 0.5
CONCENTRATION = 10.0
SAMPLES = 4


def _probit(toxicLoad, median, sigma=SIGMA):
    """The documented base-10 log-normal dose response."""
    return norm.cdf(numpy.log10(numpy.asarray(toxicLoad, dtype=float) / median) / sigma)


def _descriptor(*levels):
    parameters = {"Severe": {"TL_50": SEVERE_TL50, "sigma": SIGMA},
                  "Light": {"TL_50": LIGHT_TL50, "sigma": SIGMA}}
    return {
        "type": "Lognormal10",
        "calculator": {"TenBerge": {}},
        "parameters": {
            "type": "Lognormal10DoseResponse",
            "levels": list(levels),
            "parameters": {name: parameters[name] for name in levels},
        },
    }


@pytest.fixture()
def injury():
    """Two log-normal levels, ordered most severe first, ten Berge n = 1."""
    return injuryfactory.getInjury("inhalation", _descriptor("Severe", "Light"),
                                   tenbergeCoefficient=1)


@pytest.fixture()
def single_level_injury():
    return injuryfactory.getInjury("inhalation", _descriptor("Severe"),
                                   tenbergeCoefficient=1)


def _timeseries(samples=SAMPLES, **devices):
    """A minute-spaced concentration series per device.

    ``attrs[None]`` carries the concentration units, which is the key
    ``calculateToxicLoads`` ends up asking the calculator for (see the
    module docstring).
    """
    devices = devices or {"P1": CONCENTRATION}
    index = pandas.to_datetime(
        ["2020-01-01 00:%02d" % minute for minute in range(samples)])
    index.name = "datetime"
    frame = pandas.DataFrame(
        {name: [float(value)] * samples for name, value in devices.items()},
        index=index)
    frame.attrs[None] = 1 * ureg.mg / ureg.m ** 3
    return frame


def _cumulative_loads(concentration=CONCENTRATION, samples=SAMPLES):
    """Haber's law on a minute-spaced constant series: C, 2C, 3C, ..."""
    return [concentration * step for step in range(1, samples)]


# ---------------------------------------------------------------------------
# Shape
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPointWiseShape:
    def test_the_documented_columns_are_produced(self, injury):
        result = injury.calculatePointWiseFractionInjured(_timeseries())

        assert list(result.columns) == ["injuryPercent", "deviceName", "level"]

    def test_the_time_stays_in_the_index(self, injury):
        """Documented input shape: "the time is an index (with the name
        'datetime')" -- and it survives into the output."""
        result = injury.calculatePointWiseFractionInjured(_timeseries())

        assert result.index.name == "datetime"

    def test_one_row_per_level_per_device_per_integrated_step(self, injury):
        result = injury.calculatePointWiseFractionInjured(
            _timeseries(P1=CONCENTRATION, P2=1.0))

        assert len(result) == 2 * 2 * len(_cumulative_loads())

    def test_the_last_sample_is_consumed_by_the_differencing(self, injury):
        """Haber's law needs a dt per step, so N samples yield N-1 loads."""
        result = injury.calculatePointWiseFractionInjured(_timeseries(samples=5))

        assert len(result) == 2 * 4

    def test_every_device_column_appears(self, injury):
        result = injury.calculatePointWiseFractionInjured(
            _timeseries(P1=CONCENTRATION, P2=1.0, P3=0.5))

        assert sorted(result["deviceName"].unique()) == ["P1", "P2", "P3"]

    def test_every_level_appears_under_its_own_name(self, injury):
        result = injury.calculatePointWiseFractionInjured(_timeseries())

        assert sorted(result["level"].unique()) == ["Light", "Severe"]

    def test_the_levels_are_emitted_in_the_declared_order(self, injury):
        result = injury.calculatePointWiseFractionInjured(_timeseries())

        assert result["level"].tolist()[0] == "Severe"
        assert result["level"].tolist()[-1] == "Light"

    def test_a_single_level_injury_gives_one_block(self, single_level_injury):
        result = single_level_injury.calculatePointWiseFractionInjured(_timeseries())

        assert result["level"].unique().tolist() == ["Severe"]
        assert len(result) == len(_cumulative_loads())

    def test_the_timestamps_are_the_leading_samples(self, injury):
        frame = _timeseries()
        result = injury.calculatePointWiseFractionInjured(frame)

        severe = result[result["level"] == "Severe"]
        assert severe.index.tolist() == list(frame.index[:-1])


# ---------------------------------------------------------------------------
# Values
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPointWiseValues:
    def test_the_most_severe_level_reports_the_plain_probit(self, injury):
        """Nothing is more severe than the first level, so no subtraction."""
        result = injury.calculatePointWiseFractionInjured(_timeseries())

        severe = result[result["level"] == "Severe"]["injuryPercent"].tolist()
        assert severe == pytest.approx(_probit(_cumulative_loads(), SEVERE_TL50))

    def test_a_milder_level_reports_the_difference_from_the_level_above(self, injury):
        """The subtraction is what stops a casualty being counted twice."""
        result = injury.calculatePointWiseFractionInjured(_timeseries())

        light = result[result["level"] == "Light"]["injuryPercent"].tolist()
        expected = (_probit(_cumulative_loads(), LIGHT_TL50)
                    - _probit(_cumulative_loads(), SEVERE_TL50))
        assert light == pytest.approx(expected)

    def test_the_fractions_are_probabilities(self, injury):
        result = injury.calculatePointWiseFractionInjured(
            _timeseries(P1=CONCENTRATION, P2=1000.0))

        assert result["injuryPercent"].min() >= 0.0
        assert result["injuryPercent"].max() <= 1.0

    def test_the_load_accumulates_so_the_fraction_rises_with_time(self, injury):
        result = injury.calculatePointWiseFractionInjured(_timeseries(samples=6))

        severe = result[result["level"] == "Severe"]["injuryPercent"].tolist()
        assert numpy.all(numpy.diff(severe) > 0)

    def test_a_higher_concentration_injures_more_at_the_same_time(self, injury):
        result = injury.calculatePointWiseFractionInjured(
            _timeseries(P1=1.0, P2=100.0))

        severe = result[result["level"] == "Severe"].set_index("deviceName")
        assert severe.loc["P2", "injuryPercent"].max() > \
            severe.loc["P1", "injuryPercent"].max()

    def test_a_device_is_evaluated_independently_of_its_neighbours(self, injury):
        """Each column is its own point, so adding a second device must not
        change the first one's answer."""
        alone = injury.calculatePointWiseFractionInjured(_timeseries(P1=CONCENTRATION))
        together = injury.calculatePointWiseFractionInjured(
            _timeseries(P1=CONCENTRATION, P2=1000.0))

        wanted = together[(together["deviceName"] == "P1")
                          & (together["level"] == "Severe")]["injuryPercent"]
        assert wanted.tolist() == pytest.approx(
            alone[alone["level"] == "Severe"]["injuryPercent"].tolist())

    def test_a_zero_concentration_injures_nobody(self, injury):
        result = injury.calculatePointWiseFractionInjured(_timeseries(P1=0.0))

        assert result["injuryPercent"].tolist() == pytest.approx([0.0] * len(result))

    def test_half_the_population_is_severely_injured_at_the_median_load(self, injury):
        """The anchor of the whole model: at TL_50 the probit is one half.
        A concentration of TL_50 mg/m^3 for one minute is exactly that load.
        """
        result = injury.calculatePointWiseFractionInjured(
            _timeseries(samples=2, P1=SEVERE_TL50))

        severe = result[result["level"] == "Severe"]["injuryPercent"].tolist()
        assert severe == pytest.approx([0.5])

    def test_a_doubled_breathing_rate_doubles_the_load(self, injury):
        """The calculator scales the dose by breathingRate / injury rate, and
        the injury levels were calibrated at 10 L/min."""
        atRest = injury.calculatePointWiseFractionInjured(
            _timeseries(samples=2, P1=SEVERE_TL50 / 2))
        working = injury.calculatePointWiseFractionInjured(
            _timeseries(samples=2, P1=SEVERE_TL50 / 2),
            breathingRate=20 * ureg.L / ureg.min)

        assert atRest[atRest["level"] == "Severe"]["injuryPercent"].tolist() == \
            pytest.approx(_probit([SEVERE_TL50 / 2], SEVERE_TL50))
        assert working[working["level"] == "Severe"]["injuryPercent"].tolist() == \
            pytest.approx([0.5])

    def test_the_supplied_frame_is_not_mutated(self, injury):
        frame = _timeseries()
        before = frame.copy(deep=True)

        injury.calculatePointWiseFractionInjured(frame)

        pandas.testing.assert_frame_equal(frame, before)

    def test_a_custom_time_index_name_is_honoured(self, injury):
        frame = _timeseries()
        frame.index.name = "tt"

        result = injury.calculatePointWiseFractionInjured(frame, time="tt")

        assert len(result) == 2 * len(_cumulative_loads())
