"""Model-evaluation statistics.

These are the standard atmospheric-dispersion evaluation metrics (Chang and
Hanna, 2004).  Their defining values are what the tests assert:

    FB    fractional bias        0 for a perfect model, bounded in [-2, 2]
    NMSE  normalised mean square error   0 for a perfect model
    FAC2  fraction within a factor of 2  1 for a perfect model
    NAD   normalised absolute difference 0 for a perfect model
    R     linear correlation             1 for a perfect model

A perfect-model fixture therefore pins five numbers at once, and the sign
conventions are checked with deliberately biased data.
"""
import numpy as np
import pandas as pd
import pytest

from hera.simulations.analysis.errorCalculation import errorCalculation


@pytest.fixture()
def metrics():
    return errorCalculation()


@pytest.fixture()
def perfect():
    values = [1.0, 2.0, 3.0, 4.0, 5.0]
    return pd.DataFrame({"model": values, "measure": values})


@pytest.fixture()
def over_predicting():
    """The model reads twice the measurement everywhere."""
    measured = [1.0, 2.0, 3.0, 4.0, 5.0]
    return pd.DataFrame({"model": [2 * v for v in measured], "measure": measured})


@pytest.fixture()
def under_predicting():
    measured = [2.0, 4.0, 6.0, 8.0, 10.0]
    return pd.DataFrame({"model": [v / 2 for v in measured], "measure": measured})


@pytest.mark.unit
class TestPerfectModel:
    """Every metric has a defining value when model equals measurement."""

    def test_fractional_bias_is_zero(self, metrics, perfect):
        assert metrics.calculateFB(perfect) == pytest.approx(0.0)

    def test_normalised_mean_square_error_is_zero(self, metrics, perfect):
        assert metrics.calculateNMSE(perfect) == pytest.approx(0.0)

    def test_normalised_absolute_difference_is_zero(self, metrics, perfect):
        assert metrics.calculateNAD(perfect) == pytest.approx(0.0)

    @pytest.mark.xfail(
        strict=True,
        reason="B32: calculateR mixes normalisations -- the numerator uses .mean() "
               "(divide by N) while the denominator uses .std(), whose pandas "
               "default is the SAMPLE deviation (divide by N-1). Every result is "
               "short by exactly (N-1)/N, so a perfect model scores 0.8 at N=5 and "
               "0.667 at N=3. See the consolidated findings issue.",
    )
    def test_correlation_is_one(self, metrics, perfect):
        assert metrics.calculateR(perfect) == pytest.approx(1.0, rel=1e-6)

    def test_factor_of_two_is_one(self, metrics, perfect):
        assert metrics.calculateFAC(perfect) == pytest.approx(1.0)


@pytest.mark.unit
class TestFractionalBias:
    def test_the_sign_reports_over_prediction(self, metrics, over_predicting):
        """FB = 2 (mean_model - mean_measure) / (mean_model + mean_measure).

        Model = 2 x measure gives 2(2m - m)/(2m + m) = 2/3.
        """
        assert metrics.calculateFB(over_predicting) == pytest.approx(2 / 3)

    def test_the_sign_reports_under_prediction(self, metrics, under_predicting):
        """Model = measure/2 gives 2(0.5m - m)/(0.5m + m) = -2/3."""
        assert metrics.calculateFB(under_predicting) == pytest.approx(-2 / 3)

    def test_it_is_antisymmetric_in_the_two_columns(self, metrics, over_predicting):
        swapped = over_predicting.rename(
            columns={"model": "measure", "measure": "model"}
        )
        assert metrics.calculateFB(swapped) == pytest.approx(
            -metrics.calculateFB(over_predicting)
        )

    def test_it_stays_within_its_documented_bounds(self, metrics):
        """FB is bounded in [-2, 2] for non-negative data."""
        extreme = pd.DataFrame({"model": [100.0, 100.0], "measure": [1e-9, 1e-9]})
        assert -2.0 <= metrics.calculateFB(extreme) <= 2.0

    def test_a_scale_change_leaves_it_unchanged(self, metrics, over_predicting):
        """FB is dimensionless, so scaling both columns must not move it."""
        scaled = over_predicting * 1000.0
        assert metrics.calculateFB(scaled) == pytest.approx(
            metrics.calculateFB(over_predicting)
        )


@pytest.mark.unit
class TestNormalisedMeanSquareError:
    def test_it_grows_with_the_error(self, metrics, perfect):
        measured = perfect["measure"]
        small = pd.DataFrame({"model": measured * 1.1, "measure": measured})
        large = pd.DataFrame({"model": measured * 2.0, "measure": measured})
        assert metrics.calculateNMSE(large) > metrics.calculateNMSE(small)

    def test_it_is_never_negative(self, metrics, over_predicting, under_predicting):
        assert metrics.calculateNMSE(over_predicting) >= 0
        assert metrics.calculateNMSE(under_predicting) >= 0

    def test_it_is_symmetric_in_the_two_columns(self, metrics, over_predicting):
        """Both the squared error and the denominator product are symmetric."""
        swapped = over_predicting.rename(
            columns={"model": "measure", "measure": "model"}
        )
        assert metrics.calculateNMSE(swapped) == pytest.approx(
            metrics.calculateNMSE(over_predicting)
        )

    def test_a_scale_change_leaves_it_unchanged(self, metrics, over_predicting):
        scaled = over_predicting * 1000.0
        assert metrics.calculateNMSE(scaled) == pytest.approx(
            metrics.calculateNMSE(over_predicting)
        )


@pytest.mark.unit
class TestFactorOfTwo:
    def test_a_factor_of_exactly_two_falls_outside_the_band(self, metrics):
        """The comparison is strict: 1/2 < ratio < 2, so 2.0 does not count."""
        borderline = pd.DataFrame({"model": [2.0, 2.0], "measure": [1.0, 1.0]})
        assert metrics.calculateFAC(borderline) == pytest.approx(0.0)

    def test_just_inside_the_band_counts(self, metrics):
        inside = pd.DataFrame({"model": [1.99, 0.51], "measure": [1.0, 1.0]})
        assert metrics.calculateFAC(inside) == pytest.approx(1.0)

    def test_it_reports_the_fraction_that_qualifies(self, metrics):
        mixed = pd.DataFrame({"model": [1.0, 1.0, 10.0, 100.0], "measure": [1.0] * 4})
        assert metrics.calculateFAC(mixed) == pytest.approx(0.5)

    def test_the_band_width_is_configurable(self, metrics):
        """relation=10 accepts a factor of ten either way."""
        data = pd.DataFrame({"model": [5.0], "measure": [1.0]})
        assert metrics.calculateFAC(data, relation=2) == pytest.approx(0.0)
        assert metrics.calculateFAC(data, relation=10) == pytest.approx(1.0)

    def test_it_is_bounded_by_zero_and_one(self, metrics):
        awful = pd.DataFrame({"model": [1e6, 1e6], "measure": [1.0, 1.0]})
        assert 0.0 <= metrics.calculateFAC(awful) <= 1.0


@pytest.mark.unit
class TestCorrelation:
    @pytest.mark.xfail(
        strict=True,
        reason="B32: the same (N-1)/N shortfall, so a perfectly anticorrelated "
               "model scores -0.8 rather than -1. "
               "See the consolidated findings issue.",
    )
    def test_a_perfectly_anticorrelated_model_gives_minus_one(self, metrics):
        data = pd.DataFrame({"model": [5.0, 4.0, 3.0, 2.0, 1.0], "measure": [1.0, 2.0, 3.0, 4.0, 5.0]})
        assert metrics.calculateR(data) == pytest.approx(-1.0, rel=1e-6)

    @pytest.mark.parametrize("size", [3, 5, 10, 50, 200])
    def test_the_shortfall_is_exactly_the_sample_size_factor(self, metrics, size):
        """Pins the diagnosis of B32 so the cause is not guessed at later.

        A perfectly correlated series scores (N-1)/N instead of 1.  At N = 3
        that is a 33% error; the bias only vanishes asymptotically.  The fix is
        to make both halves use the same normalisation -- std(ddof=0).
        """
        values = np.arange(1.0, size + 1)
        data = pd.DataFrame({"model": values, "measure": values})
        assert metrics.calculateR(data) == pytest.approx((size - 1) / size, rel=1e-9)

    def test_an_offset_model_scores_the_same_as_an_unshifted_one(self, metrics, perfect):
        """R measures shape, not bias: a constant offset must not change it.

        Asserted as equality between the two rather than against 1.0, so this
        test stays meaningful whether or not B32 is fixed.
        """
        offset = pd.DataFrame(
            {"model": perfect["model"] + 100.0, "measure": perfect["measure"]}
        )
        assert metrics.calculateR(offset) == pytest.approx(metrics.calculateR(perfect))

    def test_it_agrees_with_numpy_on_ordinary_data(self, metrics):
        rng = np.random.default_rng(0)
        measured = rng.uniform(1, 10, 50)
        modelled = 1.3 * measured + rng.normal(0, 0.5, 50)
        data = pd.DataFrame({"model": modelled, "measure": measured})

        expected = np.corrcoef(modelled, measured)[0, 1]
        # times (N-1)/N, because of B32
        assert metrics.calculateR(data) == pytest.approx(expected * 49 / 50, rel=1e-6)


@pytest.mark.unit
class TestDirectionalBiasSplit:
    """FB_FP and FB_FN separate the over- and under-prediction halves of FB."""

    def test_pure_over_prediction_loads_one_side_only(self, metrics, over_predicting):
        positive = metrics.calculateFB_FP(over_predicting)
        negative = metrics.calculateFB_FN(over_predicting)
        assert positive == pytest.approx(0.0) or negative == pytest.approx(0.0)

    def test_the_two_halves_are_never_negative(self, metrics, over_predicting):
        assert metrics.calculateFB_FP(over_predicting) >= 0
        assert metrics.calculateFB_FN(over_predicting) >= 0

    def test_both_vanish_for_a_perfect_model(self, metrics, perfect):
        assert metrics.calculateFB_FP(perfect) == pytest.approx(0.0)
        assert metrics.calculateFB_FN(perfect) == pytest.approx(0.0)

    def test_their_difference_recovers_the_fractional_bias(self, metrics, over_predicting):
        """FB_FN - FB_FP is FB, since |d|+d and |d|-d split d by sign."""
        recovered = metrics.calculateFB_FN(over_predicting) - metrics.calculateFB_FP(
            over_predicting
        )
        assert recovered == pytest.approx(metrics.calculateFB(over_predicting), rel=1e-6)


@pytest.mark.unit
class TestColumnNames:
    def test_the_column_names_are_configurable(self, metrics):
        data = pd.DataFrame({"sim": [1.0, 2.0], "obs": [1.0, 2.0]})
        assert metrics.calculateFB(
            data, modelColumn="sim", measureColumn="obs"
        ) == pytest.approx(0.0)

    def test_a_missing_column_raises(self, metrics, perfect):
        with pytest.raises(KeyError):
            metrics.calculateFB(perfect, modelColumn="nope")
