r"""``riskassessment/docs/python/HaberToBerge.py``: the helper that converts
a Haber-law dose-response curve into a ten Berge one.

Haber's law integrates concentration linearly in time; ten Berge integrates
``C^n``.  Thresholds calibrated as two-minute Haber dosages therefore have
to be re-expressed before they can be used with a ten Berge calculator, and
the module's docstring gives the transform:

    TB-threshold = 2 * (Dosage / 2)^n

``EstimateTBProbit`` wraps that transform and adds ``plotGuess``, a visual
check that overlays the *correct* transformed curve on a *guessed*
log-normal fit, printing the guessed LD50 and probit alongside.  The class
is a working aid (its module body is a block of commented-out invocations
for the agents in the docs), so the tests below establish the transform
algebraically and assert that the plot draws the two curves it promises.

The transform's properties are derived from the formula, not from hera:

* it fixes the calibration point: ``transform(2) == 2`` for every ``n``;
* ``n == 1`` is the identity (Haber is ten Berge with n = 1);
* ``n > 1`` stretches dosages above 2 and compresses those below it.

``plotGuess`` calls ``plt.show()``, which the conftest's Agg backend turns
into a warning-and-no-op, and the assertions are on the artists it leaves
on the current axes.

Deliberately not covered here: the commented-out ``GBestimator.plotGuess``
calls in the module body -- they are documentation of the agents' published
LD50/probit pairs, not code.  Nothing verifies that the *guess* is a good
one; that is the human judgement the plot exists to support.
"""
import numpy
import pytest

from hera.riskassessment.docs.python.HaberToBerge import (GBestimator,
                                                          EstimateTBProbit)


@pytest.mark.unit
class TestConstruction:
    def test_the_exponent_is_stored(self):
        assert EstimateTBProbit(1.5).n == 1.5

    def test_the_module_level_estimator_uses_the_documented_exponent(self):
        """The shared instance the commented-out examples in the module body
        call into."""
        assert isinstance(GBestimator, EstimateTBProbit)
        assert GBestimator.n == 1.5


@pytest.mark.unit
class TestTransform:
    def test_the_calibration_dosage_is_a_fixed_point(self):
        """2*(2/2)^n == 2 for any n: the two-minute reference dosage is
        unchanged by construction."""
        for exponent in (0.5, 1.0, 1.5, 2.0, 8.0):
            assert EstimateTBProbit(exponent).transform(2.0) == pytest.approx(2.0)

    def test_an_exponent_of_one_is_the_identity(self):
        """Haber's law is ten Berge with n = 1."""
        identity = EstimateTBProbit(1.0)

        for dosage in (0.5, 1.0, 2.0, 18.0, 99.0):
            assert identity.transform(dosage) == pytest.approx(dosage)

    def test_it_matches_the_documented_formula(self):
        estimator = EstimateTBProbit(1.5)

        for dosage in (0.5, 1.0, 8.0, 18.0):
            assert estimator.transform(dosage) == pytest.approx(
                2.0 * (dosage / 2.0) ** 1.5)

    def test_doubling_the_reference_dosage_scales_by_two_to_the_n(self):
        """transform(4)/transform(2) == 2^n."""
        estimator = EstimateTBProbit(1.5)

        assert estimator.transform(4.0) / estimator.transform(2.0) == pytest.approx(
            2.0 ** 1.5)

    def test_an_exponent_above_one_stretches_dosages_above_the_reference(self):
        estimator = EstimateTBProbit(1.5)

        assert estimator.transform(18.0) > 18.0

    def test_an_exponent_above_one_compresses_dosages_below_the_reference(self):
        estimator = EstimateTBProbit(1.5)

        assert estimator.transform(0.5) < 0.5

    def test_it_is_monotonically_increasing(self):
        estimator = EstimateTBProbit(1.5)
        dosages = numpy.array([0.1, 0.5, 1.0, 2.0, 8.0, 18.0, 99.0])

        assert numpy.all(numpy.diff(estimator.transform(dosages)) > 0)

    def test_an_array_of_dosages_is_transformed_elementwise(self):
        estimator = EstimateTBProbit(1.5)
        dosages = numpy.array([2.0, 8.0])

        assert estimator.transform(dosages) == pytest.approx(
            [estimator.transform(2.0), estimator.transform(8.0)])

    def test_a_zero_dosage_transforms_to_zero(self):
        assert EstimateTBProbit(1.5).transform(0.0) == pytest.approx(0.0)


@pytest.mark.unit
class TestPlotGuess:
    """Documented as plotting "the transformed data (which is the correct) vs.
    the guess of the lognormal parameter"."""

    NEW_PROBIT = 2.6666
    E50 = 18.0
    PROBIT = 4.0

    @pytest.fixture()
    def drawn(self):
        import matplotlib.pyplot as plt

        estimator = EstimateTBProbit(1.5)
        with pytest.warns(UserWarning, match="non-interactive"):
            returned = estimator.plotGuess(self.NEW_PROBIT, self.E50, self.PROBIT)
        return returned, plt.gca(), estimator

    def test_it_returns_nothing(self, drawn):
        returned, _, _ = drawn
        assert returned is None

    def test_it_draws_both_curves(self, drawn):
        _, axes, _ = drawn
        assert len(axes.lines) == 2

    def test_the_curves_are_labelled_correct_and_guess(self, drawn):
        _, axes, _ = drawn
        assert [line.get_label() for line in axes.lines] == [
            "The tenberge dose response", "Guess"]

    def test_the_dose_axis_is_logarithmic(self, drawn):
        """The dosages span three decades, hence semilogx."""
        _, axes, _ = drawn
        assert axes.get_xscale() == "log"

    def test_both_curves_share_the_transformed_dose_axis(self, drawn):
        """Both are plotted against the *transformed* dose list, which is what
        makes the two comparable at a glance."""
        _, axes, estimator = drawn

        expected = estimator.transform(numpy.arange(1, 100, 0.01))
        for line in axes.lines:
            assert line.get_xdata() == pytest.approx(expected)

    def test_both_curves_are_probabilities(self, drawn):
        _, axes, _ = drawn

        for line in axes.lines:
            values = line.get_ydata()
            assert values.min() >= 0.0
            assert values.max() <= 1.0

    def test_both_curves_are_monotonically_increasing(self, drawn):
        """A dose-response curve cannot fall as the dose rises."""
        _, axes, _ = drawn

        for line in axes.lines:
            assert numpy.all(numpy.diff(line.get_ydata()) >= 0)

    def test_a_legend_is_added(self, drawn):
        _, axes, _ = drawn
        assert axes.get_legend() is not None

    def test_the_transformed_ld50_is_printed(self, capsys):
        estimator = EstimateTBProbit(1.5)
        with pytest.warns(UserWarning):
            estimator.plotGuess(self.NEW_PROBIT, self.E50, self.PROBIT)

        printed = capsys.readouterr().out
        assert "New LD50 %s" % estimator.transform(self.E50) in printed

    def test_the_transformed_probit_is_printed(self, capsys):
        """Documented in the module comments: the probit divides by n and the
        geometric standard deviation multiplies by it."""
        estimator = EstimateTBProbit(1.5)
        with pytest.warns(UserWarning):
            estimator.plotGuess(self.NEW_PROBIT, self.E50, self.PROBIT)

        printed = capsys.readouterr().out
        assert "New gStd %s probit %s" % (self.PROBIT / estimator.n,
                                          estimator.n / self.PROBIT) in printed

    def test_the_correct_curve_passes_through_one_half_at_the_transformed_ld50(self):
        """The whole point of the transform: the median stays the median.  The
        first curve is the exact transformation, so it must cross 0.5 at
        transform(E50)."""
        import matplotlib.pyplot as plt

        estimator = EstimateTBProbit(1.5)
        with pytest.warns(UserWarning):
            estimator.plotGuess(self.NEW_PROBIT, self.E50, self.PROBIT)

        correct = plt.gca().lines[0]
        crossing = numpy.interp(0.5, correct.get_ydata(), correct.get_xdata())
        assert crossing == pytest.approx(estimator.transform(self.E50), rel=1e-3)

    def test_an_exponent_of_one_makes_the_guess_agree_with_the_correct_curve(self):
        """With n = 1 nothing is transformed, so passing the unchanged probit
        as the guess has to reproduce the same curve -- a self-consistency
        check on the two lognorm.cdf calls."""
        import matplotlib.pyplot as plt

        estimator = EstimateTBProbit(1.0)
        with pytest.warns(UserWarning):
            estimator.plotGuess(self.PROBIT, self.E50, self.PROBIT)

        correct, guess = plt.gca().lines
        assert guess.get_ydata() == pytest.approx(correct.get_ydata())
