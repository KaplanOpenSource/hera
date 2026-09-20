"""measurements/meteorology/lowfreqdata/analysis.py: the parts of the
low-frequency ``analysis`` class that batch-31's file left uncovered --
``datalayer``, ``calcHourlyDist``, ``_calculateCov`` and
``resampleSecondMoments``.  ``addDatesColumns`` and the season declarations
are already covered by test_meteorology_lowfreq_analysis.py.

Unlike its twin in ``measurements/meteorology/analysis.py`` (whose methods
are declared without ``self`` -- bug B104), this copy of the class is
declared correctly, so every method here is called the ordinary way on an
instance.

The assertions come from the documented statistics rather than from the
implementation:

* ``calcHourlyDist`` bins time-of-day against a field.  The hour axis is
  pinned to the full 24-hour day and the value axis is forced to include
  zero and to reach at least one, so the returned bin centres are
  arithmetic consequences of ``bins`` alone.
* ``resampleSecondMoments`` re-aggregates second moments from a fine
  averaging window onto a coarser one.  The formula it applies is the law
  of total covariance,

      cov_W(x,y) = mean_j cov_j(x,y) + mean_j(x_j y_j) - X_W Y_W

  which -- for equal-length sub-windows j -- is *exactly* the covariance
  computed in one pass over the coarse window W.  Every expected value
  below is that one-pass covariance, computed from the raw series and never
  from the code under test.
"""
import numpy
import pandas
import pytest

from hera.measurements.meteorology.lowfreqdata.analysis import analysis

BINS = 30
HOUR_SPAN = 24.0


@pytest.fixture()
def analyser():
    """None is a fine data layer: none of these methods reach for one."""
    return analysis(None)


def _daily_frame(values, hours, day="2020-06-15", tz=None):
    """A frame whose index carries the given hours of one day."""
    index = pandas.to_datetime([f"{day} {int(h):02d}:{int(round((h % 1) * 60)):02d}:00"
                                for h in hours])
    if tz is not None:
        index = index.tz_localize(tz)
    return pandas.DataFrame({"WS": list(values)}, index=index)


def _raw_wind(n=240, seed=3):
    """A 1 Hz, four-minute sonic-like record."""
    rng = numpy.random.default_rng(seed)
    index = pandas.date_range("2020-01-01 00:00:00", periods=n, freq="1s")
    return pandas.DataFrame({
        "u": 4.0 + rng.normal(0.0, 0.7, n) + numpy.linspace(0.0, 2.0, n),
        "w": rng.normal(0.0, 0.3, n) - numpy.linspace(0.0, 0.5, n),
    }, index=index)


def _fine_moments(raw, fine="30s"):
    """First and second moments over the fine averaging window.

    This is the shape ``resampleSecondMoments`` expects: one row per fine
    window, carrying ``u_bar``/``w_bar`` and the pre-combined second-moment
    columns ``uu``, ``uw``, ``ww``.
    """
    grouped = raw.resample(fine)
    means = grouped.mean()
    fluct = raw - means.reindex(raw.index, method="ffill")
    out = pandas.DataFrame({
        "u_bar": means["u"],
        "w_bar": means["w"],
        "uu": (fluct["u"] * fluct["u"]).resample(fine).mean(),
        "uw": (fluct["u"] * fluct["w"]).resample(fine).mean(),
        "ww": (fluct["w"] * fluct["w"]).resample(fine).mean(),
    })
    return out


def _one_pass_covariance(raw, coarse, x, y):
    """cov(x, y) computed in a single pass over each coarse window."""
    means = raw.resample(coarse).mean()
    fluct = raw - means.reindex(raw.index, method="ffill")
    return (fluct[x] * fluct[y]).resample(coarse).mean()


@pytest.mark.unit
class TestDatalayerProperty:
    def test_the_constructor_argument_is_exposed_as_the_datalayer(self):
        marker = object()
        assert analysis(marker).datalayer is marker

    def test_the_datalayer_defaults_to_nothing_only_when_none_is_passed(self):
        assert analysis(None).datalayer is None

    def test_two_analysers_do_not_share_a_datalayer(self):
        first, second = analysis("A"), analysis("B")
        assert (first.datalayer, second.datalayer) == ("A", "B")


@pytest.mark.unit
class TestCalcHourlyDistAxes:
    def test_the_hour_axis_always_spans_the_whole_day(self, analyser):
        data = _daily_frame([1.0, 2.0, 3.0], [1.0, 12.0, 23.0])
        x_mid, _, _ = analyser.calcHourlyDist(data, "WS", bins=BINS)
        width = HOUR_SPAN / BINS
        assert len(x_mid) == BINS
        assert x_mid[0] == pytest.approx(width / 2)
        assert x_mid[-1] == pytest.approx(HOUR_SPAN - width / 2)

    def test_the_hour_axis_is_independent_of_the_hours_present_in_the_data(self, analyser):
        early = _daily_frame([1.0, 2.0], [0.5, 1.5])
        late = _daily_frame([1.0, 2.0], [20.5, 21.5])
        assert numpy.allclose(analyser.calcHourlyDist(early, "WS")[0],
                              analyser.calcHourlyDist(late, "WS")[0])

    def test_the_value_axis_is_anchored_at_zero_for_an_all_positive_field(self, analyser):
        data = _daily_frame([4.0, 6.0, 10.0], [2.0, 8.0, 16.0])
        _, y_mid, _ = analyser.calcHourlyDist(data, "WS", bins=BINS)
        width = 10.0 / BINS
        assert y_mid[0] == pytest.approx(width / 2)
        assert y_mid[-1] == pytest.approx(10.0 - width / 2)

    def test_the_value_axis_reaches_at_least_one_for_a_near_zero_field(self, analyser):
        data = _daily_frame([0.05, 0.1, 0.15], [3.0, 9.0, 15.0])
        _, y_mid, _ = analyser.calcHourlyDist(data, "WS", bins=BINS)
        width = 1.0 / BINS
        assert y_mid[-1] == pytest.approx(1.0 - width / 2)

    def test_a_negative_field_value_pushes_the_value_axis_below_zero(self, analyser):
        data = _daily_frame([-3.0, 1.0, 2.0], [4.0, 10.0, 20.0])
        _, y_mid, _ = analyser.calcHourlyDist(data, "WS", bins=BINS)
        width = (2.0 - (-3.0)) / BINS
        assert y_mid[0] == pytest.approx(-3.0 + width / 2)


@pytest.mark.unit
class TestCalcHourlyDistBinning:
    def test_one_time_of_day_populates_exactly_one_hour_column(self, analyser):
        data = _daily_frame([2.0, 3.0, 4.0, 5.0], [6.5, 6.5, 6.5, 6.5])
        _, _, matrix = analyser.calcHourlyDist(data, "WS", bins=BINS)
        populated = numpy.flatnonzero(matrix.sum(axis=0) > 0)
        assert populated.tolist() == [int(6.5 / (HOUR_SPAN / BINS))]

    def test_the_matrix_has_one_row_per_value_bin_and_one_column_per_hour_bin(self, analyser):
        data = _daily_frame([1.0, 2.0, 3.0], [2.0, 8.0, 20.0])
        x_mid, y_mid, matrix = analyser.calcHourlyDist(data, "WS", bins=BINS)
        assert matrix.shape == (len(y_mid), len(x_mid))

    def test_max_normalisation_makes_the_fullest_bin_exactly_one(self, analyser):
        data = _daily_frame([1.0, 1.0, 1.0, 5.0], [6.0, 6.0, 6.0, 18.0])
        _, _, matrix = analyser.calcHourlyDist(data, "WS", normalization="max_normalized")
        assert matrix.max() == pytest.approx(1.0)

    def test_density_normalisation_conserves_the_sample_count(self, analyser):
        hours = [1.0, 4.0, 7.5, 13.0, 19.0, 22.5]
        data = _daily_frame([2.0, 3.0, 4.0, 5.0, 6.0, 7.0], hours)
        x_mid, y_mid, matrix = analyser.calcHourlyDist(data, "WS", bins=BINS,
                                                       normalization="density")
        cell = (x_mid[1] - x_mid[0]) * (y_mid[1] - y_mid[0])
        assert matrix.sum() * cell == pytest.approx(len(hours))

    def test_an_unknown_normalisation_is_rejected(self, analyser):
        data = _daily_frame([1.0, 2.0], [3.0, 4.0])
        with pytest.raises(ValueError, match="normalization"):
            analyser.calcHourlyDist(data, "WS", normalization="nonsense")


@pytest.mark.unit
class TestCalcHourlyDistCleaning:
    def _sample_count(self, analyser, data, **kwargs):
        x_mid, y_mid, matrix = analyser.calcHourlyDist(data, "WS", bins=BINS,
                                                       normalization="density", **kwargs)
        cell = (x_mid[1] - x_mid[0]) * (y_mid[1] - y_mid[0])
        return matrix.sum() * cell

    def test_a_missing_field_value_is_dropped(self, analyser):
        data = _daily_frame([2.0, numpy.nan, 4.0], [3.0, 9.0, 15.0])
        assert self._sample_count(analyser, data) == pytest.approx(2)

    def test_the_minus_9999_sentinel_is_dropped(self, analyser):
        data = _daily_frame([2.0, -9999.0, 4.0, -9999.0], [3.0, 6.0, 9.0, 12.0])
        assert self._sample_count(analyser, data) == pytest.approx(2)

    def test_the_sentinel_threshold_is_exclusive_at_minus_5000(self, analyser):
        kept = _daily_frame([2.0, -4999.0], [3.0, 9.0])
        dropped = _daily_frame([2.0, -5000.0], [3.0, 9.0])
        assert self._sample_count(analyser, kept) == pytest.approx(2)
        assert self._sample_count(analyser, dropped) == pytest.approx(1)

    def test_the_caller_s_frame_is_left_untouched(self, analyser):
        data = _daily_frame([2.0, -9999.0, 4.0], [3.0, 9.0, 15.0])
        before = data.copy(deep=True)
        analyser.calcHourlyDist(data, "WS")
        pandas.testing.assert_frame_equal(data, before)

    def test_a_timezone_aware_index_is_binned_in_utc(self, analyser):
        offset_hours = 2
        local = _daily_frame([3.0, 3.0, 3.0], [6.0, 6.0, 6.0], tz=f"Etc/GMT-{offset_hours}")
        _, _, matrix = analyser.calcHourlyDist(local, "WS", bins=BINS)
        populated = numpy.flatnonzero(matrix.sum(axis=0) > 0)
        assert populated.tolist() == [int((6.0 - offset_hours) / (HOUR_SPAN / BINS))]


@pytest.mark.unit
class TestCalcHourlyDistParquetPath:
    def test_a_path_string_is_read_as_parquet(self, analyser, tmp_path):
        data = _daily_frame([2.0, 4.0, 6.0], [3.0, 9.0, 15.0])
        target = tmp_path / "lowfreq.parquet"
        data.to_parquet(target)
        from_frame = analyser.calcHourlyDist(data, "WS", bins=BINS)
        from_path = analyser.calcHourlyDist(str(target), "WS", bins=BINS)
        for expected, actual in zip(from_frame, from_path):
            assert numpy.allclose(expected, actual)


@pytest.mark.unit
class TestResampleSecondMoments:
    def test_the_first_moments_are_resampled_onto_the_coarse_window(self, analyser):
        raw = _raw_wind()
        fine = _fine_moments(raw)
        out = analyser.resampleSecondMoments(fine, "60s", ["u_bar", "w_bar"], ["u", "w"])
        expected = raw.resample("60s").mean()
        assert numpy.allclose(out["u_bar"], expected["u"])
        assert numpy.allclose(out["w_bar"], expected["w"])

    @pytest.mark.parametrize("column, x, y", [("uu", "u", "u"),
                                              ("uw", "u", "w"),
                                              ("ww", "w", "w")])
    def test_each_second_moment_matches_the_one_pass_coarse_covariance(self, analyser,
                                                                       column, x, y):
        raw = _raw_wind()
        out = analyser.resampleSecondMoments(_fine_moments(raw), "60s",
                                             ["u_bar", "w_bar"], ["u", "w"])
        assert numpy.allclose(out[column], _one_pass_covariance(raw, "60s", x, y))

    def test_the_variances_it_produces_are_non_negative(self, analyser):
        out = analyser.resampleSecondMoments(_fine_moments(_raw_wind()), "60s",
                                             ["u_bar", "w_bar"], ["u", "w"])
        assert (out["uu"] >= 0).all()
        assert (out["ww"] >= 0).all()

    def test_the_cross_moment_obeys_the_cauchy_schwarz_bound(self, analyser):
        out = analyser.resampleSecondMoments(_fine_moments(_raw_wind()), "60s",
                                             ["u_bar", "w_bar"], ["u", "w"])
        assert (out["uw"].abs() <= numpy.sqrt(out["uu"] * out["ww"]) + 1e-12).all()

    def test_only_the_upper_triangle_of_the_moment_tensor_is_produced(self, analyser):
        out = analyser.resampleSecondMoments(_fine_moments(_raw_wind()), "60s",
                                             ["u_bar", "w_bar"], ["u", "w"])
        assert list(out.columns) == ["u_bar", "w_bar", "uu", "uw", "ww"]
        assert "wu" not in out.columns

    def test_resampling_onto_the_same_window_reproduces_the_input_moments(self, analyser):
        fine = _fine_moments(_raw_wind())
        out = analyser.resampleSecondMoments(fine, "30s", ["u_bar", "w_bar"], ["u", "w"])
        for column in ("uu", "uw", "ww"):
            assert numpy.allclose(out[column], fine[column])

    def test_a_missing_pre_combined_moment_column_is_reported(self, analyser):
        fine = _fine_moments(_raw_wind()).drop(columns=["uw"])
        with pytest.raises(KeyError, match="uw"):
            analyser.resampleSecondMoments(fine, "60s", ["u_bar", "w_bar"], ["u", "w"])

    def test_no_second_moment_fields_leaves_only_the_first_moments(self, analyser):
        fine = _fine_moments(_raw_wind())
        out = analyser.resampleSecondMoments(fine, "60s", ["u_bar", "w_bar"], [])
        assert list(out.columns) == ["u_bar", "w_bar"]


@pytest.mark.unit
class TestCalculateCov:
    def test_it_writes_the_moment_into_the_resampled_frame_in_place(self, analyser):
        raw = _raw_wind()
        fine = _fine_moments(raw)
        resampled = fine[["u_bar", "w_bar"]].resample("60s").mean()
        analyser._calculateCov(fine, resampled, "u", "w", "60s")
        assert "uw" in resampled.columns
        assert numpy.allclose(resampled["uw"], _one_pass_covariance(raw, "60s", "u", "w"))

    def test_it_returns_the_analyser_rather_than_the_frame(self, analyser):
        fine = _fine_moments(_raw_wind())
        resampled = fine[["u_bar", "w_bar"]].resample("60s").mean()
        assert analyser._calculateCov(fine, resampled, "u", "w", "60s") is analyser

    def test_the_moment_column_name_is_the_two_field_names_concatenated(self, analyser):
        fine = _fine_moments(_raw_wind())
        resampled = fine[["u_bar", "w_bar"]].resample("60s").mean()
        analyser._calculateCov(fine, resampled, "w", "w", "60s")
        assert "ww" in resampled.columns

    def test_it_needs_the_bar_columns_in_the_resampled_frame(self, analyser):
        fine = _fine_moments(_raw_wind())
        resampled = fine[["u_bar"]].resample("60s").mean()
        with pytest.raises(KeyError, match="w_bar"):
            analyser._calculateCov(fine, resampled, "u", "w", "60s")
