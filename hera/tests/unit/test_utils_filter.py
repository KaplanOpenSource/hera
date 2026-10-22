"""Filter: threshold and interval filtering over a pandas DataFrame."""
import pandas as pd
import pytest

from hera.utils.filter_immediate import Filter


@pytest.fixture()
def frame():
    """Six rows, index 0..5, with a signed column to exercise the abs_* branches."""
    return pd.DataFrame({"v": [-3.0, -1.0, 0.0, 1.0, 3.0, 5.0]})


@pytest.mark.unit
class TestThresholdOperators:
    """The ten documented prepositions, checked against the values they keep."""

    @pytest.mark.parametrize(
        "preposition, bound, expected",
        [
            ("lt", 1.0, [-3.0, -1.0, 0.0]),
            ("lte", 1.0, [-3.0, -1.0, 0.0, 1.0]),
            ("gt", 1.0, [3.0, 5.0]),
            ("gte", 1.0, [1.0, 3.0, 5.0]),
            ("eq", 3.0, [3.0]),
            ("neq", 3.0, [-3.0, -1.0, 0.0, 1.0, 5.0]),
            ("abs_lt", 1.0, [0.0]),
            ("abs_lte", 1.0, [-1.0, 0.0, 1.0]),
            ("abs_gt", 3.0, [5.0]),
            ("abs_gte", 3.0, [-3.0, 3.0, 5.0]),
        ],
    )
    def test_operator_keeps_the_right_rows(self, frame, preposition, bound, expected):
        result = Filter(frame).threshold(preposition, bound, column="v")
        assert result.data["v"].tolist() == expected

    def test_unknown_preposition_raises(self, frame):
        with pytest.raises(ValueError, match="preposition must be one of"):
            Filter(frame).threshold("approximately", 1.0, column="v")

    def test_column_none_filters_on_the_index(self, frame):
        """Documented: column=None uses the index instead of a column."""
        result = Filter(frame).threshold("lt", 2, column=None)
        assert result.data.index.tolist() == [0, 1]

    def test_no_row_matches_gives_an_empty_frame(self, frame):
        result = Filter(frame).threshold("gt", 1000.0, column="v")
        assert len(result.data) == 0
        assert list(result.data.columns) == ["v"]


@pytest.mark.unit
class TestInplaceSemantics:
    """inplace=False must leave the original untouched; True must mutate."""

    def test_default_returns_a_new_filter_and_leaves_the_source_alone(self, frame):
        original = Filter(frame)
        returned = original.threshold("gt", 0.0, column="v")

        assert returned is not original
        assert len(original.data) == 6, "the source Filter was mutated"
        assert len(returned.data) == 3

    def test_inplace_mutates_and_returns_self(self, frame):
        original = Filter(frame, inplace=True)
        returned = original.threshold("gt", 0.0, column="v")

        assert returned is original
        assert len(original.data) == 3

    def test_inplace_flag_is_carried_into_the_new_filter(self, frame):
        """A chained call must keep behaving the way the first one did."""
        chained = Filter(frame).threshold("gt", -2.0, column="v")
        assert chained.inplace is False

    def test_chaining_narrows_progressively(self, frame):
        result = (
            Filter(frame)
            .threshold("gt", -2.0, column="v")
            .threshold("lt", 4.0, column="v")
        )
        assert result.data["v"].tolist() == [-1.0, 0.0, 1.0, 3.0]


@pytest.mark.unit
class TestOutsideInterval:
    """Keeps rows outside [lower, upper) -- lower inclusive, upper exclusive."""

    def test_keeps_values_beyond_both_ends(self, frame):
        result = Filter(frame).outsideInterval(-1.0, 3.0, column="v")
        assert result.data["v"].tolist() == [-3.0, 3.0, 5.0]

    def test_lower_bound_is_excluded_from_the_result(self, frame):
        """value == lower is INSIDE the interval, so it is dropped."""
        result = Filter(frame).outsideInterval(0.0, 3.0, column="v")
        assert 0.0 not in result.data["v"].tolist()

    def test_upper_bound_is_kept(self, frame):
        """value == upper is OUTSIDE the interval, so it survives."""
        result = Filter(frame).outsideInterval(-1.0, 3.0, column="v")
        assert 3.0 in result.data["v"].tolist()

    def test_interval_covering_everything_empties_the_frame(self, frame):
        result = Filter(frame).outsideInterval(-100.0, 100.0, column="v")
        assert len(result.data) == 0

    def test_inverted_interval_keeps_everything(self, frame):
        """lower > upper makes both conditions true for every row."""
        result = Filter(frame).outsideInterval(10.0, -10.0, column="v")
        assert len(result.data) == 6

    def test_column_none_filters_on_the_index(self, frame):
        result = Filter(frame).outsideInterval(1, 4, column=None)
        assert result.data.index.tolist() == [0, 4, 5]

    def test_respects_inplace(self, frame):
        original = Filter(frame, inplace=True)
        returned = original.outsideInterval(-1.0, 3.0, column="v")
        assert returned is original
        assert original.data["v"].tolist() == [-3.0, 3.0, 5.0]
