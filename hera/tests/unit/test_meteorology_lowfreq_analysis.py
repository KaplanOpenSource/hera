"""Low-frequency meteorology: date decomposition and seasonal classification.

The seasons are the standard meteorological ones, and the module declares them
itself in `seasonsdict`: DJF winter, MAM spring, JJA summer, SON autumn.  The
tests check every month against that declaration rather than against the
pd.cut bin edges, so a change to either one has to be deliberate.
"""
import numpy as np
import pandas as pd
import pytest

from hera.measurements.meteorology.lowfreqdata.analysis import (
    AUTUMN,
    SPRING,
    SUMMER,
    WINTER,
    analysis,
    seasonsdict,
)

MONTH_TO_SEASON = {
    1: WINTER, 2: WINTER, 12: WINTER,
    3: SPRING, 4: SPRING, 5: SPRING,
    6: SUMMER, 7: SUMMER, 8: SUMMER,
    9: AUTUMN, 10: AUTUMN, 11: AUTUMN,
}


@pytest.fixture()
def analyser():
    """The date helpers never touch the data layer."""
    return analysis(None)


def one_row_per_month(year=2020):
    """A frame with one timestamped row in each month of a year."""
    index = pd.to_datetime([f"{year}-{month:02d}-15 06:30:00" for month in range(1, 13)])
    return pd.DataFrame({"value": np.arange(12, dtype=float)}, index=index)


@pytest.mark.unit
class TestSeasonDeclaration:
    """seasonsdict is the module's own statement of the convention."""

    def test_every_month_appears_exactly_once(self):
        months = [m for season in seasonsdict.values() for m in season["months"]]
        assert sorted(months) == list(range(1, 13))

    @pytest.mark.parametrize("season, months", [
        (WINTER, [12, 1, 2]), (SPRING, [3, 4, 5]),
        (SUMMER, [6, 7, 8]), (AUTUMN, [9, 10, 11]),
    ])
    def test_the_declared_months_are_the_meteorological_ones(self, season, months):
        assert seasonsdict[season]["months"] == months

    @pytest.mark.parametrize("season, code", [
        (WINTER, "[DJF]"), (SPRING, "[MAM]"), (SUMMER, "[JJA]"), (AUTUMN, "[SON]"),
    ])
    def test_the_short_codes_match_the_months(self, season, code):
        assert seasonsdict[season]["strmonths"] == code


@pytest.mark.unit
class TestSeasonAssignment:
    @pytest.mark.parametrize("month", list(range(1, 13)))
    def test_each_month_lands_in_its_declared_season(self, analyser, month):
        """Checked against seasonsdict, not against the pd.cut edges."""
        frame = one_row_per_month()
        decorated = analyser.addDatesColumns(frame)
        assigned = decorated[decorated["monthonly"] == month]["season"].iloc[0]
        assert assigned == MONTH_TO_SEASON[month]

    def test_december_and_january_share_a_season(self, analyser):
        """The wrap-around case: pd.cut needs a duplicate Winter label for it."""
        decorated = analyser.addDatesColumns(one_row_per_month())
        december = decorated[decorated["monthonly"] == 12]["season"].iloc[0]
        january = decorated[decorated["monthonly"] == 1]["season"].iloc[0]
        assert december == january == WINTER

    def test_all_four_seasons_are_produced(self, analyser):
        decorated = analyser.addDatesColumns(one_row_per_month())
        assert set(decorated["season"].astype(str)) == {WINTER, SPRING, SUMMER, AUTUMN}

    def test_each_season_gets_three_months(self, analyser):
        decorated = analyser.addDatesColumns(one_row_per_month())
        counts = decorated["season"].astype(str).value_counts().to_dict()
        assert all(count == 3 for count in counts.values())


@pytest.mark.unit
class TestDateColumns:
    def test_the_year_is_extracted(self, analyser):
        decorated = analyser.addDatesColumns(one_row_per_month(year=1998))
        assert set(decorated["yearonly"]) == {1998}

    def test_the_day_is_extracted(self, analyser):
        decorated = analyser.addDatesColumns(one_row_per_month())
        assert set(decorated["dayonly"]) == {15}

    def test_the_numeric_time_column_is_hhmm(self, analyser):
        """06:30 becomes 630, so hours and minutes sort as one number."""
        decorated = analyser.addDatesColumns(one_row_per_month())
        assert set(decorated["Time"]) == {630}

    @pytest.mark.parametrize("hour, minute, expected", [
        (0, 0, 0), (0, 5, 5), (9, 45, 945), (13, 0, 1300), (23, 59, 2359),
    ])
    def test_the_hhmm_encoding_across_the_day(self, analyser, hour, minute, expected):
        index = pd.to_datetime([f"2020-06-15 {hour:02d}:{minute:02d}:00"])
        frame = pd.DataFrame({"value": [1.0]}, index=index)
        assert analyser.addDatesColumns(frame)["Time"].iloc[0] == expected

    def test_the_original_columns_survive(self, analyser):
        decorated = analyser.addDatesColumns(one_row_per_month())
        assert "value" in decorated.columns
        assert len(decorated) == 12

    def test_the_input_frame_is_not_mutated(self, analyser):
        """assign() returns a copy; a caller's frame must come back untouched."""
        frame = one_row_per_month()
        before = list(frame.columns)
        analyser.addDatesColumns(frame)
        assert list(frame.columns) == before


@pytest.mark.unit
class TestColumnSelection:
    def test_an_explicit_date_column_is_used(self, analyser):
        frame = pd.DataFrame({
            "when": pd.to_datetime(["2020-07-04 12:00:00"]),
            "value": [1.0],
        })
        decorated = analyser.addDatesColumns(frame, datecolumn="when")
        assert decorated["yearonly"].iloc[0] == 2020
        assert decorated["season"].iloc[0] == SUMMER

    def test_an_explicit_month_column_drives_the_season(self, analyser):
        """The month column can disagree with the date; the season follows it."""
        frame = pd.DataFrame({
            "when": pd.to_datetime(["2020-07-04 12:00:00"]),
            "mon": [1],
            "value": [1.0],
        })
        decorated = analyser.addDatesColumns(frame, datecolumn="when", monthcolumn="mon")
        assert decorated["season"].iloc[0] == WINTER

    def test_a_timezone_aware_index_is_made_naive(self, analyser):
        """Documented behaviour: tz_convert(None) before the field extraction."""
        index = pd.to_datetime(["2020-07-04 12:00:00"]).tz_localize("UTC")
        frame = pd.DataFrame({"value": [1.0]}, index=index)
        decorated = analyser.addDatesColumns(frame)
        assert decorated["yearonly"].iloc[0] == 2020
        assert decorated["Time"].iloc[0] == 1200


@pytest.mark.unit
class TestParquetPath:
    def test_a_path_is_read_from_disk(self, analyser, tmp_path):
        """Documented: a string argument is treated as a parquet path."""
        target = tmp_path / "data.parquet"
        one_row_per_month().to_parquet(target)

        decorated = analyser.addDatesColumns(str(target))
        assert len(decorated) == 12
        assert set(decorated["season"].astype(str)) == {
            WINTER, SPRING, SUMMER, AUTUMN,
        }

    def test_a_missing_path_raises(self, analyser, tmp_path):
        with pytest.raises(Exception):
            analyser.addDatesColumns(str(tmp_path / "absent.parquet"))
