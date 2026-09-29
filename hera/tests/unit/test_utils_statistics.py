"""calcDist2d: a normalised 2-D histogram.

The fixture below is small enough to compute by hand, which is the only way an
assertion here can be a statement about correctness rather than a snapshot.

    x = [0.5, 0.5, 1.5]      bins = 2 over [0, 2)  ->  edges [0, 1, 2]
    y = [0.5, 1.5, 1.5]

    counts M[x_bin][y_bin] = [[1, 1],
                              [0, 1]]

    returned M.T             = [[1, 0],
                                [1, 1]]
"""
import numpy as np
import pandas as pd
import pytest

from hera.utils.statistics import calcDist2d

X = [0.5, 0.5, 1.5]
Y = [0.5, 1.5, 1.5]
RANGES = dict(bins=2, x_range=(0, 2), y_range=(0, 2))


@pytest.mark.unit
class TestBinGeometry:
    def test_midpoints_are_bin_centres(self):
        x_mid, y_mid, _ = calcDist2d(X, Y, **RANGES)
        assert x_mid == pytest.approx([0.5, 1.5])
        assert y_mid == pytest.approx([0.5, 1.5])

    def test_returned_matrix_is_transposed(self):
        """Shape is (n_y, n_x): the docstring returns M.T, not M."""
        _, _, matrix = calcDist2d(X, Y, bins=(2, 3), x_range=(0, 2), y_range=(0, 3))
        assert matrix.shape == (3, 2)

    def test_midpoint_count_matches_bin_count(self):
        x_mid, y_mid, matrix = calcDist2d(X, Y, **RANGES)
        assert len(x_mid) == matrix.shape[1]
        assert len(y_mid) == matrix.shape[0]


@pytest.mark.unit
class TestNormalizations:
    def test_max_normalized_peaks_at_one(self):
        _, _, matrix = calcDist2d(X, Y, normalization="max_normalized", **RANGES)
        assert matrix.max() == pytest.approx(1.0)
        assert matrix == pytest.approx(np.array([[1.0, 0.0], [1.0, 1.0]]))

    def test_max_normalized_is_the_default(self):
        _, _, default = calcDist2d(X, Y, **RANGES)
        _, _, explicit = calcDist2d(X, Y, normalization="max_normalized", **RANGES)
        assert default == pytest.approx(explicit)

    def test_density_divides_by_bin_area(self):
        """Unit bins here, so density must equal the raw counts."""
        _, _, matrix = calcDist2d(X, Y, normalization="density", **RANGES)
        assert matrix == pytest.approx(np.array([[1.0, 0.0], [1.0, 1.0]]))

    def test_density_scales_with_bin_area(self):
        """Halving the bin width quarters the area, so density must quadruple."""
        _, _, coarse = calcDist2d(X, Y, normalization="density", bins=1,
                                  x_range=(0, 2), y_range=(0, 2))
        _, _, fine = calcDist2d(X, Y, normalization="density", bins=2,
                                x_range=(0, 2), y_range=(0, 2))
        assert coarse.sum() == pytest.approx(3 / 4.0)
        assert fine.sum() == pytest.approx(3 / 1.0)

    def test_y_normalized_makes_every_populated_row_sum_to_one(self):
        """Documented as 'normalize each row to sum to 1' -- rows of the result."""
        _, _, matrix = calcDist2d(X, Y, normalization="y_normalized", **RANGES)
        assert matrix == pytest.approx(np.array([[1.0, 0.0], [0.5, 0.5]]))
        for row in matrix:
            assert row.sum() == pytest.approx(1.0)

    def test_y_normalized_leaves_empty_rows_at_zero(self):
        """Division by a zero column sum must not produce NaN or inf."""
        _, _, matrix = calcDist2d([0.5], [0.5], normalization="y_normalized", **RANGES)
        assert np.isfinite(matrix).all()
        assert matrix[1].sum() == pytest.approx(0.0)

    def test_unknown_normalization_raises(self):
        with pytest.raises(ValueError, match="normalization must be one of"):
            calcDist2d(X, Y, normalization="whatever", **RANGES)


@pytest.mark.unit
class TestDataFrameInput:
    def test_column_names_resolve_against_data(self):
        frame = pd.DataFrame({"a": X, "b": Y})
        _, _, from_frame = calcDist2d("a", "b", data=frame, **RANGES)
        _, _, from_arrays = calcDist2d(X, Y, **RANGES)
        assert from_frame == pytest.approx(from_arrays)

    def test_missing_column_raises(self):
        frame = pd.DataFrame({"a": X, "b": Y})
        with pytest.raises(KeyError):
            calcDist2d("a", "nope", data=frame, **RANGES)


@pytest.mark.unit
class TestRangeHandling:
    def test_range_is_ignored_unless_both_axes_are_given(self):
        """The code sets a range only when x_range AND y_range are present.

        Pinned because the asymmetry is easy to trip over: passing one alone
        silently falls back to the data extent.
        """
        _, _, only_x = calcDist2d(X, Y, bins=2, x_range=(0, 10))
        _, _, neither = calcDist2d(X, Y, bins=2)
        assert only_x == pytest.approx(neither)

    def test_explicit_range_widens_the_bins(self):
        x_mid, _, _ = calcDist2d(X, Y, bins=2, x_range=(0, 4), y_range=(0, 4))
        assert x_mid == pytest.approx([1.0, 3.0])
