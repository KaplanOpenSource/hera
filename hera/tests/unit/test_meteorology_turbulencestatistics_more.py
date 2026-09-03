"""turbulencestatistics.py: the structure-function family, the Spark
``fluctuations`` override, and ``InMemoryRawData``'s HDF5 round trip -- the
parts test_meteorology_turbulencestatistics.py left uncovered.

Skips and deliberate limits, with reasons
-----------------------------------------
``StrucFunDir``/``StrucFun``/``ThirdStrucFun`` are exercised here for their
*contract only* -- what they return, which columns they declare, what they
register in ``_CalculatedParams``, how ``StrucFun`` decomposes the mean-wind
basis -- and never by materialising a structure function.  The reason is
measured, not assumed: the tau-lagged self-join these methods build
(``united_data.merge(data_tau).repartition(freq="1W")`` fed into
``.resample(...)``, on top of six chained ``ui += ...``/``uj += ...``
assignments) makes dask 2024.8.0's expression optimiser blow up.  Computing
one tau of one component took **208 s for a 120-sample series and 284 s for
an 8-sample one** on this environment -- the cost is in the expression tree,
not the data, so no amount of shrinking the fixture helps.  Dropping either
accumulator loop brings the same graph back to 0.2 s, which is what
identifies the cause.  A numeric structure-function test therefore belongs
in the integration suite, and the physics of the estimator is not asserted
here.

These methods are also dask-only by construction (``.repartition(freq=...)``
and an unconditional ``.compute()``), which is B305 below.

``SinglePointStatisticsSpark.fluctuations`` needs no pyspark at all -- the
module imports none -- so it is covered in full against a pandas frame.

``InMemoryRawData.to_hdf``/``read_hdf`` need PyTables for the real HDF5
round trip and PyTables is not installed here (``h5py`` is, but pandas' HDF
support is PyTables-only), so the one end-to-end round-trip test carries
``pytest.importorskip("tables")``.  Everything these two methods actually
own -- the extension rewriting and the JSON sidecar -- is exercised without
it, by patching ``pandas.read_hdf`` and ``pandas.DataFrame.to_hdf``.

Defects surfaced
----------------
* B303: ``StrucFun`` rejects an unknown ``mode`` with
  ``raise("mode must be either MeanDir or 3dMeanDir")`` -- raising the
  *string*, not an exception.  Python replaces it with
  ``TypeError: exceptions must derive from BaseException``, so the caller
  never sees the diagnostic that was written for them.
* B304: ``StrucFun`` and ``ThirdStrucFun`` guard the ``u_mag`` column
  with ``if "u_mag" not in self.TemporaryData.columns`` but then write and
  register ``"u_mag" + title_additions``.  With a non-empty
  ``title_additions`` the guard inspects a name that is never created, so
  every repeated call appends another ``["u_mag<additions>", {}]`` to
  ``_CalculatedParams`` -- and ``AbstractCalculator._compute`` selects
  ``_TemporaryData`` by that list, so the duplicate propagates into the
  computed result.
* B305: the whole structure-function family is dask-only although the
  class documents ``rawData`` as "pandas.DataFrame or
  dask.dataframe.DataFrame" and every other method honours that.
  ``StrucFunDir`` calls ``.repartition(freq=...)`` and ``StrucFun``/
  ``ThirdStrucFun`` call ``.compute()``, none of which pandas has, so a
  pandas-backed calculator dies with ``AttributeError``.
* B306: ``SinglePointStatisticsSpark.fluctuations`` never sets
  ``self.data``, which the base implementation does as its last act.
  ``getData()`` therefore still answers ``None`` after the override has
  computed everything -- and its own docstring says ``None`` means "no
  computations have been performed yet".
* B307: the same override puts ``wind_dir_bar`` into ``_TemporaryData``
  but -- unlike the base implementation -- leaves it out of
  ``_CalculatedParams``.  Since ``compute()`` projects ``_TemporaryData``
  onto that list, the mean wind direction is silently dropped from the
  result.
* B308: ``InMemoryRawData.to_hdf`` derives both output paths with
  ``path_or_buf.rpartition('.')[0]``.  For a path with no extension that
  prefix is the empty string, so a request to write
  ``/data/station/rawdata`` writes ``.hdf`` and ``.json`` into the *current
  working directory* instead of into ``/data/station``.
"""
import json
import os

import dask.dataframe
import numpy
import pandas
import pytest

from hera.measurements.meteorology.highfreqdata.analysis.turbulencestatistics import (
    InMemoryRawData,
    SinglePointStatisticsSpark,
    singlePointTurbulenceStatistics,
)

WINDOW = "30s"
NORTH, EAST, SOUTH, WEST = 0.0, 90.0, 180.0, 270.0


def _raw(n=120, freq="1s", start="2020-01-01 00:00:00", seed=0):
    rng = numpy.random.default_rng(seed)
    index = pandas.date_range(start, periods=n, freq=freq, name="Time")
    return pandas.DataFrame({
        "u": 3.0 + rng.normal(0, 0.5, n),
        "v": 1.0 + rng.normal(0, 0.5, n),
        "w": rng.normal(0, 0.2, n),
        "T": 20.0 + rng.normal(0, 0.1, n),
    }, index=index)


def _metadata(raw, window=WINDOW, **overrides):
    meta = dict(isMissingData=False, samplingWindow=window,
                start=raw.index[0], end=raw.index[-1] + pandas.Timedelta("1s"),
                height=10.0, buildingHeight=5.0, averagedHeight=3.0)
    meta.update(overrides)
    return meta


@pytest.fixture()
def lazy():
    """A dask-backed calculator: the structure functions accept nothing else."""
    raw = _raw()
    return singlePointTurbulenceStatistics(
        rawData=dask.dataframe.from_pandas(raw, npartitions=1),
        metadata=_metadata(raw))


@pytest.fixture()
def eager():
    """A pandas-backed calculator, for the dask-only defect."""
    raw = _raw()
    return singlePointTurbulenceStatistics(rawData=raw, metadata=_metadata(raw))


@pytest.fixture()
def direction_frame(lazy):
    """Direction vectors on the sampling-window grid, deliberately unnormalised."""
    index = pandas.date_range("2020-01-01 00:00:00", periods=4, freq=WINDOW, name="Time")
    return pandas.DataFrame({"a": [2.0, 2.0, 0.0, 0.0],
                             "b": [0.0, 0.0, 3.0, 3.0],
                             "c": [0.0, 0.0, 0.0, 0.0]}, index=index)


def _mean_wind(index=None):
    """An externally supplied mean-velocity frame with a tilted 3-D direction."""
    if index is None:
        index = pandas.date_range("2020-01-01 00:00:00", periods=4, freq=WINDOW,
                                  name="Time")
    return pandas.DataFrame({"u_bar": [3.0, 4.0, 0.0, -2.0],
                             "v_bar": [1.0, 0.0, 5.0, 2.0],
                             "w_bar": [0.5, -0.5, 0.0, 1.0]}, index=index)


def _param_names(calculator):
    return [entry[0] for entry in calculator._CalculatedParams]


# ---------------------------------------------------------------------------
# StrucFunDir -- contract only, see the module docstring
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestStrucFunDirContract:
    def test_it_returns_the_calculator_for_chaining(self, lazy, direction_frame):
        result = lazy.StrucFunDir(tau_range=[1.0], dir1_data=direction_frame,
                                  u_dir1="a", v_dir1="b", w_dir1="c")
        assert result is lazy

    def test_the_column_is_named_after_the_title_and_the_lag(self, lazy, direction_frame):
        lazy.StrucFunDir(tau_range=[1.0], dir1_data=direction_frame,
                         u_dir1="a", v_dir1="b", w_dir1="c", title="11")
        assert "D11_1.0s" in lazy.TemporaryData.columns

    def test_an_integer_lag_keeps_its_integer_spelling(self, lazy, direction_frame):
        lazy.StrucFunDir(tau_range=[2], dir1_data=direction_frame,
                         u_dir1="a", v_dir1="b", w_dir1="c")
        assert "D_2s" in lazy.TemporaryData.columns

    def test_one_column_per_requested_lag(self, lazy, direction_frame):
        lazy.StrucFunDir(tau_range=[1.0, 2.0, 5.0], dir1_data=direction_frame,
                         u_dir1="a", v_dir1="b", w_dir1="c", title="11")
        for tau in (1.0, 2.0, 5.0):
            assert "D11_%ss" % tau in lazy.TemporaryData.columns

    def test_each_lag_is_registered_as_a_calculated_parameter(self, lazy,
                                                              direction_frame):
        lazy.StrucFunDir(tau_range=[1.0, 2.0], dir1_data=direction_frame,
                         u_dir1="a", v_dir1="b", w_dir1="c", title="11")
        assert _param_names(lazy)[-2:] == ["D11_1.0s", "D11_2.0s"]

    def test_the_mean_fields_are_computed_first(self, lazy, direction_frame):
        lazy.StrucFunDir(tau_range=[1.0], dir1_data=direction_frame,
                         u_dir1="a", v_dir1="b", w_dir1="c")
        assert "u_bar" in lazy.TemporaryData.columns
        assert "up" in lazy.RawData.columns

    def test_repeating_the_same_request_is_a_no_op(self, lazy, direction_frame):
        lazy.StrucFunDir(tau_range=[1.0], dir1_data=direction_frame,
                         u_dir1="a", v_dir1="b", w_dir1="c", title="11")
        before = list(lazy._CalculatedParams)
        columns_before = list(lazy.TemporaryData.columns)
        lazy.StrucFunDir(tau_range=[1.0], dir1_data=direction_frame,
                         u_dir1="a", v_dir1="b", w_dir1="c", title="11")
        assert lazy._CalculatedParams == before
        assert list(lazy.TemporaryData.columns) == columns_before

    def test_a_new_lag_is_added_alongside_the_old_one(self, lazy, direction_frame):
        lazy.StrucFunDir(tau_range=[1.0], dir1_data=direction_frame,
                         u_dir1="a", v_dir1="b", w_dir1="c", title="11")
        lazy.StrucFunDir(tau_range=[3.0], dir1_data=direction_frame,
                         u_dir1="a", v_dir1="b", w_dir1="c", title="11")
        assert {"D11_1.0s", "D11_3.0s"} <= set(lazy.TemporaryData.columns)

    def test_a_different_title_gives_a_different_column(self, lazy, direction_frame):
        lazy.StrucFunDir(tau_range=[1.0], dir1_data=direction_frame,
                         u_dir1="a", v_dir1="b", w_dir1="c", title="11")
        lazy.StrucFunDir(tau_range=[1.0], dir1_data=direction_frame,
                         u_dir1="a", v_dir1="b", w_dir1="c", title="22")
        assert {"D11_1.0s", "D22_1.0s"} <= set(lazy.TemporaryData.columns)

    def test_the_second_direction_defaults_to_the_first(self, lazy, direction_frame):
        """Omitting dir2_data must not leave it as None."""
        lazy.StrucFunDir(tau_range=[1.0], dir1_data=direction_frame,
                         u_dir1="a", v_dir1="b", w_dir1="c", title="11")
        assert "D11_1.0s" in lazy.TemporaryData.columns

    def test_an_explicit_second_direction_is_accepted(self, lazy, direction_frame):
        second = direction_frame.rename(columns={"a": "p", "b": "q", "c": "r"})
        lazy.StrucFunDir(tau_range=[1.0], dir1_data=direction_frame,
                         u_dir1="a", v_dir1="b", w_dir1="c",
                         dir2_data=second, u_dir2="p", v_dir2="q", w_dir2="r",
                         title="12")
        assert "D12_1.0s" in lazy.TemporaryData.columns

    def test_the_supplied_direction_frame_is_not_mutated(self, lazy, direction_frame):
        before = direction_frame.copy(deep=True)
        lazy.StrucFunDir(tau_range=[1.0], dir1_data=direction_frame,
                         u_dir1="a", v_dir1="b", w_dir1="c")
        pandas.testing.assert_frame_equal(direction_frame, before)


# ---------------------------------------------------------------------------
# StrucFun -- the mean-wind basis
# ---------------------------------------------------------------------------

def _record_delegations(monkeypatch):
    """Replace StrucFunDir on the CLASS with a recorder.

    Patching the class (never an instance) keeps the substitution scoped to
    the test, and lets the basis-vector algebra be checked without touching
    the tau-lagged join at all.
    """
    calls = []

    def recorder(self, **kwargs):
        calls.append(kwargs)
        return self

    monkeypatch.setattr(singlePointTurbulenceStatistics, "StrucFunDir", recorder)
    return calls


@pytest.mark.unit
class TestStrucFunMeanDirection:
    def test_it_returns_the_calculator_for_chaining(self, lazy, monkeypatch):
        _record_delegations(monkeypatch)
        assert lazy.StrucFun(tau_range=[1.0], ubar_data=_mean_wind()) is lazy

    def test_mean_dir_asks_for_the_single_longitudinal_component(self, lazy,
                                                                 monkeypatch):
        calls = _record_delegations(monkeypatch)
        lazy.StrucFun(tau_range=[1.0], ubar_data=_mean_wind())
        assert len(calls) == 1
        assert calls[0]["title"] == "11"

    def test_the_direction_is_the_supplied_mean_velocity(self, lazy, monkeypatch):
        calls = _record_delegations(monkeypatch)
        ubar = _mean_wind()
        lazy.StrucFun(tau_range=[1.0], ubar_data=ubar)
        assert calls[0]["dir1_data"] is ubar
        assert (calls[0]["u_dir1"], calls[0]["v_dir1"], calls[0]["w_dir1"]) == \
               ("u_bar", "v_bar", "w_bar")

    def test_the_title_additions_are_appended_to_the_component_label(self, lazy,
                                                                     monkeypatch):
        calls = _record_delegations(monkeypatch)
        lazy.StrucFun(tau_range=[1.0], ubar_data=_mean_wind(), title_additions="_day")
        assert calls[0]["title"] == "11_day"

    def test_the_lag_range_is_passed_through(self, lazy, monkeypatch):
        calls = _record_delegations(monkeypatch)
        lazy.StrucFun(tau_range=[1.0, 4.0], ubar_data=_mean_wind())
        assert calls[0]["tau_range"] == [1.0, 4.0]

    def test_without_a_mean_field_it_uses_its_own_windowed_means(self, lazy,
                                                                 monkeypatch):
        calls = _record_delegations(monkeypatch)
        lazy.StrucFun(tau_range=[1.0])
        supplied = calls[0]["dir1_data"]
        assert list(supplied.columns) == ["u_bar", "v_bar", "w_bar"]
        assert isinstance(supplied, pandas.DataFrame)

    def test_the_mean_wind_magnitude_is_recorded(self, lazy, monkeypatch):
        _record_delegations(monkeypatch)
        lazy.StrucFun(tau_range=[1.0], ubar_data=_mean_wind())
        assert "u_mag" in lazy.TemporaryData.columns
        assert _param_names(lazy)[-1] == "u_mag"

    def test_the_mean_wind_magnitude_is_the_euclidean_norm(self, lazy, monkeypatch):
        _record_delegations(monkeypatch)
        ubar = _mean_wind()
        lazy.StrucFun(tau_range=[1.0], ubar_data=ubar)
        expected = numpy.sqrt(ubar["u_bar"] ** 2 + ubar["v_bar"] ** 2
                              + ubar["w_bar"] ** 2)
        computed = lazy.TemporaryData["u_mag"]
        assert numpy.allclose(numpy.asarray(computed), expected.to_numpy())


@pytest.mark.unit
class TestStrucFunThreeDimensionalBasis:
    @pytest.fixture()
    def basis(self, lazy, monkeypatch):
        calls = _record_delegations(monkeypatch)
        lazy.StrucFun(tau_range=[1.0], ubar_data=_mean_wind(), mode="3dMeanDir")
        return calls

    def test_all_six_independent_components_are_requested(self, basis):
        assert len(basis) == 6

    def test_the_components_are_the_upper_triangle_of_the_tensor(self, basis):
        assert [call["title"] for call in basis] == ["11", "22", "33", "12", "13", "23"]

    def test_the_title_additions_reach_every_component(self, lazy, monkeypatch):
        calls = _record_delegations(monkeypatch)
        lazy.StrucFun(tau_range=[1.0], ubar_data=_mean_wind(), mode="3dMeanDir",
                      title_additions="X")
        assert [call["title"] for call in calls] == \
               ["11X", "22X", "33X", "12X", "13X", "23X"]

    @staticmethod
    def _vectors(frame, prefix):
        return frame[["%s1" % prefix, "%s2" % prefix, "%s3" % prefix]].to_numpy()

    def test_the_first_axis_is_the_mean_velocity_itself(self, basis):
        frame = basis[0]["dir1_data"]
        assert numpy.allclose(self._vectors(frame, "x_hat"),
                              _mean_wind().to_numpy())

    def test_the_second_axis_is_horizontal(self, basis):
        frame = basis[1]["dir1_data"]
        assert numpy.allclose(self._vectors(frame, "y_hat")[:, 2], 0.0)

    def test_the_second_axis_is_perpendicular_to_the_mean_wind(self, basis):
        frame = basis[1]["dir1_data"]
        x_hat = self._vectors(frame, "x_hat")
        y_hat = self._vectors(frame, "y_hat")
        assert numpy.allclose((x_hat * y_hat).sum(axis=1), 0.0)

    def test_the_third_axis_is_the_cross_product_of_the_first_two(self, basis):
        frame = basis[2]["dir1_data"]
        x_hat = self._vectors(frame, "x_hat")
        y_hat = self._vectors(frame, "y_hat")
        assert numpy.allclose(self._vectors(frame, "z_hat"),
                              numpy.cross(x_hat, y_hat))

    def test_the_three_axes_are_mutually_orthogonal(self, basis):
        frame = basis[2]["dir1_data"]
        x_hat = self._vectors(frame, "x_hat")
        y_hat = self._vectors(frame, "y_hat")
        z_hat = self._vectors(frame, "z_hat")
        assert numpy.allclose((x_hat * z_hat).sum(axis=1), 0.0)
        assert numpy.allclose((y_hat * z_hat).sum(axis=1), 0.0)

    def test_the_basis_is_right_handed(self, basis):
        frame = basis[2]["dir1_data"]
        x_hat = self._vectors(frame, "x_hat")
        y_hat = self._vectors(frame, "y_hat")
        z_hat = self._vectors(frame, "z_hat")
        triple = (numpy.cross(x_hat, y_hat) * z_hat).sum(axis=1)
        assert (triple >= 0).all()

    @pytest.mark.parametrize("index, first, second", [
        (3, "x_hat", "y_hat"), (4, "x_hat", "z_hat"), (5, "y_hat", "z_hat"),
    ])
    def test_the_off_diagonal_components_pair_the_right_two_axes(self, basis,
                                                                 index, first, second):
        call = basis[index]
        assert call["u_dir1"] == "%s1" % first
        assert call["u_dir2"] == "%s1" % second


@pytest.mark.unit
class TestStrucFunRejectsAnUnknownMode:
    @pytest.mark.xfail(
        strict=True,
        reason="B303: StrucFun rejects an unknown mode with "
               "raise(\"mode must be either MeanDir or 3dMeanDir\") -- it "
               "raises the string rather than an exception, so Python "
               "substitutes TypeError('exceptions must derive from "
               "BaseException') and the diagnostic never reaches the "
               "caller. See the consolidated findings issue.",
    )
    def test_the_error_should_explain_which_modes_exist(self, lazy):
        with pytest.raises(Exception, match="mode must be either"):
            lazy.StrucFun(tau_range=[1.0], ubar_data=_mean_wind(), mode="sideways")

    def test_the_error_currently_says_nothing_about_modes(self, lazy):
        """Characterisation of B303."""
        with pytest.raises(TypeError, match="must derive from BaseException"):
            lazy.StrucFun(tau_range=[1.0], ubar_data=_mean_wind(), mode="sideways")


@pytest.mark.unit
class TestUMagGuard:
    def test_without_title_additions_the_magnitude_is_registered_once(self, lazy,
                                                                      monkeypatch):
        _record_delegations(monkeypatch)
        ubar = _mean_wind()
        lazy.StrucFun(tau_range=[1.0], ubar_data=ubar)
        lazy.StrucFun(tau_range=[1.0], ubar_data=ubar)
        assert _param_names(lazy).count("u_mag") == 1

    @pytest.mark.xfail(
        strict=True,
        reason="B304: StrucFun guards the magnitude column with "
               "`if \"u_mag\" not in self.TemporaryData.columns` but writes "
               "and registers `\"u_mag\" + title_additions`, so with a "
               "non-empty title_additions the guard inspects a name that is "
               "never created and every repeated call appends another "
               "duplicate entry to _CalculatedParams. See the consolidated "
               "findings issue.",
    )
    def test_with_title_additions_the_magnitude_should_also_be_registered_once(
            self, lazy, monkeypatch):
        _record_delegations(monkeypatch)
        ubar = _mean_wind()
        lazy.StrucFun(tau_range=[1.0], ubar_data=ubar, title_additions="X")
        lazy.StrucFun(tau_range=[1.0], ubar_data=ubar, title_additions="X")
        assert _param_names(lazy).count("u_magX") == 1

    def test_with_title_additions_it_is_currently_registered_twice(self, lazy,
                                                                   monkeypatch):
        """Characterisation of B304."""
        _record_delegations(monkeypatch)
        ubar = _mean_wind()
        lazy.StrucFun(tau_range=[1.0], ubar_data=ubar, title_additions="X")
        lazy.StrucFun(tau_range=[1.0], ubar_data=ubar, title_additions="X")
        assert _param_names(lazy).count("u_magX") == 2
        assert list(lazy.TemporaryData.columns).count("u_magX") == 1


# ---------------------------------------------------------------------------
# ThirdStrucFun
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestThirdStrucFunContract:
    def test_it_returns_the_calculator_for_chaining(self, lazy):
        assert lazy.ThirdStrucFun(tau_range=[1.0], ubar_data=_mean_wind()) is lazy

    def test_the_column_names_the_longitudinal_triple_product_and_the_lag(self, lazy):
        lazy.ThirdStrucFun(tau_range=[1.0], ubar_data=_mean_wind())
        assert "D111_1.0s" in lazy.TemporaryData.columns

    def test_one_column_per_requested_lag(self, lazy):
        lazy.ThirdStrucFun(tau_range=[1.0, 2.0], ubar_data=_mean_wind())
        assert {"D111_1.0s", "D111_2.0s"} <= set(lazy.TemporaryData.columns)

    def test_each_lag_is_registered_as_a_calculated_parameter(self, lazy):
        lazy.ThirdStrucFun(tau_range=[1.0], ubar_data=_mean_wind())
        assert "D111_1.0s" in _param_names(lazy)

    def test_the_title_additions_go_between_the_component_and_the_lag(self, lazy):
        lazy.ThirdStrucFun(tau_range=[1.0], ubar_data=_mean_wind(),
                           title_additions="_night")
        assert "D111_night_1.0s" in lazy.TemporaryData.columns

    def test_it_records_the_mean_wind_magnitude(self, lazy):
        lazy.ThirdStrucFun(tau_range=[1.0], ubar_data=_mean_wind())
        assert "u_mag" in lazy.TemporaryData.columns
        assert "u_mag" in _param_names(lazy)

    def test_the_mean_fields_are_computed_first(self, lazy):
        lazy.ThirdStrucFun(tau_range=[1.0], ubar_data=_mean_wind())
        assert "up" in lazy.RawData.columns

    def test_repeating_the_same_request_is_a_no_op(self, lazy):
        lazy.ThirdStrucFun(tau_range=[1.0], ubar_data=_mean_wind())
        before = list(lazy._CalculatedParams)
        lazy.ThirdStrucFun(tau_range=[1.0], ubar_data=_mean_wind())
        assert lazy._CalculatedParams == before

    def test_without_a_mean_field_it_uses_its_own_windowed_means(self, lazy):
        lazy.ThirdStrucFun(tau_range=[1.0])
        assert "D111_1.0s" in lazy.TemporaryData.columns


@pytest.mark.unit
class TestTheStructureFunctionsAreDaskOnly:
    @pytest.mark.xfail(
        strict=True,
        reason="B305: singlePointTurbulenceStatistics documents rawData as "
               "pandas *or* dask and every other method honours that, but "
               "StrucFunDir calls .repartition(freq=...) and StrucFun/"
               "ThirdStrucFun call .compute(), so a pandas-backed "
               "calculator cannot compute a structure function at all. See "
               "the consolidated findings issue.",
    )
    def test_a_pandas_backed_calculator_should_accept_a_structure_function(self, eager):
        eager.ThirdStrucFun(tau_range=[1.0], ubar_data=_mean_wind())

    def test_a_pandas_backed_third_structure_function_currently_raises(self, eager):
        """Characterisation of B305."""
        with pytest.raises(AttributeError, match="repartition"):
            eager.ThirdStrucFun(tau_range=[1.0], ubar_data=_mean_wind())

    def test_a_pandas_backed_default_mean_field_currently_raises(self, eager):
        """Characterisation of B305: the unconditional .compute()."""
        with pytest.raises(AttributeError, match="compute"):
            eager.StrucFun(tau_range=[1.0])

    def test_a_pandas_backed_direction_structure_function_currently_raises(self, eager,
                                                                          direction_frame):
        """Characterisation of B305."""
        with pytest.raises(AttributeError, match="repartition"):
            eager.StrucFunDir(tau_range=[1.0], dir1_data=direction_frame,
                              u_dir1="a", v_dir1="b", w_dir1="c")


# ---------------------------------------------------------------------------
# SinglePointStatisticsSpark.fluctuations
# ---------------------------------------------------------------------------

def _spark(raw, window=WINDOW):
    return SinglePointStatisticsSpark(rawData=raw, metadata=_metadata(raw, window))


def _uniform(u, v, w=0.0, T=20.0, n=60, start="2020-01-01 00:00:00"):
    index = pandas.date_range(start, periods=n, freq="1s", name="Time")
    return pandas.DataFrame({"u": numpy.full(n, float(u)),
                             "v": numpy.full(n, float(v)),
                             "w": numpy.full(n, float(w)),
                             "T": numpy.full(n, float(T))}, index=index)


@pytest.mark.unit
class TestSparkFluctuations:
    def test_it_returns_the_calculator_for_chaining(self):
        calculator = _spark(_raw())
        assert calculator.fluctuations() is calculator

    def test_it_publishes_the_windowed_means(self):
        calculator = _spark(_raw()).fluctuations()
        assert set(calculator.TemporaryData.columns) == \
               {"u_bar", "v_bar", "w_bar", "T_bar", "wind_dir_bar"}

    def test_the_windowed_means_are_the_plain_resampled_means(self):
        raw = _raw()
        calculator = _spark(raw).fluctuations()
        expected = raw.resample(WINDOW).mean()
        assert numpy.allclose(calculator.TemporaryData["u_bar"], expected["u"])
        assert numpy.allclose(calculator.TemporaryData["T_bar"], expected["T"])

    def test_it_adds_the_reynolds_fluctuations_to_the_raw_data(self):
        calculator = _spark(_raw()).fluctuations()
        assert {"up", "vp", "wp", "Tp", "wind_dir", "wind_dir_p"} <= \
               set(calculator.RawData.columns)

    def test_each_fluctuation_is_the_deviation_from_its_window_mean(self):
        data = _spark(_raw()).fluctuations().RawData
        for component, mean in (("u", "u_bar"), ("v", "v_bar"),
                                ("w", "w_bar"), ("T", "T_bar")):
            assert numpy.allclose(data["%sp" % component],
                                  data[component] - data[mean])

    def test_the_first_sample_is_snapped_onto_the_averaging_grid(self):
        raw = _raw(start="2020-01-01 00:00:07")
        calculator = _spark(raw).fluctuations()
        assert calculator.RawData.index[0] == pandas.Timestamp("2020-01-01 00:00:00")

    def test_snapping_lets_the_first_sample_pick_up_its_window_mean(self):
        raw = _raw(start="2020-01-01 00:00:07")
        calculator = _spark(raw).fluctuations()
        assert not numpy.isnan(calculator.RawData["u_bar"].iloc[0])

    @pytest.mark.parametrize("u, v, expected", [
        (1.0, 0.0, WEST),
        (-1.0, 0.0, EAST),
        (0.0, 1.0, SOUTH),
        (0.0, -1.0, NORTH),
    ])
    def test_the_wind_direction_is_the_meteorological_bearing_it_blows_from(self, u, v,
                                                                           expected):
        data = _spark(_uniform(u, v)).fluctuations().RawData
        assert numpy.allclose(data["wind_dir"] % 360.0, expected)

    def test_the_mean_wind_direction_follows_the_same_convention(self):
        calculator = _spark(_uniform(1.0, 0.0)).fluctuations()
        assert numpy.allclose(calculator.TemporaryData["wind_dir_bar"] % 360.0, WEST)

    def test_a_south_westerly_flow_reads_between_west_and_south(self):
        data = _spark(_uniform(1.0, 1.0)).fluctuations().RawData
        assert numpy.allclose(data["wind_dir"], 225.0)

    def test_the_direction_fluctuation_is_a_non_negative_angular_distance(self):
        data = _spark(_raw()).fluctuations().RawData
        assert (data["wind_dir_p"] >= 0).all()
        assert (data["wind_dir_p"] <= 180.0).all()

    def test_a_steady_wind_has_no_direction_fluctuation(self):
        data = _spark(_uniform(2.0, 1.0)).fluctuations().RawData
        assert numpy.allclose(data["wind_dir_p"], 0.0)

    def test_calling_it_twice_does_not_recompute(self):
        calculator = _spark(_raw()).fluctuations()
        first = calculator.RawData
        calculator.fluctuations()
        assert calculator.RawData is first

    def test_with_no_sampling_window_a_single_global_mean_is_used(self):
        raw = _raw()
        calculator = SinglePointStatisticsSpark(rawData=raw,
                                                metadata=_metadata(raw, window=None))
        calculator.fluctuations()
        assert len(calculator.TemporaryData) == 1
        assert numpy.allclose(calculator.RawData["u_bar"], raw["u"].mean())

    @pytest.mark.xfail(
        strict=True,
        reason="B306: SinglePointStatisticsSpark.fluctuations never sets "
               "self.data, which the base implementation does as its last "
               "act, so getData() still answers None after the override has "
               "computed everything -- and getData's docstring says None "
               "means no computation has happened yet. See the consolidated "
               "findings issue.",
    )
    def test_the_computed_frame_should_be_reachable_through_get_data(self):
        calculator = _spark(_raw()).fluctuations()
        assert calculator.getData() is not None

    def test_get_data_currently_stays_empty(self):
        """Characterisation of B306."""
        calculator = _spark(_raw()).fluctuations()
        assert calculator.getData() is None
        assert "up" in calculator.RawData.columns

    @pytest.mark.xfail(
        strict=True,
        reason="B307: the Spark override adds wind_dir_bar to "
               "_TemporaryData but, unlike the base implementation, leaves "
               "it out of _CalculatedParams. compute() projects "
               "_TemporaryData onto that list, so the mean wind direction is "
               "silently dropped from the result. See the consolidated "
               "findings issue.",
    )
    def test_the_mean_wind_direction_should_survive_compute(self):
        calculator = _spark(_raw()).fluctuations()
        assert "wind_dir_bar" in calculator.compute().columns

    def test_the_mean_wind_direction_is_currently_dropped_by_compute(self):
        """Characterisation of B307."""
        calculator = _spark(_raw()).fluctuations()
        assert "wind_dir_bar" in calculator.TemporaryData.columns
        assert _param_names(calculator) == ["u_bar", "v_bar", "w_bar", "T_bar"]
        assert "wind_dir_bar" not in calculator.compute().columns

    def test_the_base_implementation_does_register_it(self):
        raw = _raw()
        base = singlePointTurbulenceStatistics(rawData=raw, metadata=_metadata(raw))
        base.fluctuations()
        assert "wind_dir_bar" in _param_names(base)


# ---------------------------------------------------------------------------
# InMemoryRawData: HDF5 plus JSON sidecar
# ---------------------------------------------------------------------------

@pytest.fixture()
def recorded_hdf_writes(monkeypatch):
    """Patch pandas' HDF writer on the CLASS, so PyTables is not needed.

    ``InMemoryRawData.to_hdf`` copies itself first, and the copy is a plain
    ``pandas.DataFrame``, so this is the call it makes.
    """
    writes = []

    def fake_to_hdf(self, path_or_buf, key=None, **kwargs):
        writes.append(dict(path=path_or_buf, key=key, kwargs=kwargs,
                           columns=list(self.columns)))

    monkeypatch.setattr(pandas.DataFrame, "to_hdf", fake_to_hdf)
    return writes


@pytest.mark.unit
class TestToHdf:
    def test_the_extension_is_replaced_with_hdf(self, recorded_hdf_writes, tmp_path):
        data = InMemoryRawData(pandas.DataFrame({"u": [1.0, 2.0]}))
        data.to_hdf(str(tmp_path / "sonic.h5"), "raw")
        assert recorded_hdf_writes[0]["path"] == str(tmp_path / "sonic.hdf")

    def test_the_key_and_the_extra_arguments_are_forwarded(self, recorded_hdf_writes,
                                                           tmp_path):
        data = InMemoryRawData(pandas.DataFrame({"u": [1.0]}))
        data.to_hdf(str(tmp_path / "sonic.hdf"), "raw", mode="w", complevel=9)
        assert recorded_hdf_writes[0]["key"] == "raw"
        assert recorded_hdf_writes[0]["kwargs"] == dict(mode="w", complevel=9)

    def test_the_data_itself_is_handed_to_pandas(self, recorded_hdf_writes, tmp_path):
        data = InMemoryRawData(pandas.DataFrame({"u": [1.0], "v": [2.0]}))
        data.to_hdf(str(tmp_path / "sonic.hdf"), "raw")
        assert recorded_hdf_writes[0]["columns"] == ["u", "v"]

    def test_attributes_are_written_to_a_json_sidecar_beside_the_data(
            self, recorded_hdf_writes, tmp_path):
        data = InMemoryRawData(pandas.DataFrame({"u": [1.0]}))
        data._Attrs = {"height": 10, "station": "YAVNEEL"}
        data.to_hdf(str(tmp_path / "sonic.hdf"), "raw")
        sidecar = tmp_path / "sonic.json"
        assert json.loads(sidecar.read_text()) == {"height": 10, "station": "YAVNEEL"}

    def test_no_sidecar_is_written_when_there_are_no_attributes(
            self, recorded_hdf_writes, tmp_path):
        data = InMemoryRawData(pandas.DataFrame({"u": [1.0]}))
        data.to_hdf(str(tmp_path / "sonic.hdf"), "raw")
        assert not (tmp_path / "sonic.json").exists()

    def test_a_second_write_merges_into_the_existing_sidecar(self, recorded_hdf_writes,
                                                             tmp_path):
        target = str(tmp_path / "sonic.hdf")
        first = InMemoryRawData(pandas.DataFrame({"u": [1.0]}))
        first._Attrs = {"height": 10}
        first.to_hdf(target, "raw")
        second = InMemoryRawData(pandas.DataFrame({"u": [2.0]}))
        second._Attrs = {"station": "YAVNEEL"}
        second.to_hdf(target, "raw")
        assert json.loads((tmp_path / "sonic.json").read_text()) == \
               {"height": 10, "station": "YAVNEEL"}

    def test_a_later_write_wins_on_a_shared_key(self, recorded_hdf_writes, tmp_path):
        target = str(tmp_path / "sonic.hdf")
        first = InMemoryRawData(pandas.DataFrame({"u": [1.0]}))
        first._Attrs = {"height": 10}
        first.to_hdf(target, "raw")
        second = InMemoryRawData(pandas.DataFrame({"u": [2.0]}))
        second._Attrs = {"height": 20}
        second.to_hdf(target, "raw")
        assert json.loads((tmp_path / "sonic.json").read_text()) == {"height": 20}

    @pytest.mark.xfail(
        strict=True,
        reason="B308: to_hdf derives both output paths with "
               "path_or_buf.rpartition('.')[0], which is the empty string "
               "for a path that has no extension, so a request to write "
               "<dir>/rawdata writes '.hdf' and '.json' into the current "
               "working directory instead of into <dir>. See the "
               "consolidated findings issue.",
    )
    def test_an_extensionless_path_should_still_write_beside_the_data(
            self, recorded_hdf_writes, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        target = tmp_path / "station"
        target.mkdir()
        data = InMemoryRawData(pandas.DataFrame({"u": [1.0]}))
        data._Attrs = {"height": 10}
        data.to_hdf(str(target / "rawdata"), "raw")
        assert (target / "rawdata.json").exists()

    def test_an_extensionless_path_currently_writes_into_the_working_directory(
            self, recorded_hdf_writes, tmp_path, monkeypatch):
        """Characterisation of B308."""
        monkeypatch.chdir(tmp_path)
        target = tmp_path / "station"
        target.mkdir()
        data = InMemoryRawData(pandas.DataFrame({"u": [1.0]}))
        data._Attrs = {"height": 10}
        data.to_hdf(str(target / "rawdata"), "raw")
        assert recorded_hdf_writes[0]["path"] == ".hdf"
        assert (tmp_path / ".json").exists()
        assert not (target / "rawdata.json").exists()


@pytest.fixture()
def stub_hdf_reader(monkeypatch):
    """Patch ``pandas.read_hdf`` so the sidecar logic is reachable without PyTables."""
    frame = pandas.DataFrame({"u": [1.0, 2.0, 3.0]})
    reads = []

    def fake_read_hdf(path_or_buf, key=None, **kwargs):
        reads.append(dict(path=path_or_buf, key=key, kwargs=kwargs))
        return frame

    monkeypatch.setattr(pandas, "read_hdf", fake_read_hdf)
    return reads, frame


@pytest.mark.unit
class TestReadHdf:
    def test_it_returns_the_in_memory_container_not_a_bare_frame(self, stub_hdf_reader,
                                                                 tmp_path):
        loaded = InMemoryRawData.read_hdf(str(tmp_path / "sonic.hdf"), key="raw")
        assert isinstance(loaded, InMemoryRawData)

    def test_the_data_comes_from_the_file(self, stub_hdf_reader, tmp_path):
        _, frame = stub_hdf_reader
        loaded = InMemoryRawData.read_hdf(str(tmp_path / "sonic.hdf"), key="raw")
        assert loaded["u"].tolist() == frame["u"].tolist()

    def test_the_path_and_key_are_forwarded_to_pandas(self, stub_hdf_reader, tmp_path):
        reads, _ = stub_hdf_reader
        target = str(tmp_path / "sonic.hdf")
        InMemoryRawData.read_hdf(target, key="raw", columns=["u"])
        assert reads[0] == dict(path=target, key="raw", kwargs=dict(columns=["u"]))

    def test_the_key_may_be_omitted(self, stub_hdf_reader, tmp_path):
        reads, _ = stub_hdf_reader
        InMemoryRawData.read_hdf(str(tmp_path / "sonic.hdf"))
        assert reads[0]["key"] is None

    def test_a_json_sidecar_beside_the_file_becomes_the_attributes(self,
                                                                   stub_hdf_reader,
                                                                   tmp_path):
        (tmp_path / "sonic.json").write_text(json.dumps({"height": 10}))
        loaded = InMemoryRawData.read_hdf(str(tmp_path / "sonic.hdf"), key="raw")
        assert loaded._Attrs == {"height": 10}

    def test_without_a_sidecar_the_attributes_stay_empty(self, stub_hdf_reader,
                                                          tmp_path):
        loaded = InMemoryRawData.read_hdf(str(tmp_path / "sonic.hdf"), key="raw")
        assert loaded._Attrs == {}

    def test_a_sidecar_belonging_to_another_file_is_not_picked_up(self,
                                                                   stub_hdf_reader,
                                                                   tmp_path):
        (tmp_path / "other.json").write_text(json.dumps({"height": 10}))
        loaded = InMemoryRawData.read_hdf(str(tmp_path / "sonic.hdf"), key="raw")
        assert loaded._Attrs == {}


@pytest.mark.unit
class TestHdfRoundTrip:
    def test_the_data_and_the_attributes_survive_a_round_trip(self, tmp_path):
        pytest.importorskip("tables",
                            reason="pandas' HDF5 support is PyTables-only and "
                                   "PyTables is not installed in this environment")
        data = InMemoryRawData(pandas.DataFrame({"u": [1.0, 2.0], "v": [3.0, 4.0]}))
        data._Attrs = {"height": 10, "station": "YAVNEEL"}
        target = str(tmp_path / "sonic.hdf")
        data.to_hdf(target, "raw")

        loaded = InMemoryRawData.read_hdf(target, key="raw")
        assert loaded["u"].tolist() == [1.0, 2.0]
        assert loaded._Attrs == {"height": 10, "station": "YAVNEEL"}
        assert os.path.isfile(str(tmp_path / "sonic.json"))
