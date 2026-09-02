"""ProtectionPolicy: chained protection actions applied to a concentration field.

Pure computation over xarray/pandas/numpy/pint -- no DB, no toolkit, no
stubbed dependency.  One live defect surfaced while probing this module:

* B54: ``ActionIndoor`` constructed with ``alpha=`` (rather than
  ``turnover=``) calls ``.ureg(...)`` on a pint ``Quantity``, which has no
  such method. Only the ``turnover=`` keyword -- not ``alpha=``, despite
  being the class's own documented and first-listed example -- can actually
  build an indoor action.

B53, resolved independently: ``addActions`` given a string used to check
``os.JSONpath.exists`` -- an attribute ``os`` does not have -- so the
file-loading branch could never run. This was a real, accurate finding at
the time; PR #1010 (commit 14508043) fixed the typo back to
``os.path.exists`` on master, unrelated to this test-expansion effort.
See test_addactions_from_a_json_string_loads_the_file.

B58, retracted: an earlier pass claimed calling ``.indoor(...)``/
``.masks(...)`` with no begin/end/enter/stay crashed, because
``data[datetimename].to_series()[0]`` (a positional lookup on a
DatetimeIndex-backed Series) raises ``KeyError`` on newer pandas. That was
verified against a locally-drifted pandas (3.0.2), not the pandas actually
pinned in requirements.txt (2.2.3) -- under 2.2.3 the positional fallback
still works, with only a FutureWarning. See
test_compute_with_no_time_window_at_all, which now asserts the real
(working) result instead of xfailing. Most tests below still pass an
explicit window via ``_window()``, simply because it makes the assertions
easier to reason about, not to route around a defect.
"""
import numpy
import pandas
import pytest
import xarray

from hera.riskassessment.protectionpolicy.ProtectionPolicy import (
    ActionIndoor,
    ActionMasks,
    ProtectionPolicy,
    abstractAction,
)
from hera.utils import ureg


def _dataset(n=5, value=10.0):
    times = pandas.date_range("2020-01-01", periods=n, freq="1h")
    C = xarray.DataArray(numpy.full(n, value), dims=["datetime"], coords={"datetime": times})
    ds = xarray.Dataset({"C": C})
    ds.attrs = {}
    return ds


def _window(ds):
    """An explicit begin/end spanning the whole dataset, for easier assertions."""
    return dict(begin=ds.datetime[0].values, end=ds.datetime[-1].values)


@pytest.mark.unit
class TestProtectionPolicyConstruction:
    def test_default_axis_and_datetime_names(self):
        policy = ProtectionPolicy()
        assert policy.xname == "x"
        assert policy.yname == "y"
        assert policy.datetimename == "datetime"

    def test_custom_axis_and_datetime_names(self):
        policy = ProtectionPolicy(x="lon", y="lat", datetime="time")
        assert (policy.xname, policy.yname, policy.datetimename) == ("lon", "lat", "time")

    def test_the_final_field_name_defaults_to_c(self):
        assert ProtectionPolicy().finalname == "C"

    def test_data_is_none_before_compute(self):
        assert ProtectionPolicy().data is None

    def test_params_is_never_populated_by_the_constructor(self):
        """Characterisation: _params is a class attribute the constructor never assigns."""
        assert ProtectionPolicy().params is None

    def test_an_empty_action_list_leaves_no_actions(self):
        policy = ProtectionPolicy(actionList=[])
        assert policy.hdfkey == ""


@pytest.mark.unit
class TestActionDispatch:
    def test_indoor_returns_the_policy_for_chaining(self):
        policy = ProtectionPolicy()
        assert policy.indoor(turnover=2 * ureg.h) is policy

    def test_masks_returns_the_policy_for_chaining(self):
        policy = ProtectionPolicy()
        assert policy.masks(protectionFactor=1000) is policy

    def test_indoor_without_alpha_or_turnover_raises(self):
        with pytest.raises(ValueError, match="alpha or turnover"):
            ProtectionPolicy().indoor()

    def test_masks_without_protection_factor_raises(self):
        with pytest.raises(ValueError, match="protectionFactor"):
            ProtectionPolicy().masks()

    def test_add_action_with_an_unknown_name_raises(self):
        with pytest.raises(Exception, match="was not found"):
            ProtectionPolicy().addAction("NoSuchAction", {})

    def test_get_action_titles_the_name_before_lookup(self):
        """indoor()/masks() pass already-titled names; addAction lower-cases too."""
        policy = ProtectionPolicy()
        action = abstractAction.getAction(1, policy, "indoor", {"turnover": 2 * ureg.h})
        assert isinstance(action, ActionIndoor)

    def test_add_action_returns_the_new_action(self):
        policy = ProtectionPolicy()
        action = policy.addAction("Masks", {"protectionFactor": 500})
        assert isinstance(action, ActionMasks)
        assert action.actionid == 1

    def test_second_action_gets_the_next_sequential_id(self):
        policy = ProtectionPolicy()
        policy.addAction("Masks", {"protectionFactor": 500})
        second = policy.addAction("Masks", {"protectionFactor": 500})
        assert second.actionid == 2

    def test_addactions_from_a_dict_reads_the_actions_key(self):
        policy = ProtectionPolicy(actionList=[
            {"name": "Masks", "params": {"protectionFactor": 100}},
        ])
        assert len(policy._actionList) == 1
        assert policy._actionList[0].actiontype == "masks"

    def test_addactions_from_a_json_string_loads_the_file(self, tmp_path):
        """B53, fixed upstream: this typo'd `os.JSONpath.exists(...)` to
        `os.path.exists(...)` in PR #1010 (commit 14508043), independently
        of this test-expansion effort."""
        jsonfile = tmp_path / "actions.json"
        jsonfile.write_text('{"actions": [{"name": "Masks", "params": {"protectionFactor": 10}}]}')
        policy = ProtectionPolicy()
        policy.addActions(str(jsonfile))
        assert len(policy._actionList) == 1


@pytest.mark.unit
class TestActionIndoor:
    def test_alpha_is_the_inverse_of_turnover(self):
        action = ActionIndoor(1, ProtectionPolicy(), turnover=2 * ureg.h)
        assert action.alpha.m_as(1 / ureg.h) == pytest.approx(0.5)

    def test_turnover_property_reads_back_the_inverse_of_alpha(self):
        action = ActionIndoor(1, ProtectionPolicy(), turnover=2 * ureg.h)
        assert action.turnover.m_as(ureg.h) == pytest.approx(2.0)

    def test_missing_alpha_and_turnover_raises(self):
        with pytest.raises(ValueError, match="alpha or turnover"):
            ActionIndoor(1, ProtectionPolicy())

    @pytest.mark.xfail(
        strict=True,
        reason="B54: the alpha branch does `kwargs['alpha'].ureg(1/ureg.h)`, "
               "but pint Quantity objects have no `.ureg` method (that call "
               "raises AttributeError). Only turnover=, not the "
               "class-docstring's own alpha=... example, actually works. "
               "See the consolidated findings issue.",
    )
    def test_constructing_with_alpha_directly(self):
        action = ActionIndoor(1, ProtectionPolicy(), alpha=0.5 / ureg.h)
        assert action.alpha.m_as(1 / ureg.h) == pytest.approx(0.5)

    def test_hdfkey_encodes_the_turnover_in_minutes(self):
        ds = _dataset()
        policy = ProtectionPolicy().indoor(turnover=2 * ureg.h, **_window(ds))
        policy.compute(ds)
        assert policy._actionList[0].hdfkey.startswith("indoorT120min")

    def test_indoor_field_is_always_zero_today(self):
        """B59: the recurrence never actually writes into Cin.

        ``Cin[curstep].values = ...`` assigns into the DataArray returned by
        dict-based ``__getitem__``, which is a fresh copy on this xarray
        version rather than a view -- so every write in the loop is
        discarded and ``indoor_1`` comes out as all zeros regardless of
        alpha, dt or the outdoor concentration.
        """
        ds = _dataset(n=5, value=10.0)
        policy = ProtectionPolicy().indoor(turnover=2 * ureg.h, **_window(ds))
        data = policy.compute(ds)
        assert numpy.array_equal(data["indoor_1"].values, numpy.zeros(5))

    @pytest.mark.xfail(
        strict=True,
        reason="B59: Cin[curstep].values = ... discards its write because "
               "dict-indexing (`Cin[{'datetime': I}]`) returns a copy, not "
               "a view, on this xarray version. The indoor recurrence "
               "therefore never runs -- indoor_1 is all zeros no matter "
               "what alpha, dt or the outdoor concentration are. "
               "See the consolidated findings issue.",
    )
    def test_compute_builds_a_low_pass_filtered_field_towards_the_outdoor_value(self):
        ds = _dataset(n=5, value=10.0)
        policy = ProtectionPolicy().indoor(turnover=2 * ureg.h, **_window(ds))
        data = policy.compute(ds)
        indoor = data["indoor_1"].values
        # dt=1h, turnover=2h => alpha*dt=0.5 exactly; Cin[0]=0 by construction.
        assert indoor[0] == pytest.approx(0.0)
        assert indoor[1] == pytest.approx((0 + 0.5 * 10) / 1.5)
        assert indoor[2] == pytest.approx((indoor[1] + 0.5 * 10) / 1.5)
        assert numpy.all(numpy.diff(indoor) > 0)
        assert numpy.all(indoor <= 10.0)

    def test_compute_writes_the_final_field_from_the_indoor_action(self):
        ds = _dataset()
        policy = ProtectionPolicy().indoor(turnover=2 * ureg.h, **_window(ds))
        data = policy.compute(ds)
        assert numpy.array_equal(data[policy.finalname].values, data["indoor_1"].values)

    def test_compute_records_the_action_metadata_on_the_dataset(self):
        ds = _dataset()
        policy = ProtectionPolicy().indoor(turnover=2 * ureg.h, **_window(ds))
        data = policy.compute(ds)
        assert data.attrs["1"]["type"] == "indoor"
        assert data.attrs["1"]["outputs"] == ["indoor_1"]

    def test_compute_with_no_time_window_at_all(self):
        """Not B58 under the pinned pandas (see the module docstring):
        `.to_series()[0]`/`[-1]` still fall back to positional lookup here,
        with only a FutureWarning -- so this simplest possible call succeeds
        and picks up the dataset's actual first/last timestamps."""
        ds = _dataset()
        policy = ProtectionPolicy().indoor(turnover=2 * ureg.h)
        data = policy.compute(ds)
        assert data.attrs["1"]["params"]["begin"] == pandas.Timestamp(ds.datetime[0].values)
        assert data.attrs["1"]["params"]["end"] == pandas.Timestamp(ds.datetime[-1].values)


@pytest.mark.unit
class TestActionMasks:
    def test_protection_factor_is_stored_as_a_float(self):
        action = ActionMasks(1, ProtectionPolicy(), protectionFactor="1000")
        assert action.protectionFactor == pytest.approx(1000.0)

    def test_missing_protection_factor_raises(self):
        with pytest.raises(ValueError, match="protectionFactor"):
            ActionMasks(1, ProtectionPolicy())

    def test_hdfkey_encodes_the_protection_factor(self):
        ds = _dataset()
        policy = ProtectionPolicy().masks(protectionFactor=1000, **_window(ds))
        policy.compute(ds)
        assert policy._actionList[0].hdfkey.startswith("maskPF1000")

    def test_compute_divides_the_outdoor_field_by_the_protection_factor(self):
        ds = _dataset(value=10.0)
        policy = ProtectionPolicy().masks(protectionFactor=1000, **_window(ds))
        data = policy.compute(ds)
        assert numpy.allclose(data["masks_1"].values, 10.0 / 1000)

    def test_compute_writes_the_final_field_from_the_masking_action(self):
        ds = _dataset(value=10.0)
        policy = ProtectionPolicy().masks(protectionFactor=1000, **_window(ds))
        data = policy.compute(ds)
        assert numpy.array_equal(data[policy.finalname].values, data["masks_1"].values)


@pytest.mark.unit
class TestChainedPolicy:
    def test_indoor_then_masks_both_leave_metadata_and_the_policy_hdfkey_joins_both(self):
        ds = _dataset()
        window = _window(ds)
        policy = ProtectionPolicy().indoor(turnover=2 * ureg.h, **window).masks(protectionFactor=100, **window)
        data = policy.compute(ds)
        assert set(data.attrs) >= {"0", "1", "2"}
        indoor_key, masks_key = policy.hdfkey.split("/")
        assert indoor_key.startswith("indoorT120minEnter0 days 00:00:00Stay")
        assert masks_key.startswith("maskPF100Wear0 days 00:00:00Duration")

    def test_compute_seeds_outdoor_0_from_the_input_field(self):
        policy = ProtectionPolicy()
        data = policy.compute(_dataset(value=7.0))
        assert numpy.allclose(data["outdoor_0"].values, 7.0)
        assert data.attrs["0"] == {"type": "outdoor"}


@pytest.mark.unit
class TestAbstractActionHdfkeyIsUnimplemented:
    def test_the_base_class_hdfkey_is_a_no_op_stub(self):
        action = abstractAction(actionID=1, actiontype="x", policy=None)
        assert action.hdfkey is None
