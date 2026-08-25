"""The high-frequency cache lookup, and what the query bug actually does to it.

P1 (batch 0) established that dictToMongoQuery indexes list values.  This file
determines the real consequence for the turbulence cache in
abstractcalculator, which passes params__all a list at lines 172 and 190.

The answer is subtler than "the lookup never matches".  mongoengine recognises
`all` as an operator even in the middle of the flattened key, so

    params__all=["u", "v"]   ->   desc__params__all__0 : "u"
                                  desc__params__all__1 : "v"
                             ->   {"desc.params.0": {"$all": ["u"]},
                                   "desc.params.1": {"$all": ["v"]}}

which is a positional-PREFIX match, not set containment.  It succeeds when the
requested parameters happen to be a prefix of the stored list in the same
order, and silently misses otherwise.  Both halves are pinned below.
"""
import pytest

from hera.utils.query import dictToMongoQuery


@pytest.mark.unit
class TestTheQueryTheCacheBuilds:
    def test_a_list_valued_all_operator_is_indexed_by_position(self):
        built = dictToMongoQuery({"params__all": ["u", "v", "w"]}, prefix="desc")
        assert built == {
            "desc__params__all__0": "u",
            "desc__params__all__1": "v",
            "desc__params__all__2": "w",
        }

    def test_the_operator_no_longer_applies_to_the_whole_list(self):
        built = dictToMongoQuery({"params__all": ["u"]}, prefix="desc")
        assert "desc__params__all" not in built


@pytest.mark.unit
class TestPositionalPrefixBehaviour:
    """Stored document holds params = [u, v, w]."""

    @pytest.fixture(autouse=True)
    def _stored(self, unit_project):
        unit_project.addCacheDocument(
            resource="/cache/turbulence.parquet",
            dataFormat="parquet",
            type="TurbulenceCalculator",
            desc={"params": ["u", "v", "w"], "station": "YAVNEEL"},
        )

    @pytest.mark.parametrize("requested", [["u"], ["u", "v"], ["u", "v", "w"]])
    def test_a_prefix_in_order_matches(self, unit_project, requested):
        found = unit_project.getCacheDocuments(
            type="TurbulenceCalculator", params__all=requested
        )
        assert len(found) == 1

    @pytest.mark.xfail(
        strict=True,
        reason="P1: $all is defined as order-independent set containment, but the "
               "flattening turns it into a positional match, so asking for a "
               "parameter that is present but not at position 0 misses. "
               "See the consolidated findings issue.",
    )
    @pytest.mark.parametrize("requested", [["v"], ["w"], ["v", "u"], ["w", "u"]])
    def test_every_element_present_should_match_whatever_the_order(
        self, unit_project, requested
    ):
        found = unit_project.getCacheDocuments(
            type="TurbulenceCalculator", params__all=requested
        )
        assert len(found) == 1

    def test_the_miss_is_silent(self, unit_project):
        """No exception and no warning -- just a recomputation the caller pays for."""
        found = unit_project.getCacheDocuments(
            type="TurbulenceCalculator", params__all=["w"]
        )
        assert len(found) == 0, "if this now matches, P1 has been fixed"

    def test_an_exact_equality_query_still_works(self, unit_project):
        """Plain equality on the list is unaffected -- only the operator degrades."""
        found = unit_project.getCacheDocuments(
            type="TurbulenceCalculator", params=["u", "v", "w"]
        )
        assert len(found) == 1
