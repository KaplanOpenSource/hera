"""hermesWorkflowToolkit: the pure name-manipulation static methods.

B102: splitWorkflowName joins the base-name parts with "".join(...)
instead of "_".join(...), so any base name containing its own underscore
gets mangled -- not the inverse of getworkFlowName, which always
reconstructs as f"{baseName}_{formatted_number}".
"""
import pytest

from hera.simulations.hermesWorkflowToolkit import hermesWorkflowToolkit


@pytest.mark.unit
class TestGetWorkFlowName:
    def test_the_id_is_zero_padded_to_four_digits(self):
        assert hermesWorkflowToolkit.getworkFlowName("myflow", 7) == "myflow_0007"

    def test_a_large_id_is_not_truncated(self):
        assert hermesWorkflowToolkit.getworkFlowName("myflow", 12345) == "myflow_12345"


@pytest.mark.unit
class TestSplitWorkflowName:
    def test_a_name_without_an_underscore_has_no_id(self):
        assert hermesWorkflowToolkit.splitWorkflowName("myflow") == ("myflow", None)

    def test_a_simple_base_name_round_trips_with_get_work_flow_name(self):
        name = hermesWorkflowToolkit.getworkFlowName("myflow", 7)
        assert hermesWorkflowToolkit.splitWorkflowName(name) == ("myflow", "0007")

    @pytest.mark.xfail(
        strict=True,
        reason="B102: joins base-name parts with ''.join instead of "
               "'_'.join, so a base name containing its own underscore "
               "is not reconstructed correctly. "
               "See the consolidated findings issue.",
    )
    def test_a_base_name_with_its_own_underscore_should_round_trip(self):
        name = hermesWorkflowToolkit.getworkFlowName("my_flow", 7)
        assert hermesWorkflowToolkit.splitWorkflowName(name) == ("my_flow", "0007")

    def test_a_base_name_with_its_own_underscore_currently_loses_it(self):
        """Characterisation of B102."""
        base, flow_id = hermesWorkflowToolkit.splitWorkflowName("my_flow_0007")
        assert base == "myflow"
        assert flow_id == "0007"


@pytest.mark.unit
class TestIsInFolder:
    def test_a_subpath_is_inside_the_base(self):
        assert hermesWorkflowToolkit.isInFolder("/a/b", "/a/b/c") is True

    def test_a_sibling_path_is_not_inside_the_base(self):
        assert hermesWorkflowToolkit.isInFolder("/a/b", "/a/x") is False

    def test_the_base_path_itself_is_inside_itself(self):
        assert hermesWorkflowToolkit.isInFolder("/a/b", "/a/b") is True
