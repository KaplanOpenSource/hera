"""hermesWorkflowToolkit: the pure name-manipulation static methods.

B102, fixed: splitWorkflowName joined the base-name parts with
"".join(...) instead of "_".join(...), so any base name containing its own
underscore got mangled -- it was not the inverse of getworkFlowName, which
always reconstructs as f"{baseName}_{formatted_number}".
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

    def test_a_base_name_with_its_own_underscore_round_trips(self):
        name = hermesWorkflowToolkit.getworkFlowName("my_flow", 7)
        assert hermesWorkflowToolkit.splitWorkflowName(name) == ("my_flow", "0007")

    def test_every_underscore_of_the_base_name_is_kept(self):
        """B102: the parts used to be joined with "" instead of "_"."""
        assert hermesWorkflowToolkit.splitWorkflowName("a_b_c_0007") == (
            "a_b_c", "0007",
        )


@pytest.mark.unit
class TestIsInFolder:
    def test_a_subpath_is_inside_the_base(self):
        assert hermesWorkflowToolkit.isInFolder("/a/b", "/a/b/c") is True

    def test_a_sibling_path_is_not_inside_the_base(self):
        assert hermesWorkflowToolkit.isInFolder("/a/b", "/a/x") is False

    def test_the_base_path_itself_is_inside_itself(self):
        assert hermesWorkflowToolkit.isInFolder("/a/b", "/a/b") is True
