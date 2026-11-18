"""openFoam/lagrangian/LSM/toolkit.py: OFLSMToolkit."""
import pytest

from hera.simulations.openFoam.lagrangian.LSM.toolkit import OFLSMToolkit


@pytest.mark.unit
class TestConstructor:
    def test_it_is_constructible(self, unit_files_directory):
        from hera import toolkitHome

        toolkitHome.getToolkit(
            toolkitName="OF_LSM", projectName="OFLSM_TEST", filesDirectory=unit_files_directory,
        )

    def test_it_holds_the_vector_topography_toolkit(self, unit_files_directory):
        """The only uses of self.topography in this module call
        analysis.addHeight, which lives on the vector topography toolkit."""
        from hera import toolkitHome

        lsm = toolkitHome.getToolkit(
            toolkitName="OF_LSM", projectName="OFLSM_TEST", filesDirectory=unit_files_directory,
        )
        assert hasattr(lsm.topography.analysis, "addHeight")


@pytest.mark.unit
class TestPropertiesViaBypass:
    """Verified via __new__, independently of the datalayer."""

    @pytest.fixture()
    def t(self):
        return OFLSMToolkit.__new__(OFLSMToolkit)

    def test_doctype_is_a_fixed_string(self, t):
        assert t.doctype == "LSMRuns"

    def test_case_path_getter_setter(self, t):
        t.casePath = "/case"
        assert t.casePath == "/case"

    def test_cloud_name_is_coerced_to_str(self, t):
        t.cloudName = 123
        assert t.cloudName == "123"

    def test_parallel_case_getter_setter(self, t):
        t.parallelCase = True
        assert t.parallelCase is True
