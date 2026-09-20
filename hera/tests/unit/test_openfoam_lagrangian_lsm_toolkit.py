"""openFoam/lagrangian/LSM/toolkit.py: OFLSMToolkit.

B103: __init__ always crashes -- it references
toolkitHome.GIS_TOPOGRAPHY, a constant that does not exist (the real ones
are GIS_RASTER_TOPOGRAPHY/GIS_VECTOR_TOPOGRAPHY). No instance of
OFLSMToolkit can ever be constructed via toolkitHome.getToolkit. The
property tests below use __new__ to bypass the broken constructor and
verify the properties work correctly once B103 is fixed.
"""
import pytest

from hera.simulations.openFoam.lagrangian.LSM.toolkit import OFLSMToolkit


@pytest.mark.unit
class TestConstructorIsBroken:
    def test_gis_topography_constant_does_not_exist(self):
        """Characterisation of B103, at its exact source."""
        from hera import toolkitHome

        assert not hasattr(toolkitHome, "GIS_TOPOGRAPHY")

    @pytest.mark.xfail(
        strict=True,
        reason="B103: __init__ references toolkitHome.GIS_TOPOGRAPHY, a "
               "constant that does not exist (the real ones are "
               "GIS_RASTER_TOPOGRAPHY/GIS_VECTOR_TOPOGRAPHY). No instance "
               "can ever be constructed. See the consolidated findings issue.",
    )
    def test_it_should_be_constructible(self, unit_files_directory):
        from hera import toolkitHome

        toolkitHome.getToolkit(
            toolkitName="OF_LSM", projectName="OFLSM_TEST", filesDirectory=unit_files_directory,
        )

    def test_it_currently_raises_on_every_construction(self, unit_files_directory):
        """Characterisation of B103."""
        from hera import toolkitHome

        with pytest.raises(AttributeError, match="GIS_TOPOGRAPHY"):
            toolkitHome.getToolkit(
                toolkitName="OF_LSM", projectName="OFLSM_TEST", filesDirectory=unit_files_directory,
            )


@pytest.mark.unit
class TestPropertiesViaBypass:
    """Verified independently of the broken constructor, via __new__."""

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
