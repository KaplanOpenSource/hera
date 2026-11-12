"""riskassessment/riskToolkit.py: the analysis layer's LSM/datalayer
properties. getRiskAreas itself needs a full agent/effect/dispersion
setup and is out of scope here."""
import pytest


@pytest.mark.unit
class TestAnalysisProperties:
    def test_lsm_exposes_an_lsm_toolkit_for_the_same_project(self, unit_toolkit_factory):
        from hera import toolkitHome
        from hera.simulations.LSM.toolkit import LSMToolkit

        tk = unit_toolkit_factory(toolkitHome.RISKASSESSMENT)
        assert isinstance(tk.analysis.LSM, LSMToolkit)

    def test_datalayer_returns_the_owning_toolkit(self, unit_toolkit_factory):
        from hera import toolkitHome

        tk = unit_toolkit_factory(toolkitHome.RISKASSESSMENT)
        assert tk.analysis.datalayer is tk
