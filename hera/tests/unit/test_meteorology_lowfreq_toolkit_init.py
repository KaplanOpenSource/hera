"""lowFreqToolKit.__init__/docType. Probing this ad hoc (outside pytest)
hung -- the docstring says it "looks up 'IMSData' in the public database",
suggesting a real network attempt. Using unit_toolkit_factory here instead,
since conftest's _block_network autouse fixture turns any such attempt
into an immediate RuntimeError instead of a hang, so a construction that's
actually hermetic passes cleanly and one that isn't fails fast and loud."""
import pytest


@pytest.mark.unit
class TestLowFreqToolKitConstruction:
    def test_it_constructs_without_touching_the_network(self, unit_toolkit_factory):
        from hera import toolkitHome

        toolkit = unit_toolkit_factory(toolkitHome.METEOROLOGY_LOWFREQ)
        assert toolkit is not None

    def test_doc_type_is_set(self, unit_toolkit_factory):
        from hera import toolkitHome

        toolkit = unit_toolkit_factory(toolkitHome.METEOROLOGY_LOWFREQ)
        assert toolkit.docType is not None
