"""GIS/raster/landcover.py: getCodingMap and roughnesslength2sandgrainroughness.
Constructed via unit_toolkit_factory -- an ad hoc construction outside
pytest hung, presumably a real lookup attempt that mongomock/_block_network
intercepts cleanly inside the fixture but not in a bare script."""
import pytest


@pytest.fixture()
def lc(unit_toolkit_factory):
    from hera import toolkitHome

    return unit_toolkit_factory(toolkitHome.GIS_LANDCOVER)


@pytest.mark.unit
class TestGetCodingMap:
    def test_type_1_maps_known_codes_to_names(self, lc):
        mapping = lc.getCodingMap("Type-1")
        assert mapping[0] == "Water"
        assert mapping[13] == "Urban and built-up"

    def test_an_unknown_datasource_returns_an_empty_map(self, lc):
        assert lc.getCodingMap("NoSuchType") == {}


@pytest.mark.unit
class TestRoughnessLengthToSandGrainRoughness:
    def test_it_scales_by_30(self, lc):
        assert lc.roughnesslength2sandgrainroughness(0.1) == pytest.approx(3.0)
