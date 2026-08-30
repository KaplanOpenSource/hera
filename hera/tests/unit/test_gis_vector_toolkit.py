"""VectorToolkit: GeoJSON export and the region-cutting dispatch.

``cutRegionFromSource`` is tested against a hand-built fake datasource
document (just `.desc['desc']['crs']` and `.getData(**kwargs)`) rather than
a real registered shapefile datasource, so the CRS-matching and
bbox/mask dispatch logic can be verified without file I/O.

B80: the CRS-mismatch reprojection is silently discarded. ``dct`` --
the dict eventually passed to ``getData(**dct)`` -- is built from
``regionWithCRS`` *before* the CRS check runs. When the shape's CRS
differs from the datasource's, the mismatch branch does
``regionWithCRS = regionWithCRS.to_crs(...)``, which returns a *new*
object rather than mutating in place -- so the reassignment never
reaches ``dct``, which still holds the original, un-reprojected shape.
(The "shape has no CRS" branch happens to work, because it assigns
``regionWithCRS.crs = ...`` directly onto the object ``dct`` already
references, mutating it in place -- the one branch that actually needed
a transform is the one that's broken.)
"""
import geopandas
import pytest
from shapely.geometry import box

from hera.measurements.GIS.vector.toolkit import VectorToolkit


class _FakeDatasourceDoc:
    def __init__(self, crs):
        self.desc = {"desc": {"crs": crs}}
        self.calls = []

    def getData(self, **kwargs):
        self.calls.append(kwargs)
        return "FAKE_RESULT"


@pytest.mark.unit
class TestGeopandasToGeoJson:
    def test_a_geodataframe_becomes_a_geojson_feature_collection_string(self):
        gdf = geopandas.GeoDataFrame({"geometry": [box(0, 0, 1, 1)]}, crs=4326)
        result = VectorToolkit.geopandasToGeoJson(gdf)
        assert isinstance(result, str)
        assert '"type": "FeatureCollection"' in result

    def test_anything_else_raises(self):
        with pytest.raises(ValueError, match="only GeoDataFrame"):
            VectorToolkit.geopandasToGeoJson({"not": "a geodataframe"})

    def test_a_plain_list_also_raises(self):
        with pytest.raises(ValueError, match="only GeoDataFrame"):
            VectorToolkit.geopandasToGeoJson([1, 2, 3])


@pytest.mark.unit
class TestCutRegionFromSource:
    @pytest.fixture()
    def toolkit(self):
        return VectorToolkit.__new__(VectorToolkit)

    def test_isbounds_true_calls_getdata_with_bbox(self, toolkit):
        doc = _FakeDatasourceDoc(crs=2039)
        toolkit.cutRegionFromSource(doc, shape=[0, 0, 10, 10], isBounds=True, inputCRS=4326)
        assert list(doc.calls[0].keys()) == ["bbox"]

    def test_isbounds_false_calls_getdata_with_mask(self, toolkit):
        doc = _FakeDatasourceDoc(crs=2039)
        toolkit.cutRegionFromSource(doc, shape=[0, 0, 10, 10], isBounds=False, inputCRS=4326)
        assert list(doc.calls[0].keys()) == ["mask"]

    @pytest.mark.xfail(
        strict=True,
        reason="B80: dct (the dict handed to getData) is built from "
               "regionWithCRS before the CRS-mismatch branch runs. That "
               "branch does regionWithCRS = regionWithCRS.to_crs(...), a "
               "non-mutating reassignment, so the reprojection never "
               "reaches dct -- getData() receives the original CRS. "
               "See the consolidated findings issue.",
    )
    def test_the_shape_should_be_reprojected_to_the_datasource_crs(self, toolkit):
        doc = _FakeDatasourceDoc(crs=2039)
        toolkit.cutRegionFromSource(doc, shape=[0, 0, 10, 10], isBounds=True, inputCRS=4326)
        reprojected = doc.calls[0]["bbox"]
        assert reprojected.crs.to_epsg() == 2039

    def test_the_shape_currently_keeps_its_original_crs_despite_the_mismatch(self, toolkit):
        """Characterisation of B80."""
        doc = _FakeDatasourceDoc(crs=2039)
        toolkit.cutRegionFromSource(doc, shape=[0, 0, 10, 10], isBounds=True, inputCRS=4326)
        unreprojected = doc.calls[0]["bbox"]
        assert unreprojected.crs.to_epsg() == 4326

    def test_when_the_shape_has_no_crs_it_does_pick_up_the_datasource_crs(self, toolkit):
        """The other branch (bare .crs = ... assignment mutates in place),
        so this one actually works -- confirms the bug is specific to the
        reprojection branch, not the whole CRS-handling block."""
        doc = _FakeDatasourceDoc(crs=2039)
        toolkit.cutRegionFromSource(doc, shape=[0, 0, 10, 10], isBounds=True, inputCRS=None)
        result = doc.calls[0]["bbox"]
        assert result.crs.to_epsg() == 2039

    def test_it_returns_whatever_the_datasource_getdata_returns(self, toolkit):
        doc = _FakeDatasourceDoc(crs=2039)
        result = toolkit.cutRegionFromSource(doc, shape=[0, 0, 10, 10], isBounds=True, inputCRS=4326)
        assert result == "FAKE_RESULT"
