"""hera/measurements/GIS/vector/buildings/analysis.py::analysis.LambdaFromBuildingData
-- the one function in that module ``test_gis_buildings_analysis.py`` leaves
uncovered.

It is the cached front end to the ``Blocks`` pipeline that file already tests:
validate the input, build a cache key from the data's bounds/wind/resolution/
CRS, and either return the stored GeoJSON or compute the urban-canopy
parameters and store them.

How the inputs were produced
----------------------------
``_buildings`` places one 20 m x 20 m footprint inside each 100 m block of a
200 m domain, exactly as ``test_gis_buildings_analysis.py`` does, so the
lambda values coming back have hand-computable expectations
(``lambdaP = 400/10000``) and this file can concentrate on the caching and
validation the wrapper owns rather than re-testing the morphology maths.  The
CRS is set explicitly, because the method requires it.

Bugs pinned (each has an xfail(strict) for the intended behaviour plus a
passing characterisation of what actually happens)
--------------------------------------------------------------------------
* **B242** -- the cache is write-only.  The descriptor includes
  ``"bounds": data.total_bounds``, a ``numpy.ndarray``, so the stored document
  carries an array-valued field; ``getCacheDocuments`` then builds a query
  term from the same array and comparing it against that field evaluates an
  array in a boolean context -- *"ValueError: The truth value of an array with
  more than one element is ambiguous"*.  The first call, against an empty
  collection, has nothing to compare and therefore succeeds; from then on
  every call with ``saveCache=True`` fails, ``overwrite=True`` included,
  because the lookup runs before the overwrite check.  Two whole branches --
  the "return found data in DB" path and the "update old record" save path,
  along with its documented ``FileNotFoundError`` advice -- are consequently
  dead code.

Deliberately not tested, and why
--------------------------------
* The numeric content of ``lambdaP``/``lambdaF``/``hc`` beyond one
  hand-computed check: owned by ``Blocks``, covered in
  ``test_gis_buildings_analysis.py``.
* The behaviour of the cache-hit and record-update branches: unreachable
  while B242 stands, so the tests below assert that they are unreachable
  rather than asserting on dead code.
* ``analysis.__init__``'s config bootstrap: exercised by every fixture in this
  file (it runs on construction) and asserted directly in the first class;
  the rest is left to ``test_gis_buildings_analysis.py``.
"""
import os

import geopandas
import numpy
import pytest
from shapely.geometry import box

from hera import toolkitHome
from hera.measurements.GIS.utils import ITM
from hera.measurements.GIS.vector.buildings.analysis import (
    BUILDINGS_LAMBDA_RESOLUTION,
    BUILDINGS_LAMBDA_WIND_DIRECTION,
)

_RESOLUTION = 100
_WIND = 270


def _buildings(crs=ITM):
    """One 20x20 footprint in each 100 m block of a 200 m domain."""
    return geopandas.GeoDataFrame(
        {
            "geometry": [box(10, 10, 30, 30), box(150, 150, 170, 170)],
            "BLDG_HT": [10.0, 15.0],
            "HT_LAND": [5.0, 5.0],
            "FTYPE": [1, 1],
            "HI_PNT_Z": [15.0, 20.0],
        },
        crs=crs,
    )


@pytest.fixture()
def buildings_toolkit(unit_toolkit_factory):
    return unit_toolkit_factory(toolkitHome.GIS_BUILDINGS)


@pytest.fixture()
def analysis_layer(buildings_toolkit):
    return buildings_toolkit.analysis


def _call(layer, data=None, **kwargs):
    kwargs.setdefault("windMeteorologicalDirection", _WIND)
    kwargs.setdefault("resolution", _RESOLUTION)
    return layer.LambdaFromBuildingData(
        buildingsData=_buildings() if data is None else data, **kwargs
    )


@pytest.mark.unit
class TestConfigBootstrap:
    def test_constructing_the_layer_seeds_the_cache_counter(self, buildings_toolkit):
        assert buildings_toolkit.getConfig()["analysis_CacheCounter"] == 0

    def test_an_existing_counter_is_not_reset(self, buildings_toolkit):
        from hera.measurements.GIS.vector.buildings.analysis import analysis

        buildings_toolkit.setConfig(analysis_CacheCounter=7)
        analysis(dataLayer=buildings_toolkit)
        assert buildings_toolkit.getConfig()["analysis_CacheCounter"] == 7


@pytest.mark.unit
class TestLambdaFromBuildingDataValidation:
    """The validation runs before the cache lookup, so B242 does not mask it."""

    def test_neither_a_string_nor_a_geodataframe_is_rejected(self, analysis_layer):
        with pytest.raises(ValueError, match="must be str or geopandas"):
            _call(analysis_layer, data=12345)

    def test_the_rejection_names_the_type_it_was_given(self, analysis_layer):
        with pytest.raises(ValueError, match="int"):
            _call(analysis_layer, data=12345)

    def test_a_plain_pandas_frame_is_rejected_too(self, analysis_layer):
        import pandas

        with pytest.raises(ValueError, match="must be str or geopandas"):
            _call(analysis_layer, data=pandas.DataFrame({"a": [1]}))

    def test_a_geodataframe_without_a_crs_is_rejected(self, analysis_layer):
        with pytest.raises(ValueError, match="crs must be set"):
            _call(analysis_layer, data=_buildings(crs=None))

    def test_a_geojson_string_is_accepted_and_carries_its_own_crs(self, analysis_layer):
        """readGeoJSONString gives the frame a CRS, so validation lets it through."""
        result = _call(analysis_layer, data=_buildings().to_json())
        assert "lambdaP" in result.columns

    def test_a_string_that_is_not_geojson_is_rejected_by_the_parser(self, analysis_layer):
        with pytest.raises(ValueError):
            _call(analysis_layer, data="not geojson at all")


@pytest.mark.unit
class TestLambdaFromBuildingDataComputation:
    def test_it_returns_the_urban_canopy_parameters(self, analysis_layer):
        result = _call(analysis_layer)
        assert {"lambdaP", "lambdaF", "hc"} <= set(result.columns)

    def test_one_row_per_block_of_the_domain(self, analysis_layer):
        """A 200 m span of footprints at 100 m resolution gives a 2x2 grid."""
        assert len(_call(analysis_layer)) == 4

    def test_the_plan_area_fraction_matches_the_hand_computation(self, analysis_layer):
        """A 20x20 footprint (area 400) in a 100x100 block."""
        result = _call(analysis_layer)
        occupied = result[(result["i0"] == 0) & (result["j0"] == 0)].iloc[0]
        assert occupied["lambdaP"] == pytest.approx(400 / 10000)

    def test_a_coarser_resolution_gives_a_single_block(self, analysis_layer):
        """The resolution argument reaches the Blocks constructor."""
        assert len(_call(analysis_layer, resolution=200, saveCache=False)) == 1

    def test_the_wind_direction_reaches_the_frontal_area_calculation(self, analysis_layer):
        """A square footprint's frontal area is invariant under 90-degree
        rotations, so 270 and 360 must agree while the argument still lands."""
        at_270 = _call(analysis_layer, windMeteorologicalDirection=270, saveCache=False)
        at_360 = _call(analysis_layer, windMeteorologicalDirection=360, saveCache=False)
        assert at_270["lambdaF"].sum() == pytest.approx(at_360["lambdaF"].sum())


@pytest.mark.unit
class TestLambdaFromBuildingDataCacheWrite:
    """The first call, against an empty cache -- the only one B242 permits."""

    def test_the_result_is_stored_as_a_cache_document(self, analysis_layer, buildings_toolkit):
        _call(analysis_layer)
        assert len(buildings_toolkit.getCacheDocuments(type="Lambda_Buildings")) == 1

    def test_the_cache_key_records_bounds_wind_resolution_and_crs(self, analysis_layer, buildings_toolkit):
        _call(analysis_layer)
        desc = buildings_toolkit.getCacheDocuments(type="Lambda_Buildings")[0].desc
        assert desc[BUILDINGS_LAMBDA_WIND_DIRECTION] == _WIND
        assert desc[BUILDINGS_LAMBDA_RESOLUTION] == _RESOLUTION
        assert desc["crs"] == ITM
        assert list(desc["bounds"]) == pytest.approx(list(_buildings().total_bounds))

    def test_the_geojson_is_written_into_the_files_directory(self, analysis_layer, buildings_toolkit):
        _call(analysis_layer)
        expected = os.path.join(buildings_toolkit.filesDirectory, "LambdaGeoJson_0.geojson")
        assert os.path.exists(expected)

    def test_the_document_points_at_the_file_that_was_written(self, analysis_layer, buildings_toolkit):
        _call(analysis_layer)
        doc = buildings_toolkit.getCacheDocuments(type="Lambda_Buildings")[0]
        assert doc.resource == os.path.join(
            buildings_toolkit.filesDirectory, "LambdaGeoJson_0.geojson"
        )
        assert os.path.exists(doc.resource)

    def test_the_file_name_comes_from_the_config_counter(self, analysis_layer, buildings_toolkit):
        buildings_toolkit.setConfig(analysis_CacheCounter=41)
        _call(analysis_layer)
        assert os.listdir(buildings_toolkit.filesDirectory) == ["LambdaGeoJson_41.geojson"]

    def test_the_counter_is_advanced_for_the_next_caller(self, analysis_layer, buildings_toolkit):
        _call(analysis_layer)
        assert buildings_toolkit.getConfig()["analysis_CacheCounter"] == 1

    def test_the_stored_file_holds_the_same_blocks_that_were_returned(self, analysis_layer, buildings_toolkit):
        result = _call(analysis_layer)
        written = geopandas.read_file(
            buildings_toolkit.getCacheDocuments(type="Lambda_Buildings")[0].resource
        )
        assert len(written) == len(result)
        assert written["lambdaP"].sum() == pytest.approx(result["lambdaP"].sum())


@pytest.mark.unit
class TestLambdaFromBuildingDataCacheRead:
    """B242: the cache is write-only -- reading it back always raises."""

    @pytest.mark.xfail(
        strict=True,
        reason="B242: the cache descriptor carries \"bounds\": "
               "data.total_bounds, a numpy.ndarray, so the document is stored "
               "with an array-valued field.  getCacheDocuments then turns the "
               "same array into a query term and the comparison against that "
               "field evaluates an array in a boolean context -- "
               "'ValueError: The truth value of an array with more than one "
               "element is ambiguous'.  The lookup on line 147 runs before "
               "the overwrite check, so once anything at all has been cached "
               "under type='Lambda_Buildings' EVERY subsequent call fails, "
               "overwrite=True included.  The cache can be written but never "
               "read. See the consolidated findings issue.",
    )
    def test_a_second_identical_call_reuses_the_stored_result(self, analysis_layer):
        first = _call(analysis_layer)
        second = _call(analysis_layer)
        assert second["lambdaP"].sum() == pytest.approx(first["lambdaP"].sum())

    def test_a_second_call_raises_on_the_array_valued_cache_key(self, analysis_layer):
        """Characterisation of B242."""
        _call(analysis_layer)
        with pytest.raises(ValueError, match="truth value of an array"):
            _call(analysis_layer)

    def test_even_a_call_with_a_different_key_is_poisoned(self, analysis_layer):
        """Characterisation of B242: the query scans every Lambda_Buildings
        document, so one cached entry breaks lookups for all of them."""
        _call(analysis_layer, windMeteorologicalDirection=270)
        with pytest.raises(ValueError, match="truth value of an array"):
            _call(analysis_layer, windMeteorologicalDirection=45)

    def test_overwrite_does_not_help(self, analysis_layer):
        """Characterisation of B242: the lookup precedes the overwrite test."""
        _call(analysis_layer)
        with pytest.raises(ValueError, match="truth value of an array"):
            _call(analysis_layer, overwrite=True)

    def test_the_descriptor_field_that_causes_it_is_the_bounds_array(self):
        """Characterisation of B242: the root cause, stated directly."""
        assert isinstance(_buildings().total_bounds, numpy.ndarray)
        assert _buildings().total_bounds.size == 4

    def test_the_same_key_with_bounds_as_a_list_queries_cleanly(self, buildings_toolkit):
        """Characterisation of B242: nothing else in the key is at fault.

        The identical four descriptor fields round-trip through write and
        lookup without complaint as soon as ``bounds`` is a plain list, which
        places the defect exactly at ``total_bounds`` and nowhere else.
        """
        data = _buildings()
        desc = {
            "bounds": list(data.total_bounds),
            BUILDINGS_LAMBDA_WIND_DIRECTION: _WIND,
            BUILDINGS_LAMBDA_RESOLUTION: _RESOLUTION,
            "crs": data.crs.to_epsg(),
        }
        buildings_toolkit.addCacheDocument(
            resource="", dataFormat="geopandas", type="Lambda_Buildings", desc=desc
        )
        assert len(buildings_toolkit.getCacheDocuments(type="Lambda_Buildings", **desc)) == 1


@pytest.mark.unit
class TestLambdaFromBuildingDataCacheControl:
    def test_save_cache_false_computes_without_storing_anything(self, analysis_layer, buildings_toolkit):
        result = _call(analysis_layer, saveCache=False)
        assert "lambdaP" in result.columns
        assert len(buildings_toolkit.getCacheDocuments(type="Lambda_Buildings")) == 0
        assert os.listdir(buildings_toolkit.filesDirectory) == []

    def test_save_cache_false_leaves_the_counter_alone(self, analysis_layer, buildings_toolkit):
        _call(analysis_layer, saveCache=False)
        assert buildings_toolkit.getConfig()["analysis_CacheCounter"] == 0

    def test_save_cache_false_is_repeatable(self, analysis_layer):
        """Because it never populates the cache, it never trips B242 either --
        the only way to call this method more than once today."""
        first = _call(analysis_layer, saveCache=False)
        second = _call(analysis_layer, saveCache=False)
        assert second["lambdaP"].sum() == pytest.approx(first["lambdaP"].sum())

    def test_the_recovery_path_for_a_vanished_cache_file_is_unreachable(self, analysis_layer, buildings_toolkit):
        """Characterisation of B242: the FileNotFoundError branch is dead code.

        ``LambdaFromBuildingData`` promises a clear "the cached data ... is
        not found on the disk ... use overwrite=True" error, but reaching it
        requires a successful cache lookup, which B242 makes impossible: the
        ValueError arrives first, with no such advice.
        """
        _call(analysis_layer)
        os.remove(buildings_toolkit.getCacheDocuments(type="Lambda_Buildings")[0].resource)
        with pytest.raises(ValueError, match="truth value of an array"):
            _call(analysis_layer)

    def test_the_update_an_existing_record_path_is_unreachable_too(self, analysis_layer, buildings_toolkit):
        """Characterisation of B242: the ``len(dataDoc) > 0`` save branch --
        the one that reuses the stored resource path instead of allocating a
        new file name -- can never run, so a second document is never
        written and the counter never advances past its first value."""
        _call(analysis_layer)
        with pytest.raises(ValueError):
            _call(analysis_layer, overwrite=True)
        assert os.listdir(buildings_toolkit.filesDirectory) == ["LambdaGeoJson_0.geojson"]
        assert buildings_toolkit.getConfig()["analysis_CacheCounter"] == 1
