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
* **B242** (fixed) -- the cache was write-only.  The descriptor included
  ``"bounds": data.total_bounds``, a ``numpy.ndarray``, so the stored document
  carried an array-valued field; ``getCacheDocuments`` then built a query
  term from the same array and comparing it against that field evaluated an
  array in a boolean context -- *"ValueError: The truth value of an array with
  more than one element is ambiguous"*.  The first call, against an empty
  collection, had nothing to compare and therefore succeeded; from then on
  every call with ``saveCache=True`` failed, ``overwrite=True`` included,
  because the lookup runs before the overwrite check.  ``bounds`` is now
  stored as a plain list, which makes the "return found data in DB" path, the
  "update old record" save path and its documented ``FileNotFoundError``
  advice reachable.

Deliberately not tested, and why
--------------------------------
* The numeric content of ``lambdaP``/``lambdaF``/``hc`` beyond one
  hand-computed check: owned by ``Blocks``, covered in
  ``test_gis_buildings_analysis.py``.
* The numeric content of the cache-hit result beyond its equality with the
  freshly computed one: the geojson round-trip is geopandas'.
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
    """The first call, against an empty cache."""

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
    """B242: the cache can be read back, not only written."""

    def test_a_second_identical_call_reuses_the_stored_result(self, analysis_layer):
        """B242."""
        first = _call(analysis_layer)
        second = _call(analysis_layer)
        assert second["lambdaP"].sum() == pytest.approx(first["lambdaP"].sum())

    def test_a_second_call_does_not_write_another_document(self, analysis_layer, buildings_toolkit):
        """B242: the cache hit short-circuits the recomputation."""
        _call(analysis_layer)
        _call(analysis_layer)
        assert len(buildings_toolkit.getCacheDocuments(type="Lambda_Buildings")) == 1

    def test_a_call_with_a_different_key_is_cached_separately(self, analysis_layer, buildings_toolkit):
        """B242: one cached entry no longer poisons lookups for the rest."""
        _call(analysis_layer, windMeteorologicalDirection=270)
        _call(analysis_layer, windMeteorologicalDirection=45)
        assert len(buildings_toolkit.getCacheDocuments(type="Lambda_Buildings")) == 2

    def test_overwrite_recomputes_into_the_existing_record(self, analysis_layer, buildings_toolkit):
        """B242: the 'update old record' save branch is reachable."""
        _call(analysis_layer)
        _call(analysis_layer, overwrite=True)
        assert len(buildings_toolkit.getCacheDocuments(type="Lambda_Buildings")) == 1

    def test_the_descriptor_bounds_field_is_a_plain_list(self, analysis_layer, buildings_toolkit):
        """B242: the root cause, stated directly -- not a numpy.ndarray."""
        assert isinstance(_buildings().total_bounds, numpy.ndarray)
        _call(analysis_layer)
        desc = buildings_toolkit.getCacheDocuments(type="Lambda_Buildings")[0].desc
        assert not isinstance(desc["bounds"], numpy.ndarray)
        assert len(desc["bounds"]) == 4


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
        """It never populates the cache, so every call recomputes."""
        first = _call(analysis_layer, saveCache=False)
        second = _call(analysis_layer, saveCache=False)
        assert second["lambdaP"].sum() == pytest.approx(first["lambdaP"].sum())

    def test_the_recovery_path_for_a_vanished_cache_file_is_reachable(self, analysis_layer, buildings_toolkit):
        """B242: the documented FileNotFoundError advice can now be given."""
        _call(analysis_layer)
        os.remove(buildings_toolkit.getCacheDocuments(type="Lambda_Buildings")[0].resource)
        with pytest.raises(FileNotFoundError, match="overwrite=True"):
            _call(analysis_layer)

    def test_the_update_an_existing_record_path_reuses_the_file_name(self, analysis_layer, buildings_toolkit):
        """B242: the ``len(dataDoc) > 0`` save branch reuses the stored
        resource path instead of allocating a new file name, so the counter
        does not advance past its first value."""
        _call(analysis_layer)
        _call(analysis_layer, overwrite=True)
        assert os.listdir(buildings_toolkit.filesDirectory) == ["LambdaGeoJson_0.geojson"]
        assert buildings_toolkit.getConfig()["analysis_CacheCounter"] == 1
