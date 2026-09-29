"""datahandler: format registry, dispatch, and save/load round trips.

The registry decides which storage format a Python object gets, so it decides
whether hera's own policy in CLAUDE.md -- tabular data as parquet, never CSV --
is actually followed.  That makes dispatch a correctness concern.

Every round trip here writes to tmp_path, so nothing touches the test data set.
"""
import json

import numpy as np
import pandas as pd
import pytest

import hera.datalayer.datahandler as datahandler
from hera.datalayer.datahandler import datatypes

FORMAT_CONSTANTS = sorted(
    name for name in dir(datatypes) if name.isupper() and not name.startswith("_")
)

# The registry keys objects by their fully-qualified class path, which for
# third-party classes is an internal module path that moves between releases.
# Where the installed pandas exposes a path the map does not know, detection
# degrades to pickle -- see B20.  Computed at collection so the expectation is
# correct in both the pinned CI environment and a drifted local one.
_PANDAS_KEY_IS_KNOWN = (
    datatypes.get_obj_or_instance_fullName(pd.DataFrame())
    in datatypes.typeDatatypeMap
)


@pytest.mark.unit
class TestFormatRegistry:
    @pytest.mark.parametrize("constant", FORMAT_CONSTANTS)
    def test_every_declared_format_has_a_handler(self, constant):
        """A constant with no handler is a trap: it only fails at save time."""
        formatName = getattr(datatypes, constant)
        assert datatypes.getHandler(formatName) is not None

    @pytest.mark.parametrize("constant", FORMAT_CONSTANTS)
    def test_every_handler_exposes_both_halves(self, constant):
        handler = datatypes.getHandler(getattr(datatypes, constant))
        assert callable(getattr(handler, "getData", None))
        assert callable(getattr(handler, "saveData", None))

    def test_unknown_format_raises_and_names_the_type(self):
        with pytest.raises(ValueError, match="not known"):
            datatypes.getHandler("no_such_format")

    def test_handler_name_follows_the_convention(self):
        assert datatypes.getHandler(datatypes.STRING) is datahandler.DataHandler_string

    def test_module_level_helpers_agree_with_the_class(self):
        assert datahandler.getHandler(datatypes.STRING) is datatypes.getHandler(
            datatypes.STRING
        )

    def test_format_constants_are_unique(self):
        values = [getattr(datatypes, name) for name in FORMAT_CONSTANTS]
        assert len(values) == len(set(values))


@pytest.mark.unit
class TestFullyQualifiedNames:
    def test_builtin_types_are_unqualified(self):
        assert datatypes.get_obj_or_instance_fullName("text") == "str"
        assert datatypes.get_obj_or_instance_fullName(42) == "int"

    def test_a_class_and_an_instance_give_the_same_name(self):
        assert datatypes.get_obj_or_instance_fullName(
            pd.DataFrame
        ) == datatypes.get_obj_or_instance_fullName(pd.DataFrame())

    def test_third_party_types_are_module_qualified(self):
        assert datatypes.get_obj_or_instance_fullName(np.array([1])) == "numpy.ndarray"


@pytest.mark.unit
class TestFormatDetection:
    def test_a_string_is_stored_as_text(self):
        assert datatypes.getDataFormatName("text") == datatypes.STRING
        assert datatypes.getDataFormatExtension("text") == "txt"

    def test_a_numpy_array_is_stored_as_npy(self):
        assert datatypes.getDataFormatName(np.array([1])) == datatypes.NUMPY_ARRAY
        assert datatypes.getDataFormatExtension(np.array([1])) == "npy"

    def test_an_unknown_object_falls_back_to_pickle(self):
        """Documented fallback -- 'object' is a key in the map."""
        assert datatypes.getDataFormatName(object()) == datatypes.PICKLE

    def test_a_geodataframe_is_stored_as_geopandas(self):
        import geopandas as gpd

        assert datatypes.getDataFormatName(gpd.GeoDataFrame()) == datatypes.GEOPANDAS

    def test_guess_handler_matches_explicit_lookup(self):
        assert datatypes.guessHandler("text") is datatypes.getHandler(datatypes.STRING)

    @pytest.mark.xfail(
        condition=not _PANDAS_KEY_IS_KNOWN,
        strict=True,
        reason="B20: typeDatatypeMap is keyed on third-party internal module paths "
               "('pandas.core.frame.DataFrame'). pandas 3 reports "
               "'pandas.DataFrame', which the map does not know, so every "
               "DataFrame silently degrades to pickle instead of parquet -- "
               "against the policy in CLAUDE.md, with no error. "
               "See the consolidated findings issue.",
    )
    def test_a_dataframe_is_stored_as_parquet(self):
        """CLAUDE.md: tabular data is parquet, and parquet is preferred over CSV."""
        assert datatypes.getDataFormatName(pd.DataFrame({"a": [1]})) == datatypes.PARQUET
        assert datatypes.getDataFormatExtension(pd.DataFrame({"a": [1]})) == "parquet"

    @pytest.mark.xfail(
        condition=not _PANDAS_KEY_IS_KNOWN,
        strict=True,
        reason="B20: same stale-key problem for Series. "
               "See the consolidated findings issue.",
    )
    def test_a_series_is_stored_as_json(self):
        assert datatypes.getDataFormatName(pd.Series([1])) == datatypes.JSON_PANDAS


@pytest.mark.unit
class TestStringHandler:
    def test_get_returns_the_resource_itself(self):
        """The resource IS the data for this format -- no file is read."""
        assert datahandler.DataHandler_string.getData("hello") == "hello"

    def test_save_writes_the_text(self, tmp_path):
        target = tmp_path / "out.txt"
        datahandler.DataHandler_string.saveData("hello", str(target))
        assert target.read_text(encoding="utf-8") == "hello"

    def test_save_returns_a_dict_of_extra_metadata(self, tmp_path):
        result = datahandler.DataHandler_string.saveData(
            "hello", str(tmp_path / "out.txt")
        )
        assert result == {}


@pytest.mark.unit
class TestJsonDictHandler:
    def test_round_trip_through_a_file(self, tmp_path):
        target = tmp_path / "conf.json"
        payload = {"a": 1, "nested": {"b": [1, 2]}}

        datahandler.DataHandler_JSON_dict.saveData(payload, str(target))
        assert datahandler.DataHandler_JSON_dict.getData(str(target)) == payload

    def test_saved_file_is_valid_json(self, tmp_path):
        target = tmp_path / "conf.json"
        datahandler.DataHandler_JSON_dict.saveData({"a": 1}, str(target))
        assert json.loads(target.read_text(encoding="utf-8")) == {"a": 1}

    def test_get_accepts_a_json_string_directly(self):
        """loadJSON takes a string, a path or a dict, so all three must work."""
        assert datahandler.DataHandler_JSON_dict.getData('{"a": 1}') == {"a": 1}

    def test_get_accepts_a_dict_unchanged(self):
        assert datahandler.DataHandler_JSON_dict.getData({"a": 1}) == {"a": 1}


@pytest.mark.unit
class TestParquetHandler:
    def test_round_trip_preserves_the_frame(self, tmp_path):
        target = tmp_path / "table.parquet"
        frame = pd.DataFrame({"a": [1, 2, 3], "b": ["x", "y", "z"]})

        datahandler.DataHandler_parquet.saveData(frame, str(target))
        restored = datahandler.DataHandler_parquet.getData(str(target))
        if hasattr(restored, "compute"):
            restored = restored.compute()

        # check_dtype=False on purpose: a text column's exact dtype depends on
        # the pandas version (object under the pinned 2.2.3, str/string under
        # pandas 3), which says nothing about hera.  The dtypes that matter --
        # numeric and boolean -- get their own test below.
        pd.testing.assert_frame_equal(
            restored.reset_index(drop=True),
            frame.reset_index(drop=True),
            check_dtype=False,
        )

    def test_round_trip_preserves_numeric_dtypes(self, tmp_path):
        """The reason parquet is preferred over CSV: types survive."""
        target = tmp_path / "table.parquet"
        frame = pd.DataFrame({"i": [1, 2], "f": [1.5, 2.5], "b": [True, False]})

        datahandler.DataHandler_parquet.saveData(frame, str(target))
        restored = datahandler.DataHandler_parquet.getData(str(target))
        if hasattr(restored, "compute"):
            restored = restored.compute()

        assert restored["i"].dtype == frame["i"].dtype
        assert restored["f"].dtype == frame["f"].dtype
        assert restored["b"].dtype == frame["b"].dtype


@pytest.mark.unit
class TestNumpyArrayHandler:
    def test_round_trip_preserves_values_and_shape(self, tmp_path):
        target = tmp_path / "data.npy"
        array = np.arange(6).reshape(2, 3)

        datahandler.DataHandler_numpy_array.saveData(array, str(target))
        restored = datahandler.DataHandler_numpy_array.getData(str(target))

        assert restored.shape == array.shape
        assert np.array_equal(restored, array)

    def test_round_trip_preserves_dtype(self, tmp_path):
        target = tmp_path / "data.npy"
        array = np.arange(3, dtype=np.float32)

        datahandler.DataHandler_numpy_array.saveData(array, str(target))
        assert datahandler.DataHandler_numpy_array.getData(str(target)).dtype == np.float32


@pytest.mark.unit
class TestPickleHandler:
    def test_round_trip_of_a_plain_dict(self, tmp_path):
        target = tmp_path / "data.pckle"
        payload = {"a": 1, "b": [1, 2]}

        datahandler.DataHandler_pickle.saveData(payload, str(target))
        assert datahandler.DataHandler_pickle.getData(str(target)) == payload

    def test_round_trip_of_a_list(self, tmp_path):
        target = tmp_path / "data.pckle"
        datahandler.DataHandler_pickle.saveData([1, "two", 3.0], str(target))
        assert datahandler.DataHandler_pickle.getData(str(target)) == [1, "two", 3.0]
