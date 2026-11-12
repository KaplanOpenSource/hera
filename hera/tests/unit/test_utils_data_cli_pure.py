"""utils/data/CLI.py: the two pure helpers that don't need the lazily
imported toolkit/DB machinery -- _parse_query_value and load_project_name.
The ~26 argparse command handlers (project_list, toolkit_register, etc.)
all touch a real Project/toolkit via _setup()'s deferred imports and are
left uncovered; they'd need a much larger DB-backed harness."""
import json

import pytest

from hera.utils.data.CLI import _parse_query_value, load_project_name


@pytest.mark.unit
class TestParseQueryValue:
    @pytest.mark.parametrize("raw,expected", [
        ("none", None), ("None", None), ("NONE", None),
        ("true", True), ("True", True),
        ("false", False), ("False", False),
    ])
    def test_recognized_keywords_are_case_insensitive(self, raw, expected):
        assert _parse_query_value(raw) is expected

    def test_an_integer_string_becomes_an_int(self):
        assert _parse_query_value("42") == 42
        assert isinstance(_parse_query_value("42"), int)

    def test_a_float_string_becomes_a_float(self):
        assert _parse_query_value("3.14") == pytest.approx(3.14)
        assert isinstance(_parse_query_value("3.14"), float)

    def test_a_quoted_string_has_its_quotes_stripped(self):
        assert _parse_query_value('"hello"') == "hello"
        assert _parse_query_value("'hello'") == "hello"

    def test_mismatched_quotes_are_left_untouched(self):
        assert _parse_query_value('"hello\'') == '"hello\''

    def test_a_bare_word_falls_back_to_the_raw_string(self):
        assert _parse_query_value("myValue") == "myValue"

    def test_surrounding_whitespace_is_stripped(self):
        assert _parse_query_value("  42  ") == 42


@pytest.mark.unit
class TestLoadProjectName:
    def test_it_reads_the_project_name_from_the_config_file(self, tmp_path):
        config = tmp_path / "caseConfiguration.json"
        config.write_text(json.dumps({"projectName": "MY_PROJECT"}))
        assert load_project_name(str(config)) == "MY_PROJECT"

    def test_a_missing_file_raises_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="Could not find"):
            load_project_name(str(tmp_path / "noSuchFile.json"))

    def test_a_config_without_projectname_raises_value_error(self, tmp_path):
        config = tmp_path / "caseConfiguration.json"
        config.write_text(json.dumps({"other": "field"}))
        with pytest.raises(ValueError, match="does not contain 'projectName'"):
            load_project_name(str(config))

    def test_an_empty_projectname_also_raises(self, tmp_path):
        config = tmp_path / "caseConfiguration.json"
        config.write_text(json.dumps({"projectName": ""}))
        with pytest.raises(ValueError, match="does not contain 'projectName'"):
            load_project_name(str(config))
