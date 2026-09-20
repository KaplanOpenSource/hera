"""Foam_parser_FieldDescription: the one CLI command with no toolkit/DB
dependency -- writes a JSON field-description template to disk.
"""
import json
from argparse import Namespace

import pytest

from hera.simulations.openFoam.CLI import Foam_parser_FieldDescription


@pytest.mark.unit
class TestFoamParserFieldDescription:
    def test_with_no_fields_it_writes_a_single_example_field(self, tmp_path):
        target = tmp_path / "fields.json"
        Foam_parser_FieldDescription(Namespace(fields=None, fileName=str(target)))
        data = json.loads(target.read_text())
        assert list(data.keys()) == ["exampleField"]
        assert data["exampleField"]["dimensions"] == "[0 0 0 0 0 0 0]"

    def test_with_explicit_fields_it_writes_one_entry_each(self, tmp_path):
        target = tmp_path / "fields.json"
        Foam_parser_FieldDescription(Namespace(fields=["U", "p"], fileName=str(target)))
        data = json.loads(target.read_text())
        assert set(data.keys()) == {"U", "p"}
        assert data["U"]["componentNames"] is None

    def test_an_empty_fields_list_falls_back_to_the_example(self, tmp_path):
        target = tmp_path / "fields.json"
        Foam_parser_FieldDescription(Namespace(fields=[], fileName=str(target)))
        data = json.loads(target.read_text())
        assert list(data.keys()) == ["exampleField"]
