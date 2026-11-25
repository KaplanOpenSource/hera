"""datalayer/project.py: export/load round trip, the `all` property, and
AbstractCollection.deleteDocumentByID.

addDocumentFromJSON is left uncovered: probing it showed inconsistent
mongomock behavior (.save() reports success and a new id, but the
document is not found immediately after by that id) that needs isolating
from a possible mongomock quirk before it can be asserted as a real
hera defect -- out of scope to chase down here.
"""
import os

import pandas
import pytest

from hera.datalayer import Project


@pytest.mark.unit
class TestAllProperty:
    def test_all_exposes_the_abstract_collection(self, unit_project):
        assert unit_project.all is not None


@pytest.mark.unit
class TestDeleteDocumentByID:
    def test_it_removes_the_document(self, unit_project):
        doc = unit_project.saveMeasurementData("m1", pandas.DataFrame({"a": [1]}), desc={})
        unit_project.all.deleteDocumentByID(doc.id)
        assert list(unit_project.getMeasurementsDocuments(type="m1")) == []


@pytest.mark.unit
class TestExportAndLoad:
    def test_export_on_the_default_project_is_a_no_op(self, unit_project, tmp_path):
        default_proj = Project(projectName=unit_project.DEFAULTPROJECT, filesDirectory=str(tmp_path))
        target = tmp_path / "export.zip"
        default_proj.export(str(target))
        assert not target.exists()

    def test_export_then_load_round_trips_measurement_data(self, unit_project, tmp_path):
        unit_project.saveMeasurementData("m1", pandas.DataFrame({"a": [1]}), desc={"k": "v"})
        target = tmp_path / "export.zip"
        unit_project.export(str(target))
        assert target.exists()

        other = Project(projectName="IMPORT_TARGET_PROJECT", filesDirectory=str(tmp_path / "other"))
        Project.load(other, str(target), is_hard_import=False)
        docs = other.getMeasurementsDocuments(type="m1")
        assert len(docs) == 1
        assert docs[0].desc["k"] == "v"
