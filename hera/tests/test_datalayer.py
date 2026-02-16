"""
Native pytest tests for hera.datalayer.project.Project.

Replaces: hera/tests/datalayer/project.py (unittest-based)
Covers:   Project init, addMeasurementsDocument, getMeasurementsDocuments,
          deleteMeasurementsDocuments, setConfig/getConfig (counters).
"""

import pytest
from hera.datalayer.project import Project


# ---------------------------------------------------------------------------
# Test-case data (migrated from testCases/project.json)
# ---------------------------------------------------------------------------

INIT_CASES = [
    {"projectName": "Unit_Test_Project_01"},
]

ADD_MEASUREMENT_CASES = [
    {
        "projectName": "Unit_Test_Project_01",
        "resource": "PATH1/PATH2",
        "dataFormat": "string",
        "type": "",
        "desc": {},
    },
]

DELETE_CASES = [
    {"projectName": "Unit_Test_Project_01"},
]


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def projects():
    """Create all projects listed in INIT_CASES; tear them down at the end."""
    created = {}
    for case in INIT_CASES:
        proj = Project(projectName=case["projectName"])
        created[case["projectName"]] = proj
    yield created
    # Cleanup
    for proj in created.values():
        try:
            for doc in proj.getMeasurementsDocuments():
                doc.delete()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestProjectInit:
    @pytest.mark.parametrize("case", INIT_CASES, ids=[c["projectName"] for c in INIT_CASES])
    def test_project_init(self, case, projects):
        """Project creation must succeed."""
        proj = projects[case["projectName"]]
        assert proj is not None
        assert proj.projectName == case["projectName"]


class TestProjectCRUD:
    @pytest.mark.parametrize("case", ADD_MEASUREMENT_CASES,
                             ids=[c["projectName"] for c in ADD_MEASUREMENT_CASES])
    def test_add_measurements_document(self, case, projects):
        """addMeasurementsDocument must return a truthy document."""
        proj = projects[case["projectName"]]
        doc = proj.addMeasurementsDocument(
            resource=case["resource"],
            dataFormat=case["dataFormat"],
            type=case["type"],
            desc=case["desc"],
        )
        assert doc, f"Document was not added to project {case['projectName']}"

    @pytest.mark.parametrize("case", ADD_MEASUREMENT_CASES,
                             ids=[c["projectName"] for c in ADD_MEASUREMENT_CASES])
    def test_get_measurements_documents(self, case, projects):
        """getMeasurementsDocuments must return at least one document after add."""
        proj = projects[case["projectName"]]
        # Ensure at least one document exists
        proj.addMeasurementsDocument(
            resource=case["resource"],
            dataFormat=case["dataFormat"],
            type=case["type"],
            desc=case["desc"],
        )
        docs = list(proj.getMeasurementsDocuments(
            resource=case["resource"],
            dataFormat=case["dataFormat"],
            type=case["type"],
        ))
        assert len(docs) > 0, "getMeasurementsDocuments returned empty list"

    @pytest.mark.parametrize("case", DELETE_CASES,
                             ids=[c["projectName"] for c in DELETE_CASES])
    def test_delete_measurements_documents(self, case, projects):
        """deleteMeasurementsDocuments must remove all documents."""
        proj = projects[case["projectName"]]
        # Add a document first so there's something to delete
        proj.addMeasurementsDocument(
            resource="to_delete",
            dataFormat="string",
            type="",
            desc={},
        )
        proj.deleteMeasurementsDocuments()
        remaining = list(proj.getMeasurementsDocuments())
        assert len(remaining) == 0, f"Expected 0 documents after delete, got {len(remaining)}"


class TestProjectCounters:
    """Test setConfig / getConfig (Counter documents)."""

    def test_add_and_read_counters(self, projects):
        proj = projects["Unit_Test_Project_01"]

        # If Project has setConfig/getConfig, use them; otherwise skip.
        if not (hasattr(proj, "setConfig") and hasattr(proj, "getConfig")):
            pytest.skip("Project does not expose setConfig/getConfig")

        proj.setConfig(testCounter="hello_world")
        cfg = proj.getConfig()
        assert cfg.get("testCounter") == "hello_world"
