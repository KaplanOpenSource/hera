"""RiskToolkit: the agent-loading and protection-policy entry point.

Constructing it also constructs an ``LSMToolkit``, so it needs the same
mongomock seam as every other toolkit test. ``analysis.getRiskAreas`` chains
into a live LSM simulation and a full agent dose-response calculation --
out of scope here; it needs LSM/S3 data and belongs with an integration
test. ``loadData`` turned out to have three real defects around its
``saveMode`` argument, found by exercising the documented modes directly:

* B55: ``TOOLKIT_SAVEMODE_NOSAVE`` only has an effect when a matching
  document already exists (the ``elif`` branch). The very first load of a
  new agent goes through the ``if agentDoc is None`` branch, which calls
  ``addDataSource`` unconditionally -- so "no save" still saves.
* B56: ``TOOLKIT_SAVEMODE_FILEANDDB_REPLACE`` looks up the existing document
  by the *new* version being loaded, not by name alone. Replacing with a
  different version therefore never finds the old entry, falls into the
  "create new" branch, and leaves a duplicate datasource under the same
  name instead of replacing anything.
* B57: even when the version matches and the replace path is taken, only
  ``desc['version']`` is updated -- the rest of ``desc`` (e.g.
  ``effectParameters``) is left stale while ``resource`` (the raw JSON) is
  fully overwritit, so the document's metadata and payload diverge.
"""
import pytest

from hera.riskassessment.agents.Agents import Agent
from hera.riskassessment.presentation.casualtiesFigs import casualtiesPlot
from hera.riskassessment.protectionpolicy.ProtectionPolicy import ProtectionPolicy
from hera.riskassessment.riskToolkit import RiskToolkit
from hera.toolkit import (
    TOOLKIT_SAVEMODE_FILEANDDB,
    TOOLKIT_SAVEMODE_FILEANDDB_REPLACE,
    TOOLKIT_SAVEMODE_NOSAVE,
)


@pytest.fixture()
def risk(unit_project):
    return RiskToolkit(projectName=unit_project.projectName)


def _desc(name="TestAgent", version=(0, 0, 1), **extra):
    d = {"name": name, "version": list(version), "effectParameters": {}, "effects": {}}
    d.update(extra)
    return d


@pytest.mark.unit
class TestConstruction:
    def test_the_protection_policy_property_is_the_class_itself(self, risk):
        assert risk.ProtectionPolicy is ProtectionPolicy

    def test_the_presentation_layer_is_a_casualties_plot(self, risk):
        assert isinstance(risk.presentation, casualtiesPlot)

    def test_the_analysis_layer_points_back_at_the_toolkit(self, risk):
        assert risk.analysis.datalayer is risk

    def test_a_fresh_project_has_no_agents(self, risk):
        assert risk.listAgentsNames() == []


@pytest.mark.unit
class TestLoadData:
    def test_a_new_agent_is_added_to_the_database(self, risk):
        risk.loadData(_desc(), saveMode=TOOLKIT_SAVEMODE_FILEANDDB)
        assert risk.listAgentsNames() == ["TestAgent"]

    def test_loading_the_same_name_and_version_twice_raises(self, risk):
        risk.loadData(_desc(), saveMode=TOOLKIT_SAVEMODE_FILEANDDB)
        with pytest.raises(ValueError, match="in the database"):
            risk.loadData(_desc(), saveMode=TOOLKIT_SAVEMODE_FILEANDDB)

    def test_a_json_string_payload_is_accepted(self, risk):
        import json

        risk.loadData(json.dumps(_desc(name="FromString")), saveMode=TOOLKIT_SAVEMODE_FILEANDDB)
        assert "FromString" in risk.listAgentsNames()

    def test_neither_a_dict_a_json_string_nor_a_file_path_raises(self, risk):
        with pytest.raises(ValueError, match="must be a file"):
            risk.loadData(12345, saveMode=TOOLKIT_SAVEMODE_FILEANDDB)

    @pytest.mark.xfail(
        strict=True,
        reason="B55: saveMode is only consulted in the `elif agentDoc is "
               "not None` branch. When no matching document exists yet -- "
               "true for the very first load of any agent -- the `if "
               "agentDoc is None` branch runs unconditionally and calls "
               "addDataSource regardless of saveMode, so "
               "TOOLKIT_SAVEMODE_NOSAVE still persists to the database. "
               "See the consolidated findings issue.",
    )
    def test_nosave_does_not_persist_a_new_agent(self, risk):
        risk.loadData(_desc(), saveMode=TOOLKIT_SAVEMODE_NOSAVE)
        assert risk.listAgentsNames() == []

    def test_nosave_currently_does_persist_a_new_agent(self, risk):
        """Characterisation of B55, so the defect above is a number, not an adjective."""
        risk.loadData(_desc(), saveMode=TOOLKIT_SAVEMODE_NOSAVE)
        assert risk.listAgentsNames() == ["TestAgent"]

    def test_replacing_the_same_name_and_version_updates_in_place(self, risk):
        risk.loadData(_desc(), saveMode=TOOLKIT_SAVEMODE_FILEANDDB)
        risk.loadData(_desc(), saveMode=TOOLKIT_SAVEMODE_FILEANDDB_REPLACE)
        assert risk.listAgentsNames() == ["TestAgent"]

    @pytest.mark.xfail(
        strict=True,
        reason="B56: the replace path looks up the existing document with "
               "getDataSourceDocument(datasourceName=name, version=<the "
               "NEW version>), so replacing with a different version never "
               "finds the old entry, falls into the 'create' branch, and "
               "leaves two datasources named TestAgent instead of one. "
               "See the consolidated findings issue.",
    )
    def test_replacing_with_a_different_version_updates_the_single_entry(self, risk):
        risk.loadData(_desc(version=(0, 0, 1)), saveMode=TOOLKIT_SAVEMODE_FILEANDDB)
        risk.loadData(_desc(version=(0, 0, 2)), saveMode=TOOLKIT_SAVEMODE_FILEANDDB_REPLACE)
        assert risk.listAgentsNames() == ["TestAgent"]

    def test_replacing_with_a_different_version_currently_duplicates_the_entry(self, risk):
        """Characterisation of B56."""
        risk.loadData(_desc(version=(0, 0, 1)), saveMode=TOOLKIT_SAVEMODE_FILEANDDB)
        risk.loadData(_desc(version=(0, 0, 2)), saveMode=TOOLKIT_SAVEMODE_FILEANDDB_REPLACE)
        assert risk.listAgentsNames() == ["TestAgent", "TestAgent"]

    @pytest.mark.xfail(
        strict=True,
        reason="B57: the replace path only assigns agentDoc.desc['version'] "
               "before saving -- the rest of desc (e.g. effectParameters) "
               "is never refreshed, even though `resource` (the raw JSON) "
               "is fully overwritten with the new description. "
               "See the consolidated findings issue.",
    )
    def test_replacing_with_the_same_version_refreshes_the_metadata(self, risk):
        risk.loadData(_desc(), saveMode=TOOLKIT_SAVEMODE_FILEANDDB)
        doc = risk.loadData(
            _desc(effectParameters={"tenbergeCoefficient": 2}),
            saveMode=TOOLKIT_SAVEMODE_FILEANDDB_REPLACE,
        )
        assert doc.desc["effectParameters"] == {"tenbergeCoefficient": 2}


@pytest.mark.unit
class TestGetAgent:
    def test_getting_an_agent_by_name_looks_it_up_in_the_database(self, risk):
        risk.loadData(_desc(), saveMode=TOOLKIT_SAVEMODE_FILEANDDB)
        agent = risk.getAgent("TestAgent")
        assert isinstance(agent, Agent)

    def test_an_unknown_agent_name_raises(self, risk):
        with pytest.raises(ValueError, match="not found"):
            risk.getAgent("NoSuchAgent")

    def test_getting_an_agent_from_a_descriptor_dict_skips_the_database(self, risk):
        agent = risk.getAgent(_desc(name="Inline"))
        assert isinstance(agent, Agent)
        assert "Inline" not in risk.listAgentsNames()

    def test_a_descriptor_that_is_neither_a_name_nor_a_dict_raises(self, risk):
        with pytest.raises(ValueError, match="must be the agent name"):
            risk.getAgent(12345)

    def test_loadagent_stamps_the_name_and_version_onto_the_description(self, risk):
        risk.loadAgent("Agent2", {"effectParameters": {}, "effects": {}}, version=[1, 0, 0])
        assert risk.listAgentsNames() == ["Agent2"]
