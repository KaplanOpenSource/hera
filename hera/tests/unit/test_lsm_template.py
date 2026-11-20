"""LSM/template.py: LSMTemplate's properties, parameter preparation, and
simulation lookups. run()/_toNetcdf() need a real LSM binary, real
input-file templates and subprocess execution -- out of scope for a
hermetic unit test.

The document is constructed by hand rather than via LSMToolkit.loadData(),
because loadData carries B92 (already documented in batch27): it never
persists the parsed template JSON's own desc (params, modelFolder, etc.)
into the datasource document it registers, so `params`/`modelFolder`
always KeyError on a template loaded that way. Building the document
directly here isolates LSMTemplate's own logic from that separate,
already-pinned defect.

B99: ``getSimulationByID`` calls ``self.getDocumentByID(id)``, but
``getDocumentByID`` is a ``Project``/toolkit method, not one
``LSMTemplate`` has itself -- it lives on ``self.Toolkit``, not ``self``.
Every call raises ``AttributeError``.
"""
import pytest

import hera.simulations.LSM.template as template_mod
from hera.simulations.LSM.template import LSMTemplate
from hera.utils import ureg


def _doc(**desc_overrides):
    desc = dict(params={"a": 1}, version=[1, 0, 0], datasourceName="myTemplate", modelFolder="sub")
    desc.update(desc_overrides)
    return {"resource": "/some/dir", "desc": desc}


@pytest.fixture()
def tpl(unit_toolkit_factory):
    from hera import toolkitHome

    toolkit = unit_toolkit_factory(toolkitHome.LSM)
    return LSMTemplate(_doc(), toolkit)


@pytest.mark.unit
class TestProperties:
    def test_toolkit_is_the_constructor_argument(self, tpl, unit_toolkit_factory):
        assert tpl.Toolkit is not None

    def test_doctype_simulation_is_a_fixed_string(self, tpl):
        assert tpl.doctype_simulation == "LSM_run"

    def test_dir_path_reads_the_document_resource(self, tpl):
        assert tpl.dirPath == "/some/dir"

    def test_params_and_version_and_template_name(self, tpl):
        assert tpl.params == {"a": 1}
        assert tpl.version == [1, 0, 0]
        assert tpl.templateName == "myTemplate"

    def test_model_folder_stays_absolute_if_already_absolute(self, unit_toolkit_factory):
        from hera import toolkitHome

        toolkit = unit_toolkit_factory(toolkitHome.LSM)
        tpl = LSMTemplate(_doc(modelFolder="/abs/sub"), toolkit)
        assert tpl.modelFolder == "/abs/sub"

    def test_model_folder_is_joined_to_files_directory_if_relative(self, tpl, unit_toolkit_factory):
        assert tpl.modelFolder == tpl.Toolkit.filesDirectory + "sub"


@pytest.mark.unit
class TestPrepareParams:
    def test_with_no_template_desc_params_pass_through_unchanged(self):
        result = LSMTemplate.prepareParams(template_desc=None, paramsToPrepare={"a": 1})
        assert result["a"] == 1

    def test_units_are_converted_to_the_templates_declared_unit(self):
        template_desc = {"params": {}, "units": {"speed": "m/s"}}
        result = LSMTemplate.prepareParams(
            template_desc=template_desc, paramsToPrepare={"speed": "36 km/hour"}
        )
        assert result["speed"] == pytest.approx(10.0)

    def test_duration_is_treated_as_minutes_even_without_units_declared(self):
        template_desc = {"params": {}, "units": {"duration": "min"}}
        result = LSMTemplate.prepareParams(template_desc=template_desc, paramsToPrepare={"duration": 5})
        assert result["duration"] == pytest.approx(5.0)

    def test_topo_xn_and_yn_are_cast_to_int(self):
        result = LSMTemplate.prepareParams(template_desc=None, paramsToPrepare={"TopoXn": 5.0, "TopoYn": 3.0})
        assert (result["TopoXn"], result["TopoYn"]) == (5, 3)
        assert isinstance(result["TopoXn"], int)


@pytest.mark.unit
class TestGetSimulationByIDIsBroken:
    @pytest.mark.xfail(
        strict=True,
        reason="B99: getSimulationByID calls self.getDocumentByID(id), but "
               "getDocumentByID is a Project/toolkit method LSMTemplate "
               "does not have itself -- it lives on self.Toolkit. "
               "See the consolidated findings issue.",
    )
    def test_it_should_find_a_registered_simulation(self, tpl):
        doc = tpl.Toolkit.addSimulationsDocument(
            type="LSM_run", resource="res", dataFormat="string",
            desc=dict(version=[1, 0, 0], templateName="myTemplate", simulationName="sim1", params={"a": 1}),
        )
        tpl.getSimulationByID(doc.id)

    def test_it_currently_raises_attributeerror(self, tpl):
        """Characterisation of B99."""
        doc = tpl.Toolkit.addSimulationsDocument(
            type="LSM_run", resource="res", dataFormat="string",
            desc=dict(version=[1, 0, 0], templateName="myTemplate", simulationName="sim1", params={"a": 1}),
        )
        with pytest.raises(AttributeError, match="getDocumentByID"):
            tpl.getSimulationByID(doc.id)


@pytest.mark.unit
class TestSimulationLookups:
    """SingleSimulation is monkeypatched to a passthrough, since building a
    real one needs actual netcdf files on disk -- these tests check the
    query/dispatch logic in LSMTemplate, not SingleSimulation itself."""

    @pytest.fixture(autouse=True)
    def _fake_single_simulation(self, monkeypatch):
        class _FakeSim:
            def __init__(self, doc):
                self.doc = doc

        monkeypatch.setattr(template_mod, "SingleSimulation", _FakeSim)

    def test_get_simulation_by_name_finds_the_matching_document(self, tpl):
        tpl.Toolkit.addSimulationsDocument(
            type="LSM_run", resource="res", dataFormat="string",
            desc=dict(version=[1, 0, 0], templateName="myTemplate", simulationName="sim1", params={"a": 1}),
        )
        found = tpl.getSimulationByName("sim1")
        assert len(found) == 1
        assert found[0].doc.desc["simulationName"] == "sim1"

    def test_get_simulations_defaults_to_the_templates_own_params(self, tpl):
        """With no query kwargs, getSimulations filters by the template's
        own default params (params={"a": 1} here) -- not "everything"."""
        tpl.Toolkit.addSimulationsDocument(
            type="LSM_run", resource="res1", dataFormat="string",
            desc=dict(version=[1, 0, 0], templateName="myTemplate", simulationName="sim1", params={"a": 1}),
        )
        tpl.Toolkit.addSimulationsDocument(
            type="LSM_run", resource="res2", dataFormat="string",
            desc=dict(version=[1, 0, 0], templateName="myTemplate", simulationName="sim2", params={"a": 2}),
        )
        found = tpl.getSimulations()
        assert len(found) == 1
        assert found[0].doc.desc["simulationName"] == "sim1"

    def test_get_simulations_query_overrides_the_default_params(self, tpl):
        tpl.Toolkit.addSimulationsDocument(
            type="LSM_run", resource="res2", dataFormat="string",
            desc=dict(version=[1, 0, 0], templateName="myTemplate", simulationName="sim2", params={"a": 2}),
        )
        found = tpl.getSimulations(a=2)
        assert len(found) == 1
        assert found[0].doc.desc["simulationName"] == "sim2"
