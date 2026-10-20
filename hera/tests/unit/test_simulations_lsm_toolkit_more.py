"""LSMToolkit: the three members test_simulations_lsm_toolkit.py left out --
``getSimulations``, ``getSimulationsList`` and ``prepareSlurmLSMExecution``.

The datasource documents are registered with ``addDataSource`` directly
rather than through ``loadData``, so the ``params``/``units`` each lookup
needs can be stated here instead of coming from a template JSON file.

``SingleSimulation`` is replaced with a pass-through stand-in for the
lookup tests: building a real one needs netcdf files on disk, and what
LSMToolkit owns is the query it issues, not the reader.
``slurm.prepareSlurmScriptExecution`` is replaced with a recorder so that
nothing is submitted to a scheduler.

Bugs pinned here, each with an xfail(strict) for the intended behaviour
and a passing characterisation of what actually happens:

  * B285: the same method logs ``logger.error(...)`` for a non-dict
    ``baseParameters`` and for an unusable ``jsonVariations`` without
    returning, so execution runs on into an unrelated failure.
"""
import json
import os

import pytest

import hera.simulations.LSM.toolkit as lsmToolkit_mod
from hera.utils import slurm
from hera.utils.unitHandler import ureg

TEMPLATE_NAME = "v4-general"
TEMPLATE_PARAMS = dict(TopoXmin=0, TopoXmax=100)


@pytest.fixture()
def lsm(unit_toolkit_factory):
    from hera import toolkitHome

    return unit_toolkit_factory(toolkitHome.LSM)


@pytest.fixture()
def registerTemplate(lsm):
    def _register(name=TEMPLATE_NAME, units=None, params=None):
        extra = dict(params=dict(TEMPLATE_PARAMS if params is None else params),
                     modelFolder="sub")
        if units is not None:
            extra["units"] = units
        return lsm.addDataSource(
            dataSourceName=name, resource="/templates", dataFormat="string",
            version=name, **extra
        )

    return _register


@pytest.fixture()
def registerSimulation(lsm):
    def _register(simulationName, **params):
        return lsm.addSimulationsDocument(
            type="LSM_run",
            resource=f"resource_{simulationName}",
            dataFormat="string",
            desc=dict(simulationName=simulationName, params=dict(params)),
        )

    return _register


@pytest.fixture()
def fakeSingleSimulation(monkeypatch):
    class _FakeSimulation:
        def __init__(self, document):
            self.document = document

    monkeypatch.setattr(lsmToolkit_mod, "SingleSimulation", _FakeSimulation)
    return _FakeSimulation


@pytest.fixture()
def slurmRecorder(monkeypatch):
    calls = []
    monkeypatch.setattr(
        slurm, "prepareSlurmScriptExecution", lambda **kwargs: calls.append(kwargs)
    )
    return calls


# ---------------------------------------------------------------------------
# getSimulations
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetSimulationsWithoutAUnitsTemplate:
    """getTemplateByName returns None when the template is absent, so a
    project that never registered one must simply fall back to the path for
    a template that declares no units."""

    def test_a_project_without_the_units_template_still_lists_simulations(
        self, lsm, registerSimulation, fakeSingleSimulation
    ):
        registerSimulation("sim1", TopoXmax=100)
        assert len(lsm.getSimulations(TopoXmax=100)) == 1

    def test_an_empty_project_without_the_units_template_lists_nothing(self, lsm):
        assert lsm.getSimulations() == []

    def test_naming_a_template_that_does_not_exist_behaves_the_same_way(
        self, lsm, registerTemplate, registerSimulation, fakeSingleSimulation
    ):
        registerTemplate(name=TEMPLATE_NAME)
        registerSimulation("sim1", TopoXmax=100)
        found = lsm.getSimulations(unitsTemplateVersion="noSuchTemplate", TopoXmax=100)
        assert [item.document.desc["simulationName"] for item in found] == ["sim1"]


@pytest.mark.unit
class TestGetSimulationsWithoutDeclaredUnits:
    def test_a_matching_simulation_is_returned(
        self, lsm, registerTemplate, registerSimulation, fakeSingleSimulation
    ):
        registerTemplate()
        registerSimulation("sim1", TopoXmax=100)
        found = lsm.getSimulations(TopoXmax=100)
        assert len(found) == 1
        assert found[0].document.desc["simulationName"] == "sim1"

    def test_a_non_matching_simulation_is_filtered_out(
        self, lsm, registerTemplate, registerSimulation, fakeSingleSimulation
    ):
        registerTemplate()
        registerSimulation("sim1", TopoXmax=100)
        registerSimulation("sim2", TopoXmax=999)
        found = lsm.getSimulations(TopoXmax=999)
        assert [item.document.desc["simulationName"] for item in found] == ["sim2"]

    def test_no_query_returns_every_simulation_of_the_type(
        self, lsm, registerTemplate, registerSimulation, fakeSingleSimulation
    ):
        registerTemplate()
        registerSimulation("sim1", TopoXmax=100)
        registerSimulation("sim2", TopoXmax=999)
        assert len(lsm.getSimulations()) == 2

    def test_the_simulation_name_narrows_the_lookup(
        self, lsm, registerTemplate, registerSimulation, fakeSingleSimulation
    ):
        registerTemplate()
        registerSimulation("sim1", TopoXmax=100)
        registerSimulation("sim2", TopoXmax=100)
        found = lsm.getSimulations(simulationName="sim2")
        assert [item.document.desc["simulationName"] for item in found] == ["sim2"]

    def test_an_unknown_simulation_name_returns_nothing(
        self, lsm, registerTemplate, registerSimulation, fakeSingleSimulation
    ):
        registerTemplate()
        registerSimulation("sim1", TopoXmax=100)
        assert lsm.getSimulations(simulationName="nothingLikeIt") == []

    def test_an_empty_project_returns_nothing(
        self, lsm, registerTemplate, fakeSingleSimulation
    ):
        registerTemplate()
        assert lsm.getSimulations() == []


@pytest.mark.unit
class TestGetSimulationsWithDeclaredUnits:
    def test_a_query_in_the_declared_unit_matches(
        self, lsm, registerTemplate, registerSimulation, fakeSingleSimulation
    ):
        registerTemplate(units=dict(TopoXmax="m"))
        registerSimulation("sim1", TopoXmax=100)
        found = lsm.getSimulations(TopoXmax=100 * ureg.m)
        assert [item.document.desc["simulationName"] for item in found] == ["sim1"]

    def test_a_query_in_another_length_unit_is_converted_first(
        self, lsm, registerTemplate, registerSimulation, fakeSingleSimulation
    ):
        registerTemplate(units=dict(TopoXmax="m"))
        registerSimulation("sim1", TopoXmax=100)
        found = lsm.getSimulations(TopoXmax=0.1 * ureg.km)
        assert [item.document.desc["simulationName"] for item in found] == ["sim1"]

    def test_a_query_in_the_wrong_unit_does_not_match(
        self, lsm, registerTemplate, registerSimulation, fakeSingleSimulation
    ):
        registerTemplate(units=dict(TopoXmax="m"))
        registerSimulation("sim1", TopoXmax=100)
        assert lsm.getSimulations(TopoXmax=100 * ureg.km) == []

    def test_a_bare_number_where_units_are_declared_is_refused(
        self, lsm, registerTemplate, registerSimulation, fakeSingleSimulation
    ):
        registerTemplate(units=dict(TopoXmax="m"))
        registerSimulation("sim1", TopoXmax=100)
        with pytest.raises(ValueError, match="must use either pint or unum"):
            lsm.getSimulations(TopoXmax=100)

    def test_a_key_outside_the_units_block_is_passed_through_unconverted(
        self, lsm, registerTemplate, registerSimulation, fakeSingleSimulation
    ):
        registerTemplate(units=dict(TopoXmax="m"))
        registerSimulation("sim1", TopoXmin=7)
        found = lsm.getSimulations(TopoXmin=7)
        assert [item.document.desc["simulationName"] for item in found] == ["sim1"]


@pytest.mark.unit
class TestGetSimulationsReaderFailures:
    def test_an_unreadable_simulation_is_skipped_without_printing(
        self, lsm, registerTemplate, registerSimulation, capsys
    ):
        """A document whose data cannot be read is still skipped, but the
        diagnostic goes to the toolkit's logger, not to stdout."""
        registerTemplate()
        registerSimulation("sim1", TopoXmax=100)
        found = lsm.getSimulations(TopoXmax=100)
        assert found == []
        assert capsys.readouterr().out == ""

    def test_a_keyboard_interrupt_is_not_swallowed(
        self, lsm, registerTemplate, registerSimulation, monkeypatch
    ):
        registerTemplate()
        registerSimulation("sim1", TopoXmax=100)

        def _interrupt(document):
            raise KeyboardInterrupt

        monkeypatch.setattr(lsmToolkit_mod, "SingleSimulation", _interrupt)
        with pytest.raises(KeyboardInterrupt):
            lsm.getSimulations(TopoXmax=100)


# ---------------------------------------------------------------------------
# getSimulationsList
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetSimulationsListLongFormat:
    def test_an_empty_project_is_reported_as_a_missing_file(self, lsm):
        with pytest.raises(FileNotFoundError, match="No simulations.old found"):
            lsm.getSimulationsList()

    def test_one_simulation_is_melted_into_one_row_per_field(
        self, lsm, registerSimulation
    ):
        document = registerSimulation("sim1", TopoXmax=100)
        table = lsm.getSimulationsList()
        assert set(table.columns) == {"variable", "value"}
        assert set(table.index) == {document.id}

    def test_the_metadata_and_the_parameters_both_appear(self, lsm, registerSimulation):
        registerSimulation("sim1", TopoXmax=100)
        table = lsm.getSimulationsList()
        pairs = dict(zip(table["variable"], table["value"]))
        assert pairs["simulationName"] == "sim1"
        assert pairs["params__TopoXmax"] == 100

    def test_the_id_becomes_the_index_rather_than_a_field(
        self, lsm, registerSimulation
    ):
        registerSimulation("sim1", TopoXmax=100)
        table = lsm.getSimulationsList()
        assert "id" not in set(table["variable"])

    def test_a_query_narrows_the_listing(self, lsm, registerSimulation):
        registerSimulation("sim1", TopoXmax=100)
        registerSimulation("sim2", TopoXmax=999)
        table = lsm.getSimulationsList(params__TopoXmax=999)
        pairs = dict(zip(table["variable"], table["value"]))
        assert pairs["simulationName"] == "sim2"

    def test_two_simulations_each_melt_to_their_own_fields(
        self, lsm, registerSimulation
    ):
        """Two documents with four fields each melt to 2*4 = 8 rows."""
        registerSimulation("sim1", TopoXmax=100, extra=1)
        registerSimulation("sim2", TopoXmax=100, extra=2)
        table = lsm.getSimulationsList()
        assert len(table) == 8

    def test_each_run_keeps_its_own_parameters(self, lsm, registerSimulation):
        registerSimulation("sim1", TopoXmax=100, extra=1)
        registerSimulation("sim2", TopoXmax=100, extra=2)
        table = lsm.getSimulationsList()
        extras = table[table["variable"] == "params__extra"]
        byId = {}
        for identifier, value in zip(extras.index, extras["value"]):
            byId.setdefault(identifier, set()).add(value)
        assert sorted(byId.values(), key=sorted) == [{1}, {2}]


@pytest.mark.unit
class TestGetSimulationsListWideFormat:
    def test_one_simulation_becomes_a_single_row_of_fields(
        self, lsm, registerSimulation
    ):
        registerSimulation("sim1", TopoXmax=100)
        table = lsm.getSimulationsList(wideFormat=True)
        assert len(table) == 1
        assert table["simulationName"].iloc[0] == "sim1"
        assert table["params__TopoXmax"].iloc[0] == 100

    def test_an_empty_project_is_reported_as_a_missing_file(self, lsm):
        with pytest.raises(FileNotFoundError, match="No simulations.old found"):
            lsm.getSimulationsList(wideFormat=True)

    def test_a_pivot_failure_is_not_reported_as_a_missing_file(
        self, lsm, registerSimulation, monkeypatch
    ):
        """Only emptiness may become FileNotFoundError; a reshape failure on
        a project that does hold simulations must surface as itself."""
        registerSimulation("sim1", TopoXmax=100)

        def boom(*args, **kwargs):
            raise ValueError("Index contains duplicate entries, cannot reshape")

        monkeypatch.setattr(lsmToolkit_mod.pandas.DataFrame, "pivot", boom)
        with pytest.raises(ValueError, match="cannot reshape"):
            lsm.getSimulationsList(wideFormat=True)

    def test_two_simulations_should_give_two_rows(self, lsm, registerSimulation):
        registerSimulation("sim1", TopoXmax=100, extra=1)
        registerSimulation("sim2", TopoXmax=100, extra=2)
        assert len(lsm.getSimulationsList(wideFormat=True)) == 2

    def test_each_row_carries_its_own_runs_parameters(self, lsm, registerSimulation):
        registerSimulation("sim1", TopoXmax=100, extra=1)
        registerSimulation("sim2", TopoXmax=100, extra=2)
        table = lsm.getSimulationsList(wideFormat=True)
        pairs = set(zip(table["simulationName"], table["params__extra"]))
        assert pairs == {("sim1", 1), ("sim2", 2)}


# ---------------------------------------------------------------------------
# prepareSlurmLSMExecution
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPrepareSlurmLSMExecution:
    VARIATIONS = [{"TopoXmax": [100, 200]}]

    def _prepare(self, lsm, **overrides):
        arguments = dict(
            baseParameters=dict(TEMPLATE_PARAMS),
            jsonVariations=list(self.VARIATIONS),
            templateName=TEMPLATE_NAME,
        )
        arguments.update(overrides)
        return lsm.prepareSlurmLSMExecution(**arguments)

    def _scriptsDirectory(self, lsm):
        return os.path.join(lsm.filesDirectory, "simulationsScripts")

    def test_it_writes_one_run_script_per_variation(self, lsm, slurmRecorder):
        self._prepare(lsm)
        assert sorted(
            name for name in os.listdir(self._scriptsDirectory(lsm))
            if name.startswith("LSM_")
        ) == ["LSM_Simulation_0", "LSM_Simulation_1"]

    def test_the_case_list_names_every_simulation(self, lsm, slurmRecorder):
        self._prepare(lsm)
        with open(os.path.join(self._scriptsDirectory(lsm), "cases.txt")) as cases:
            assert cases.read().split() == ["LSM_Simulation_0", "LSM_Simulation_1"]

    def test_the_scheduler_is_handed_the_case_list(self, lsm, slurmRecorder):
        self._prepare(lsm)
        assert len(slurmRecorder) == 1
        assert str(slurmRecorder[0]["jobDirListFilePath"]).endswith("cases.txt")

    def test_stations_are_written_out_for_the_run_scripts(
        self, lsm, slurmRecorder, tmp_path
    ):
        import pandas

        stations = pandas.DataFrame(dict(x=[0.0], y=[0.0], u=[3.0]))
        self._prepare(lsm, stations=stations)
        assert os.path.isfile(
            os.path.join(self._scriptsDirectory(lsm), "stations.parquet")
        )

    def test_a_variations_file_on_disk_is_read(self, lsm, slurmRecorder, tmp_path):
        """The str branch loads the same variations the list branch takes."""
        path = tmp_path / "variations.json"
        path.write_text(json.dumps(self.VARIATIONS))
        self._prepare(lsm, jsonVariations=str(path))
        assert sorted(
            name for name in os.listdir(self._scriptsDirectory(lsm))
            if name.startswith("LSM_")
        ) == ["LSM_Simulation_0", "LSM_Simulation_1"]

    def test_an_unreadable_variations_path_is_refused(self, lsm, slurmRecorder):
        with pytest.raises(FileNotFoundError):
            self._prepare(lsm, jsonVariations="/no/such/variations.json")


@pytest.mark.unit
class TestPrepareSlurmLSMExecutionInputValidation:
    VARIATIONS = [{"TopoXmax": [100]}]

    @pytest.mark.xfail(
        strict=True,
        reason="B285: prepareSlurmLSMExecution reports bad input with "
               "`logger.error('Slurm preparation can only handle dictionary "
               "of parameters to run variations')` and then falls straight "
               "through into JSONVariations(baseParameters, ...), which "
               "fails on `dict(base)` with 'dictionary update sequence "
               "element #0 has length 1; 2 is required'. The caller gets an "
               "opaque pandas/builtins error instead of the diagnostic that "
               "was written for them; the same shape of mistake is repeated "
               "for jsonVariations a few lines below. Both need a `return` "
               "or a raise. See the consolidated findings issue.",
    )
    def test_a_non_dict_baseParameters_should_be_refused_clearly(self, lsm):
        with pytest.raises((ValueError, TypeError), match="dictionary of parameters"):
            lsm.prepareSlurmLSMExecution(
                baseParameters="notADictionary",
                jsonVariations=list(self.VARIATIONS),
                templateName=TEMPLATE_NAME,
            )

    def test_a_non_dict_baseParameters_currently_fails_deeper_in(self, lsm):
        """Characterisation of B285."""
        with pytest.raises(ValueError, match="dictionary update sequence"):
            lsm.prepareSlurmLSMExecution(
                baseParameters="notADictionary",
                jsonVariations=list(self.VARIATIONS),
                templateName=TEMPLATE_NAME,
            )

    @pytest.mark.xfail(
        strict=True,
        reason="B285: the jsonVariations check has the same shape -- "
               "`elif not isinstance(jsonVariations, list): logger.error("
               "'Slurm preparation only supports json variation input as "
               "path or dict')` with no return -- so a dict (which the "
               "message itself invites) falls through into JSONVariations "
               "and fails on iteration instead. See the consolidated "
               "findings issue.",
    )
    def test_an_unsupported_variations_type_should_be_refused_clearly(self, lsm):
        with pytest.raises(
            (ValueError, TypeError), match="only supports json variation input"
        ):
            lsm.prepareSlurmLSMExecution(
                baseParameters=dict(TEMPLATE_PARAMS),
                jsonVariations={"TopoXmax": [100]},
                templateName=TEMPLATE_NAME,
            )

    def test_an_unsupported_variations_type_currently_fails_deeper_in(self, lsm):
        """Characterisation of B285: the error message says 'path or
        dict', yet a dict is exactly what the isinstance check rejects."""
        with pytest.raises(Exception) as raised:
            lsm.prepareSlurmLSMExecution(
                baseParameters=dict(TEMPLATE_PARAMS),
                jsonVariations={"TopoXmax": [100]},
                templateName=TEMPLATE_NAME,
            )
        assert "only supports json variation input" not in str(raised.value)
