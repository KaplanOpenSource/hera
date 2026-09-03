"""LSMToolkit: the three members test_simulations_lsm_toolkit.py left out --
``getSimulations``, ``getSimulationsList`` and ``prepareSlurmLSMExecution``.

The datasource documents are registered with ``addDataSource`` directly
rather than through ``loadData``, because ``loadData`` never persists the
parsed template JSON's own desc (B92, already pinned in
test_simulations_lsm_toolkit.py) -- a template registered the normal way
has no ``params`` and no ``units`` to look up, so it could not exercise
either lookup here.

``SingleSimulation`` is replaced with a pass-through stand-in for the
lookup tests: building a real one needs netcdf files on disk, and what
LSMToolkit owns is the query it issues, not the reader.
``slurm.prepareSlurmScriptExecution`` is replaced with a recorder so that
nothing is submitted to a scheduler.

Bugs pinned here, each with an xfail(strict) for the intended behaviour
and a passing characterisation of what actually happens:

  * B281: ``getSimulations`` dereferences
    ``getTemplateByName(unitsTemplateVersion)`` without checking for None,
    so it dies with AttributeError in any project that does not happen to
    hold a template called "v4-general" (its default).
  * B282: ``getSimulations`` wraps ``SingleSimulation(doc)`` in a bare
    ``except:`` that prints a warning and drops the document, so callers
    silently receive a short list.
  * B283: ``getSimulationsList(wideFormat=True)`` reports
    FileNotFoundError("No simulations.old found") for any project holding
    more than one matching simulation -- the pivot's "Index contains
    duplicate entries" ValueError is swallowed by the ``except ValueError``
    that exists to report emptiness.
  * B284: ``prepareSlurmLSMExecution`` calls
    ``LSMTemplate.prepareParams(desc=None, ...)``, but that parameter is
    named ``template_desc`` -- every call raises TypeError.
  * B285: the same method logs ``logger.error(...)`` for a non-dict
    ``baseParameters`` and for an unusable ``jsonVariations`` without
    returning, so execution runs on into an unrelated failure.

B276 (the ``product()`` cross-join, pinned against
LSMTemplate._getSimulationsList in test_lsm_template_more.py) is repeated
verbatim in ``getSimulationsList``; its consequences here are
characterised rather than re-pinned as a separate finding.
"""
import json
import os

import pytest

import hera.simulations.LSM.toolkit as lsmToolkit_mod
from hera.simulations.LSM.template import LSMTemplate
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
class TestGetSimulationsNeedsAUnitsTemplate:
    @pytest.mark.xfail(
        strict=True,
        reason="B281: getSimulations looks the units template up with "
               "`template = self.getTemplateByName(unitsTemplateVersion)` "
               "and immediately reads `template._document['desc']`. "
               "getTemplateByName is documented to return None when the "
               "template is absent (and does), so any project that has not "
               "registered a template literally named 'v4-general' -- the "
               "default -- cannot list its simulations at all: every call "
               "raises AttributeError: 'NoneType' object has no attribute "
               "'_document'. See the consolidated findings issue.",
    )
    def test_a_project_without_the_units_template_should_still_list_simulations(
        self, lsm, registerSimulation, fakeSingleSimulation
    ):
        registerSimulation("sim1", TopoXmax=100)
        assert len(lsm.getSimulations(TopoXmax=100)) == 1

    def test_a_project_without_the_units_template_currently_raises(self, lsm):
        """Characterisation of B281."""
        with pytest.raises(AttributeError, match="'NoneType' object has no attribute '_document'"):
            lsm.getSimulations()

    def test_naming_a_template_that_does_not_exist_raises_the_same_way(
        self, lsm, registerTemplate
    ):
        """Characterisation of B281: the default name is not the problem."""
        registerTemplate(name=TEMPLATE_NAME)
        with pytest.raises(AttributeError, match="_document"):
            lsm.getSimulations(unitsTemplateVersion="noSuchTemplate")


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
class TestGetSimulationsSwallowsReaderFailures:
    @pytest.mark.xfail(
        strict=True,
        reason="B282: getSimulations wraps SingleSimulation(doc) in a "
               "bare `except:` that prints 'Warning: could not find data "
               "with the following document: ...' and moves on. A document "
               "whose data cannot be read is therefore dropped from the "
               "result without the caller being able to tell, the diagnostic "
               "goes to stdout instead of the toolkit's logger, and the bare "
               "except also catches KeyboardInterrupt and SystemExit. See "
               "the consolidated findings issue.",
    )
    def test_an_unreadable_simulation_should_not_vanish_silently(
        self, lsm, registerTemplate, registerSimulation, capsys
    ):
        registerTemplate()
        registerSimulation("sim1", TopoXmax=100)
        found = lsm.getSimulations(TopoXmax=100)
        assert capsys.readouterr().out == ""
        assert len(found) == 1

    def test_an_unreadable_simulation_is_currently_dropped_with_a_printed_warning(
        self, lsm, registerTemplate, registerSimulation, capsys
    ):
        """Characterisation of B282."""
        registerTemplate()
        registerSimulation("sim1", TopoXmax=100)
        found = lsm.getSimulations(TopoXmax=100)
        assert found == []
        assert "could not find data with the following document" in capsys.readouterr().out


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

    def test_two_simulations_currently_produce_the_cross_join(
        self, lsm, registerSimulation
    ):
        """Characterisation of B276, repeated verbatim in this method:
        two documents with three fields each melt to 2*2*4 = 16 rows."""
        registerSimulation("sim1", TopoXmax=100, extra=1)
        registerSimulation("sim2", TopoXmax=100, extra=2)
        table = lsm.getSimulationsList()
        assert len(table) == 16

    def test_the_cross_join_currently_mixes_one_run_with_anothers_parameters(
        self, lsm, registerSimulation
    ):
        """Characterisation of B276."""
        registerSimulation("sim1", TopoXmax=100, extra=1)
        registerSimulation("sim2", TopoXmax=100, extra=2)
        table = lsm.getSimulationsList()
        extras = table[table["variable"] == "params__extra"]
        byId = {}
        for identifier, value in zip(extras.index, extras["value"]):
            byId.setdefault(identifier, set()).add(value)
        assert all(values == {1, 2} for values in byId.values())


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

    @pytest.mark.xfail(
        strict=True,
        reason="B283: with two or more matching simulations the wide "
               "format cannot be built at all. The product() cross-join "
               "(B276) puts each document id into the melted index twice "
               "with the same 'variable', so pandas.pivot raises "
               "ValueError('Index contains duplicate entries, cannot "
               "reshape') -- and that lands in the `except ValueError` "
               "whose only purpose is to report an empty project, so the "
               "caller is told FileNotFoundError('No simulations.old "
               "found') about a project that holds two. See the "
               "consolidated findings issue.",
    )
    def test_two_simulations_should_give_two_rows(self, lsm, registerSimulation):
        registerSimulation("sim1", TopoXmax=100, extra=1)
        registerSimulation("sim2", TopoXmax=100, extra=2)
        assert len(lsm.getSimulationsList(wideFormat=True)) == 2

    def test_two_simulations_are_currently_reported_as_none_found(
        self, lsm, registerSimulation
    ):
        """Characterisation of B283."""
        registerSimulation("sim1", TopoXmax=100, extra=1)
        registerSimulation("sim2", TopoXmax=100, extra=2)
        with pytest.raises(FileNotFoundError, match="No simulations.old found"):
            lsm.getSimulationsList(wideFormat=True)

    def test_the_long_format_of_the_same_project_still_returns_data(
        self, lsm, registerSimulation
    ):
        """Characterisation of B283: the data is plainly there."""
        registerSimulation("sim1", TopoXmax=100, extra=1)
        registerSimulation("sim2", TopoXmax=100, extra=2)
        assert len(lsm.getSimulationsList(wideFormat=False)) > 0


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

    @pytest.mark.xfail(
        strict=True,
        reason="B284: prepareSlurmLSMExecution calls "
               "`LSMTemplate.prepareParams(desc=None, "
               "paramsToPrepare=jsonConfig)`, but that staticmethod's first "
               "parameter is named `template_desc`. Every call therefore "
               "raises TypeError: prepareParams() got an unexpected keyword "
               "argument 'desc' on the first variation, before a single "
               "run script is written. The whole Slurm entry point is "
               "unusable. See the consolidated findings issue.",
    )
    def test_it_should_write_one_run_script_per_variation(self, lsm, slurmRecorder):
        self._prepare(lsm)
        scriptsDirectory = os.path.join(lsm.filesDirectory, "simulationsScripts")
        assert sorted(
            name for name in os.listdir(scriptsDirectory) if name.startswith("LSM_")
        ) == ["LSM_Simulation_0", "LSM_Simulation_1"]

    def test_it_currently_raises_on_the_first_variation(self, lsm, slurmRecorder):
        """Characterisation of B284."""
        with pytest.raises(TypeError, match="unexpected keyword argument 'desc'"):
            self._prepare(lsm)

    def test_the_parameter_it_passes_is_not_the_one_declared(self):
        """Characterisation of B284: `template_desc`, not `desc`."""
        import inspect

        parameters = list(
            inspect.signature(LSMTemplate.prepareParams).parameters
        )
        assert parameters == ["template_desc", "paramsToPrepare"]

    def test_nothing_reaches_the_scheduler(self, lsm, slurmRecorder):
        """Characterisation of B284."""
        with pytest.raises(TypeError):
            self._prepare(lsm)
        assert slurmRecorder == []

    def test_the_scripts_directory_is_created_before_the_failure(
        self, lsm, slurmRecorder
    ):
        """Characterisation of B284: a half-built layout is left behind."""
        with pytest.raises(TypeError):
            self._prepare(lsm)
        assert os.path.isdir(
            os.path.join(lsm.filesDirectory, "simulationsScripts")
        )

    def test_stations_are_written_out_before_the_failure(
        self, lsm, slurmRecorder, tmp_path
    ):
        import pandas

        stations = pandas.DataFrame(dict(x=[0.0], y=[0.0], u=[3.0]))
        with pytest.raises(TypeError):
            self._prepare(lsm, stations=stations)
        assert os.path.isfile(
            os.path.join(lsm.filesDirectory, "simulationsScripts", "stations.parquet")
        )

    def test_a_variations_file_on_disk_is_read(self, lsm, slurmRecorder, tmp_path):
        """The str branch is taken before the B284 failure."""
        path = tmp_path / "variations.json"
        path.write_text(json.dumps(self.VARIATIONS))
        with pytest.raises(TypeError, match="unexpected keyword argument 'desc'"):
            self._prepare(lsm, jsonVariations=str(path))

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
