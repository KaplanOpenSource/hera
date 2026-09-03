"""measurements/experiment/CLI.py: the argparse command handlers and the
three private scaffolding helpers.

Covered: ``_resolve_project_name``, ``experiments_list``,
``experiments_table``, ``get_experiment_data``, ``create_experiment``,
``_create_empty_class``, ``_create_repository``,
``_make_runtimeExperimentData`` and ``load_experiment_to_project``.

The experiment handlers all reach for an ``experimentHome`` toolkit, the
``dataToolkit`` or the ``hera.utils.data`` project CLI.  Each of those is
replaced with a light recorder, so what is under test is the CLI plumbing:
how the project name is resolved, which toolkit call a subcommand maps to,
what is printed, what is written to disk and which errors bad input
raises.  Nothing here touches MongoDB, an Argos zip file or the network,
and every path written to lives under ``tmp_path``.

Deliberately not tested: the generated experiment class file is checked for
structure only, not imported or executed -- it subclasses
``experimentSetupWithData``, which needs a project and a real datasource
tree.  The argparse wiring in ``hera/bin/hera-experiment`` also stays
untested: it sits under ``if __name__ == "__main__"`` and cannot be
imported.

Two bugs are pinned here:

* B129 -- ``_create_repository`` writes into ``runtimeExperimentData/``
  before that directory is created, so ``create_experiment`` always dies
  with FileNotFoundError.
* B130 -- ``create_experiment`` writes ``createNodeRedDeviceMap.sh`` into
  the process cwd instead of the experiment directory.
"""
import json
import os
from argparse import Namespace
from types import SimpleNamespace

import pytest

from hera.measurements.experiment import CLI


# ---------------------------------------------------------------------------
# fakes
# ---------------------------------------------------------------------------

class _FakeExperiment:
    def __init__(self, name):
        self.name = name
        self.dataCalls = []

    def getExperimentData(self):
        return self

    def getData(self, deviceType, deviceName=None, perDevice=None):
        self.dataCalls.append(dict(deviceType=deviceType, deviceName=deviceName,
                                   perDevice=perDevice))
        return {"value": [1, 2]}


class _FakeExperimentHome:
    """Recorder standing in for experimentHome; the class collects instances."""

    instances = []
    experimentNames = ("bigExperiment", "smallExperiment")
    table = "<experiments table>"

    def __init__(self, projectName, **kwargs):
        self.projectName = projectName
        self.kwargs = kwargs
        self.requestedExperiments = []
        type(self).instances.append(self)

    def getExperimentsMap(self):
        return {name: f"<document {name}>" for name in type(self).experimentNames}

    def getExperimentsTable(self):
        return type(self).table

    def getExperiment(self, experimentName):
        self.requestedExperiments.append(experimentName)
        return _FakeExperiment(experimentName)


class _TablelessExperimentHome:
    """An experimentHome without getExperimentsTable, for the fallback branch."""

    instances = []
    experimentNames = ("alpha", "beta")

    def __init__(self, projectName, **kwargs):
        self.projectName = projectName
        type(self).instances.append(self)

    def getExperimentsMap(self):
        return {name: None for name in type(self).experimentNames}


class _FakeDataToolkit:
    instances = []

    def __init__(self, *args, **kwargs):
        self.calls = []
        self.loadError = None
        type(self).instances.append(self)

    def addRepository(self, **kwargs):
        self.calls.append(("addRepository", kwargs))

    def loadAllDatasourcesInRepositoryToProject(self, *args, **kwargs):
        self.calls.append(("loadAllDatasourcesInRepositoryToProject", args, kwargs))
        if self.loadError is not None:
            raise self.loadError

    def deleteDataSource(self, **kwargs):
        self.calls.append(("deleteDataSource", kwargs))

    def names(self):
        return [entry[0] for entry in self.calls]


class _FakeProjectCLI:
    def __init__(self):
        self.calls = []

    def project_create(self, arguments):
        self.calls.append(("project_create", arguments))

    def repository_load(self, arguments):
        self.calls.append(("repository_load", arguments))

    def names(self):
        return [entry[0] for entry in self.calls]


class _FakeEntity:
    def __init__(self, perDevice):
        self.properties = {"StoreDataPerDevice": perDevice}


class _FakeZipFile:
    """Stands in for argos' ExperimentZipFile."""

    def __init__(self, entities):
        self._entities = entities
        self.entityType = {}
        for entity in entities:
            byName = self.entityType.setdefault(entity["entityTypeName"], {})
            byName[entity["entityName"]] = _FakeEntity(entity["perDevice"])

    def getExperimentEntities(self):
        return [dict(entityTypeName=entity["entityTypeName"],
                     entityName=entity["entityName"]) for entity in self._entities]


@pytest.fixture()
def experimentHome(monkeypatch):
    """Inject the recorder over the module the handlers import from."""
    _FakeExperimentHome.instances = []
    _FakeExperimentHome.experimentNames = ("bigExperiment", "smallExperiment")
    monkeypatch.setattr("hera.measurements.experiment.experiment.experimentHome",
                        _FakeExperimentHome)
    return _FakeExperimentHome


@pytest.fixture()
def projectList(monkeypatch):
    """A controllable datalayer.getProjectList()."""
    projects = ["EXISTING_PROJECT"]
    monkeypatch.setattr(CLI, "datalayer",
                        SimpleNamespace(getProjectList=lambda: list(projects)))
    return projects


@pytest.fixture()
def dataToolkit(monkeypatch):
    _FakeDataToolkit.instances = []
    monkeypatch.setattr(CLI, "dataToolkit", _FakeDataToolkit)
    return _FakeDataToolkit


@pytest.fixture()
def projectCLI(monkeypatch):
    fake = _FakeProjectCLI()
    monkeypatch.setattr(CLI, "projectCLI", fake)
    return fake


def _experimentDirectory(tmp_path, name="myExperiment", withRuntime=True):
    """An experiment directory as create_experiment would have it mid-run.

    ``runtimeExperimentData`` is pre-created because of B129: without it,
    _create_repository cannot write its configuration file at all.
    """
    path = tmp_path / name
    (path / "code").mkdir(parents=True)
    (path / "data").mkdir()
    if withRuntime:
        (path / "runtimeExperimentData").mkdir()
    return path


# ---------------------------------------------------------------------------
# project name resolution
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestResolveProjectName:
    def test_an_explicit_projectname_wins(self, tmp_path):
        arguments = Namespace(projectName="FROM_FLAG",
                              configurationFile=str(tmp_path / "missing.json"))
        assert CLI._resolve_project_name(arguments) == "FROM_FLAG"

    def test_without_a_projectname_it_falls_back_to_the_configuration_file(self, tmp_path):
        config = tmp_path / "caseConfiguration.json"
        config.write_text(json.dumps({"projectName": "FROM_FILE"}))
        arguments = Namespace(projectName=None, configurationFile=str(config))
        assert CLI._resolve_project_name(arguments) == "FROM_FILE"

    def test_an_empty_projectname_also_falls_back(self, tmp_path):
        config = tmp_path / "caseConfiguration.json"
        config.write_text(json.dumps({"projectName": "FROM_FILE"}))
        arguments = Namespace(projectName="", configurationFile=str(config))
        assert CLI._resolve_project_name(arguments) == "FROM_FILE"

    def test_with_no_configurationfile_argument_it_reads_the_default_filename(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "caseConfiguration.json").write_text(json.dumps({"projectName": "CWD_PROJECT"}))
        assert CLI._resolve_project_name(Namespace()) == "CWD_PROJECT"

    def test_a_configuration_without_a_projectname_raises_keyerror(self, tmp_path):
        config = tmp_path / "caseConfiguration.json"
        config.write_text(json.dumps({"other": "field"}))
        with pytest.raises(KeyError, match="projectName"):
            CLI._resolve_project_name(Namespace(projectName=None, configurationFile=str(config)))

    def test_a_missing_configuration_file_raises_value_error(self, tmp_path):
        arguments = Namespace(projectName=None,
                              configurationFile=str(tmp_path / "nothingHere.json"))
        with pytest.raises(ValueError):
            CLI._resolve_project_name(arguments)


# ---------------------------------------------------------------------------
# list / table
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestExperimentsList:
    def test_it_prints_the_experiment_names_of_the_project(self, experimentHome, capsys):
        CLI.experiments_list(Namespace(projectName="MY_PROJECT"))

        assert experimentHome.instances[0].projectName == "MY_PROJECT"
        assert capsys.readouterr().out.strip() == "['bigExperiment', 'smallExperiment']"

    def test_it_prints_an_empty_list_when_the_project_has_no_experiments(self, experimentHome, capsys):
        experimentHome.experimentNames = ()
        CLI.experiments_list(Namespace(projectName="EMPTY"))
        assert capsys.readouterr().out.strip() == "[]"

    def test_the_project_name_may_come_from_the_configuration_file(self, experimentHome, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "caseConfiguration.json").write_text(json.dumps({"projectName": "CONFIG_PROJECT"}))
        CLI.experiments_list(Namespace(projectName=None))
        assert experimentHome.instances[0].projectName == "CONFIG_PROJECT"

    def test_it_builds_the_home_directly_rather_than_through_the_toolkit_registry(self, experimentHome, monkeypatch):
        """The handler deliberately bypasses toolkitHome.getToolkit, which
        would inject projectName twice."""
        def _refuse(*args, **kwargs):
            raise AssertionError("toolkitHome.getToolkit must not be used here")

        # Patch the class, not the singleton instance: `getToolkit` lives on
        # ToolkitHome, so patching the instance makes monkeypatch "restore" it
        # as a permanent *instance* attribute on teardown, which then shadows
        # the class for the rest of the session and silently defeats any later
        # class-level patch of it. That broke all 64 tests in
        # test_simulations_cli.py whenever this file happened to run first.
        monkeypatch.setattr(type(CLI.toolkitHome), "getToolkit", _refuse)
        CLI.experiments_list(Namespace(projectName="MY_PROJECT"))
        assert len(experimentHome.instances) == 1


@pytest.mark.unit
class TestExperimentsTable:
    def test_it_prints_the_table_when_the_home_can_build_one(self, experimentHome, capsys):
        CLI.experiments_table(Namespace(projectName="MY_PROJECT"))
        assert capsys.readouterr().out.strip() == "<experiments table>"

    def test_without_a_table_method_it_prints_a_numbered_fallback(self, monkeypatch, capsys):
        _TablelessExperimentHome.instances = []
        _TablelessExperimentHome.experimentNames = ("alpha", "beta")
        monkeypatch.setattr("hera.measurements.experiment.experiment.experimentHome",
                            _TablelessExperimentHome)

        CLI.experiments_table(Namespace(projectName="MY_PROJECT"))

        printed = [line for line in capsys.readouterr().out.splitlines() if line.strip()]
        assert printed[0] == "Experiments in project 'MY_PROJECT':"
        assert printed[1].strip() == "1. alpha"
        assert printed[2].strip() == "2. beta"

    def test_the_fallback_reports_an_empty_project_as_none(self, monkeypatch, capsys):
        _TablelessExperimentHome.instances = []
        _TablelessExperimentHome.experimentNames = ()
        monkeypatch.setattr("hera.measurements.experiment.experiment.experimentHome",
                            _TablelessExperimentHome)

        CLI.experiments_table(Namespace(projectName="EMPTY"))

        assert "(none)" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetExperimentData:
    def _arguments(self, **overrides):
        arguments = Namespace(projectName="EXISTING_PROJECT",
                              experiment="bigExperiment",
                              deviceType="Sonic")
        for key, value in overrides.items():
            setattr(arguments, key, value)
        return arguments

    def test_an_unknown_project_is_rejected_before_any_toolkit_is_built(self, experimentHome, projectList):
        with pytest.raises(ValueError, match="does not exists"):
            CLI.get_experiment_data(self._arguments(projectName="NO_SUCH_PROJECT"))
        assert experimentHome.instances == []

    def test_it_prints_the_device_data_as_a_dataframe(self, experimentHome, projectList, capsys):
        CLI.get_experiment_data(self._arguments())

        printed = capsys.readouterr().out
        assert "value" in printed
        assert experimentHome.instances[0].requestedExperiments == ["bigExperiment"]

    def test_the_device_type_is_forwarded_and_the_optionals_default_to_none(self, experimentHome, projectList, monkeypatch):
        recorded = []
        monkeypatch.setattr(_FakeExperiment, "getData",
                            lambda self, *args, **kwargs: recorded.append((args, kwargs)) or {})
        CLI.get_experiment_data(self._arguments())
        assert recorded == [(("Sonic",), {"deviceName": None, "perDevice": None})]

    def test_an_explicit_devicename_and_perdevice_are_forwarded(self, experimentHome, projectList, monkeypatch):
        recorded = []
        monkeypatch.setattr(_FakeExperiment, "getData",
                            lambda self, *args, **kwargs: recorded.append((args, kwargs)) or {})
        CLI.get_experiment_data(self._arguments(deviceName="sonic_01", perDevice=True))
        assert recorded == [(("Sonic",), {"deviceName": "sonic_01", "perDevice": True})]


# ---------------------------------------------------------------------------
# scaffolding helpers
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestCreateEmptyClass:
    def test_it_writes_the_three_boilerplate_classes(self, tmp_path):
        path = _experimentDirectory(tmp_path, "windTunnel")
        CLI._create_empty_class(str(path), "windTunnel")

        source = (path / "code" / "windTunnel.py").read_text()
        assert "class windTunnel(experimentSetupWithData):" in source
        assert "class windTunnelAnalysis(experimentSetupWithData):" in source
        assert "class windTunnelPresentation(experimentSetupWithData):" in source

    def test_the_generated_module_is_syntactically_valid_python(self, tmp_path):
        path = _experimentDirectory(tmp_path, "windTunnel")
        CLI._create_empty_class(str(path), "windTunnel")
        compile((path / "code" / "windTunnel.py").read_text(), "windTunnel.py", "exec")

    def test_it_overwrites_a_previous_version_of_the_file(self, tmp_path):
        path = _experimentDirectory(tmp_path, "windTunnel")
        target = path / "code" / "windTunnel.py"
        target.write_text("# stale content\n")
        CLI._create_empty_class(str(path), "windTunnel")
        assert "stale content" not in target.read_text()


@pytest.mark.unit
class TestMakeRuntimeExperimentData:
    def test_it_creates_the_directory_and_the_configuration(self, tmp_path):
        path = tmp_path / "exp"
        path.mkdir()

        CLI._make_runtimeExperimentData(None, str(path), "exp")

        config = path / "runtimeExperimentData" / "Datasources_Configurations.json"
        assert json.loads(config.read_text()) == {"experimentName": "exp"}

    def test_it_is_safe_to_run_twice(self, tmp_path):
        path = tmp_path / "exp"
        path.mkdir()
        CLI._make_runtimeExperimentData(None, str(path), "exp")
        CLI._make_runtimeExperimentData(None, str(path), "exp")
        assert (path / "runtimeExperimentData").is_dir()

    def test_a_supplied_zip_is_copied_in_under_the_experiment_name(self, tmp_path):
        path = tmp_path / "exp"
        path.mkdir()
        zipFile = tmp_path / "argosDownload.zip"
        zipFile.write_bytes(b"PK\x03\x04not-really-a-zip")

        CLI._make_runtimeExperimentData(str(zipFile), str(path), "exp")

        copied = path / "runtimeExperimentData" / "exp.zip"
        assert copied.read_bytes() == b"PK\x03\x04not-really-a-zip"


@pytest.mark.unit
class TestCreateRepository:
    def test_a_relative_repository_records_a_dot_resource(self, tmp_path):
        path = _experimentDirectory(tmp_path, "exp")
        CLI._create_repository(None, str(path), "exp", relative=True)

        repo = json.loads((path / "exp_repository.json").read_text())
        item = repo["experiment"]["DataSource"]["exp"]
        assert item["isRelativePath"] == "True"
        assert item["item"]["resource"] == "."
        assert item["item"]["dataFormat"] == "string"

    def test_an_absolute_repository_records_the_experiment_path(self, tmp_path):
        path = _experimentDirectory(tmp_path, "exp")
        CLI._create_repository(None, str(path), "exp", relative=False)

        repo = json.loads((path / "exp_repository.json").read_text())
        item = repo["experiment"]["DataSource"]["exp"]
        assert item["isRelativePath"] == "False"
        assert item["item"]["resource"] == str(path)

    def test_without_a_zip_there_are_no_measurements(self, tmp_path):
        path = _experimentDirectory(tmp_path, "exp")
        CLI._create_repository(None, str(path), "exp", relative=True)
        repo = json.loads((path / "exp_repository.json").read_text())
        assert repo["experiment"]["Measurements"] == {}

    def test_a_per_device_entity_gets_a_measurement_keyed_by_device_name(self, tmp_path, monkeypatch):
        path = _experimentDirectory(tmp_path, "exp")
        monkeypatch.setattr(CLI, "ExperimentZipFile", lambda _zip: _FakeZipFile([
            dict(entityTypeName="Sonic", entityName="sonic_01", perDevice=True),
        ]))

        CLI._create_repository("argos.zip", str(path), "exp", relative=True)

        measurements = json.loads((path / "exp_repository.json").read_text())["experiment"]["Measurements"]
        assert list(measurements) == ["sonic_01"]
        item = measurements["sonic_01"]["item"]
        assert item["resource"] == os.path.join("data", "sonic_01.parquet")
        assert item["dataFormat"] == "parquet"
        assert item["desc"] == {"deviceType": "Sonic", "experimentName": "exp",
                                "deviceName": "sonic_01"}

    def test_entities_stored_together_share_one_measurement_per_device_type(self, tmp_path, monkeypatch):
        path = _experimentDirectory(tmp_path, "exp")
        monkeypatch.setattr(CLI, "ExperimentZipFile", lambda _zip: _FakeZipFile([
            dict(entityTypeName="Sonic", entityName="sonic_01", perDevice=False),
            dict(entityTypeName="Sonic", entityName="sonic_02", perDevice=False),
        ]))

        CLI._create_repository("argos.zip", str(path), "exp", relative=True)

        measurements = json.loads((path / "exp_repository.json").read_text())["experiment"]["Measurements"]
        assert list(measurements) == ["Sonic"]
        # the first entity wins the desc; the second never overwrites it
        assert measurements["Sonic"]["item"]["desc"]["deviceName"] == "sonic_01"

    def test_it_also_writes_the_datasources_configuration(self, tmp_path):
        path = _experimentDirectory(tmp_path, "exp")
        CLI._create_repository(None, str(path), "exp", relative=True)
        config = path / "runtimeExperimentData" / "Datasources_Configurations.json"
        assert json.loads(config.read_text()) == {"experimentName": "exp"}

    def test_it_fails_when_the_runtime_directory_does_not_exist_yet(self, tmp_path):
        """Characterisation of B129: the repository JSON is written, but the
        configuration write into the not-yet-created runtimeExperimentData
        directory raises."""
        path = _experimentDirectory(tmp_path, "exp", withRuntime=False)

        with pytest.raises(FileNotFoundError):
            CLI._create_repository(None, str(path), "exp", relative=True)

        assert (path / "exp_repository.json").exists()
        assert not (path / "runtimeExperimentData").exists()


# ---------------------------------------------------------------------------
# create
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestCreateExperiment:
    def _arguments(self, path, **overrides):
        arguments = Namespace(experimentName="windTunnel", path=path, zip=None,
                              relative=False)
        for key, value in overrides.items():
            setattr(arguments, key, value)
        return arguments

    def test_it_creates_the_code_and_data_directories(self, tmp_path, projectCLI):
        target = tmp_path / "windTunnel"
        with pytest.raises(FileNotFoundError):  # B129, pinned below
            CLI.create_experiment(self._arguments(str(target)))

        assert (target / "code").is_dir()
        assert (target / "data").is_dir()
        assert (target / "code" / "windTunnel.py").exists()

    def test_without_a_path_the_experiment_goes_under_the_current_directory(self, tmp_path, projectCLI, monkeypatch):
        monkeypatch.chdir(tmp_path)
        with pytest.raises(FileNotFoundError):  # B129
            CLI.create_experiment(self._arguments(None))
        assert (tmp_path / "windTunnel" / "code").is_dir()

    def test_it_asks_the_project_cli_to_create_a_project_named_after_the_experiment(self, tmp_path, projectCLI):
        target = tmp_path / "windTunnel"
        arguments = self._arguments(str(target))
        with pytest.raises(FileNotFoundError):  # B129
            CLI.create_experiment(arguments)

        assert projectCLI.names() == ["project_create"]
        assert arguments.projectName == "windTunnel"
        assert arguments.directory == str(target)
        assert arguments.overwrite is True
        assert arguments.loadRepositories is True

    def test_an_existing_class_file_is_not_regenerated(self, tmp_path, projectCLI):
        target = tmp_path / "windTunnel"
        (target / "code").mkdir(parents=True)
        (target / "code" / "windTunnel.py").write_text("# hand written\n")

        with pytest.raises(FileNotFoundError):  # B129
            CLI.create_experiment(self._arguments(str(target)))

        assert (target / "code" / "windTunnel.py").read_text() == "# hand written\n"

    @pytest.mark.xfail(
        strict=True,
        reason="B129: create_experiment calls _create_repository before "
               "_make_runtimeExperimentData, but _create_repository writes "
               "runtimeExperimentData/Datasources_Configurations.json -- a "
               "directory only _make_runtimeExperimentData creates. Nothing "
               "earlier in the flow creates it (createProjectDirectory only "
               "writes caseConfiguration.json), so `hera-experiment create` "
               "always aborts with FileNotFoundError, leaving a half-built "
               "experiment and no loaded repository. See the consolidated "
               "findings issue.",
    )
    def test_creating_an_experiment_should_finish(self, tmp_path, projectCLI):
        target = tmp_path / "windTunnel"
        CLI.create_experiment(self._arguments(str(target)))
        assert (target / "windTunnel_repository.json").exists()

    def test_it_currently_aborts_before_loading_the_repository(self, tmp_path, projectCLI):
        """Characterisation of B129."""
        target = tmp_path / "windTunnel"
        with pytest.raises(FileNotFoundError, match="Datasources_Configurations"):
            CLI.create_experiment(self._arguments(str(target)))

        assert projectCLI.names() == ["project_create"]  # repository_load never runs
        assert not (target / "runtimeExperimentData").exists()

    # -- the second bug, reachable only past the first ---------------------

    def _createPastTheFirstBug(self, tmp_path, arguments, monkeypatch):
        """Run create_experiment with runtimeExperimentData pre-created, so
        B129 does not mask what comes after it.

        The chdir is not cosmetic: past B129 the run reaches the
        createNodeRedDeviceMap.sh write, which lands in the cwd (B130) --
        without it the test litters the repository root.
        """
        (tmp_path / "windTunnel" / "runtimeExperimentData").mkdir(parents=True)
        cwd = tmp_path / "elsewhere"
        cwd.mkdir(exist_ok=True)
        monkeypatch.chdir(cwd)
        CLI.create_experiment(arguments)
        return cwd

    def test_it_loads_the_generated_repository_into_the_project(self, tmp_path, projectCLI, monkeypatch):
        target = tmp_path / "windTunnel"
        arguments = self._arguments(str(target))
        self._createPastTheFirstBug(tmp_path, arguments, monkeypatch)

        assert projectCLI.names() == ["project_create", "repository_load"]
        assert arguments.repositoryName == str(target / "windTunnel_repository.json")
        assert (target / "windTunnel_repository.json").exists()

    @pytest.mark.xfail(
        strict=True,
        reason="B130: create_experiment writes createNodeRedDeviceMap.sh "
               "with a bare filename instead of "
               "os.path.join(experiment_path, ...), so with --path the "
               "script lands in the process cwd rather than in the new "
               "experiment directory. See the consolidated findings issue.",
    )
    def test_the_nodered_script_should_land_in_the_experiment_directory(self, tmp_path, projectCLI, monkeypatch):
        target = tmp_path / "windTunnel"
        self._createPastTheFirstBug(tmp_path, self._arguments(str(target)), monkeypatch)

        assert (target / "createNodeRedDeviceMap.sh").exists()

    def test_the_nodered_script_currently_lands_in_the_current_directory(self, tmp_path, projectCLI, monkeypatch):
        """Characterisation of B130."""
        target = tmp_path / "windTunnel"
        cwd = self._createPastTheFirstBug(tmp_path, self._arguments(str(target)), monkeypatch)

        script = cwd / "createNodeRedDeviceMap.sh"
        assert script.read_text() == "argos-experiment-manager nodered createDeviceMap"
        assert os.access(script, os.X_OK)
        assert not (target / "createNodeRedDeviceMap.sh").exists()


# ---------------------------------------------------------------------------
# load
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestLoadExperimentToProject:
    def _experimentWithRepository(self, tmp_path, name="exp"):
        path = tmp_path / name
        path.mkdir()
        (path / f"{name}_repository.json").write_text(json.dumps({"experiment": {}}))
        return path

    def test_it_registers_loads_and_then_removes_the_repository_datasource(self, tmp_path, dataToolkit, projectList):
        path = self._experimentWithRepository(tmp_path)

        CLI.load_experiment_to_project(
            Namespace(experiment=str(path), projectName="EXISTING_PROJECT", overwrite=True)
        )

        toolkit = dataToolkit.instances[0]
        assert toolkit.names() == ["addRepository",
                                   "loadAllDatasourcesInRepositoryToProject",
                                   "deleteDataSource"]
        assert toolkit.calls[0][1] == dict(repositoryName="exp_repository.json",
                                           repositoryPath=str(path / "exp_repository.json"),
                                           overwrite=False)
        assert toolkit.calls[1][1] == ("EXISTING_PROJECT",)
        assert toolkit.calls[1][2] == dict(repositoryName="exp_repository.json",
                                           overwrite=True)
        assert toolkit.calls[2][1] == dict(datasourceName="exp_repository.json")

    def test_the_overwrite_flag_reaches_only_the_load_step(self, tmp_path, dataToolkit, projectList):
        path = self._experimentWithRepository(tmp_path)
        CLI.load_experiment_to_project(
            Namespace(experiment=str(path), projectName="EXISTING_PROJECT", overwrite=False)
        )
        toolkit = dataToolkit.instances[0]
        assert toolkit.calls[0][1]["overwrite"] is False
        assert toolkit.calls[1][2]["overwrite"] is False

    def test_without_an_experiment_path_the_current_directory_is_scanned(self, tmp_path, dataToolkit, projectList, monkeypatch):
        path = self._experimentWithRepository(tmp_path)
        monkeypatch.chdir(path)

        CLI.load_experiment_to_project(
            Namespace(experiment=None, projectName="EXISTING_PROJECT", overwrite=False)
        )

        assert dataToolkit.instances[0].calls[0][1]["repositoryName"] == "exp_repository.json"

    def test_a_directory_without_a_repository_is_rejected(self, tmp_path, dataToolkit, projectList):
        empty = tmp_path / "notAnExperiment"
        empty.mkdir()
        with pytest.raises(ValueError, match="Can't find repository file"):
            CLI.load_experiment_to_project(
                Namespace(experiment=str(empty), projectName="EXISTING_PROJECT", overwrite=False)
            )
        assert dataToolkit.instances == []

    def test_two_repositories_in_one_directory_are_rejected(self, tmp_path, dataToolkit, projectList):
        path = self._experimentWithRepository(tmp_path)
        (path / "other_repository.json").write_text("{}")
        with pytest.raises(ValueError, match="More than 1 repositories"):
            CLI.load_experiment_to_project(
                Namespace(experiment=str(path), projectName="EXISTING_PROJECT", overwrite=False)
            )
        assert dataToolkit.instances == []

    def test_an_unknown_project_is_only_reported_not_refused(self, tmp_path, dataToolkit, projectList):
        """Unlike `data`, `load` is allowed to create the project."""
        path = self._experimentWithRepository(tmp_path)

        CLI.load_experiment_to_project(
            Namespace(experiment=str(path), projectName="BRAND_NEW", overwrite=False)
        )

        assert "BRAND_NEW" not in projectList
        assert dataToolkit.instances[0].calls[1][1] == ("BRAND_NEW",)

    def test_a_failing_load_still_removes_the_temporary_datasource(self, tmp_path, dataToolkit, projectList, monkeypatch):
        path = self._experimentWithRepository(tmp_path)

        def _failingLoad(self, *args, **kwargs):
            self.calls.append(("loadAllDatasourcesInRepositoryToProject", args, kwargs))
            raise ValueError("datasource already exists")

        monkeypatch.setattr(_FakeDataToolkit, "loadAllDatasourcesInRepositoryToProject",
                            _failingLoad)

        CLI.load_experiment_to_project(
            Namespace(experiment=str(path), projectName="EXISTING_PROJECT", overwrite=False)
        )

        assert dataToolkit.instances[0].names()[-1] == "deleteDataSource"
