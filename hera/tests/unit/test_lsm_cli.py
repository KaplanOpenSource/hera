"""hera/simulations/LSM/CLI.py -- the `hera-LSM` command handlers and the
deferred-import machinery in front of them.

The module keeps four names (`ToolkitHome`, `toolkitHome`, `loadJSON`,
`_get_logger`) at None until `_setup()` runs, and wraps every public
handler in `_lazy_setup` so the first call performs the imports.  `logging`
is a `_LoggingShim` instance that proxies the standard library so that
pre-existing `logging.get_logger(...)` calls still resolve.

Covered: `_setup`, `_lazy_setup` (and its `wrapper`),
`_LoggingShim.get_logger`, `_LoggingShim.__getattr__`,
`_confirm_project_name`, `list_templates`, `setup_template`.

The Namespace each handler receives is built by hand, mirroring the dests
declared in `hera/bin/hera-LSM` (templateName, codeDir, templateVersion,
projectName).  The toolkit is replaced by a recording stand-in bound to
the module's own `toolkitHome` global -- with `_initialized` pinned True,
`_lazy_setup`'s `_setup()` early-returns and leaves the stand-in in place.
Nothing here touches MongoDB, the network, or the real home directory;
`setup_template` runs entirely inside `tmp_path`.

Not covered: the argparse wiring in `hera/bin/hera-LSM`, which is a
script, not an importable module.  (Worth noting even so: its line 41 is
`parsed = parser.parse_args()()` -- it *calls* the Namespace, so
`hera-LSM` cannot dispatch any subcommand at all.  That is outside this
file's scope and is reported separately.)

Bugs pinned here, each with an xfail(strict) for the intended behaviour
and a passing characterisation of what actually happens:

  * B267: `_LoggingShim.get_logger(name)` forwards a bare logger-*name*
    string into hera's `get_logger(instance, name=None)`, which reads its
    first argument as an *instance*.  Every handler's logger therefore
    ends up named "builtins.str".
  * B268: `setup_template`'s "a symlink already exists ... pointing to a
    different file" guard is dead code -- both destinations are
    unconditionally `os.unlink`ed earlier in the same function, so a stale
    link into a different code directory is silently replaced.
  * B269: the sibling "already exists and isn't a symlink" guard is dead
    for the same reason; a real directory in that position escapes as a
    raw `IsADirectoryError` from the earlier `os.unlink`.
  * B270: a clean, successful first-time setup logs two `logger.error`
    records ("file ... not found") for the entirely expected absence of
    the links it is about to create.
"""
import functools
import json
import logging as stdlib_logging
import os
from argparse import Namespace

import pytest

from hera import ToolkitHome
from hera.simulations.LSM import CLI


# ---------------------------------------------------------------------------
# recording stand-ins
# ---------------------------------------------------------------------------

class RecordingLogger:
    """Collects (level, message) pairs instead of configuring logging."""

    def __init__(self):
        self.records = []

    def debug(self, message, *args, **kwargs):
        self.records.append(("debug", str(message)))

    def info(self, message, *args, **kwargs):
        self.records.append(("info", str(message)))

    def warning(self, message, *args, **kwargs):
        self.records.append(("warning", str(message)))

    def error(self, message, *args, **kwargs):
        self.records.append(("error", str(message)))

    def levels(self):
        return [level for level, _ in self.records]

    def messages(self, level=None):
        return [
            message
            for recordLevel, message in self.records
            if level is None or recordLevel == level
        ]


class RecordingShim:
    """Stands in for CLI.logging, handing out one RecordingLogger."""

    def __init__(self, logger):
        self.logger = logger
        self.requestedNames = []

    def get_logger(self, *args, **kwargs):
        self.requestedNames.append(args[0] if args else None)
        return self.logger


class FakeTemplate:
    def __init__(self, templateName, version, dirPath, modelFolder):
        self.templateName = templateName
        self.version = version
        self.dirPath = dirPath
        self.modelFolder = modelFolder


class FakeLSMToolkit:
    def __init__(self, templates=(), byName=None):
        self._templates = list(templates)
        self._byName = dict(byName or {})
        self.getTemplatesCalls = []
        self.getTemplateByNameCalls = []

    def getTemplates(self):
        self.getTemplatesCalls.append({})
        return self._templates

    def getTemplateByName(self, templateName, templateVersion=None):
        self.getTemplateByNameCalls.append((templateName, templateVersion))
        return self._byName.get(templateName)


class RecordingToolkitHome:
    def __init__(self, toolkit):
        self.toolkit = toolkit
        self.calls = []

    def getToolkit(self, **kwargs):
        self.calls.append(kwargs)
        return self.toolkit


def _refuseLoadJSON(*args, **kwargs):
    raise AssertionError(f"loadJSON must not be called; got {args!r}")


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def wiring(monkeypatch):
    """Pin the deferred globals to recording stand-ins.

    `_initialized` is pinned True so that `_lazy_setup`'s `_setup()` call
    early-returns instead of overwriting the stand-ins with the real
    objects.  Every patch is a module-attribute rebind, restored by
    monkeypatch at teardown.
    """
    logger = RecordingLogger()
    shim = RecordingShim(logger)

    monkeypatch.setattr(CLI, "_initialized", True)
    monkeypatch.setattr(CLI, "ToolkitHome", ToolkitHome)
    monkeypatch.setattr(CLI, "logging", shim)
    monkeypatch.setattr(CLI, "loadJSON", _refuseLoadJSON)

    class _Wiring:
        def __init__(self):
            self.logger = logger
            self.shim = shim
            self.home = None

        def install(self, toolkit):
            self.home = RecordingToolkitHome(toolkit)
            monkeypatch.setattr(CLI, "toolkitHome", self.home)
            return self.home

        def loadJSON(self, function):
            monkeypatch.setattr(CLI, "loadJSON", function)

    return _Wiring()


@pytest.fixture()
def modelFolder(tmp_path):
    """An unpacked template directory, as `hera-project` would leave it."""
    folder = tmp_path / "templateModel"
    folder.mkdir()
    return folder


@pytest.fixture()
def codeDir(tmp_path):
    """A code directory holding the compiled model and its met data."""
    root = tmp_path / "code"
    (root / "LSM" / "tozaot" / "Meteorology").mkdir(parents=True)
    (root / "LSM" / "a.out").write_text("#!/bin/false\n")
    return root


def _setupArguments(codeDir, templateName="LSM-v1", templateVersion="1",
                    projectName="UNIT_TEST_PROJECT"):
    return Namespace(
        templateName=templateName,
        codeDir=str(codeDir),
        templateVersion=templateVersion,
        projectName=projectName,
    )


# ---------------------------------------------------------------------------
# _setup
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestSetup:
    def test_it_binds_every_deferred_name_to_the_real_object(self, monkeypatch):
        import hera
        from hera.utils.jsonutils import loadJSON
        from hera.utils.logging import get_logger

        monkeypatch.setattr(CLI, "_initialized", False)
        monkeypatch.setattr(CLI, "ToolkitHome", None)
        monkeypatch.setattr(CLI, "toolkitHome", None)
        monkeypatch.setattr(CLI, "loadJSON", None)
        monkeypatch.setattr(CLI, "_get_logger", None)

        CLI._setup()

        assert CLI.ToolkitHome is hera.ToolkitHome
        assert CLI.toolkitHome is hera.toolkitHome
        assert CLI.loadJSON is loadJSON
        assert CLI._get_logger is get_logger

    def test_it_records_that_it_ran(self, monkeypatch):
        monkeypatch.setattr(CLI, "_initialized", False)
        CLI._setup()
        assert CLI._initialized is True

    def test_a_second_call_rebinds_nothing(self, monkeypatch):
        """The guard is what makes _lazy_setup cheap on every later call."""
        marker = object()
        monkeypatch.setattr(CLI, "_initialized", True)
        monkeypatch.setattr(CLI, "loadJSON", marker)
        monkeypatch.setattr(CLI, "toolkitHome", marker)

        CLI._setup()

        assert CLI.loadJSON is marker
        assert CLI.toolkitHome is marker

    def test_it_is_idempotent_when_called_repeatedly(self, monkeypatch):
        monkeypatch.setattr(CLI, "_initialized", False)
        CLI._setup()
        first = (CLI.ToolkitHome, CLI.toolkitHome, CLI.loadJSON, CLI._get_logger)
        CLI._setup()
        assert (CLI.ToolkitHome, CLI.toolkitHome, CLI.loadJSON, CLI._get_logger) == first


# ---------------------------------------------------------------------------
# _lazy_setup / wrapper
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestLazySetup:
    def test_the_wrapper_runs_setup_before_the_function(self, monkeypatch):
        order = []
        monkeypatch.setattr(CLI, "_setup", lambda: order.append("setup"))

        @CLI._lazy_setup
        def probe():
            order.append("probe")

        probe()
        assert order == ["setup", "probe"]

    def test_setup_runs_on_every_call_not_just_the_first(self, monkeypatch):
        """The idempotence guard lives in _setup, not in the wrapper."""
        calls = []
        monkeypatch.setattr(CLI, "_setup", lambda: calls.append(1))

        @CLI._lazy_setup
        def probe():
            return None

        probe()
        probe()
        assert len(calls) == 2

    def test_positional_and_keyword_arguments_are_forwarded(self, monkeypatch):
        monkeypatch.setattr(CLI, "_setup", lambda: None)
        seen = {}

        @CLI._lazy_setup
        def probe(first, second, third=None):
            seen.update(first=first, second=second, third=third)
            return "done"

        assert probe(1, second=2, third=3) == "done"
        assert seen == {"first": 1, "second": 2, "third": 3}

    def test_the_return_value_is_passed_straight_through(self, monkeypatch):
        monkeypatch.setattr(CLI, "_setup", lambda: None)
        sentinel = object()

        @CLI._lazy_setup
        def probe():
            return sentinel

        assert probe() is sentinel

    def test_exceptions_are_not_swallowed(self, monkeypatch):
        monkeypatch.setattr(CLI, "_setup", lambda: None)

        @CLI._lazy_setup
        def probe():
            raise KeyError("boom")

        with pytest.raises(KeyError, match="boom"):
            probe()

    def test_the_identity_of_the_wrapped_function_survives(self, monkeypatch):
        monkeypatch.setattr(CLI, "_setup", lambda: None)

        @CLI._lazy_setup
        def probe():
            """A docstring worth keeping."""

        assert probe.__name__ == "probe"
        assert probe.__doc__ == "A docstring worth keeping."
        assert probe.__wrapped__.__name__ == "probe"

    def test_the_public_handlers_are_wrapped(self):
        for handler in (CLI._confirm_project_name, CLI.list_templates, CLI.setup_template):
            assert hasattr(handler, "__wrapped__"), handler


# ---------------------------------------------------------------------------
# _LoggingShim
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestLoggingShimGetLogger:
    def test_before_setup_it_falls_back_to_the_standard_library(self, monkeypatch):
        monkeypatch.setattr(CLI, "_get_logger", None)
        shim = CLI._LoggingShim()
        assert shim.get_logger("hera.bin.probe") is stdlib_logging.getLogger("hera.bin.probe")

    def test_the_fallback_logger_carries_the_requested_name(self, monkeypatch):
        monkeypatch.setattr(CLI, "_get_logger", None)
        assert CLI._LoggingShim().get_logger("hera.bin.probe").name == "hera.bin.probe"

    def test_after_setup_it_delegates_to_heras_own_get_logger(self, monkeypatch):
        calls = []
        sentinel = object()

        def fake(*args, **kwargs):
            calls.append((args, kwargs))
            return sentinel

        monkeypatch.setattr(CLI, "_get_logger", fake)
        assert CLI._LoggingShim().get_logger("a.b", name="c") is sentinel
        assert calls == [(("a.b",), {"name": "c"})]

    def test_the_module_level_logging_name_is_a_shim_instance(self):
        assert isinstance(CLI.logging, CLI._LoggingShim)

    @pytest.mark.xfail(
        strict=True,
        reason="B267: _LoggingShim.get_logger forwards its arguments "
               "verbatim to hera's get_logger(instance, name=None), but the "
               "callers in this module pass a logger *name* as the single "
               "positional argument. hera reads it as an instance and "
               "returns getClassLogger(type(name)), so every hera-LSM "
               "handler logs to a logger named 'builtins.str' instead of "
               "'hera.bin.hera_lsm.load_template'. The shim exists to keep "
               "name-based calls working, which is exactly what it breaks. "
               "See the consolidated findings issue.",
    )
    def test_a_name_based_call_should_produce_a_logger_of_that_name(self):
        CLI._setup()
        logger = CLI.logging.get_logger("hera.bin.hera_lsm.load_template")
        assert logger.name == "hera.bin.hera_lsm.load_template"

    def test_a_name_based_call_currently_produces_a_logger_named_builtins_str(self):
        """Characterisation of B267."""
        CLI._setup()
        assert CLI.logging.get_logger("hera.bin.hera_lsm.load_template").name == "builtins.str"


@pytest.mark.unit
class TestLoggingShimGetattr:
    def test_a_standard_library_constant_is_proxied(self):
        shim = CLI._LoggingShim()
        assert shim.DEBUG == stdlib_logging.DEBUG
        assert shim.ERROR == stdlib_logging.ERROR

    def test_a_standard_library_function_is_proxied_by_identity(self):
        assert CLI._LoggingShim().getLogger is stdlib_logging.getLogger

    def test_an_unknown_attribute_raises_attributeerror(self):
        with pytest.raises(AttributeError, match="noSuchAttribute"):
            CLI._LoggingShim().noSuchAttribute

    def test_get_logger_is_the_shims_own_and_not_proxied(self):
        """The standard library has no get_logger, so the override is the point."""
        assert not hasattr(stdlib_logging, "get_logger")
        assert callable(CLI._LoggingShim().get_logger)


# ---------------------------------------------------------------------------
# _confirm_project_name
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestConfirmProjectName:
    def test_an_explicit_project_name_is_left_alone(self, wiring):
        arguments = Namespace(projectName="EXPLICIT_PROJECT")
        CLI._confirm_project_name(arguments, wiring.logger)
        assert arguments.projectName == "EXPLICIT_PROJECT"

    def test_an_explicit_project_name_reads_no_file(self, wiring):
        """loadJSON is wired to refuse, so any read fails the test."""
        CLI._confirm_project_name(Namespace(projectName="EXPLICIT_PROJECT"), wiring.logger)

    def test_an_explicit_project_name_logs_nothing(self, wiring):
        CLI._confirm_project_name(Namespace(projectName="EXPLICIT_PROJECT"), wiring.logger)
        assert wiring.logger.records == []

    def test_a_missing_project_name_is_read_from_the_case_configuration(
        self, wiring, tmp_path, monkeypatch
    ):
        from hera.utils.jsonutils import loadJSON

        wiring.loadJSON(loadJSON)
        monkeypatch.chdir(tmp_path)
        (tmp_path / "caseConfiguration.json").write_text(
            json.dumps({"projectName": "FROM_THE_CASE_FILE", "other": 1})
        )

        arguments = Namespace(projectName=None)
        CLI._confirm_project_name(arguments, wiring.logger)

        assert arguments.projectName == "FROM_THE_CASE_FILE"

    def test_the_file_it_reads_is_named_caseConfiguration_json(self, wiring):
        seen = []

        def fake(path):
            seen.append(path)
            return {"projectName": "P"}

        wiring.loadJSON(fake)
        CLI._confirm_project_name(Namespace(projectName=None), wiring.logger)
        assert seen == ["caseConfiguration.json"]

    def test_the_fallback_is_announced_at_debug_level(self, wiring):
        wiring.loadJSON(lambda path: {"projectName": "P"})
        CLI._confirm_project_name(Namespace(projectName=None), wiring.logger)
        assert wiring.logger.levels() == ["debug"]
        assert "caseConfiguration.json" in wiring.logger.messages("debug")[0]

    def test_a_case_configuration_without_a_project_name_raises(self, wiring):
        wiring.loadJSON(lambda path: {"somethingElse": 1})
        with pytest.raises(KeyError, match="projectName"):
            CLI._confirm_project_name(Namespace(projectName=None), wiring.logger)

    def test_a_missing_case_configuration_file_raises(self, wiring, tmp_path, monkeypatch):
        from hera.utils.jsonutils import loadJSON

        wiring.loadJSON(loadJSON)
        monkeypatch.chdir(tmp_path)
        with pytest.raises(ValueError, match="caseConfiguration.json"):
            CLI._confirm_project_name(Namespace(projectName=None), wiring.logger)


# ---------------------------------------------------------------------------
# list_templates
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestListTemplates:
    def test_an_empty_project_says_so_and_stops(self, wiring, capsys):
        wiring.install(FakeLSMToolkit(templates=[]))
        CLI.list_templates(Namespace(projectName="P"))
        assert capsys.readouterr().out.strip() == "There are no templates"

    def test_the_toolkit_is_requested_by_the_lsm_constant(self, wiring):
        home = wiring.install(FakeLSMToolkit(templates=[]))
        CLI.list_templates(Namespace(projectName="P"))
        assert home.calls == [{"toolkitName": ToolkitHome.LSM, "projectName": "P"}]

    def test_the_project_name_is_forwarded_to_the_toolkit(self, wiring):
        home = wiring.install(FakeLSMToolkit(templates=[]))
        CLI.list_templates(Namespace(projectName="SOME_PROJECT"))
        assert home.calls[0]["projectName"] == "SOME_PROJECT"

    def test_a_missing_project_name_is_taken_from_the_case_configuration(self, wiring):
        home = wiring.install(FakeLSMToolkit(templates=[]))
        wiring.loadJSON(lambda path: {"projectName": "FROM_THE_CASE_FILE"})
        CLI.list_templates(Namespace(projectName=None))
        assert home.calls[0]["projectName"] == "FROM_THE_CASE_FILE"

    def test_every_template_is_listed_with_its_four_fields(self, wiring, capsys):
        wiring.install(
            FakeLSMToolkit(
                templates=[
                    FakeTemplate("LSM-v1", [1, 0, 0], "/dir/one", "/model/one"),
                    FakeTemplate("LSM-v2", [2, 0, 0], "/dir/two", "/model/two"),
                ]
            )
        )
        CLI.list_templates(Namespace(projectName="P"))
        out = capsys.readouterr().out

        for expected in (
            "LSM-v1", "[1, 0, 0]", "/dir/one", "/model/one",
            "LSM-v2", "[2, 0, 0]", "/dir/two", "/model/two",
        ):
            assert expected in out

    def test_a_non_empty_listing_prints_a_header(self, wiring, capsys):
        wiring.install(
            FakeLSMToolkit(templates=[FakeTemplate("LSM-v1", "1", "/d", "/m")])
        )
        CLI.list_templates(Namespace(projectName="P"))
        out = capsys.readouterr().out
        assert "List of template" in out
        assert "There are no templates" not in out

    def test_the_listing_labels_each_field(self, wiring, capsys):
        wiring.install(
            FakeLSMToolkit(templates=[FakeTemplate("LSM-v1", "1", "/d", "/m")])
        )
        CLI.list_templates(Namespace(projectName="P"))
        out = capsys.readouterr().out
        assert "version:" in out
        assert "path:" in out
        assert "model folder:" in out

    def test_it_asks_the_toolkit_for_templates_exactly_once(self, wiring):
        toolkit = FakeLSMToolkit(templates=[])
        wiring.install(toolkit)
        CLI.list_templates(Namespace(projectName="P"))
        assert len(toolkit.getTemplatesCalls) == 1


# ---------------------------------------------------------------------------
# setup_template
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestSetupTemplateRefusals:
    def test_an_unknown_template_is_refused_with_repository_advice(
        self, wiring, codeDir
    ):
        wiring.install(FakeLSMToolkit(byName={}))
        with pytest.raises(RuntimeError, match="updateRepositories"):
            CLI.setup_template(_setupArguments(codeDir))

    def test_the_refusal_names_the_template_that_was_asked_for(self, wiring, codeDir):
        wiring.install(FakeLSMToolkit(byName={}))
        with pytest.raises(RuntimeError, match="LSM-v1"):
            CLI.setup_template(_setupArguments(codeDir, templateName="LSM-v1"))

    def test_the_requested_name_and_version_reach_the_toolkit(self, wiring, codeDir):
        toolkit = FakeLSMToolkit(byName={})
        wiring.install(toolkit)
        with pytest.raises(RuntimeError):
            CLI.setup_template(
                _setupArguments(codeDir, templateName="LSM-v1", templateVersion="7")
            )
        assert toolkit.getTemplateByNameCalls == [("LSM-v1", "7")]

    def test_nothing_is_created_when_the_template_is_unknown(
        self, wiring, codeDir, modelFolder
    ):
        wiring.install(FakeLSMToolkit(byName={}))
        with pytest.raises(RuntimeError):
            CLI.setup_template(_setupArguments(codeDir))
        assert list(modelFolder.iterdir()) == []


@pytest.mark.unit
class TestSetupTemplateHappyPath:
    @pytest.fixture()
    def run(self, wiring, codeDir, modelFolder):
        def _run(templateName="LSM-v1"):
            wiring.install(
                FakeLSMToolkit(
                    byName={
                        templateName: FakeTemplate(
                            templateName, "1", "/dir", str(modelFolder)
                        )
                    }
                )
            )
            CLI.setup_template(_setupArguments(codeDir, templateName=templateName))

        return _run

    def test_the_output_directories_are_created(self, run, modelFolder):
        run()
        assert (modelFolder / "tozaot").is_dir()
        assert (modelFolder / "tozaot" / "machsan").is_dir()

    def test_the_executable_is_symlinked_from_the_code_directory(
        self, run, modelFolder, codeDir
    ):
        run()
        link = modelFolder / "a.out"
        assert link.is_symlink()
        assert os.readlink(link) == str(codeDir / "LSM" / "a.out")

    def test_the_meteorology_directory_is_symlinked_from_the_code_directory(
        self, run, modelFolder, codeDir
    ):
        run()
        link = modelFolder / "tozaot" / "Meteorology"
        assert link.is_symlink()
        assert os.readlink(link) == str(codeDir / "LSM" / "tozaot" / "Meteorology")

    def test_the_code_subdirectory_is_the_template_name_before_the_first_dash(
        self, wiring, codeDir, modelFolder
    ):
        """'LSM-v1' and 'LSM-v4-general' must both resolve to <codeDir>/LSM."""
        wiring.install(
            FakeLSMToolkit(
                byName={
                    "LSM-v4-general": FakeTemplate(
                        "LSM-v4-general", "4", "/dir", str(modelFolder)
                    )
                }
            )
        )
        CLI.setup_template(_setupArguments(codeDir, templateName="LSM-v4-general"))
        assert os.readlink(modelFolder / "a.out") == str(codeDir / "LSM" / "a.out")

    def test_the_links_resolve_to_the_real_files(self, run, modelFolder):
        run()
        assert (modelFolder / "a.out").resolve().read_text() == "#!/bin/false\n"
        assert (modelFolder / "tozaot" / "Meteorology").resolve().is_dir()

    def test_running_it_twice_leaves_the_same_links(self, run, modelFolder):
        run()
        first = (
            os.readlink(modelFolder / "a.out"),
            os.readlink(modelFolder / "tozaot" / "Meteorology"),
        )
        run()
        assert (
            os.readlink(modelFolder / "a.out"),
            os.readlink(modelFolder / "tozaot" / "Meteorology"),
        ) == first

    def test_the_successful_setup_is_announced(self, run, wiring):
        run()
        messages = " ".join(wiring.logger.messages("info"))
        assert "Found template LSM-v1" in messages
        assert "setup directories successfully" in messages


@pytest.mark.unit
class TestSetupTemplateGuardsAreDead:
    """Both symlink destinations are os.unlink'ed before the guards run."""

    @pytest.fixture()
    def install(self, wiring, modelFolder):
        def _install(templateName="LSM-v1"):
            wiring.install(
                FakeLSMToolkit(
                    byName={
                        templateName: FakeTemplate(
                            templateName, "1", "/dir", str(modelFolder)
                        )
                    }
                )
            )

        return _install

    @pytest.mark.xfail(
        strict=True,
        reason="B268: setup_template raises 'a symlink already exists at "
               "... pointing to a different file' only if os.path.exists() "
               "is true at that point -- but the same function has already "
               "run os.unlink() on that exact path a few lines earlier, so "
               "the destination is always gone and the guard is dead code. "
               "A template still linked against an old code directory is "
               "silently repointed instead of being reported. See the "
               "consolidated findings issue.",
    )
    def test_a_stale_link_to_a_different_code_directory_should_be_reported(
        self, install, codeDir, modelFolder, tmp_path
    ):
        install()
        other = tmp_path / "otherCode"
        (other / "LSM").mkdir(parents=True)
        (other / "LSM" / "a.out").write_text("stale\n")
        (modelFolder / "a.out").symlink_to(other / "LSM" / "a.out")

        with pytest.raises(Exception, match="pointing to a different file"):
            CLI.setup_template(_setupArguments(codeDir))

    def test_a_stale_link_is_currently_replaced_without_a_word(
        self, install, codeDir, modelFolder, tmp_path, wiring
    ):
        """Characterisation of B268."""
        install()
        other = tmp_path / "otherCode"
        (other / "LSM").mkdir(parents=True)
        (other / "LSM" / "a.out").write_text("stale\n")
        (modelFolder / "a.out").symlink_to(other / "LSM" / "a.out")

        CLI.setup_template(_setupArguments(codeDir))

        assert os.readlink(modelFolder / "a.out") == str(codeDir / "LSM" / "a.out")
        assert "pointing to a different file" not in " ".join(
            wiring.logger.messages()
        )

    @pytest.mark.xfail(
        strict=True,
        reason="B269: the sibling guard 'Couldn't create symlink, ... "
               "already exists and isn't a symlink' is dead for the same "
               "reason -- os.unlink() runs first. The one thing os.unlink "
               "cannot remove is a directory, so the case the message "
               "describes escapes as a bare IsADirectoryError from the "
               "unlink instead of the intended diagnostic. See the "
               "consolidated findings issue.",
    )
    def test_a_directory_in_the_executables_place_should_be_reported(
        self, install, codeDir, modelFolder
    ):
        install()
        (modelFolder / "a.out").mkdir()
        with pytest.raises(Exception, match="isn't a symlink"):
            CLI.setup_template(_setupArguments(codeDir))

    def test_a_directory_in_the_executables_place_currently_raises_isadirectoryerror(
        self, install, codeDir, modelFolder
    ):
        """Characterisation of B269."""
        install()
        (modelFolder / "a.out").mkdir()
        with pytest.raises(IsADirectoryError):
            CLI.setup_template(_setupArguments(codeDir))

    def test_a_plain_file_in_the_executables_place_is_silently_removed(
        self, install, codeDir, modelFolder
    ):
        """Characterisation of B269: a regular file unlinks cleanly."""
        install()
        (modelFolder / "a.out").write_text("a stale build\n")
        CLI.setup_template(_setupArguments(codeDir))
        assert (modelFolder / "a.out").is_symlink()


@pytest.mark.unit
class TestSetupTemplateLogsErrorsOnSuccess:
    @pytest.fixture()
    def install(self, wiring, modelFolder):
        wiring.install(
            FakeLSMToolkit(
                byName={"LSM-v1": FakeTemplate("LSM-v1", "1", "/dir", str(modelFolder))}
            )
        )
        return wiring

    @pytest.mark.xfail(
        strict=True,
        reason="B270: the first-ever setup of a template necessarily has "
               "no tozaot/Meteorology link and no a.out to unlink, yet both "
               "FileNotFoundError handlers call logger.error('file ... not "
               "found'). A completely successful run reports two errors, "
               "which makes the log useless for spotting real ones. See the "
               "consolidated findings issue.",
    )
    def test_a_clean_successful_setup_should_log_no_errors(self, install, codeDir):
        CLI.setup_template(_setupArguments(codeDir))
        assert install.logger.messages("error") == []

    def test_a_clean_successful_setup_currently_logs_two_errors(self, install, codeDir):
        """Characterisation of B270."""
        CLI.setup_template(_setupArguments(codeDir))
        errors = install.logger.messages("error")
        assert len(errors) == 2
        assert all("not found" in message for message in errors)
