"""openFoam/CLI.py: the three branches test_openfoam_cli.py left uncovered.

That file drives every handler in the module and pins thirteen defects in
them; what its 98% left behind are three conditional branches nothing entered:

* ``stochasticLagrangian_source_makeEscapedMassFile`` with ``--dt`` given --
  the resampling path, which joins a regular time axis onto the log and
  interpolates the mass differences onto it;
* ``stochasticLagrangian_postProcess_toVTK`` with no ``projectName``
  attribute -- the caseConfiguration.json fallback;
* ``stochasticLagrangian_dispersionFlow_list`` with ``--file`` omitted -- the
  print-only path.

The handlers are driven with hand-built argparse Namespaces, mirroring the
dests ``hera/bin/hera-openFoam`` declares, and the toolkit is a recording
stand-in installed by patching ``getToolkit`` on the ToolkitHome *class*
(patching the singleton instance would leave a permanent instance attribute
behind on teardown).  Nothing here touches MongoDB, OpenFOAM, Luigi or
paraview; everything written lands under ``tmp_path``, which is also the cwd
for the handlers that write relative to it.

``toolkitHome.OF_LSM`` has to be supplied for the mass-file handler to run at
all -- that is B177, already pinned in the sibling file; it is patched onto
the class here so the branch behind it can be reached, exactly as the sibling
file does.

The ``--format json`` branch of the list handler is deliberately left to the
sibling file: ``DataFrame.to_json()`` goes through pandas' C encoder, which
counts against the interpreter's recursion limit, and with a tracer installed
(``pytest --cov``) it raises ``OverflowError: Maximum recursion level
reached`` before returning anything.  Every existing assertion on that branch
already fails under coverage for that reason, and adding another buys nothing.
The ``--file`` omitted branch this file targets is reached through the csv and
latex formats instead.

Bug pinned here (with a strict xfail for the intended behaviour and a passing
characterisation of what happens today):

* B295 ``stochasticLagrangian_source_makeEscapedMassFile --dt`` writes an
  OpenFOAM ``scalarField`` whose declared record count does not match the
  number of records it then writes.  The count is ``len(data)`` -- the frame
  *after* the regular axis has been outer-joined onto the log, so it holds the
  union of both time sets -- while the values are written by looping over
  ``timesteps``, which is only ``numpy.arange(min, max, dt)`` and excludes the
  final log time.  A log at t=1,2,3 with ``--dt 0.5`` produces "5" followed by
  four values, which OpenFOAM refuses to parse ("size 5 is not equal to the
  given value of 4").  Without ``--dt`` the two are the same list, which is
  why the covered path never showed it.
"""
import json
import os
from argparse import Namespace

import pandas
import pytest

from hera.simulations.openFoam.CLI import (
    stochasticLagrangian_dispersionFlow_list,
    stochasticLagrangian_postProcess_toVTK,
    stochasticLagrangian_source_makeEscapedMassFile,
)


# ---------------------------------------------------------------------------
# stand-ins
# ---------------------------------------------------------------------------

class _Recorder:
    def __init__(self):
        self.calls = []

    def call(self, name):
        return [entry for entry in self.calls if entry[0] == name]


class _FakeMassAnalysis:
    def __init__(self, calls, data):
        self._calls = calls
        self._data = data

    def getMassFromLog(self, logFile, solver="StochasticLagrangianSolver"):
        self._calls.append(("getMassFromLog", dict(logFile=logFile, solver=solver)))
        return self._data.copy()


class _FakeLSMToolkit(_Recorder):
    def __init__(self, data):
        super().__init__()
        self.analysis = _FakeMassAnalysis(self.calls, data)


class _FakeVTKToolkit(_Recorder):
    """The shape stochasticLagrangian_postProcess_toVTK drives."""

    def __init__(self, documents=()):
        super().__init__()
        self._documents = list(documents)
        self.stochasticLagrangian = self
        self.presentation = self

    def getWorkflowDocumentFromDB(self, name):
        self.calls.append(("getWorkflowDocumentFromDB", dict(name=name)))
        return self._documents

    def getCaseResults(self, **kwargs):
        self.calls.append(("getCaseResults", kwargs))
        return pandas.DataFrame(dict(x=[0.0], y=[1.0], z=[2.0]))

    def toUnstructuredVTK(self, **kwargs):
        self.calls.append(("toUnstructuredVTK", kwargs))


class _FakeCompareToolkit(_Recorder):
    def __init__(self, frame):
        super().__init__()
        self._frame = frame

    def workflowCompare(self, **kwargs):
        self.calls.append(("workflowCompare", kwargs))
        return self._frame

    def compareWorkflows(self, **kwargs):  # the name the module should be using (B174)
        self.calls.append(("compareWorkflows", kwargs))
        return self._frame


@pytest.fixture()
def installToolkit(monkeypatch):
    """Make getToolkit hand the CLI our stand-in; return the request log."""
    from hera import toolkitHome

    def install(toolkit):
        asked = []

        def getToolkit(self, toolkitName=None, filesDirectory=None, **kwargs):
            asked.append(dict(toolkitName=toolkitName, filesDirectory=filesDirectory,
                              kwargs=kwargs))
            return toolkit

        monkeypatch.setattr(type(toolkitHome), "getToolkit", getToolkit)
        return asked

    return install


@pytest.fixture()
def ofLSMConstant(monkeypatch):
    """Supply the constant B177 says is missing, so the handler can run."""
    from hera import toolkitHome

    monkeypatch.setattr(type(toolkitHome), "OF_LSM", "OF_LSM", raising=False)


# ---------------------------------------------------------------------------
# stochasticLagrangian_source_makeEscapedMassFile --dt
# ---------------------------------------------------------------------------

def _massLog():
    """Three escaped-mass records one second apart: diffs 0, 10, 15."""
    return pandas.DataFrame({"time": [1.0, 2.0, 3.0],
                             "mass": [0.0, 10.0, 25.0],
                             "filterType": ["outlet"] * 3,
                             "action": ["escaped"] * 3})


def _massArguments(casePath, **overrides):
    arguments = dict(casePath=casePath, patch="outlet", massFileName=None, dt=None,
                     logFile="log.solver", solver="myFoam", action="escaped")
    arguments.update(overrides)
    return Namespace(**arguments)


def _massValues(text):
    """The value lines of the emitted scalarField, in order."""
    body = text.split("//")[-1]
    return [line.strip() for line in body.splitlines()
            if line.strip() and line.strip() not in ("(", ")")]


@pytest.fixture()
def massCase(tmp_path):
    (tmp_path / "constant").mkdir()
    return tmp_path


@pytest.mark.unit
class TestMakeEscapedMassFileWithADeltaT:
    def test_the_values_are_resampled_onto_the_regular_axis(self, massCase, installToolkit,
                                                            ofLSMConstant):
        installToolkit(_FakeLSMToolkit(_massLog()))
        stochasticLagrangian_source_makeEscapedMassFile(
            _massArguments(str(massCase), dt=0.5))
        # arange(1, 3, 0.5) -> 1.0, 1.5, 2.0, 2.5; the diffs 0, 10, 15 sit at
        # 1, 2, 3, so the midpoints interpolate to 5 and 12.5
        assert [float(value) for value in _massValues(
            (massCase / "constant" / "outletMass").read_text())[1:]] == \
            [0.0, 5.0, 10.0, 12.5]

    def test_the_final_log_time_is_outside_the_half_open_axis(self, massCase,
                                                              installToolkit,
                                                              ofLSMConstant):
        """numpy.arange stops before the maximum, so t=3 (diff 15) is dropped."""
        installToolkit(_FakeLSMToolkit(_massLog()))
        stochasticLagrangian_source_makeEscapedMassFile(
            _massArguments(str(massCase), dt=0.5))
        values = [float(value) for value in _massValues(
            (massCase / "constant" / "outletMass").read_text())[1:]]
        assert 15.0 not in values

    def test_a_delta_t_larger_than_the_log_span_gives_a_single_sample(self, massCase,
                                                                      installToolkit,
                                                                      ofLSMConstant):
        installToolkit(_FakeLSMToolkit(_massLog()))
        stochasticLagrangian_source_makeEscapedMassFile(
            _massArguments(str(massCase), dt=5))
        values = _massValues((massCase / "constant" / "outletMass").read_text())[1:]
        assert len(values) == 1

    def test_a_delta_t_given_as_a_string_is_accepted(self, massCase, installToolkit,
                                                     ofLSMConstant):
        """argparse hands strings through unless the parser says otherwise."""
        installToolkit(_FakeLSMToolkit(_massLog()))
        stochasticLagrangian_source_makeEscapedMassFile(
            _massArguments(str(massCase), dt="1.0"))
        values = [float(value) for value in _massValues(
            (massCase / "constant" / "outletMass").read_text())[1:]]
        assert values == [0.0, 10.0]

    def test_the_file_still_lands_in_the_constant_directory(self, massCase, installToolkit,
                                                            ofLSMConstant):
        installToolkit(_FakeLSMToolkit(_massLog()))
        stochasticLagrangian_source_makeEscapedMassFile(
            _massArguments(str(massCase), dt=0.5))
        assert (massCase / "constant" / "outletMass").is_file()

    def test_the_openfoam_header_is_unchanged_by_the_resampling(self, massCase,
                                                                installToolkit,
                                                                ofLSMConstant):
        installToolkit(_FakeLSMToolkit(_massLog()))
        stochasticLagrangian_source_makeEscapedMassFile(
            _massArguments(str(massCase), dt=0.5))
        written = (massCase / "constant" / "outletMass").read_text()
        assert "class       scalarField;" in written
        assert written.rstrip().endswith(")")

    def test_the_log_is_still_filtered_by_patch_and_action_first(self, massCase,
                                                                 installToolkit,
                                                                 ofLSMConstant):
        log = pandas.concat([_massLog(),
                             pandas.DataFrame({"time": [1.0, 2.0, 3.0],
                                               "mass": [900.0, 900.0, 900.0],
                                               "filterType": ["inlet"] * 3,
                                               "action": ["escaped"] * 3})],
                            ignore_index=True)
        installToolkit(_FakeLSMToolkit(log))
        stochasticLagrangian_source_makeEscapedMassFile(
            _massArguments(str(massCase), dt=0.5))
        values = [float(value) for value in _massValues(
            (massCase / "constant" / "outletMass").read_text())[1:]]
        assert all(value < 100.0 for value in values)

    @pytest.mark.xfail(
        strict=True,
        reason="B295: with --dt the declared record count is len(data) -- "
               "the frame after the regular axis has been outer-joined onto "
               "the log times, i.e. the union of both -- while the values are "
               "written by looping over `timesteps`, which is only "
               "numpy.arange(min, max, dt) and excludes the last log time.  A "
               "log at t=1,2,3 with --dt 0.5 therefore emits '5' followed by "
               "four values, and OpenFOAM refuses to read the field.  Without "
               "--dt the two are the same list, which is why the covered path "
               "never showed it.  See the consolidated findings issue.",
    )
    def test_the_declared_count_should_match_the_number_of_values(self, massCase,
                                                                  installToolkit,
                                                                  ofLSMConstant):
        installToolkit(_FakeLSMToolkit(_massLog()))
        stochasticLagrangian_source_makeEscapedMassFile(
            _massArguments(str(massCase), dt=0.5))
        values = _massValues((massCase / "constant" / "outletMass").read_text())
        assert int(values[0]) == len(values) - 1

    def test_the_count_currently_overstates_the_values_by_the_dropped_time(
            self, massCase, installToolkit, ofLSMConstant):
        """Characterisation of B295."""
        installToolkit(_FakeLSMToolkit(_massLog()))
        stochasticLagrangian_source_makeEscapedMassFile(
            _massArguments(str(massCase), dt=0.5))
        values = _massValues((massCase / "constant" / "outletMass").read_text())
        assert int(values[0]) == 5
        assert len(values) - 1 == 4

    def test_without_a_delta_t_the_count_and_the_values_agree(self, massCase,
                                                              installToolkit,
                                                              ofLSMConstant):
        """Characterisation of B295: --dt is the trigger."""
        installToolkit(_FakeLSMToolkit(_massLog()))
        stochasticLagrangian_source_makeEscapedMassFile(_massArguments(str(massCase)))
        values = _massValues((massCase / "constant" / "outletMass").read_text())
        assert int(values[0]) == len(values) - 1 == 3


# ---------------------------------------------------------------------------
# stochasticLagrangian_postProcess_toVTK without a projectName attribute
# ---------------------------------------------------------------------------

def _vtkArguments(**overrides):
    """The parser always declares --projectName, so a Namespace without it is
    what a *programmatic* caller produces -- which is the only way to reach
    the caseConfiguration.json fallback (see B167 for the parser side)."""
    arguments = dict(dispersionName="myDispersion", cloudName="kinematicCloud",
                     overwrite=False)
    arguments.update(overrides)
    return Namespace(**arguments)


@pytest.fixture()
def vtkCwd(tmp_path, monkeypatch):
    """The handler writes VTK/ relative to the process cwd."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "myDispersion").mkdir()
    return tmp_path


@pytest.mark.unit
class TestPostProcessToVTKProjectNameFallback:
    def test_the_project_name_is_read_from_caseconfiguration_json(self, vtkCwd,
                                                                  installToolkit):
        (vtkCwd / "caseConfiguration.json").write_text(json.dumps(dict(projectName="FROM_FILE")))
        asked = installToolkit(_FakeVTKToolkit())
        stochasticLagrangian_postProcess_toVTK(_vtkArguments())
        assert asked[0]["kwargs"]["projectName"] == "FROM_FILE"

    def test_the_configuration_file_name_is_overridable(self, vtkCwd, installToolkit):
        (vtkCwd / "other.json").write_text(json.dumps(dict(projectName="FROM_OTHER")))
        asked = installToolkit(_FakeVTKToolkit())
        stochasticLagrangian_postProcess_toVTK(
            _vtkArguments(configurationFile="other.json"))
        assert asked[0]["kwargs"]["projectName"] == "FROM_OTHER"

    def test_a_missing_configuration_file_is_not_silently_ignored(self, vtkCwd,
                                                                  installToolkit):
        installToolkit(_FakeVTKToolkit())
        with pytest.raises(Exception):
            stochasticLagrangian_postProcess_toVTK(_vtkArguments())

    def test_the_openfoam_toolkit_is_the_one_asked_for(self, vtkCwd, installToolkit):
        (vtkCwd / "caseConfiguration.json").write_text(json.dumps(dict(projectName="P")))
        asked = installToolkit(_FakeVTKToolkit())
        stochasticLagrangian_postProcess_toVTK(_vtkArguments())
        from hera import toolkitHome

        assert asked[0]["toolkitName"] == toolkitHome.SIMULATIONS_OPENFOAM

    def test_the_dispersion_name_is_used_as_a_directory_when_the_db_has_nothing(
            self, vtkCwd, installToolkit):
        (vtkCwd / "caseConfiguration.json").write_text(json.dumps(dict(projectName="P")))
        toolkit = _FakeVTKToolkit()
        installToolkit(toolkit)
        stochasticLagrangian_postProcess_toVTK(_vtkArguments())
        assert toolkit.call("getCaseResults")[0][1]["caseDescriptor"] == "myDispersion"

    def test_the_vtk_directory_is_created_under_the_process_cwd(self, vtkCwd,
                                                                installToolkit):
        (vtkCwd / "caseConfiguration.json").write_text(json.dumps(dict(projectName="P")))
        installToolkit(_FakeVTKToolkit())
        stochasticLagrangian_postProcess_toVTK(_vtkArguments())
        assert os.path.isdir(vtkCwd / "VTK" / "myDispersion")

    def test_the_results_are_not_cached_on_this_route(self, vtkCwd, installToolkit):
        (vtkCwd / "caseConfiguration.json").write_text(json.dumps(dict(projectName="P")))
        toolkit = _FakeVTKToolkit()
        installToolkit(toolkit)
        stochasticLagrangian_postProcess_toVTK(_vtkArguments())
        assert toolkit.call("getCaseResults")[0][1]["cache"] is False

    def test_the_cloud_name_becomes_the_vtk_file_name(self, vtkCwd, installToolkit):
        (vtkCwd / "caseConfiguration.json").write_text(json.dumps(dict(projectName="P")))
        toolkit = _FakeVTKToolkit()
        installToolkit(toolkit)
        stochasticLagrangian_postProcess_toVTK(_vtkArguments(cloudName="myCloud"))
        assert toolkit.call("toUnstructuredVTK")[0][1]["filename"] == "myCloud"

    def test_the_frame_that_came_back_is_the_one_written(self, vtkCwd, installToolkit):
        (vtkCwd / "caseConfiguration.json").write_text(json.dumps(dict(projectName="P")))
        toolkit = _FakeVTKToolkit()
        installToolkit(toolkit)
        stochasticLagrangian_postProcess_toVTK(_vtkArguments())
        written = toolkit.call("toUnstructuredVTK")[0][1]
        assert list(written["data"].columns) == ["x", "y", "z"]
        assert (written["xcoord"], written["ycoord"], written["zcoord"]) == ("x", "y", "z")


# ---------------------------------------------------------------------------
# stochasticLagrangian_dispersionFlow_list without --file
# ---------------------------------------------------------------------------

def _listArguments(**overrides):
    arguments = dict(projectName="P", format="csv", file=None)
    arguments.update(overrides)
    return Namespace(**arguments)


@pytest.fixture()
def comparisonFrame():
    return pandas.DataFrame(dict(groupName=["flow", "flow"], alpha=[1.0, 2.0]))


@pytest.mark.unit
class TestDispersionFlowListWithoutAFile:
    def test_nothing_is_written_when_no_file_was_asked_for(self, installToolkit, tmp_path,
                                                           monkeypatch, comparisonFrame):
        monkeypatch.chdir(tmp_path)
        installToolkit(_FakeCompareToolkit(comparisonFrame))
        stochasticLagrangian_dispersionFlow_list(_listArguments())
        assert os.listdir(tmp_path) == []

    def test_the_comparison_is_printed_instead(self, installToolkit, capsys,
                                               comparisonFrame):
        installToolkit(_FakeCompareToolkit(comparisonFrame))
        stochasticLagrangian_dispersionFlow_list(_listArguments())
        printed = capsys.readouterr().out
        assert "groupName" in printed and "alpha" in printed

    def test_the_csv_format_prints_comma_separated_values(self, installToolkit, capsys,
                                                          comparisonFrame):
        installToolkit(_FakeCompareToolkit(comparisonFrame))
        stochasticLagrangian_dispersionFlow_list(_listArguments(format="csv"))
        assert "groupName,alpha" in capsys.readouterr().out

    def test_the_latex_format_prints_a_tabular_environment(self, installToolkit, capsys,
                                                           comparisonFrame):
        installToolkit(_FakeCompareToolkit(comparisonFrame))
        stochasticLagrangian_dispersionFlow_list(_listArguments(format="latex"))
        assert "tabular" in capsys.readouterr().out

    def test_an_empty_comparison_says_so_instead_of_printing_nothing(self, installToolkit,
                                                                     capsys):
        installToolkit(_FakeCompareToolkit(pandas.DataFrame(columns=["groupName"])))
        stochasticLagrangian_dispersionFlow_list(_listArguments())
        assert "Could not found any workflows to compare" in capsys.readouterr().out

    def test_the_project_name_is_forwarded_to_gettoolkit(self, installToolkit,
                                                         comparisonFrame):
        asked = installToolkit(_FakeCompareToolkit(comparisonFrame))
        stochasticLagrangian_dispersionFlow_list(_listArguments(projectName="MY_PROJECT"))
        assert asked[0]["kwargs"]["projectName"] == "MY_PROJECT"

    def test_only_the_stochastic_lagrangian_workflows_are_compared(self, installToolkit,
                                                                   comparisonFrame):
        toolkit = _FakeCompareToolkit(comparisonFrame)
        installToolkit(toolkit)
        stochasticLagrangian_dispersionFlow_list(_listArguments())
        assert toolkit.calls[0][1] == dict(workflowsType="stochasticLagrangianSolver")
