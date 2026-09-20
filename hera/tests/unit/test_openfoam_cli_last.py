"""openFoam/CLI.py: the ``--format json`` tail of the two listing handlers.

``test_openfoam_cli.py`` and ``test_openfoam_cli_more.py`` between them
drive every handler in the module and leave exactly two statements
unentered -- the ``ext = "json"`` assignment at the end of the json branch
of ``foam_solver_simulations_list`` and of
``stochasticLagrangian_dispersionFlow_list``.  Their docstrings explain
why: those existing tests hand the handler a real ``pandas.DataFrame``,
and ``DataFrame.to_json()`` goes through pandas' C encoder, which under a
coverage tracer raises ``OverflowError: Maximum recursion level reached``
before the branch finishes.

The branch itself has nothing to do with pandas.  All the handler asks of
the object ``compareWorkflows``/``workflowCompare`` returns is
``to_json()``, ``len()`` and ``str()``, so this file hands over a small
stand-in that answers those three -- which both reaches the branch and
isolates it from pandas entirely.  What is then asserted is the handler's
own contract: that ``--format json`` re-indents the payload, prints it,
and picks ``.json`` as the suffix for a ``--file`` given without one.

The toolkit stand-in is installed by patching ``getToolkit`` on the
ToolkitHome *class* (patching the singleton instance would leave a
permanent instance attribute behind after teardown).  Nothing here
touches MongoDB, OpenFOAM, hermes or Luigi; the handlers that write
relative to the cwd are run with the cwd moved to ``tmp_path``.

Deliberately not covered: everything else in the module, which the two
sibling files already drive, and the thirteen-plus defects they pin.  No
new defect surfaced in this branch -- the json path behaves as its three
siblings (pandas/latex/csv) do.
"""
import json
import os
from argparse import Namespace

import pytest

from hera.simulations.openFoam.CLI import (foam_solver_simulations_list,
                                           stochasticLagrangian_dispersionFlow_list)

PAYLOAD = {"groupName": {"0": "runA"}, "cells": {"0": 1000}}


class _Comparison:
    """The three things the handlers ask of a workflow comparison."""

    def __init__(self, payload=None, length=1):
        self._payload = PAYLOAD if payload is None else payload
        self._length = length

    def to_json(self):
        return json.dumps(self._payload)

    def __len__(self):
        return self._length

    def __str__(self):
        return "<comparison table>"


class _Document:
    def __init__(self, groupName):
        self.desc = dict(groupName=groupName)


class _Toolkit:
    projectName = "UNIT_OPENFOAM"

    def __init__(self, comparison):
        self._comparison = comparison
        self.calls = []

    def getWorkflowListOfSolvers(self, solver):
        self.calls.append(("getWorkflowListOfSolvers", solver))
        return [_Document("runA")]

    def compareWorkflows(self, groupNames, **kwargs):
        self.calls.append(("compareWorkflows", tuple(groupNames), kwargs))
        return self._comparison

    def workflowCompare(self, **kwargs):
        self.calls.append(("workflowCompare", kwargs))
        return self._comparison


@pytest.fixture()
def installToolkit(monkeypatch):
    """Hand the CLI a stand-in toolkit; return the log of getToolkit calls."""
    from hera import toolkitHome

    def install(toolkit):
        asked = []

        def getToolkit(self, toolkitName=None, **kwargs):
            asked.append(dict(toolkitName=toolkitName, kwargs=kwargs))
            return toolkit

        monkeypatch.setattr(type(toolkitHome), "getToolkit", getToolkit)
        return asked

    return install


def _simulationsArguments(**overrides):
    arguments = dict(projectName="UNIT_OPENFOAM", solver="simpleFoam",
                     longFormat=False, transpose=False, format="json", file=None)
    arguments.update(overrides)
    return Namespace(**arguments)


def _dispersionArguments(**overrides):
    arguments = dict(projectName="UNIT_OPENFOAM", format="json", file=None)
    arguments.update(overrides)
    return Namespace(**arguments)


@pytest.mark.unit
class TestSolverSimulationsListAsJson:
    def test_the_comparison_is_reprinted_as_indented_json(
        self, installToolkit, capsys
    ):
        installToolkit(_Toolkit(_Comparison()))
        foam_solver_simulations_list(_simulationsArguments())
        printed = capsys.readouterr().out
        assert json.dumps(PAYLOAD, indent=4) in printed

    def test_the_raw_single_line_json_is_not_what_is_printed(
        self, installToolkit, capsys
    ):
        installToolkit(_Toolkit(_Comparison()))
        foam_solver_simulations_list(_simulationsArguments())
        printed = capsys.readouterr().out
        assert json.dumps(PAYLOAD) not in printed

    def test_a_file_without_a_suffix_is_given_the_json_one(
        self, installToolkit, tmp_path, monkeypatch
    ):
        installToolkit(_Toolkit(_Comparison()))
        monkeypatch.chdir(tmp_path)
        foam_solver_simulations_list(_simulationsArguments(file="comparison"))
        assert os.path.isfile(tmp_path / "comparison.json")

    def test_the_written_file_holds_the_same_payload(
        self, installToolkit, tmp_path, monkeypatch
    ):
        installToolkit(_Toolkit(_Comparison()))
        monkeypatch.chdir(tmp_path)
        foam_solver_simulations_list(_simulationsArguments(file="comparison"))
        assert json.loads((tmp_path / "comparison.json").read_text()) == PAYLOAD

    def test_a_file_that_already_has_a_suffix_keeps_it(
        self, installToolkit, tmp_path, monkeypatch
    ):
        installToolkit(_Toolkit(_Comparison()))
        monkeypatch.chdir(tmp_path)
        foam_solver_simulations_list(_simulationsArguments(file="table.out"))
        assert os.path.isfile(tmp_path / "table.out")
        assert not os.path.exists(tmp_path / "table.out.json")

    def test_an_empty_comparison_is_reported_rather_than_written(
        self, installToolkit, tmp_path, monkeypatch, capsys
    ):
        installToolkit(_Toolkit(_Comparison(payload={}, length=0)))
        monkeypatch.chdir(tmp_path)
        foam_solver_simulations_list(_simulationsArguments(file="comparison"))
        assert "Could not found any workflows" in capsys.readouterr().out
        assert not os.path.exists(tmp_path / "comparison.json")

    def test_the_group_name_read_off_the_documents_reaches_compareworkflows(
        self, installToolkit
    ):
        toolkit = _Toolkit(_Comparison())
        installToolkit(toolkit)
        foam_solver_simulations_list(_simulationsArguments())
        assert ("compareWorkflows", ("runA",),
                dict(longFormat=False, transpose=False)) in toolkit.calls


@pytest.mark.unit
class TestDispersionFlowListAsJson:
    def test_the_comparison_is_reprinted_as_indented_json(
        self, installToolkit, capsys
    ):
        installToolkit(_Toolkit(_Comparison()))
        stochasticLagrangian_dispersionFlow_list(_dispersionArguments())
        assert json.dumps(PAYLOAD, indent=4) in capsys.readouterr().out

    def test_a_file_without_a_suffix_is_given_the_json_one(
        self, installToolkit, tmp_path, monkeypatch
    ):
        installToolkit(_Toolkit(_Comparison()))
        monkeypatch.chdir(tmp_path)
        stochasticLagrangian_dispersionFlow_list(_dispersionArguments(file="flows"))
        assert json.loads((tmp_path / "flows.json").read_text()) == PAYLOAD

    def test_the_workflow_type_asked_for_is_the_stochastic_lagrangian_solver(
        self, installToolkit
    ):
        toolkit = _Toolkit(_Comparison())
        installToolkit(toolkit)
        stochasticLagrangian_dispersionFlow_list(_dispersionArguments())
        assert ("workflowCompare",
                dict(workflowsType="stochasticLagrangianSolver")) in toolkit.calls

    def test_an_empty_comparison_is_reported_rather_than_written(
        self, installToolkit, tmp_path, monkeypatch, capsys
    ):
        installToolkit(_Toolkit(_Comparison(payload={}, length=0)))
        monkeypatch.chdir(tmp_path)
        stochasticLagrangian_dispersionFlow_list(_dispersionArguments(file="flows"))
        assert "Could not found any workflows" in capsys.readouterr().out
        assert not os.path.exists(tmp_path / "flows.json")
