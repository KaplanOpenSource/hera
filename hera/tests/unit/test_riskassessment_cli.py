"""riskassessment/CLI.py: createRepository and _check_if_agent.

B109: createRepository always produces an empty repository. It scans
`os.listdir(arguments.path)`, which yields bare filenames, and passes
those straight to `_check_if_agent(file_path)` without joining them with
`arguments.path` -- so `_check_if_agent` tries to open the file relative
to the current working directory, not the directory being scanned. Unless
the process happens to already be running with `arguments.path` as its
cwd, every file lookup fails and DataSource ends up empty.
"""
import json
import os

import pytest

from hera.riskassessment.CLI import _check_if_agent, createRepository


def _write_agent(path, **extra):
    desc = {"name": "agent1", "effectParameters": {"tenbergeCoefficient": 2}}
    desc.update(extra)
    path.write_text(json.dumps(desc))
    return desc


@pytest.mark.unit
class TestCheckIfAgent:
    def test_a_valid_agent_descriptor_is_returned(self, tmp_path):
        path = tmp_path / "agent1.json"
        desc = _write_agent(path)
        assert _check_if_agent(str(path)) == desc

    def test_json_without_effectparameters_returns_none(self, tmp_path):
        path = tmp_path / "plain.json"
        path.write_text(json.dumps({"foo": "bar"}))
        assert _check_if_agent(str(path)) is None

    def test_a_file_that_is_not_json_returns_false(self, tmp_path):
        path = tmp_path / "notjson.txt"
        path.write_text("not json at all")
        assert _check_if_agent(str(path)) is False

    def test_a_nonexistent_file_returns_false(self, tmp_path):
        assert _check_if_agent(str(tmp_path / "nope.json")) is False


class _Args:
    def __init__(self, path, repository_name):
        self.path = path
        self.repository_name = repository_name


@pytest.mark.unit
class TestCreateRepositoryIsBroken:
    @pytest.mark.xfail(
        strict=True,
        reason="B109: createRepository passes bare os.listdir() filenames "
               "to _check_if_agent without joining them with "
               "arguments.path, so the lookup is relative to the process "
               "cwd instead of the scanned directory -- the DataSource "
               "dict is always empty unless cwd happens to equal path. "
               "See the consolidated findings issue.",
    )
    def test_a_valid_agent_file_should_appear_in_the_repository(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path.parent)  # ensure cwd != tmp_path
        _write_agent(tmp_path / "agent1.json")
        createRepository(_Args(str(tmp_path), "repo.json"))
        repo = json.loads((tmp_path / "repo.json").read_text())
        assert "agent1" in repo["RiskAssessment"]["DataSource"]

    def test_it_currently_produces_an_empty_repository(self, tmp_path, monkeypatch):
        """Characterisation of B109."""
        monkeypatch.chdir(tmp_path.parent)
        _write_agent(tmp_path / "agent1.json")
        createRepository(_Args(str(tmp_path), "repo.json"))
        repo = json.loads((tmp_path / "repo.json").read_text())
        assert repo["RiskAssessment"]["DataSource"] == {}

    def test_a_missing_path_defaults_to_the_current_working_directory(self, tmp_path, monkeypatch):
        """When arguments.path is falsy, it defaults to os.getcwd() -- and
        then the lookup accidentally works, since listdir() and
        _check_if_agent() now agree on the same (implicit) directory."""
        monkeypatch.chdir(tmp_path)
        _write_agent(tmp_path / "agent1.json")
        args = _Args(path=None, repository_name="repo.json")
        createRepository(args)
        assert args.path == str(tmp_path)
        repo = json.loads((tmp_path / "repo.json").read_text())
        assert "agent1" in repo["RiskAssessment"]["DataSource"]
