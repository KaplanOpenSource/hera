"""riskassessment/CLI.py: createRepository and _check_if_agent."""
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
class TestCreateRepository:
    def test_a_valid_agent_file_should_appear_in_the_repository(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path.parent)  # ensure cwd != tmp_path
        _write_agent(tmp_path / "agent1.json")
        createRepository(_Args(str(tmp_path), "repo.json"))
        repo = json.loads((tmp_path / "repo.json").read_text())
        assert "agent1" in repo["RiskAssessment"]["DataSource"]

    def test_a_non_agent_file_is_left_out(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path.parent)
        (tmp_path / "plain.json").write_text(json.dumps({"foo": "bar"}))
        createRepository(_Args(str(tmp_path), "repo.json"))
        repo = json.loads((tmp_path / "repo.json").read_text())
        assert repo["RiskAssessment"]["DataSource"] == {}

    def test_a_missing_path_defaults_to_the_current_working_directory(self, tmp_path, monkeypatch):
        """When arguments.path is falsy, it defaults to os.getcwd()."""
        monkeypatch.chdir(tmp_path)
        _write_agent(tmp_path / "agent1.json")
        args = _Args(path=None, repository_name="repo.json")
        createRepository(args)
        assert args.path == str(tmp_path)
        repo = json.loads((tmp_path / "repo.json").read_text())
        assert "agent1" in repo["RiskAssessment"]["DataSource"]
