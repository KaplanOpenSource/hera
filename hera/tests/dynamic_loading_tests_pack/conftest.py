import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
import json
import pytest
from pathlib import Path


def _ensure_pyargos_in_path():
    """
    We don't want hard-coded absolute paths inside the tests logic,
    but we can allow user to provide PYTHONPATH externally.
    """
    # If user already exported PYTHONPATH, we rely on it.
    pass


def _write_minimal_argos_zip(zip_path: Path):
    """
    Create a minimal Argos-like zip that contains data.json with keys
    that Hera experiment loader expects (trialSets/entityTypes etc.).
    """
    payload = {
        "trialSets": [],
        "entityTypes": [],
        # Keep extra keys harmless; some parsers expect them to exist
        "entities": [],
        "experimentName": zip_path.stem,
    }
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("data.json", json.dumps(payload, indent=2))


def _write_minimal_experiment_folder(base_dir: Path, name: str) -> Path:
    """
    Creates a minimal experiment directory structure that matches what
    hera-experiment load expects.
    """
    exp_dir = base_dir / name
    (exp_dir / "code").mkdir(parents=True, exist_ok=True)
    (exp_dir / "data").mkdir(parents=True, exist_ok=True)
    (exp_dir / "runtimeExperimentData").mkdir(parents=True, exist_ok=True)

    # 1) minimal experiment class file
    #    IMPORTANT: name must match class name and file name
    code = f"""
from hera.measurements.experiment.experiment import experimentSetupWithData

class {name}(experimentSetupWithData):
    def __init__(self, projectName, pathToExperiment, filesDirectory):
        super().__init__(projectName=projectName, pathToExperiment=pathToExperiment, filesDirectory=filesDirectory)
"""
    (exp_dir / "code" / f"{name}.py").write_text(code)

    # 2) runtimeExperimentData config
    (exp_dir / "runtimeExperimentData" / "Datasources_Configurations.json").write_text(
        json.dumps({"experimentName": name}, indent=2)
    )

    # 3) minimal zip (Argos style)
    _write_minimal_argos_zip(exp_dir / "runtimeExperimentData" / f"{name}.zip")

    # 4) repository json that hera-experiment load searches for: "*_repository.json"
    repo = {
        "experiment": {
            "DataSource": {
                name: {
                    "isRelativePath": "False",
                    "item": {
                        "dataSourceName": name,
                        "resource": str(exp_dir),  # absolute or relative is fine; we use absolute here
                        "dataFormat": "string",
                        "overwrite": "True",
                    },
                }
            },
            "Measurements": {},
        }
    }
    (exp_dir / f"{name}_repository.json").write_text(json.dumps(repo, indent=2))

    return exp_dir


def _run_cmd(cmd, cwd=None, env=None):
    res = subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    return res.returncode, res.stdout


@pytest.fixture(scope="session")
def hera_repo_root() -> Path:
    # assumes tests executed from repo root or within it
    # find by walking up until we see "hera" package dir
    here = Path.cwd().resolve()
    for p in [here] + list(here.parents):
        if (p / "hera").is_dir() and (p / "hera" / "__init__.py").exists():
            return p
    return here


@pytest.fixture()
def temp_project_name() -> str:
    # avoid collisions; use pid
    return f"unittest_project_dynamic_{os.getpid()}"


@pytest.fixture()
def dummy_experiment_dir(tmp_path: Path) -> Path:
    # creates temp dummy experiment
    return _write_minimal_experiment_folder(tmp_path, "DummyExperiment")


@pytest.fixture()
def env_with_pyargos():
    env = os.environ.copy()
    # user already exports PYTHONPATH typically; keep it.
    # but keep it safe if empty:
    env["PYTHONPATH"] = env.get("PYTHONPATH", "")
    return env


@pytest.fixture()
def load_dummy_experiment_to_project(hera_repo_root, dummy_experiment_dir, temp_project_name, env_with_pyargos):
    """
    Loads dummy experiment to temp project using CLI.
    """
    cmd = [
        "hera-experiment",
        "load",
        temp_project_name,
        "--experiment",
        str(dummy_experiment_dir),
        "--overwrite",
    ]
    code, out = _run_cmd(cmd, cwd=str(hera_repo_root), env=env_with_pyargos)
    if code != 0:
        pytest.fail(f"hera-experiment load failed.\nCMD: {' '.join(cmd)}\nOUTPUT:\n{out}")

    yield {
        "project": temp_project_name,
        "experiment": "DummyExperiment",
        "output": out,
        "experiment_dir": str(dummy_experiment_dir),
    }

    # Teardown: remove ALL DB documents for the project (incl. the hidden
    # ``<name>__config__`` doc, so it no longer appears in getProjectList),
    # then delete its files directory (~/.hera/<project_name>).
    try:
        from hera.datalayer.collection import AbstractCollection
        AbstractCollection().deleteDocuments(projectName=temp_project_name)
    except Exception:
        pass
    project_files_dir = os.path.join(os.path.expanduser("~"), ".hera", temp_project_name)
    shutil.rmtree(project_files_dir, ignore_errors=True)
