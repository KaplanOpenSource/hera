import os
import subprocess
from pathlib import Path

import pytest


# -----------------------------
# Helpers
# -----------------------------
def _run(cmd, cwd, env):
    """
    Run a CLI command and return (returncode, combined_stdout_stderr).
    """
    p = subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    return p.returncode, p.stdout


def _find_repo_root(start: Path) -> Path:
    """
    Find repo root by walking up until we find:
    - hera/bin/hera-experiment
    - and hera/ folder
    """
    cur = start.resolve()
    for _ in range(10):
        if (cur / "hera" / "bin" / "hera-experiment").exists() and (cur / "hera").is_dir():
            return cur
        if cur.parent == cur:
            break
        cur = cur.parent
    return Path(os.getcwd()).resolve()


# -----------------------------
# Pytest fixtures
# -----------------------------
@pytest.fixture(scope="session")
def hera_repo_root():
    # this file is under: <repo>/hera/tests/...
    return _find_repo_root(Path(__file__))


@pytest.fixture(scope="session")
def env_with_pyargos():
    """
    Ensure PYTHONPATH includes pyargos-master (needed by experiment toolkit).
    """
    env = dict(os.environ)

    pyargos = os.environ.get("PYARGOS_PATH", "")
    if pyargos and os.path.isdir(pyargos):
        env["PYTHONPATH"] = f"{pyargos}:{env.get('PYTHONPATH', '')}"

    return env


@pytest.fixture()
def load_dummy_experiment_to_project(tmp_path, hera_repo_root, env_with_pyargos):
    """
    Creates a unique project name and loads the existing dummy experiment dir into it via CLI.
    Returns dict with project and experiment metadata.
    """
    project = f"unittest_project_dynamic_{os.getpid()}"
    experiment_dir = hera_repo_root / "hera" / "tests" / "dynamic_loading" / "dummy_experiment"

    if not experiment_dir.exists():
        pytest.skip(f"Dummy experiment directory not found at: {experiment_dir}")

    # Load experiment into project (overwrite OK for repeatable runs)
    cmd = [
        "hera-experiment",
        "load",
        project,
        "--experiment",
        str(experiment_dir),
        "--overwrite",
    ]
    code, out = _run(cmd, str(hera_repo_root), env_with_pyargos)
    if code != 0:
        pytest.fail(
            f"Failed to load dummy experiment into project.\n"
            f"CMD: {' '.join(cmd)}\nOUTPUT:\n{out}"
        )

    return {
        "project": project,
        "experiment": "DummyExperiment",
        "experiment_dir": str(experiment_dir),
        "load_output": out,
    }


# -----------------------------
# Tests
# -----------------------------
@pytest.mark.integration
def test_cli_list_and_table_do_not_crash(
    load_dummy_experiment_to_project,
    hera_repo_root,
    env_with_pyargos,
):
    """
    Integration test:
    - hera-experiment list --projectName <proj>  -> must include DummyExperiment
    - hera-experiment table --projectName <proj> -> must not crash, and should look like a table.
    """
    proj = load_dummy_experiment_to_project["project"]

    # 1) list must include DummyExperiment
    code, out = _run(
        ["hera-experiment", "list", "--projectName", proj],
        str(hera_repo_root),
        env_with_pyargos,
    )
    assert code == 0, f"list failed:\n{out}"
    assert "DummyExperiment" in out, f"DummyExperiment not found in list output:\n{out}"

    # 2) table must not crash and must look like a dataframe
    code, out = _run(
        ["hera-experiment", "table", "--projectName", proj],
        str(hera_repo_root),
        env_with_pyargos,
    )
    assert code == 0, f"table failed:\n{out}"

    # Stable assertions (not sensitive to pandas truncation)
    assert "dataFormat" in out, f"Expected 'dataFormat' header in table output:\n{out}"
    assert "experimentToolKit" in out, f"Expected 'experimentToolKit' in table output:\n{out}"
    # typical pandas summary line: "[1 rows x 5 columns]"
    assert "rows x" in out or "columns" in out, f"Expected dataframe-like summary in output:\n{out}"


@pytest.mark.integration
def test_cli_data_returns_sonic_parquet(
    load_dummy_experiment_to_project,
    hera_repo_root,
    env_with_pyargos,
):
    """
    Integration test:
    - hera-experiment data --projectName <proj> DummyExperiment Sonic --deviceName Sonic_1
      should not crash and should print the small Sonic parquet (5 rows).
    """
    proj = load_dummy_experiment_to_project["project"]
    experiment_name = load_dummy_experiment_to_project["experiment"]

    cmd = [
        "hera-experiment",
        "data",
        "--projectName",
        proj,
        experiment_name,
        "Sonic",
        "--deviceName",
        "Sonic_1",
    ]
    code, out = _run(cmd, str(hera_repo_root), env_with_pyargos)

    assert code == 0, f"data command failed:\nCMD: {' '.join(cmd)}\nOUTPUT:\n{out}"

    # בדיקות יציבות בסיסיות על הפלט:
    # אנחנו יודעים שה-parquet כולל deviceName וערכים 0..4 לסוניק_1.
    assert "Sonic_1" in out, f"Expected 'Sonic_1' in data output:\n{out}"
    # בדיקה שהערכים המספריים נמצאים
    for v in ["0.0", "1.0", "2.0", "3.0", "4.0"]:
        assert v in out, f"Expected value {v} in data output:\n{out}"


