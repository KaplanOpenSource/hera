"""
Pytest configuration and shared fixtures for the hera test suite.

Architecture
------------
A single **session-scoped** Hera ``Project`` is created at the start of the
test session.  The test repository JSON (``test_repository.json`` located next
to ``data_config.json`` inside ``TEST_HERA``) is loaded into this project via
``dataToolkit.loadAllDatasourcesInRepositoryJSONToProject``.  Each toolkit test
module receives a **real toolkit instance** backed by this project – no
monkey-patching of ``getConfig`` / ``getDataSourceData`` is needed.

Environment
-----------
TEST_HERA : str (env var)
    Path to the test data repository root (default: ~/hera_unittest_data).
    Must contain ``data_config.json``, ``test_repository.json``, and the
    ``measurements/`` + ``expected/`` trees.

RESULT_SET : str (env var or --result-set CLI)
    Name of the expected-output result set under ``expected/`` (default: BASELINE).

PREPARE_EXPECTED_OUTPUT : str (env var)
    Set to "1" to generate expected output files instead of comparing.
"""

import json
import math
import os
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import xarray as xr

# Project name used for ALL toolkit tests in this session
PYTEST_PROJECT_NAME = "PYTEST_HERA_PROJECT"


# ---------------------------------------------------------------------------
# CLI option: --result-set
# ---------------------------------------------------------------------------

def pytest_addoption(parser):
    parser.addoption(
        "--result-set",
        action="store",
        default=None,
        help="Name of the expected-output result set (overrides RESULT_SET env var).",
    )


# ---------------------------------------------------------------------------
# Session-scoped: test data root and configuration
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def test_hera_root():
    """Return the validated root path of the test data repository."""
    root = os.environ.get("TEST_HERA", os.path.expanduser("~/hera_unittest_data"))
    root = Path(root)
    if not root.is_dir():
        pytest.skip(f"TEST_HERA directory not found: {root}")
    return root


@pytest.fixture(scope="session")
def data_config(test_hera_root):
    """Parsed ``data_config.json`` dict."""
    cfg_path = test_hera_root / "data_config.json"
    if not cfg_path.is_file():
        pytest.skip(f"data_config.json not found at {cfg_path}")
    with open(cfg_path, encoding="utf-8") as fh:
        return json.load(fh)


@pytest.fixture(scope="session")
def result_set(request, data_config):
    """Active result-set name (CLI > env var > config default)."""
    rs = request.config.getoption("--result-set")
    if rs is None:
        rs = os.environ.get("RESULT_SET", data_config.get("default_result_set", "BASELINE"))
    return rs


@pytest.fixture(scope="session")
def expected_dir(test_hera_root, result_set):
    """``Path('<test_hera_root>/expected/<result_set>')``."""
    d = test_hera_root / "expected" / result_set
    if not d.is_dir():
        pytest.skip(f"Expected result set directory not found: {d}")
    return d


# ---------------------------------------------------------------------------
# Session-scoped Project fixture  (THE source of truth for toolkit tests)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def hera_test_project(test_hera_root):
    """
    Create a Hera Project and load the test repository into it.

    The project is populated once for the entire test session via
    ``dataToolkit.loadAllDatasourcesInRepositoryJSONToProject``.
    Every toolkit test receives its data from this project.
    """
    from hera.datalayer.project import Project
    from hera.utils.data.toolkit import dataToolkit

    repo_json_path = test_hera_root / "test_repository.json"
    if not repo_json_path.is_file():
        pytest.skip(f"test_repository.json not found at {repo_json_path}")

    # Read the repository JSON
    with open(repo_json_path, encoding="utf-8") as fh:
        repo_json = json.load(fh)

    basedir = str(test_hera_root)

    # Create the project
    proj = Project(projectName=PYTEST_PROJECT_NAME)

    # Load all datasources + configs into the project
    dt = dataToolkit()
    dt.loadAllDatasourcesInRepositoryJSONToProject(
        projectName=PYTEST_PROJECT_NAME,
        repositoryJSON=repo_json,
        basedir=basedir,
        overwrite=True,
    )

    yield proj

    # Teardown: remove all documents created during the session
    try:
        for doc in proj.getMeasurementsDocuments():
            doc.delete()
    except Exception:
        pass


@pytest.fixture(scope="session")
def hera_project_name():
    """The project name string for toolkit constructors."""
    return PYTEST_PROJECT_NAME


# ---------------------------------------------------------------------------
# Per-toolkit session-scoped fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def topo_toolkit(hera_test_project):
    """Session-scoped TopographyToolkit backed by the test project."""
    from hera.measurements.GIS.raster.topography import TopographyToolkit
    return TopographyToolkit(projectName=PYTEST_PROJECT_NAME)


@pytest.fixture(scope="session")
def lc_toolkit(hera_test_project):
    """Session-scoped LandCoverToolkit backed by the test project."""
    from hera.measurements.GIS.raster.landcover import LandCoverToolkit
    return LandCoverToolkit(projectName=PYTEST_PROJECT_NAME)


@pytest.fixture(scope="session")
def demo_toolkit(hera_test_project):
    """Session-scoped DemographyToolkit backed by the test project."""
    from hera.measurements.GIS.vector.demography import DemographyToolkit
    return DemographyToolkit(projectName=PYTEST_PROJECT_NAME)


@pytest.fixture(scope="session")
def lf_toolkit(hera_test_project):
    """Session-scoped lowFreqToolKit backed by the test project."""
    from hera.measurements.meteorology.lowfreqdata.toolkit import lowFreqToolKit
    return lowFreqToolKit(projectName=PYTEST_PROJECT_NAME)


@pytest.fixture(scope="session")
def hf_toolkit(hera_test_project):
    """Session-scoped HighFreqToolKit backed by the test project."""
    from hera.measurements.meteorology.highfreqdata.toolkit import HighFreqToolKit
    return HighFreqToolKit(projectName=PYTEST_PROJECT_NAME)


# ---------------------------------------------------------------------------
# Lightweight fixtures (function-scoped) for DB-specific tests
# ---------------------------------------------------------------------------

@pytest.fixture(scope="function")
def project_fixture():
    """Create a temporary Project and clean up on teardown (for test_datalayer)."""
    from hera.datalayer.project import Project

    project_name = "pytest_temp_project"
    proj = Project(projectName=project_name)
    yield proj
    # Cleanup: remove all documents created during the test
    try:
        for doc in proj.getMeasurementsDocuments():
            doc.delete()
    except Exception:
        pass


@pytest.fixture(scope="session")
def data_toolkit_fixture():
    """Session-scoped ``dataToolkit`` instance."""
    from hera.utils.data.toolkit import dataToolkit
    return dataToolkit()


# ---------------------------------------------------------------------------
# Comparison helpers (ported from run_all_definitions.py)
# ---------------------------------------------------------------------------

def compare_dataframes(df1, df2, rtol=1e-6, atol=1e-6):
    """Deep DataFrame comparison with tolerance for numerics, datetimes, and geometry."""
    try:
        df1 = df1.sort_index(axis=1).reset_index(drop=True)
        df2 = df2.sort_index(axis=1).reset_index(drop=True)

        if list(df1.columns) != list(df2.columns):
            return False

        has_geom = "geometry" in df1.columns and "geometry" in df2.columns
        if has_geom:
            pref = ["areaFraction", "total_pop", "age_0_14", "age_15_19",
                     "age_20_29", "age_30_64", "age_65_up"]
            sort_keys = [c for c in pref if c in df1.columns]
            if not sort_keys:
                sort_keys = [c for c in df1.columns
                             if c != "geometry" and pd.api.types.is_numeric_dtype(df1[c])]
            if sort_keys:
                df1 = df1.sort_values(by=sort_keys).reset_index(drop=True)
                df2 = df2.sort_values(by=sort_keys).reset_index(drop=True)

        for col in df1.columns:
            s1, s2 = df1[col], df2[col]

            if pd.api.types.is_datetime64_any_dtype(s1) and pd.api.types.is_datetime64_any_dtype(s2):
                s1 = pd.to_datetime(s1).dt.tz_localize(None)
                s2 = pd.to_datetime(s2).dt.tz_localize(None)
                if not s1.equals(s2):
                    return False
                continue

            if has_geom and col == "geometry":
                tol_area = float(os.environ.get("GDF_TOL_AREA", "1e-7"))
                for g1, g2 in zip(s1, s2):
                    try:
                        if not g1.equals(g2) and g1.symmetric_difference(g2).area >= tol_area:
                            return False
                    except Exception:
                        return False
                continue

            if pd.api.types.is_numeric_dtype(s1) and pd.api.types.is_numeric_dtype(s2):
                if not np.allclose(s1.fillna(0), s2.fillna(0), rtol=rtol, atol=atol, equal_nan=True):
                    return False
                continue

            if not s1.equals(s2):
                return False

        return True

    except Exception:
        return False


def compare_dataarrays(da1, da2, rtol=1e-6, atol=1e-6):
    """Compare two xarray DataArrays with tolerance."""
    try:
        if not isinstance(da1, xr.DataArray) or not isinstance(da2, xr.DataArray):
            return False
        return np.allclose(da1.values, da2.values, rtol=rtol, atol=atol, equal_nan=True)
    except Exception:
        return False


def deep_compare_with_tolerance(obj1, obj2, rel_tol=1e-6, abs_tol=1e-6):
    """Recursive comparison of nested structures (dict, list, tuple, ndarray, float)."""
    if isinstance(obj1, pd.DataFrame) and isinstance(obj2, pd.DataFrame):
        return compare_dataframes(obj1, obj2, rel_tol, abs_tol)
    if isinstance(obj1, float) and isinstance(obj2, float):
        return math.isclose(obj1, obj2, rel_tol=rel_tol, abs_tol=abs_tol)
    if isinstance(obj1, list) and isinstance(obj2, list):
        if len(obj1) != len(obj2):
            return False
        return all(deep_compare_with_tolerance(i, j, rel_tol, abs_tol) for i, j in zip(obj1, obj2))
    if isinstance(obj1, dict) and isinstance(obj2, dict):
        if set(obj1.keys()) != set(obj2.keys()):
            return False
        return all(deep_compare_with_tolerance(obj1[k], obj2[k], rel_tol, abs_tol) for k in obj1)
    if isinstance(obj1, tuple) and isinstance(obj2, tuple):
        if len(obj1) != len(obj2):
            return False
        return all(deep_compare_with_tolerance(i, j, rel_tol, abs_tol) for i, j in zip(obj1, obj2))
    if isinstance(obj1, np.ndarray) and isinstance(obj2, np.ndarray):
        return np.allclose(obj1, obj2, rtol=rel_tol, atol=abs_tol, equal_nan=True)
    return obj1 == obj2


def compare_outputs(result, expected, output_type):
    """
    Type-based comparison dispatcher.

    Parameters
    ----------
    result : any
        Actual test result.
    expected : any
        Expected reference value.
    output_type : str
        One of: dataframe, geodataframe, ndarray, npz, xarray, dataarray,
        float, int, dict, list, tuple, str, string, bytes, nondbmetadataframe.
    """
    output_type = output_type.lower()

    if output_type == "nondbmetadataframe":
        if hasattr(result, "getData"):
            result = result.getData()
        if hasattr(expected, "getData"):
            expected = expected.getData()
        output_type = "dataframe"

    if hasattr(result, "df") and isinstance(result.df, (pd.DataFrame, gpd.GeoDataFrame)):
        result = result.df
    if hasattr(expected, "df") and isinstance(expected.df, (pd.DataFrame, gpd.GeoDataFrame)):
        expected = expected.df

    if result is None and expected == "null":
        return True
    if expected is None and result == "null":
        return True

    funcs = {
        "dataframe": lambda: isinstance(result, pd.DataFrame) and compare_dataframes(result, expected),
        "metadataframe": lambda: isinstance(result, pd.DataFrame) and compare_dataframes(result, expected),
        "geodataframe": lambda: isinstance(result, (gpd.GeoDataFrame, pd.DataFrame)) and compare_dataframes(result, expected),
        "ndarray": lambda: isinstance(result, np.ndarray) and isinstance(expected, np.ndarray) and np.allclose(result, expected, rtol=1e-6, atol=1e-6, equal_nan=True),
        "npz": lambda: isinstance(result, tuple) and isinstance(expected, tuple) and all(np.allclose(r, e, rtol=1e-6, atol=1e-6, equal_nan=True) for r, e in zip(result, expected)),
        "xarray": lambda: result.equals(expected),
        "dataarray": lambda: compare_dataarrays(result, expected),
        "float": lambda: math.isclose(result, expected, rel_tol=1e-6, abs_tol=1e-6),
        "int": lambda: result == expected,
        "dict": lambda: deep_compare_with_tolerance(result, expected),
        "list": lambda: deep_compare_with_tolerance(result, expected),
        "tuple": lambda: deep_compare_with_tolerance(result, expected),
        "str": lambda: result == expected,
        "string": lambda: result == expected,
        "bytes": lambda: result == expected,
    }

    compare = funcs.get(output_type)
    if compare:
        try:
            ok = compare()
            return bool(ok)
        except Exception:
            return False
    return False


# ---------------------------------------------------------------------------
# Save / Load helpers for expected outputs
# ---------------------------------------------------------------------------

def save_expected_output(filename, data, output_type, expected_dir):
    """Save test output to the expected-results directory."""
    dest = str(expected_dir / filename) if not os.path.isabs(filename) else filename
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    output_type = output_type.lower()

    if data is None:
        with open(dest, "w") as fh:
            fh.write("null")
        return

    handlers = {
        "dataframe": lambda: data.to_parquet(dest, index=False) if dest.endswith(".parquet") else data.to_json(dest, orient="records", indent=2),
        "geodataframe": lambda: data.to_file(dest, driver="GeoJSON"),
        "xarray": lambda: data.to_netcdf(dest),
        "dataarray": lambda: data.to_netcdf(dest),
        "dict": lambda: _write_json(dest, data),
        "list": lambda: _write_json(dest, data),
        "float": lambda: _write_json(dest, float(data)),
        "int": lambda: _write_json(dest, int(data)),
        "str": lambda: _write_text(dest, str(data)),
        "string": lambda: _write_text(dest, str(data)),
        "npz": lambda: np.savez(dest, **({f"arr{i}": a for i, a in enumerate(data)} if isinstance(data, tuple) else {"data": data})),
        "ndarray": lambda: np.savez(dest, **({f"arr{i}": a for i, a in enumerate(data)} if isinstance(data, tuple) else {"data": data})),
    }
    handler = handlers.get(output_type)
    if handler:
        handler()


def load_expected_output(filename, output_type, expected_dir):
    """Load expected output from the expected-results directory."""
    output_type = output_type.lower()
    candidate = str(expected_dir / filename) if not os.path.isabs(filename) else filename

    if not os.path.exists(candidate):
        raise FileNotFoundError(f"Expected output not found: {candidate}")

    loaders = {
        "float": lambda: float(_read_json(candidate)),
        "int": lambda: int(_read_json(candidate)),
        "dict": lambda: _read_json(candidate),
        "list": lambda: _read_json(candidate),
        "dataframe": lambda: pd.read_json(candidate) if candidate.endswith(".json") else pd.read_parquet(candidate),
        "geodataframe": lambda: gpd.read_file(candidate),
        "xarray": lambda: xr.open_dataset(candidate),
        "dataarray": lambda: xr.open_dataarray(candidate),
        "npz": lambda: _read_ndarray(candidate),
        "ndarray": lambda: _read_ndarray(candidate),
        "str": lambda: _read_text(candidate),
        "string": lambda: _read_text(candidate),
    }
    loader = loaders.get(output_type)
    if loader:
        return loader()
    raise ValueError(f"Unknown output_type: {output_type}")


def _write_json(path, obj):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, default=str)


def _write_text(path, text):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


def _read_json(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def _read_text(path):
    with open(path, "r", encoding="utf-8") as fh:
        return fh.read()


def _read_ndarray(path):
    npz = np.load(path)
    if "data" in npz:
        return npz["data"]
    return tuple(npz[f] for f in npz.files)
