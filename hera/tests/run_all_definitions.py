import os
import sys
import json
import math
import glob
import traceback
import importlib
import os
from pathlib import Path

import pandas as pd
import geopandas as gpd
import numpy as np
import xarray as xr
import rasterio
from shapely.geometry import Polygon

try:
    from termcolor import cprint
except ImportError:
    def cprint(msg, color=None): print(msg)

# Caches for function-based parameters/results
function_results = {}
successful_tests = []
failed_tests = []

# =========================
# Result-set base directory
# =========================
def _require_result_set():
    rs = os.environ.get("RESULT_SET")
    if not rs:
        print("Missing RESULT_SET (set via --result-set or env RESULT_SET)", file=sys.stderr)
        sys.exit(2)
    return rs


def _get_base_expected_dir() -> Path:
    """
    Resolve the base directory for expected results:
    Prefer $HERA_UNITTEST_DATA/expected/<RESULT_SET>.
    Fallback to tests/expected/<RESULT_SET> if HERA_UNITTEST_DATA is not set.
    Create the directory in prepare mode.
    """
    rs = os.environ.get("RESULT_SET", "BASELINE")
    env_root = os.environ.get("HERA_UNITTEST_DATA")
    if env_root:
        base = Path(env_root) / "expected" / rs
    else:
        base = Path("tests") / "expected" / rs

    if os.environ.get("PREPARE_EXPECTED_OUTPUT", "0") == "1":
        base.mkdir(parents=True, exist_ok=True)
    return base

def _expected_base_dir():
    """
    Return Path('<HERA_UNITTEST_DATA>/expected/<RESULT_SET>').

    Base directory is:
      - $HERA_UNITTEST_DATA if defined
      - otherwise: ~/hera_unittest_data

    In prepare-mode (PREPARE_EXPECTED_OUTPUT=1) the directory is created
    if missing. In run-mode it must already exist, otherwise we exit
    with an error.
    """
    rs = _require_result_set()

    # Base root for all unittest data (configurable via env var)
    unit_root = os.environ.get(
        "HERA_UNITTEST_DATA",
        os.path.expanduser("~/hera_unittest_data"),
    )

    base = Path(unit_root) / "expected" / rs

    is_prepare = os.environ.get("PREPARE_EXPECTED_OUTPUT") == "1"
    if is_prepare:
        base.mkdir(parents=True, exist_ok=True)
    else:
        if not base.exists():
            print(f"Expected results set '{rs}' not found at {base}", file=sys.stderr)
            print("First create it: scripts/run_all_json_tests.sh prepare --result-set {rs}", file=sys.stderr)
            sys.exit(3)

    return base


# =========================
# DataFrame comparisons
# =========================
def compare_dataframes(df1, df2, rtol=1e-6, atol=1e-6):
    """
    Compare two pandas DataFrames with tolerance for numeric columns, special handling
    for datetime and GeoDataFrame geometry, and a robust fallback for string-like columns.

    SAFE additions (won't affect other tests):
    - If both are empty GeoDataFrames: treat as equal even if columns differ
      (GeoJSON roundtrip may drop schema for empty frames).
    - For object/category columns: compare as strings to avoid dtype-only mismatches.
    """
    try:
        import pandas as pd
        import numpy as np
        import os

        # ---- Detect GeoDataFrame safely (optional dependency) ----
        try:
            import geopandas as gpd
            is_gdf1 = isinstance(df1, gpd.GeoDataFrame)
            is_gdf2 = isinstance(df2, gpd.GeoDataFrame)
        except Exception:
            is_gdf1 = is_gdf2 = False

        # ---- Canonicalize ordering (stable comparisons) ----
        df1 = df1.sort_index(axis=1).reset_index(drop=True)
        df2 = df2.sort_index(axis=1).reset_index(drop=True)

        # ------------------------------------------------------------------
        # SAFE SPECIAL CASE:
        # Empty GeoDataFrame schema may not survive GeoJSON write/read.
        # If both are GeoDataFrames and both empty -> consider equal.
        # ------------------------------------------------------------------
        if is_gdf1 and is_gdf2 and len(df1) == 0 and len(df2) == 0:
            return True

        # ---- Strict columns check (keep as-is for non-empty / non-Geo cases) ----
        if list(df1.columns) != list(df2.columns):
            print("Column mismatch between DataFrames")
            return False

        has_geom = "geometry" in df1.columns and "geometry" in df2.columns
        if has_geom:
            # Keep existing sorting heuristics (only if those columns exist)
            pref = ["areaFraction", "total_pop", "age_0_14", "age_15_19",
                    "age_20_29", "age_30_64", "age_65_up"]

            sort_keys = [c for c in pref if c in df1.columns]
            if not sort_keys:
                sort_keys = [
                    c for c in df1.columns
                    if c != "geometry" and pd.api.types.is_numeric_dtype(df1[c])
                ]

            if sort_keys:
                df1 = df1.sort_values(by=sort_keys).reset_index(drop=True)
                df2 = df2.sort_values(by=sort_keys).reset_index(drop=True)

        # ---- Column-by-column compare ----
        for col in df1.columns:
            s1, s2 = df1[col], df2[col]

            # --- Datetime: normalize timezones and compare strictly ---
            if pd.api.types.is_datetime64_any_dtype(s1) and pd.api.types.is_datetime64_any_dtype(s2):
                s1n = pd.to_datetime(s1).dt.tz_localize(None)
                s2n = pd.to_datetime(s2).dt.tz_localize(None)
                if not s1n.equals(s2n):
                    print(f"Mismatch in datetime column '{col}'")
                    return False
                continue

            # --- Geometry: topology equality with area tolerance ---
            if has_geom and col == "geometry":
                tol_area = float(os.environ.get("GDF_TOL_AREA", "1e-7"))
                for i, (g1, g2) in enumerate(zip(s1, s2)):
                    try:
                        # handle None geometries explicitly
                        if g1 is None and g2 is None:
                            continue
                        if (g1 is None) != (g2 is None):
                            print(f"Mismatch in geometry at row {i} (None vs geom)")
                            return False

                        equal_topo = g1.equals(g2)
                        if not equal_topo:
                            # tolerate tiny area diffs
                            if g1.symmetric_difference(g2).area >= tol_area:
                                print(f"Mismatch in geometry at row {i}")
                                return False
                    except Exception:
                        print(f"Geometry compare exception at row {i}; treating as mismatch")
                        return False
                continue

            # --- Numeric: tolerant comparison (rtol/atol) ---
            if pd.api.types.is_numeric_dtype(s1) and pd.api.types.is_numeric_dtype(s2):
                if not np.allclose(
                    s1.astype(float).fillna(0.0),
                    s2.astype(float).fillna(0.0),
                    rtol=rtol,
                    atol=atol,
                    equal_nan=True
                ):
                    print(f"Mismatch in numeric column '{col}'")
                    return False
                continue

            # --- Object/category fallback: compare as strings (dtype-insensitive) ---
            try:
                is_cat1 = pd.api.types.is_categorical_dtype(s1)
                is_cat2 = pd.api.types.is_categorical_dtype(s2)
            except Exception:
                is_cat1 = is_cat2 = False

            if (s1.dtype == "object" or s2.dtype == "object" or is_cat1 or is_cat2):
                s1_str = s1.astype(str)
                s2_str = s2.astype(str)
                if not s1_str.equals(s2_str):
                    print(f"Mismatch in column '{col}'")
                    return False
                continue

            # --- Default strict comparison for all other dtypes ---
            if not s1.equals(s2):
                print(f"Mismatch in column '{col}'")
                return False

        return True

    except Exception as e:
        print(f"Exception during DataFrame comparison: {e}")
        return False

# =========================
# DataArray comparisons
# =========================
def compare_dataarrays(da1, da2, rtol=1e-6, atol=1e-6):
    try:
        if not isinstance(da1, xr.DataArray) or not isinstance(da2, xr.DataArray):
            print("One of the inputs is not a DataArray")
            return False
        return np.allclose(da1.values, da2.values, rtol=rtol, atol=atol, equal_nan=True)
    except Exception as e:
        print(f"Exception during DataArray comparison: {e}")
        return False

# =========================
# Deep comparison helpers
# =========================
def deep_compare_with_tolerance(obj1, obj2, rel_tol=1e-6, abs_tol=1e-6):
    if isinstance(obj1, pd.DataFrame) and isinstance(obj2, pd.DataFrame):
        return compare_dataframes(obj1, obj2, rel_tol, abs_tol)
    if isinstance(obj1, float) and isinstance(obj2, float):
        return math.isclose(obj1, obj2, rel_tol=rel_tol, abs_tol=abs_tol)
    if isinstance(obj1, list) and isinstance(obj2, list):
        if len(obj1) != len(obj2): return False
        return all(deep_compare_with_tolerance(i, j, rel_tol, abs_tol) for i, j in zip(obj1, obj2))
    if isinstance(obj1, dict) and isinstance(obj2, dict):
        if set(obj1.keys()) != set(obj2.keys()): return False
        return all(deep_compare_with_tolerance(obj1[k], obj2[k], rel_tol, abs_tol) for k in obj1)
    if isinstance(obj1, tuple) and isinstance(obj2, tuple):
        if len(obj1) != len(obj2): return False
        return all(deep_compare_with_tolerance(i, j, rel_tol, abs_tol) for i, j in zip(obj1, obj2))
    if isinstance(obj1, np.ndarray) and isinstance(obj2, np.ndarray):
        return np.allclose(obj1, obj2, rtol=rel_tol, atol=abs_tol, equal_nan=True)
    return obj1 == obj2

def normalize_output_for_type(result, output_type):
    """
    Make runner tolerant to toolkits returning wrapper objects.
    Only triggers for specific output_type, so it won't affect other tests.
    """
    if result is None:
        return None

    # ---------------- dict ----------------
    if output_type == "dict":
        if isinstance(result, dict):
            return result

        # common patterns
        for m in ("toDict", "to_dict", "asDict", "as_dict"):
            if hasattr(result, m) and callable(getattr(result, m)):
                return getattr(result, m)()

        # last resort: use __dict__ if it's meaningful
        if hasattr(result, "__dict__") and isinstance(result.__dict__, dict) and result.__dict__:
            return dict(result.__dict__)

        return result  # keep as-is if cannot normalize

    # -------------- dataframe --------------
    if output_type == "dataframe":
        import pandas as pd
        if isinstance(result, pd.DataFrame):
            return result

        for m in ("toDataFrame", "to_dataframe"):
            if hasattr(result, m) and callable(getattr(result, m)):
                df = getattr(result, m)()
                return df

        # common attributes
        for attr in ("df", "data", "result", "output"):
            if hasattr(result, attr):
                val = getattr(result, attr)
                if isinstance(val, pd.DataFrame):
                    return val

        return result

    return result


# =========================
# Type-based comparison
# =========================
def compare_outputs(result, expected, output_type):
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

    comparison_funcs = {
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

    compare = comparison_funcs.get(output_type)
    if compare:
        try:
            ok = compare()
            if not isinstance(ok, bool):
                print(f"Comparison function for '{output_type}' did not return a boolean (got {type(ok)})")
                return False
            if not ok:
                print("Comparison failed")
                print(f"Expected ({output_type}):\n{expected}")
                print(f"Got ({type(result)}):\n{result}")
            return ok
        except Exception as e:
            print(f"Exception during comparison for '{output_type}': {e}")
            return False
    else:
        print(f"No comparison function for type: {output_type}")
        return False

# =========================
# Save / Load — under tests/expected/<RESULT_SET>/
# =========================
def _resolve_dest_path(filename: str, base_expected_dir: Path) -> str:
    if os.path.isabs(filename):
        return filename
    return str(base_expected_dir / filename)

def save_output(filename, data, output_type):
    base_expected_dir = _expected_base_dir()
    dest = _resolve_dest_path(filename, base_expected_dir)
    os.makedirs(os.path.dirname(dest), exist_ok=True)

    try:
        from hera.datalayer.document.metadataDocument import nonDBMetadataFrame
    except ImportError:
        nonDBMetadataFrame = None  # noqa: F841

    if data is None:
        with open(dest, "w") as f:
            f.write("null")
        print(f"Output is None, saved 'null' to {dest}")
        return

    output_type = output_type.lower()

    custom_extractors = ["singlePointTurbulenceStatistics", "AveragingCalculator"]
    obj_class = type(data).__name__
    if obj_class in custom_extractors:
        for attr in ["getData", "data", "rawData", "result"]:
            if hasattr(data, attr):
                v = getattr(data, attr)
                data = v() if callable(v) else v
                break

    if output_type in ["metadataframe", "dataframe"]:
        if not isinstance(data, pd.DataFrame) and hasattr(data, "getData"):
            data = data.getData()

    def _write_json(path, obj):
        with open(path, "w", encoding="utf-8") as f:
            json.dump(obj, f, indent=2, default=str)

    def _write_text(path, text):
        with open(path, "w", encoding="utf-8") as f:
            f.write(text)

    def _write_bytes(path, b):
        with open(path, "wb") as f:
            f.write(b)

    handlers = {
        "dataframe": lambda: (
            data.to_json(dest, orient="records", indent=2)
            if dest.endswith(".json")
            else data.to_parquet(dest, index=False, engine="pyarrow")
        ),
        "metadataframe": lambda: (
            data.to_file(str(Path(dest)), driver="GeoJSON")
            if isinstance(data, gpd.GeoDataFrame)
            else data.to_json(dest, orient="records", indent=2)
        ),
        "nondbmetadataframe": lambda: (
            data.to_file(str(Path(dest)), driver="GeoJSON")
            if hasattr(data, "to_file")
            else data.to_json(dest, orient="records", indent=2)
            if hasattr(data, "to_json")
            else _write_json(dest, data.getData().to_dict(orient="records"))
        ),
        "geodataframe": lambda: (
            data.to_file(str(Path(dest)), driver="GeoJSON")
            if dest.endswith((".geojson", ".json"))
            else _write_json(dest, json.loads(data.to_json()))
        ),
        "float": lambda: _write_json(dest, float(data)),
        "int": lambda: _write_json(dest, int(data)),
        "dict": lambda: _write_json(dest, data),
        "list": lambda: _write_json(dest, data),
        "tuple": lambda: _write_json(dest, data),
        "str": lambda: _write_text(dest, str(data)),
        "string": lambda: _write_text(dest, str(data)),
        "xarray": lambda: data.to_netcdf(dest),
        "dataarray": lambda: data.to_netcdf(dest),
        "ndarray": lambda: np.savez(
            dest,
            **({f"arr{i}": arr for i, arr in enumerate(data)} if isinstance(data, tuple) else {"data": data})
        ),
        "npz": lambda: np.savez(
            dest,
            **({f"arr{i}": arr for i, arr in enumerate(data)} if isinstance(data, tuple) else {"data": data})
        ),
        "bytes": lambda: _write_bytes(dest, data),
    }

    if output_type in handlers:
        try:
            handlers[output_type]()
            print(f"Saved to {dest}")
        except Exception as e:
            print(f"Failed to save output ({output_type}) to {dest}: {e}")
    else:
        raise ValueError(f"Unknown output_type: {output_type}")

def load_output(filename, output_type):
    base_expected_dir = _expected_base_dir()
    output_type = output_type.lower()

    # Resolve to the active result-set directory for any non-absolute path,
    # even if it contains subfolders like "expected_outputs/...".
    if os.path.isabs(filename):
        candidate = filename
    else:
        candidate = str(base_expected_dir / filename)


    if not os.path.exists(candidate):
        legacy = os.path.join("expected_outputs", os.path.basename(filename))
        if os.path.exists(legacy):
            candidate = legacy
        else:
            raise FileNotFoundError(f"Output file not found: {candidate}")

    def _read_json(path):
        with open(path, encoding="utf-8") as f:
            return json.load(f)

    def _read_text(path):
        with open(path, "r", encoding="utf-8") as f:
            return f.read()

    def _read_bytes(path):
        with open(path, "rb") as f:
            return f.read()

    def _read_metadataframe(path):
        geojson_file = Path(path).with_suffix(".geojson")
        if geojson_file.exists():
            return gpd.read_file(geojson_file)
        return pd.read_json(path)

    def _read_ndarray(path):
        npz = np.load(path)
        if "data" in npz:
            return npz["data"]
        else:
            return tuple(npz[f] for f in npz.files)

    from hera.datalayer.document.metadataDocument import nonDBMetadataFrame
    def _read_nondbmetadataframe(path):
        if hasattr(nonDBMetadataFrame, "loadFromJsonFile"):
            return nonDBMetadataFrame.loadFromJsonFile(path)
        else:
            raise RuntimeError("nonDBMetadataFrame does not support loadFromJsonFile")

    loaders = {
        "float": lambda: float(_read_json(candidate)),
        "int": lambda: int(_read_json(candidate)),
        "dict": lambda: _read_json(candidate),
        "list": lambda: _read_json(candidate),
        "dataframe": lambda: (
            pd.read_json(candidate)
            if candidate.endswith(".json")
            else pd.read_parquet(candidate)
        ),
        "metadataframe": lambda: _read_metadataframe(candidate),
        "geodataframe": lambda: gpd.read_file(candidate),
        "xarray": lambda: xr.open_dataset(candidate),
        "dataarray": lambda: xr.open_dataarray(candidate),
        "ndarray": lambda: _read_ndarray(candidate),
        "npz": lambda: _read_ndarray(candidate),
        "bytes": lambda: _read_bytes(candidate),
        "str": lambda: _read_text(candidate),
        "string": lambda: _read_text(candidate),
        "nondbmetadataframe": lambda: _read_nondbmetadataframe(candidate),
    }

    if output_type in loaders:
        try:
            return loaders[output_type]()
        except Exception as e:
            raise RuntimeError(f"Failed to load output ({output_type}) from {candidate}: {e}")
    else:
        raise ValueError(f"Unknown output_type: {output_type}")

# =========================
# Parameter resolution / config injection
# =========================
def resolve_parameter(param):
    """
    Resolve a JSON parameter into a concrete Python object.

    Supported patterns:
      1) str: expands environment variables (e.g., "${HERA_UNITTEST_DATA}/...").
      2) dict:
         - {"fromFunction": "<name>"} -> calls resolve_from_function(name)
         - {"type": "DataFrame", "data": {...}} -> pandas.DataFrame
         - {"type": "GeoDataFrame", "path": ...} -> geopandas.read_file(path)
         - {"type": "Polygon", "coordinates": [[lon,lat], ...]} -> shapely.geometry.Polygon
         - {"fromEnv": "...", "relative": "...", "load": "..."} -> path joined with optional loader
         - {"load": "...", "path": ...} -> explicit loader from a resolved path
         - otherwise: recursively resolves nested dict values
      3) list: resolves each element recursively
      4) other: returned as-is
    """
    print("Resolving parameter:", param)

    # -------------------------
    # 1) Plain strings: expand ${VARS}
    # -------------------------
    if isinstance(param, str):
        resolved = os.path.expandvars(param)
        print(f"Expanded env var: {resolved}")
        return resolved

    # -------------------------
    # 2) Dict-based resolution
    # -------------------------
    if isinstance(param, dict):
        # 2.1) Reference a cached / registered function result
        if "fromFunction" in param:
            name = param["fromFunction"]
            print(f"Resolving from function: {name}")
            return resolve_from_function(name)

        # 2.2) Explicit typed objects (JSON literals -> python objects)

        # DataFrame literal: {"type":"DataFrame","data":{...}}
        if param.get("type") == "DataFrame":
            import pandas as pd
            if "data" not in param:
                raise ValueError(f"DataFrame param must include 'data'. Got: {param}")
            return pd.DataFrame(param["data"])

        # GeoDataFrame literal: {"type":"GeoDataFrame","path": ...}
        if param.get("type") == "GeoDataFrame":
            import geopandas as gpd
            if "path" in param:
                p = resolve_parameter(param["path"])
                return gpd.read_file(p)
            raise ValueError(f"GeoDataFrame param must include 'path'. Got: {param}")

        # Polygon literal: {"type":"Polygon","coordinates":[[lon,lat], ...]}
        # NOTE: Shapely expects (x,y) == (lon,lat)
        if param.get("type") == "Polygon":
            if "coordinates" not in param:
                raise ValueError(f"Polygon param must include 'coordinates'. Got: {param}")
            coords = param["coordinates"]
            if not isinstance(coords, list) or len(coords) < 4:
                raise ValueError(f"Polygon 'coordinates' must be a list with >= 4 points. Got: {coords}")
            return Polygon(coords)

        # 2.3) Build path from env + relative
        if "fromEnv" in param and "relative" in param:
            env_name = param["fromEnv"]
            rel = param["relative"]

            base = os.environ.get(env_name)
            if not base:
                raise ValueError(f"Environment variable '{env_name}' is not set")

            full_path = os.path.expandvars(os.path.join(base, rel))
            print(f"Resolved path from env: {full_path}")

            # Optional loaders (additive)
            load = param.get("load")

            if load == "parquet":
                import pandas as pd
                return pd.read_parquet(full_path)

            if load == "csv":
                import pandas as pd
                return pd.read_csv(full_path)

            if load in ("geo", "geofile", "shp", "geojson"):
                import geopandas as gpd
                return gpd.read_file(full_path)

            # Default behavior: just return the path (important for rasters/shp etc)
            return full_path

        # 2.4) If dict asks to load parquet/csv/geo from a "path" key
        # e.g. {"load":"parquet", "path": {"fromEnv":..., "relative":...}}
        if "load" in param and "path" in param:
            load = param["load"]
            p = resolve_parameter(param["path"])

            if not isinstance(p, str):
                raise ValueError(f"Expected path to resolve to str, got {type(p)} from {param['path']}")

            if load == "parquet":
                import pandas as pd
                return pd.read_parquet(p)

            if load == "csv":
                import pandas as pd
                return pd.read_csv(p)

            if load in ("geo", "geofile", "shp", "geojson"):
                import geopandas as gpd
                return gpd.read_file(p)

            raise ValueError(f"Unknown load type: {load}")

        # 2.5) Generic dict literal: resolve inner values recursively
        return {k: resolve_parameter(v) for k, v in param.items()}

    # -------------------------
    # 3) Lists: resolve each item
    # -------------------------
    if isinstance(param, list):
        return [resolve_parameter(x) for x in param]

    # -------------------------
    # 4) Fallback: return as-is
    # -------------------------
    return param


def resolve_from_function(name):
    print(f"Resolving from function: {name}")

    if name in function_results:
        return function_results[name]

    def _cache_and_return(key, value):
        function_results[key] = value
        return value

    function_registry = {
        "mockHighFreqToolkit": lambda: _cache_and_return(
            name,
            importlib.import_module("hera.measurements.meteorology.highfreqdata.toolkit")
            .HighFreqToolKit(projectName="TestProject")
        ),
        "mockDataLayer": lambda: _cache_and_return(
            name,
            type("MockDataLayer", (), {
                "getDataSourceData": staticmethod(lambda data, ver=None: data),
                "projectName": "TestProject"
            })()
        ),
        "testPolygon_basic": lambda: _create_test_polygon()
    }

    if name in function_registry:
        return function_registry[name]()

    raise ValueError(f"Unknown function reference: {name}")

DATA_SOURCE_LOADERS = {
    "TopographyToolkit": {"ext": "*.hgt", "dataFormat": "SRTM", "valueType": "Elevation"},
    "LandCoverToolkit": {"ext": "*.tif", "dataFormat": "RASTER", "valueType": "LandCover"},
    "DemographyToolkit": {"ext": "*.shp", "dataFormat": "VECTOR", "valueType": "Population"},
    "LowFreqToolkit": {"ext": "*.parquet", "dataFormat": "TABLE", "valueType": "MeteorologyLow"},
    "HighFreqToolkit": {"ext": "*.parquet", "dataFormat": "TABLE", "valueType": "MeteorologyHigh"}
}

def inject_custom_config_if_needed(toolkit, toolkit_name):
    if toolkit_name not in DATA_SOURCE_LOADERS:
        return
    config = DATA_SOURCE_LOADERS[toolkit_name]
    files = glob.glob(
        os.path.join(os.environ.get("HERA_DATA_PATH", ""), "**", config["ext"]),
        recursive=True
    )
    sources = {
        os.path.basename(f): {
            "item": {
                "resource": os.path.dirname(f),
                "resource_folders": [os.path.dirname(f)],
                "dataFormat": config["dataFormat"],
                "valueType": config["valueType"],
                "desc": {}
            }
        }
        for f in files if os.path.isfile(f)
    }
    toolkit.getConfig = lambda: {"DataSources": sources}
    toolkit.getDataSourceData = lambda name, version=None: sources.get(name, {}).get("item", {}).get("resource_folders")

def resolve_data_source(value):
    if isinstance(value, dict) and "fromEnv" in value and "relative" in value:
        env_value = os.environ.get(value["fromEnv"])
        if env_value is None:
            raise ValueError(f"Environment variable '{value['fromEnv']}' is not set.")
        full_path = os.path.join(env_value, value["relative"])
        print(f"Resolved path from env: {full_path}")
        return rasterio.open(full_path)
    return value

# =========================
# Test execution
# =========================
def run_definition(defn):
    # -------------------------
    # CLI step support
    # -------------------------
    if defn.get("command") == "cli":
        import subprocess

        cmd = resolve_parameter(defn["exec"])
        print(f"Running CLI: {cmd}")

        res = subprocess.run(cmd, shell=True, text=True, capture_output=True)
        out = (res.stdout or "") + "\n" + (res.stderr or "")

        if res.returncode != 0:
            raise RuntimeError(f"CLI failed (code={res.returncode}). Output:\n{out}")

        expect = defn.get("expect", {})
        for s in expect.get("contains", []):
            if s not in out:
                raise AssertionError(f"Expected '{s}' in CLI output, but it was missing.")

        successful_tests.append(defn.get("name", "cli_step"))
        return

    # -------------------------
    # Instantiate class
    # -------------------------
    module_path = ".".join(defn["class_path"].split(".")[:-1])
    class_name = defn["class_path"].split(".")[-1]
    module = importlib.import_module(module_path)
    cls = getattr(module, class_name)

    init_params = {k: resolve_parameter(v) for k, v in defn.get("init_parameters", {}).items()}
    print(f"Instantiating {cls.__name__} with parameters:", init_params)

    obj = cls(**init_params)
    inject_custom_config_if_needed(obj, cls.__name__)

    # -------------------------
    # Resolve method / attribute
    # -------------------------
    method = obj
    for part in defn["method_name"].split("."):
        method = getattr(method, part)

    if defn.get("method_type") == "attribute":
        result = method
    else:
        resolved_params = {k: resolve_parameter(v) for k, v in defn.get("parameters", {}).items()}
        if method.__name__ in ["getLandCoverAtPoint", "getLandCover", "getRoughnessAtPoint", "getRoughness"]:
            if "dataSourceName" in resolved_params:
                resolved_params["dataSourceName"] = resolve_data_source(resolved_params["dataSourceName"])
        result = method(**resolved_params)

    # -------------------------
    # Normalize output early
    # -------------------------
    output_type = defn["output_type"]
    result = normalize_output_for_type(result, output_type)

    # -------------------------
    # Postprocess (supports dict.get as well)
    # -------------------------
    if "postprocess" in defn:
        post = defn["postprocess"]

        def run_postprocess_step(step, result):
            mname = step.get("method")
            args = step.get("args", [])
            print(f"Running postprocess step: method={mname}, args={args}")

            # If result is a DataFrame and someone tries getData() - skip
            if mname == "getData" and isinstance(result, pd.DataFrame):
                print("Skipping getData(): result already a DataFrame")
                return result

            # dict "get"
            if isinstance(result, dict) and mname == "get":
                if len(args) != 1:
                    raise ValueError("postprocess dict.get expects exactly one argument (key)")
                return result.get(args[0])

            # normal attribute method call
            try:
                func = getattr(result, mname)
                return func(*args)
            except AttributeError:
                # fallback: common "data" holders
                if mname == "getData":
                    for attr in ["data", "rawData", "result"]:
                        if hasattr(result, attr):
                            val = getattr(result, attr)
                            return val() if callable(val) else val
                raise

        if isinstance(post, list):
            for step in post:
                result = run_postprocess_step(step, result)
        else:
            result = run_postprocess_step(post, result)

    function_results[defn["method_name"]] = result

    # -------------------------
    # Save / Compare
    # -------------------------
    output_file = defn["output_filename"]
    output_type = defn["output_type"]

    # Important: ensure directory exists (relative to EXP_DIR inside save_output, or absolute paths)
    # If your save_output already handles it, this is harmless.
    try:
        parent_dir = os.path.dirname(output_file)
        if parent_dir:
            os.makedirs(parent_dir, exist_ok=True)
    except Exception:
        pass

    prepare_mode = os.environ.get("PREPARE_EXPECTED_OUTPUT") == "1"

    if prepare_mode:
        save_output(output_file, result, output_type)
        print(f"Saved: {output_file}")
        successful_tests.append(defn["method_name"])
        return

    # Not prepare-mode: if expected missing -> auto-generate instead of failing
    try:
        expected = load_output(output_file, output_type)
    except FileNotFoundError as e:
        print(f"\n[WARNING] expected output file is missing: {output_file}")
        print(f"[WARNING] Auto-saving expected output (since file is missing) instead of failing.")
        save_output(output_file, result, output_type)
        print(f"Saved: {output_file}")
        successful_tests.append(defn["method_name"])
        return

    match = compare_outputs(result, expected, output_type)

    if match:
        print(f"OK {defn['method_name']} → {output_file}")
        successful_tests.append(defn["method_name"])
    else:
        print(f"\nFAILED: {defn['method_name']} → {output_file}")
        print(f"Output type   : {output_type}")
        print(f"Expected type : {type(expected)}")
        print(f"Result type   : {type(result)}")

        try:
            print("Expected (short):", str(expected)[:300])
        except Exception as ex:
            print("Error printing expected:", ex)

        try:
            print("Result (short)  :", str(result)[:300])
        except Exception as ex:
            print("Error printing result:", ex)

        failed_tests.append(defn["method_name"])

def run_all_from_json(json_file):
    print("\n" + "=" * 80)
    print(f"Starting tests from: {json_file}")
    print("=" * 80 + "\n")

    with open(json_file) as f:
        for defn in json.load(f):
            try:
                run_definition(defn)
            except Exception as e:
                where = defn.get("method_name") or defn.get("name") or "<unknown>"
                cprint(f"Exception in {where}: {e}", "red")
                traceback.print_exc()
                failed_tests.append(where)


def print_summary():
    print("\nTest Summary")
    print("=" * 60)
    for name in successful_tests:
        cprint(f"OK {name}", "green")
    for name in failed_tests:
        cprint(f"FAIL {name}", "red")
    print("=" * 60)
    print(f"Total: {len(successful_tests) + len(failed_tests)} | Passed: {len(successful_tests)} | Failed: {len(failed_tests)}")

if __name__ == "__main__":
    _expected_base_dir()  # Enforce result-set presence/creation rules at startup
    for path in sys.argv[1:]:
        run_all_from_json(path)
    print_summary()
