import os
import json
import importlib
import pandas as pd
import geopandas as gpd
import numpy as np
import xarray as xr
from pathlib import Path
import glob
import rasterio

function_results = {}

# -----------------------------------------------------------------------------------
# Compares outputs based on the given type using a dynamic registry of comparison functions
# -----------------------------------------------------------------------------------
def compare_outputs(result, expected, output_type):
    output_type = output_type.lower()

    comparison_funcs = {
        "dataframe": lambda: isinstance(result, pd.DataFrame) and result.equals(expected),
        "metadataframe": lambda: isinstance(result, pd.DataFrame) and result.equals(expected),
        "geodataframe": lambda: isinstance(result, (gpd.GeoDataFrame, pd.DataFrame)) and result.equals(expected),
        "ndarray": lambda: isinstance(result, np.ndarray) and np.allclose(result, expected),
        "xarray": lambda: isinstance(result, xr.Dataset) and result.equals(expected),
        "float": lambda: result == expected,
        "int": lambda: result == expected,
        "dict": lambda: result == expected,
        "str": lambda: result == expected,
        "string": lambda: result == expected,
        "bytes": lambda: result == expected,
    }

    compare = comparison_funcs.get(output_type)
    if compare:
        return compare()
    else:
        print(f"⚠️ No comparison function defined for type: {output_type}")
        return False

# -----------------------------------------------------------------------------------
# Saves output to disk in a format based on its type
# -----------------------------------------------------------------------------------
def save_output(filename, data, output_type):
    output_type = output_type.lower()

    if output_type in ["metadataframe", "dataframe"] and hasattr(data, "getData"):
        print(f"📌 {output_type} type: {type(data)} — using getData()")
        data = data.getData()

    def _write_json(path, obj):
        with open(path, "w") as f:
            json.dump(obj, f, indent=2, default=str)

    def _write_text(path, text):
        with open(path, "w") as f:
            f.write(text)

    def _write_bytes(path, b):
        with open(path, "wb") as f:
            f.write(b)

    handlers = {
        "dataframe": lambda: pd.DataFrame(data).to_json(filename, orient="records", indent=2),
        "metadataframe": lambda: (
            data.to_file(str(Path(filename).with_suffix(".geojson")), driver="GeoJSON")
            if isinstance(data, gpd.GeoDataFrame)
            else pd.DataFrame(data).to_json(filename, orient="records", indent=2)
        ),
        "geodataframe": lambda: data.to_file(str(Path(filename).with_suffix(".geojson")), driver="GeoJSON"),
        "float": lambda: _write_json(filename, data),
        "dict": lambda: _write_json(filename, data),
        "int": lambda: _write_json(filename, data),
        "str": lambda: _write_text(filename, str(data)),
        "string": lambda: _write_text(filename, str(data)),
        "xarray": lambda: data.to_netcdf(filename),
        "ndarray": lambda: np.savez(filename, **{f"arr{i}": arr for i, arr in enumerate(data)} if isinstance(data, tuple) else {"data": data}),
        "bytes": lambda: _write_bytes(filename, data),
    }

    if output_type in handlers:
        handlers[output_type]()
    else:
        raise ValueError(f"Unknown output_type: {output_type}")

# -----------------------------------------------------------------------------------
# Loads expected output from disk based on its type
# -----------------------------------------------------------------------------------
def load_output(filename, output_type):
    output_type = output_type.lower()

    def _read_json(path):
        with open(path) as f:
            return json.load(f)

    def _read_text(path):
        with open(path, "r") as f:
            return f.read()

    def _read_bytes(path):
        with open(path, "rb") as f:
            return f.read()

    def _read_metadataframe(path):
        geojson_file = Path(path).with_suffix(".geojson")
        return gpd.read_file(geojson_file) if geojson_file.exists() else pd.read_json(path)

    loaders = {
        "float": lambda: _read_json(filename),
        "dict": lambda: _read_json(filename),
        "int": lambda: _read_json(filename),
        "dataframe": lambda: pd.read_json(filename),
        "metadataframe": lambda: _read_metadataframe(filename),
        "geodataframe": lambda: gpd.read_file(filename),
        "xarray": lambda: xr.open_dataset(filename),
        "ndarray": lambda: np.load(filename)["data"],
        "bytes": lambda: _read_bytes(filename),
        "str": lambda: _read_text(filename),
        "string": lambda: _read_text(filename),
    }

    if output_type in loaders:
        return loaders[output_type]()
    else:
        raise ValueError(f"Unknown output_type: {output_type}")

# -----------------------------------------------------------------------------------
# Resolves dynamic parameters: function calls, env paths, DataFrame creation
# -----------------------------------------------------------------------------------
def resolve_parameter(param):
    print("🔍 Resolving parameter:", param)

    if isinstance(param, str):
        resolved = os.path.expandvars(param)
        print(f"🔄 Expanded env var: {resolved}")
        return resolved

    def resolve_from_function(name):
        print(f"⚙️ Resolving from function: {name}")

        if name in function_results:
            return function_results[name]

        function_registry = {
            "mockHighFreqToolkit": lambda: _cache_and_return(
                name,
                importlib.import_module("hera.measurements.meteorology.highfreqdata.toolkit").HighFreqToolKit(
                    projectName="TestProject"
                ),
            ),
            "mockDataLayer": lambda: _cache_and_return(
                name,
                type("MockDataLayer", (), {
                    "getDataSourceData": staticmethod(lambda data, ver=None: data),
                    "projectName": "TestProject",
                })(),
            ),
            "testPolygon_basic": lambda: _create_test_polygon()
        }

        if name in function_registry:
            return function_registry[name]()
        else:
            return function_results.get(name)

    def _cache_and_return(key, value):
        function_results[key] = value
        return value

    def _create_test_polygon():
        print("📐 Creating basic test polygon...")
        from shapely.geometry import Polygon
        pop_path = os.path.join(os.environ["HERA_DATA_PATH"], "measurements", "GIS", "vector", "population_lamas.shp")
        gdf = gpd.read_file(pop_path)
        geom = gdf.iloc[0].geometry
        minx, miny, maxx, maxy = geom.bounds
        buffer = 10
        return Polygon([
            (minx - buffer, miny - buffer),
            (maxx + buffer, miny - buffer),
            (maxx + buffer, maxy + buffer),
            (minx - buffer, maxy + buffer),
            (minx - buffer, miny - buffer)
        ])

    if isinstance(param, dict):
        if "fromFunction" in param:
            return resolve_from_function(param["fromFunction"])

        if param.get("type") == "DataFrame":
            return pd.DataFrame(param["data"])

        if "fromEnv" in param and "relative" in param:
            full_path = os.path.join(os.environ.get(param["fromEnv"], ""), param["relative"])
            print(f"📁 Resolved path from env: {full_path}")
            if full_path.endswith(".parquet"):
                return pd.read_parquet(full_path)
            if full_path.endswith(".csv"):
                return pd.read_csv(full_path)
            return full_path

    return param

# -----------------------------------------------------------------------------------
# Injects data sources into toolkit if needed based on toolkit name and file patterns
# -----------------------------------------------------------------------------------
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
    files = glob.glob(os.path.join(os.environ.get("HERA_DATA_PATH", ""), "**", config["ext"]), recursive=True)

    sources = {
        os.path.basename(f): {
            "item": {
                "resource": os.path.dirname(f),
                "resource_folders": [os.path.dirname(f)],
                "dataFormat": config["dataFormat"],
                "valueType": config["valueType"],
                "desc": {}
            }
        } for f in files if os.path.isfile(f)
    }

    toolkit.getConfig = lambda: {"DataSources": sources}
    toolkit.getDataSourceData = lambda name, version=None: sources.get(name, {}).get("item", {}).get("resource_folders")

# -----------------------------------------------------------------------------------
# Resolves raster source from relative path using rasterio
# -----------------------------------------------------------------------------------
def resolve_data_source(value):
    if isinstance(value, dict) and "fromEnv" in value and "relative" in value:
        env_value = os.environ.get(value["fromEnv"])
        if env_value is None:
            raise ValueError(f"Environment variable '{value['fromEnv']}' is not set.")
        full_path = os.path.join(env_value, value["relative"])
        print(f"📁 Resolved path from env: {full_path}")
        return rasterio.open(full_path)
    return value

# -----------------------------------------------------------------------------------
# Runs a single test case definition from a JSON description
# -----------------------------------------------------------------------------------
def run_definition(defn):
    module = importlib.import_module(".".join(defn["class_path"].split(".")[:-1]))
    cls = getattr(module, defn["class_path"].split(".")[-1])
    init_params = {k: resolve_parameter(v) for k, v in defn.get("init_parameters", {}).items()}
    print(f"⚙️ Instantiating {cls.__name__} with parameters:", init_params)

    obj = cls(**init_params)
    inject_custom_config_if_needed(obj, cls.__name__)
    method = obj
    for part in defn["method_name"].split("."):
        method = getattr(method, part)

    if defn.get("method_type") == "attribute":
        result = method
    else:
        resolved_params = {k: resolve_parameter(v) for k, v in defn["parameters"].items()}
        if method.__name__ in ["getLandCoverAtPoint", "getLandCover", "getRoughnessAtPoint", "getRoughness"]:
            if "dataSourceName" in resolved_params:
                resolved_params["dataSourceName"] = resolve_data_source(resolved_params["dataSourceName"])
        result = method(**resolved_params)

    if "postprocess" in defn:
        post = defn["postprocess"]
        if isinstance(post, list):
            for step in post:
                result = getattr(result, step["method"])(*step.get("args", []))
        else:
            result = getattr(result, post["method"])(*post.get("args", []))

    function_results[defn["method_name"]] = result

    output_file = defn["output_filename"]
    output_type = defn["output_type"]

    if os.environ.get("PREPARE_EXPECTED_OUTPUT") == "1":
        save_output(output_file, result, output_type)
        print(f"✅ Saved: {output_file}")
    else:
        expected = load_output(output_file, output_type)
        match = compare_outputs(result, expected, output_type)
        print(f"{'✅' if match else '❌'} {defn['method_name']} → {output_file}")
        if not match:
            print("Expected:", expected)
            print("Got     :", result)

# -----------------------------------------------------------------------------------
# Loads and runs all test definitions from a given JSON file
# -----------------------------------------------------------------------------------
def run_all_from_json(json_file):
    print("\n" + "=" * 80)
    print(f"📂 Starting tests from: {json_file}")
    print("=" * 80 + "\n")
    with open(json_file) as f:
        for defn in json.load(f):
            run_definition(defn)

# -----------------------------------------------------------------------------------
# CLI entry point for executing test JSON files
# -----------------------------------------------------------------------------------
if __name__ == "__main__":
    import sys
    for path in sys.argv[1:]:
        run_all_from_json(path)
