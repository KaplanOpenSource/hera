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

# ---------- COMPARISON FUNCTIONS ----------
def compare_pandas_dataframe(df1, df2):
    return df1.equals(df2)

def compare_numpy_array(arr1, arr2):
    return np.allclose(arr1, arr2)

def compare_xarray(xr1, xr2):
    return xr1.equals(xr2)

# ---------- SAVE / LOAD HANDLERS ----------
def save_output(filename, data, output_type):
    output_type = output_type.lower()

    if output_type in ["metadataframe", "dataframe"]:
        if hasattr(data, "getData"):
            print(f"📌 {output_type} type: {type(data)} — using getData()")
            data = data.getData()

        if isinstance(data, gpd.GeoDataFrame):
            filename = str(Path(filename).with_suffix(".geojson"))
            data.to_file(filename, driver="GeoJSON")

        elif isinstance(data, pd.DataFrame):
            data.to_json(filename, orient="records", indent=2)

        else:
            raise TypeError(f"Output type '{output_type}' expects DataFrame or GeoDataFrame, got: {type(data)}")

    elif output_type in ["float", "dict", "int"]:
        with open(filename, "w") as f:
            json.dump(data, f, indent=2, default=str)

    elif output_type in ["str", "string"]:
        with open(filename, "w") as f:
            f.write(str(data))

    elif output_type == "geodataframe":
        filename = str(Path(filename).with_suffix(".geojson"))
        data.to_file(filename, driver="GeoJSON")

    elif output_type == "xarray":
        data.to_netcdf(filename)

    elif output_type == "ndarray":
        np.savez(filename, **{f"arr{i}": arr for i, arr in enumerate(data)} if isinstance(data, tuple) else {"data": data})

    elif output_type == "bytes":
        with open(filename, "wb") as f:
            f.write(data)

    else:
        raise ValueError(f"Unknown output_type: {output_type}")

def load_output(filename, output_type):
    output_type = output_type.lower()

    if output_type in ["float", "dict", "int"]:
        with open(filename) as f:
            return json.load(f)

    elif output_type == "dataframe":
        return pd.read_json(filename)

    elif output_type == "metadataframe":
        geojson_file = Path(filename).with_suffix(".geojson")
        if geojson_file.exists():
            return gpd.read_file(geojson_file)
        else:
            return pd.read_json(filename)

    elif output_type == "geodataframe":
        return gpd.read_file(filename)

    elif output_type == "xarray":
        return xr.open_dataset(filename)

    elif output_type == "ndarray":
        return np.load(filename)["data"]

    elif output_type == "bytes":
        with open(filename, "rb") as f:
            return f.read()

    elif output_type in ["str", "string"]:
        with open(filename, "r") as f:
            return f.read()

    else:
        raise ValueError(f"Unknown output_type: {output_type}")

# ---------- PARAMETER RESOLUTION ----------
def resolve_parameter(param):
    print("🔍 Resolving parameter:", param)

    if isinstance(param, dict):
        if "fromFunction" in param:
            name = param["fromFunction"]
            print(f"⚙️ Resolving from function: {name}")

            if name == "mockHighFreqToolkit":
                if name in function_results:
                    print(f"✅ Using cached instance of {name}")
                    return function_results[name]

                from hera.measurements.meteorology.highfreqdata.toolkit import HighFreqToolKit
                instance = HighFreqToolKit(projectName="TestProject")
                function_results[name] = instance
                print(f"✅ Created new instance of HighFreqToolKit")
                return instance

            if name == "mockDataLayer":
                if name in function_results:
                    print(f"✅ Using cached instance of {name}")
                    return function_results[name]

                class MockDataLayer:
                    def getDataSourceData(self, dataSourceOrData, dataSourceVersion=None):
                        print("📦 MockDataLayer.getDataSourceData called")
                        return dataSourceOrData
                    projectName = "TestProject"

                mock = MockDataLayer()
                function_results[name] = mock
                print(f"✅ Created new instance of MockDataLayer")
                return mock

            if name == "testPolygon_basic":
                print("📐 Creating basic test polygon...")
                from shapely.geometry import Polygon
                pop_path = os.path.join(os.environ["HERA_DATA_PATH"], "measurements", "GIS", "vector", "population_lamas.shp")
                gdf = gpd.read_file(pop_path)
                geom = gdf.iloc[0].geometry
                minx, miny, maxx, maxy = geom.bounds
                buffer = 10
                polygon = Polygon([
                    (minx - buffer, miny - buffer),
                    (maxx + buffer, miny - buffer),
                    (maxx + buffer, maxy + buffer),
                    (minx - buffer, maxy + buffer),
                    (minx - buffer, miny - buffer)
                ])
                print("✅ Polygon created")
                return polygon

            return function_results.get(name)

        if param.get("type") == "DataFrame":
            print("📊 Creating DataFrame from JSON data...")
            df = pd.DataFrame(param["data"])
            print("✅ DataFrame created:", df.head())
            return df

        if "fromEnv" in param and "relative" in param:
            full_path = os.path.join(os.environ.get(param["fromEnv"], ""), param["relative"])
            print(f"📁 Resolved path from env: {full_path}")

            if full_path.endswith(".parquet"):
                return pd.read_parquet(full_path)

            if full_path.endswith(".csv"):
                return pd.read_csv(full_path)

            return full_path

    if isinstance(param, str):
        resolved = os.path.expandvars(param)
        print(f"🔄 Expanded env var: {resolved}")
        return resolved

    print("🔸 Returning raw param")
    return param

# ---------- DATA SOURCE INJECTION ----------
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

# ---------- RESOLVE RASTERIO WRAPPER ----------
def resolve_data_source(value):
    if isinstance(value, dict) and "fromEnv" in value and "relative" in value:
        env_value = os.environ.get(value["fromEnv"])
        if env_value is None:
            raise ValueError(f"Environment variable '{value['fromEnv']}' is not set.")
        full_path = os.path.join(env_value, value["relative"])
        print(f"📁 Resolved path from env: {full_path}")
        return rasterio.open(full_path)
    return value

# ---------- MAIN EXECUTION LOGIC ----------
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
                print(f"🔁 Postprocessing: calling '{step['method']}' with args: {step.get('args', [])}")
                result = getattr(result, step["method"])(*step.get("args", []))
        else:
            print(f"🔁 Postprocessing: calling '{post['method']}' with args: {post.get('args', [])}")
            result = getattr(result, post["method"])(*post.get("args", []))

    function_results[defn["method_name"]] = result

    output_file = defn["output_filename"]
    output_type = defn["output_type"]
    if os.environ.get("PREPARE_EXPECTED_OUTPUT") == "1":
        save_output(output_file, result, output_type)
        print(f"✅ Saved: {output_file}")
    else:
        expected = load_output(output_file, output_type)
        match = {
            "float": lambda: result == expected,
            "dict": lambda: result == expected,
            "dataframe": lambda: compare_pandas_dataframe(result, expected),
            "metadataframe": lambda: compare_pandas_dataframe(result, expected),
            "geodataframe": lambda: compare_pandas_dataframe(result, expected),
            "xarray": lambda: compare_xarray(result, expected),
            "ndarray": lambda: compare_numpy_array(result, expected),
            "bytes": lambda: result == expected
        }.get(output_type.lower(), lambda: False)()
        print(f"{'✅' if match else '❌'} {defn['method_name']} → {output_file}")
        if not match:
            print("Expected:", expected)
            print("Got     :", result)

def run_all_from_json(json_file):
    print("\n" + "=" * 80)
    print(f"📂 Starting tests from: {json_file}")
    print("=" * 80 + "\n")
    with open(json_file) as f:
        for defn in json.load(f):
            run_definition(defn)

if __name__ == "__main__":
    import sys
    for path in sys.argv[1:]:
        run_all_from_json(path)
