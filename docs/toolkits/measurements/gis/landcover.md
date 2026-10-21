# Land Cover

**Toolkit name:** `GIS_LandCover` (constant `toolkitHome.GIS_LANDCOVER`)

**Class:** `hera.measurements.GIS.raster.landcover.LandCoverToolkit` (internal toolkit name `LandCoverToolkit`)

Reads a land-cover classification raster (e.g. MODIS MCD12Q1, IGBP classes) and converts land-cover classes to aerodynamic roughness length (z0) for atmospheric and CFD boundary conditions. Urban cells can instead get z0 from building morphology through the [Buildings](buildings.md) toolkit.

```python
from hera import toolkitHome

# Tip: if you created the project with `hera-project project create`, you can omit projectName
lc = toolkitHome.getToolkit(toolkitHome.GIS_LANDCOVER, projectName="MY_PROJECT")

cls = lc.getLandCoverAtPoint(lon=35.0, lat=32.0)    # int class
z0 = lc.getRoughnessAtPoint(lon=35.0, lat=32.0)     # float, m
```

For the full API, see the [API Reference](../../../developer_guide/api/measurements.md).

---

## Data source format

| Property | Value |
|----------|-------|
| File format | Single-band GeoTIFF (band 1), UINT8 class values 0–16 |
| Classes (IGBP / MCD12Q1 Type 1) | 0 Water, 1–5 Forests, 6–7 Shrublands, 8–9 Savannas, 10 Grasslands, 11 Permanent wetlands, 12 Croplands, 13 Urban and built-up, 14 Cropland/natural vegetation mosaic, 15 Snow and ice, 16 Barren |
| CRS | Read from the GeoTIFF; point coordinates are reprojected to it. The grid methods assume a WGS84 raster. |
| `dataFormat` | `"string"` (the path to the `.tif`, opened with `rasterio`) is recommended. `"geotiff"` also works but needs the GDAL Python bindings. |
| `desc.type` | The classification scheme, used to pick the roughness table. Set it to `1` (see [Roughness tables](#roughness-tables)). |
| Config key | `defaultLandCover` — the datasource used when `dataSourceName` is not given. It must be set, otherwise the calls raise `KeyError`. |

Repository JSON example:

```json
{
    "GIS_LandCover": {
        "Config": { "defaultLandCover": "Type-1" },
        "DataSource": {
            "Type-1": {
                "isRelativePath": "True",
                "item": {
                    "resource": "data/landcover/lc_mcd12q1.tif",
                    "dataFormat": "string",
                    "desc": { "year": 2021, "type": 1 }
                }
            }
        }
    }
}
```

---

## Coordinate convention

`getLandCover` and `getRoughness` take a bounding box `minx, miny, maxx, maxy` and a grid spacing `dxdy` **in meters** (the grid is built in ITM and converted to WGS84).

- With `inputCRS=WSG84` (the default), `minx`/`maxx` are **latitudes** and `miny`/`maxy` are **longitudes**.
- With `inputCRS=ITM`, they are eastings and northings.

The point methods take `lon, lat` in that order.

---

## Land cover

```python
from hera.measurements.GIS import WSG84

# Class at a point (raises IndexError outside the raster)
cls = lc.getLandCoverAtPoint(lon=35.0, lat=32.0)

# Class grid: minx/maxx = latitude, miny/maxy = longitude
grid = lc.getLandCover(minx=31.9, miny=34.9, maxx=32.1, maxy=35.1, dxdy=500, inputCRS=WSG84)
grid["landcover"]    # (i, j) coordinate, alongside 'lat', 'lon' and the scalar 'dxdy'

# Class names
names = lc.getCodingMap("Type-1")    # {0: "Water", 1: "Evergreen needleleaf forest", ...}
```

| Method | Returns |
|--------|---------|
| `getLandCoverAtPoint(lon, lat, inputCRS=WSG84, dataSourceName=None)` | `int` class value |
| `getLandCover(minx, miny, maxx, maxy, dxdy=30, inputCRS=WSG84, dataSourceName=None)` | `xarray.DataArray` over `(i, j)` with `lat`, `lon`, `dxdy` and `landcover` coordinates |
| `getCodingMap(datasourceName)` | Class-to-name dict. Only the literal argument `"Type-1"` is recognised; anything else returns `{}`. |

`dataSourceName` may also be a path to a `.tif` file instead of a registered datasource name.

---

## Roughness

```python
# z0 grid from a bounding box (getLandCover + getRoughnessFromLandcover)
rough = lc.getRoughness(minx=31.9, miny=34.9, maxx=32.1, maxy=35.1, dxdy=500)
rough["z0"]

# Or add z0 to a grid you already have
rough = lc.getRoughnessFromLandcover(grid)

# Equivalent sand-grain roughness Ks = 30 * z0 (for OpenFOAM nutkRoughWallFunction)
ks = lc.roughnesslength2sandgrainroughness(0.1)
```

| Method | Description |
|--------|-------------|
| `getRoughnessAtPoint(lon, lat, inputCRS=WSG84, dataSourceName=None)` | z0 (`float`) at a point |
| `getRoughness(minx, miny, maxx, maxy, dxdy=30, inputCRS=WSG84, dataSourceName=None, isBuilding=False, windMeteorologicalDirection=None, resolution=None, GIS_BUILDINGS_dataSourceName=None)` | Land-cover grid with a `z0` coordinate |
| `getRoughnessFromLandcover(landcover, windMeteorologicalDirection=None, resolution=None, isBuilding=False, dataSourceName=None, GIS_BUILDINGS_dataSourceName=None)` | Adds `z0` to a `getLandCover` result |
| `roughnesslength2sandgrainroughness(rl)` | Static: returns `rl * 30` |

### Roughness tables

The table is chosen from the datasource `desc.type`:

- `type: 1` uses the MCD12Q1 Type-1 table from Floors et al. (2021), e.g. water 0.0001 m, forests 1 m, croplands 0.15 m, urban 0.8 m. Both the point and the grid methods use it.
- Without a `type`, `getRoughnessAtPoint` falls back to an illustrative placeholder table (0.01–1.4 m, increasing with the class number), while `getRoughnessFromLandcover`/`getRoughness` raise `AttributeError`. Always set `type: 1` for MCD12Q1 data.

### Urban areas from buildings

With `isBuilding=True`, z0 is derived from building morphology instead of the class table. The toolkit queries the [Buildings](buildings.md) toolkit of the same project over the grid area, runs `LambdaFromBuildingData(windMeteorologicalDirection, resolution, buildings)` and computes z0 and displacement height from λp, λf and hc.

```python
rough = lc.getRoughness(
    minx=32.05, miny=34.75, maxx=32.10, maxy=34.80, dxdy=100,
    isBuilding=True,
    windMeteorologicalDirection=270,   # meteorological degrees
    resolution=200,                    # m, block size for the lambda calculation
    GIS_BUILDINGS_dataSourceName="BNTL",
)
# adds 'z0', 'dd', 'lambdaF', 'lambdaP', 'hc', 'll'
```

`windMeteorologicalDirection` and `resolution` are required with `isBuilding=True` (otherwise `ValueError`), and an empty building query raises `ValueError`.

!!! note "Known issue"
    For cells with no buildings, the urban path stores the land-cover **class** (from `getLandCoverAtPoint`) in `z0` instead of a roughness (`hera/measurements/GIS/raster/landcover.py:673`).

---

## Presentation layer

`lc.presentation` draws the grid cells as coloured squares over a background image. `plot` must be a matplotlib `AxesImage` in ITM coordinates (for example `ax.get_images()[0]` after plotting a map from the [Tiles](tiles.md) toolkit). Each method opens a new figure and calls `plt.show()`; nothing is returned.

| Method | Description |
|--------|-------------|
| `plotLandcover(plot, landcover, alpha=0.2, figsize=(28, 28))` | Land-cover classes in fixed colours |
| `plotRoughness(plot, landcover, alpha=0.5, figsize=(28, 28))` | `z0` with a viridis colour bar |
| `plotLambdas(field, plot, landcover, alpha=0.2, figsize=(28, 28))` | One of the urban fields (`"lambdaF"`, `"lambdaP"`, `"hc"`, ...) |

---

## Command line

```bash
hera-GIS landcover getLandcover --minx 31.9 --miny 34.9 --maxx 32.1 --maxy 35.1 \
    --dxdy 500 --filePath landcover.csv --projectName MY_PROJECT
```

Writes a CSV with `lat`, `lon`, `landcover` and (by default) `z0`.

!!! warning "Known issue"
    In master the bounding-box options are passed to the toolkit as strings, `--roughness`/`--isBuilding` are `type=bool` (any value, even `False`, means True), and `--windDirection` hits a typo (`arguments.windDirectionis`) in `hera/measurements/GIS/CLI.py`. Prefer the Python API.

---

## Complete example

```python
from hera import toolkitHome

lc = toolkitHome.getToolkit(toolkitHome.GIS_LANDCOVER, projectName="MY_PROJECT")

# Register the raster once
lc.addDataSource("Type-1", resource="/data/landcover/lc_mcd12q1.tif", dataFormat="string",
                 version=(0, 0, 1), desc={"type": 1})
lc.setDataSourceDefaultVersion("Type-1", (0, 0, 1))
lc.setConfig(defaultLandCover="Type-1")

# Roughness for a 20 km x 20 km area at 500 m resolution (lat first for WGS84)
rough = lc.getRoughness(minx=31.9, miny=34.9, maxx=32.1, maxy=35.1, dxdy=500)

import pandas as pd
df = pd.DataFrame({f: rough[f].values.ravel() for f in ["lat", "lon", "landcover", "z0"]})
df.to_parquet("roughness.parquet")
```
