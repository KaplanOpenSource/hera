# Topography (Raster)

**Toolkit name:** `GIS_Raster_Topography` (constant `toolkitHome.GIS_RASTER_TOPOGRAPHY`)

**Class:** `hera.measurements.GIS.raster.topography.TopographyToolkit` (internal toolkit name `TopographyToolkit`)

Reads SRTM elevation tiles (`.hgt` files) and provides terrain height for single points, point lists, regular grids, and STL surfaces for CFD meshing.

```python
from hera import toolkitHome

# Tip: if you created the project with `hera-project project create`, you can omit projectName
topo = toolkitHome.getToolkit(toolkitHome.GIS_RASTER_TOPOGRAPHY, projectName="MY_PROJECT")

# Elevation (m above sea level) at a single WGS84 point
elev = topo.getPointElevation(lat=33.459, long=35.755)
```

For the full API, see the [API Reference](../../../developer_guide/api/measurements.md).

---

## Data source format

| Property | Value |
|----------|-------|
| File format | SRTM `.hgt` tiles (1 arc-second SRTMGL1 or 3 arc-second SRTM3), read with `rasterio` |
| File naming | `N{lat}E{lon:03d}.hgt`, e.g. `N33E035.hgt`. The toolkit builds this name from the integer part of the coordinates, so only northern/eastern-hemisphere tiles are found. |
| CRS | WGS84 (EPSG:4326) |
| Datasource resource | A **directory** that contains the `.hgt` files (`dataFormat: "string"`). A list of directories is also accepted; they are searched in order. |
| Config key | `defaultSRTM` — the datasource used when `dataSourceName` is not given. It must be set, otherwise the calls raise `KeyError`. |

Repository JSON example:

```json
{
    "GIS_Raster_Topography": {
        "Config": { "defaultSRTM": "SRTMGL1" },
        "DataSource": {
            "SRTMGL1": {
                "isRelativePath": "True",
                "item": {
                    "resource": "data/raster/SRTM",
                    "dataFormat": "string",
                    "desc": {}
                }
            }
        }
    }
}
```

Or register it directly from Python:

```python
topo.addDataSource("SRTMGL1", resource="/data/SRTM", dataFormat="string", version=(0, 0, 1))
topo.setDataSourceDefaultVersion("SRTMGL1", (0, 0, 1))
topo.setConfig(defaultSRTM="SRTMGL1")
```

---

## Coordinate convention

All grid methods take a bounding box `minx, miny, maxx, maxy` and a spacing `dxdy` **in meters**. The grid itself is built in ITM (EPSG:2039) and converted back to WGS84.

- With `inputCRS=WSG84` (the default), `minx`/`maxx` are **latitudes** and `miny`/`maxy` are **longitudes**.
- With `inputCRS=ITM`, `minx`/`maxx` are eastings and `miny`/`maxy` are northings.

Grids larger than 1,000,000 points raise `MemoryError`; increase `dxdy` or shrink the area.

---

## Point elevations

```python
# Single point (WGS84); bilinear interpolation inside the tile
elev = topo.getPointElevation(lat=33.85, long=35.15)

# Explicit datasource
elev = topo.getPointElevation(lat=33.85, long=35.15, dataSourceName="SRTMGL1")
```

```python
import pandas as pd

# A DataFrame with 'lat' and 'lon' columns (a GeoSeries of points, or a
# GeoDataFrame with a 'point' column, is also accepted)
points = pd.DataFrame({"lat": [33.85, 33.90], "lon": [35.15, 35.20]})
result = topo.getPointListElevation(points)
# result has the added columns 'filename' and 'elevation'
```

| Method | Returns |
|--------|---------|
| `getPointElevation(lat, long, dataSourceName=None)` | `float` — elevation in m. Note the argument is `long`, not `lon`. |
| `getPointListElevation(pointList, dataSourceName=None, inputCRS=WSG84)` | The input as a DataFrame with `filename` and `elevation` columns. Points must be WGS84. |

---

## Gridded elevation

```python
from hera.measurements.GIS import WSG84

# minx/maxx = latitude, miny/maxy = longitude (WGS84 convention)
ds = topo.getElevation(minx=33.85, miny=35.10, maxx=33.87, maxy=35.12, dxdy=30, inputCRS=WSG84)

ds["lat"]          # data variable over (i, j)
ds["lon"]          # data variable over (i, j)
ds["elevation"]    # coordinate over (i, j)
```

To add elevation to an existing grid, pass an `xarray.Dataset` that has `lat` and `lon` variables over `(i, j)`:

```python
grid = topo.create_xarray(minx=33.85, miny=35.10, maxx=33.87, maxy=35.12, dxdy=30)
grid = topo.getElevationOfXarray(grid)
```

| Method | Returns |
|--------|---------|
| `create_xarray(minx, miny, maxx, maxy, dxdy=30, inputCRS=WSG84)` | `xarray.Dataset` with `lat`/`lon` over `(i, j)` |
| `getElevation(minx, miny, maxx, maxy, dxdy=30, inputCRS=WSG84, dataSourceName=None)` | The same grid with an `elevation` coordinate |
| `getElevationOfXarray(xarray_dataset, dataSourceName=None)` | The input dataset with an `elevation` coordinate |

---

## STL generation

Both methods return the STL as a **string**; write it to a file yourself. The STL vertices are in ITM meters.

```python
# From a bounding box (same convention as getElevation)
stl = topo.createElevationSTL(minx=33.85, miny=35.10, maxx=33.87, maxy=35.12,
                              dxdy=30, solidName="Terrain")

# From a grid you already have (needs 'lat', 'lon' and 'elevation')
stl = topo.getElevationSTL(ds, solidName="Terrain")

with open("terrain.stl", "w") as f:
    f.write(stl)
```

| Method | Description |
|--------|-------------|
| `createElevationSTL(minx, miny, maxx, maxy, dxdy=30, shiftx=0, shifty=0, inputCRS=WSG84, dataSourceName=None, solidName="Topography")` | `getElevation` + `getElevationSTL` |
| `getElevationSTL(elevation, shiftx=0, shifty=0, solidName="Topography")` | STL string from an elevation grid. `shiftx`/`shifty` (m) are subtracted from the ITM coordinates to move the origin. |

The same operation is available from the command line:

```bash
hera-GIS topography_raster list
hera-GIS topography_raster toSTL --minx 33.85 --miny 35.10 --maxx 33.87 --maxy 35.12 \
    --dxdy 30 --fileName terrain.stl --projectName MY_PROJECT
```

!!! warning "Known issue"
    In master the `toSTL` CLI passes `--minx/--miny/--maxx/--maxy` (and `--dxdy`) to the toolkit as strings (`hera/bin/hera-GIS`), so the call fails. Use the Python API until this is fixed. The `list` subcommand takes no `--projectName` option and uses the default project.

---

## Coordinate conversion

```python
from hera.measurements.GIS import WSG84, ITM

# list of (x, y) tuples, an (N, 2) numpy array, or a DataFrame with 'x'/'y' columns
gdf = topo.convertPointsCRS([(35.2, 32.5)], inputCRS=WSG84, outputCRS=ITM)
gdf.geometry.iloc[0].x, gdf.geometry.iloc[0].y
```

`convertPointsCRS(points, inputCRS, outputCRS, **kwargs)` returns a `GeoDataFrame` with a `geometry` column. For WGS84 input, `x` is the longitude. For a DataFrame, pass `x=`/`y=` to name other columns.

---

## Analysis layer

`topo.analysis.calculateStastics(elevation)` (the method name is spelled this way in the code) returns a dict with `xmin`, `xmax`, `ymin`, `ymax`, `size`, `mean`, `std`, `domainmax`, `domainmaxlocation`, `domainmin`, `domainminlocation`.

!!! note
    It expects a dataset with `X`, `Y` and `Elevation` variables, not the `lat`/`lon`/`elevation` grid that `getElevation` returns. Rename the fields (and use metric coordinates) before calling it.

There is no presentation layer for this toolkit.
