# Topography (Vector)

**Toolkit name:** `GIS_Vector_Topography` (constant `toolkitHome.GIS_VECTOR_TOPOGRAPHY`)

**Class:** `hera.measurements.GIS.vector.topography.TopographyToolkit` (internal toolkit name `Topography`; inherits `VectorToolkit`)

Works with vector topography — height contour lines — and converts it to a regular height grid, an STL surface for CFD meshing, or a DEM text file.

```python
from hera import toolkitHome

# Tip: if you created the project with `hera-project project create`, you can omit projectName
vtopo = toolkitHome.getToolkit(toolkitHome.GIS_VECTOR_TOPOGRAPHY, projectName="MY_PROJECT")

contours = vtopo.getDataSourceData("contours")     # GeoDataFrame
stl = vtopo.geoPandasToSTL(contours, dxdy=50)       # STL string
```

For the full API, see the [API Reference](../../../developer_guide/api/measurements.md).

---

## Data source format

| Property | Value |
|----------|-------|
| File format | Shapefile, GeoJSON, or GeoPackage (`dataFormat: "geopandas"`) |
| Geometry type | `LineString` (contour lines). `MultiLineString` rows fail with the pinned shapely 2.x; run `gdf.explode(index_parts=False)` first. |
| Required columns | `geometry` and `HEIGHT` (contour elevation in m). The column name is fixed by `vtopo.stlFactory.heightColumnsNames`, which you can reassign. |
| CRS | A projected CRS in meters (e.g. ITM), since `dxdy` is applied in the data's own units. Store the EPSG code in the datasource `desc` as `crs`. |

Repository JSON example:

```json
{
    "GIS_Vector_Topography": {
        "DataSource": {
            "contours": {
                "isRelativePath": "True",
                "item": {
                    "resource": "data/GIS/contours.shp",
                    "dataFormat": "geopandas",
                    "desc": { "crs": 2039 }
                }
            }
        }
    }
}
```

---

## STL and DEM generation

```python
contours = vtopo.getDataSourceData("contours")

# Optionally clip to the area of interest with plain geopandas
from shapely.geometry import box
area = contours.clip(box(176000, 660000, 180000, 664000))

# STL surface (string). The grid spacing dxdy is in the data's units (m).
stl = vtopo.geoPandasToSTL(area, dxdy=50)
with open("terrain.stl", "w") as f:
    f.write(stl)

# DEM text: "xmin xmax ymin ymax\nNx Ny\n" followed by the height rows
dem = vtopo.toDEM(area, dxdy=50)
```

| Method | Description |
|--------|-------------|
| `geoPandasToSTL(gpandas, dxdy=50, solidName="Topography")` | Interpolates the contours to a regular grid and returns the STL string. `solidName` is currently ignored; the solid is always named `Topography`. |
| `toDEM(regionNameOrData, dxdy=50, **filters)` | Returns the DEM as a string. Missing grid cells are filled with nearest-neighbour values. Pass a `GeoDataFrame` (see the warning below for the string form). |
| `stlFactory` | The underlying `hera.measurements.GIS.utils.stlFactory` (`vectorToSTL`, `rasterizeGeopandas`, `rasterToSTL`). |
| `geopandasToGeoJson(geoData)` | Static helper (from `VectorToolkit`) that converts a GeoDataFrame to a GeoJSON string. |

!!! warning "Broken in master"
    - `cutRegionFromSource(shapeDataOrName, datasourceName, isBounds=False, crs=None)` and `regionToSTL(shapeDataOrName, dxdy, datasourceName, crs=None)` always raise `TypeError`: the override calls `VectorToolkit.cutRegionFromSource` with keyword names that the base class does not accept (`hera/measurements/GIS/vector/topography.py:77`, `vector/toolkit.py:98`). Clip the GeoDataFrame yourself, as above.
    - `toDEM` with a string calls `self.getDatasourceData`, which does not exist, so only the `GeoDataFrame` form works.

---

## Analysis layer

`vtopo.analysis.addHeight(data, groundData, coord1="x", coord2="y", coord3="z", resolution=10, saveMode=TOOLKIT_SAVEMODE_ONLYFILE, file=None, fillna=0, **kwargs)` adds `ground` and `height` (height above ground, clipped at 0) columns to a mesh DataFrame `data`, by interpolating the ground points in `groundData` onto a regular grid of cell size `resolution`.

- `saveMode` takes the `TOOLKIT_SAVEMODE_*` constants from `hera.toolkit`. With the default, the result is written to `<filesDirectory>/cellData.parquet` (or `file`), and `FileExistsError` is raised if the file already exists. Use `TOOLKIT_SAVEMODE_NOSAVE` to only return the DataFrame.
- With the `FILEANDDB` modes the file is also registered as a Cache document of type `cellData`.

!!! warning "Broken in master"
    `addHeight` raises `AttributeError` before doing any work: it calls `coordinateHandler.regularizeTimeSteps` on the `hera.simulations.utils.coordinateHandler` *module* instead of the class, and it reads `toolkit.TOOLKIT_SAVEMODE_*` from `hera.measurements.GIS.vector.toolkit`, which does not define them (`hera/measurements/GIS/vector/topography.py:262, 282`).

There is no presentation layer for this toolkit.

---

## Command line

```bash
hera-GIS topography_vector list      # prints the datasource table of the default project
```
