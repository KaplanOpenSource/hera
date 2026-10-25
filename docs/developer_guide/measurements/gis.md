# GIS Toolkits Implementation

Implementation details for the GIS toolkit package (`hera/measurements/GIS/`).

For user-facing documentation, see the individual toolkit pages under [User Guide > Toolkits > Measurements > GIS](../../toolkits/measurements/gis/topography_raster.md).

---

## Package structure

```
hera/measurements/GIS/
    __init__.py              # exports WSG84, ITM, convertCRS
    utils.py                 # CRS constants, coordinate conversion, create_xarray, STL factory
    CLI.py                   # hera-GIS CLI entry points
    raster/
        topography.py        # TopographyToolkit — SRTM elevation
        landcover.py         # LandCoverToolkit — MODIS land cover + roughness
        tiles.py             # TilesToolkit — tile server map images
    vector/
        toolkit.py           # VectorToolkit — base class for vector GIS
        topography.py        # TopographyToolkit — contour lines
        demography.py        # DemographyToolkit — census population + presentation
        _io_utils.py         # readGeoJSONString — engine-agnostic GeoJSON parsing
        buildings/
            toolkit.py       # BuildingsToolkit — footprints + 3D STL
            analysis.py      # Building analysis (lambda calculations, block iteration)
            presentation.py  # Building plots (footprints, heights, morphology, map overlay)
```

---

## Inheritance hierarchy

All GIS toolkits inherit from `abstractToolkit`. Vector-based toolkits share a common `VectorToolkit` base class.

| Class | Module | Inherits from | Layers |
|-------|--------|--------------|--------|
| `TopographyToolkit` (raster) | `raster/topography.py` | `abstractToolkit` | analysis |
| `LandCoverToolkit` | `raster/landcover.py` | `abstractToolkit` | presentation |
| `TilesToolkit` | `raster/tiles.py` | `abstractToolkit` | presentation |
| `VectorToolkit` | `vector/toolkit.py` | `abstractToolkit` | (base class) |
| `TopographyToolkit` (vector) | `vector/topography.py` | `VectorToolkit` | analysis |
| `BuildingsToolkit` | `vector/buildings/toolkit.py` | `VectorToolkit` | analysis, presentation |
| `DemographyToolkit` | `vector/demography.py` | `VectorToolkit` | analysis, presentation |

---

## Coordinate utilities (`GIS/utils.py`)

Key constants and functions used across all GIS toolkits:

| Item | Value | Description |
|------|-------|-------------|
| `WSG84` | `4326` | WGS84 CRS identifier (EPSG:4326); legacy misspelling kept for backward compatibility |
| `WGS84` | `4326` | Correctly spelled alias of `WSG84` |
| `ITM` | `2039` | Israeli Transverse Mercator CRS identifier (EPSG:2039) |
| `ED50_ZONE36N` | `23036` | ED50 / UTM zone 36N |
| `convertCRS(points, inputCRS, outputCRS)` | — | Transform coordinates between CRS |

These are exported from `hera.measurements.GIS`:
```python
from hera.measurements.GIS import WSG84, ITM, convertCRS
```

### STL factory (`stlFactory` in `utils.py`)

The `stlFactory` class generates STL binary meshes from elevation grids. Used by:
- `TopographyToolkit.createElevationSTL()` / `getElevationSTL()` (raster) — terrain mesh via `rasterToSTL`
- `TopographyToolkit.geoPandasToSTL()` / `regionToSTL()` (vector) — contour-based terrain mesh via `vectorToSTL`

Building STL export does not use `stlFactory`; it goes through FreeCAD in `BuildingsToolkit.buildingsGeopandasToSTLRasterTopography()`.

Key properties:
- `heightColumnsNames` — column names for height data in GeoDataFrames

---

## Raster toolkits

### TopographyToolkit (raster)

**Source:** `raster/topography.py`

Reads SRTM HGT files (30m resolution binary elevation data). Key internals:
- Uses `rasterio` (`_open_raster`) to read the HGT band and geotransform; GDAL is only imported optionally
- Data source is a folder path (or list of folders) containing `.hgt` files, resolved by `findElevationFile()`
- `topographyAnalysis` class provides `calculateStastics(elevation)` (sic)

### LandCoverToolkit

**Source:** `raster/landcover.py`

Reads MODIS MCD12Q1 GeoTIFF land cover classification. Key internals:
- Uses `rasterio` for raster access
- Maps IGBP class codes (0–16) to roughness lengths via lookup tables
- `presentation` class provides rectangle-based plotting methods (`plotRoughness`, `plotLandcover`, `plotLambdas`)
- `_getRoughnessFromBuildingsDataFrame()` integrates building data for urban roughness

### TilesToolkit

**Source:** `raster/tiles.py`

Fetches and assembles map tile images from XYZ tile servers. Key internals:
- `deg2tile` / `tile2deg` — coordinate ↔ tile index conversion
- Downloads individual 256×256 PNG tiles via HTTP
- Assembles into a single image with extent metadata
- `presentation.plot()` renders with CRS conversion (WGS84 ↔ ITM)

---

## Vector toolkits

### VectorToolkit (base class)

**Source:** `vector/toolkit.py`

Base class for all vector GIS toolkits. Adds:
- `geopandasToGeoJson()` (static) — convert GeoDataFrame to GeoJSON string
- `cutRegionFromSource(datasourceDocument, shape, isBounds=False, inputCRS=WGS84)` — clip a vector datasource to a shape or its bounding box
- Common patterns for data source management of vector data

### BuildingsToolkit

**Source:** `vector/buildings/toolkit.py` + `vector/buildings/analysis.py` + `vector/buildings/presentation.py`

**Toolkit methods:**

| Method | Purpose |
|--------|---------|
| `getBuildingsFromRectangle(minx, miny, maxx, maxy, dataSourceName, inputCRS, withElevation)` | Query footprints by bounding box, optionally add raster topography elevation |
| `getBuildingHeightFromRasterTopographyToolkit(buildingData, topographyDataSource)` | Add ground elevation column from SRTM at each building centroid |
| `get_buildings_height(gdf)` | Static. Build a `name`/`geometry`/`height` GeoDataFrame from OSM-style data: `height` column, fallback `building:levels * 3` |
| `filter_buildings_in_area(buildings_data, min_longitude, min_latitude, max_longitude, max_latitude)` | Filter GeoJSON features by bounding box (static method) |
| `buildingsGeopandasToSTLRasterTopography(buildingData, buildingHeightColumn, buildingElevationColumn, outputFileName, flatTerrain=False, referenceTopography=0, nonFlatTopographyShift=10)` | Convert GeoDataFrame to 3D STL mesh via FreeCAD |

There is no high-level `regionToSTL` on `BuildingsToolkit`; the `hera-GIS buildings toSTL` CLI chains `getBuildingsFromRectangle(..., withElevation=True)` and `buildingsGeopandasToSTLRasterTopography(buildings, "BLDG_HT", "elevation", fileName)`.

**Height resolution logic** (in `Blocks._LambdaP` / `Blocks._LambdaF`):

1. Use the `BLDG_HT` column
2. If `HT_LAND` is 0 but `BLDG_HT` > 0, take `HT_LAND` of the nearest building with valid ground height and use `max(BLDG_HT - HT_LAND_nearest, 0)`
3. In `_LambdaF`, if `BLDG_HT` < 2 m, use `HI_PNT_Z - HT_LAND`; if that is still < 2 m the building is skipped
4. Buildings with `FTYPE` 14 or 16, or with zero height, are excluded from morphology calculations

**Analysis layer (`analysis.py`):**

| Method | Purpose |
|--------|---------|
| `LambdaFromBuildingData(windMeteorologicalDirection, resolution, buildingsData, overwrite=False, saveCache=True)` | Compute block-averaged λp, λf, hc for a domain; caches results in DB |
| `ConvexPolygons(regionNameOrData, buffer=100)` | Group nearby buildings into convex hulls (recursive merge of overlapping hulls) |

**Blocks class** — domain decomposition engine used by `LambdaFromBuildingData`:

The `Blocks` class divides the domain bounding box into a grid of square cells and computes morphology parameters per cell:

| Method | Purpose |
|--------|---------|
| `__init__(level, df, exteriorBlock, **kwargs)` | Create a block level with division by `size`, `npxy`, `width`/`height`, or `npx`/`npy` |
| `iterBlocks(**kwargs)` | Create a nested block level one level deeper |
| `Lambda(buildings, windDirection)` | Iterate all blocks, compute λp, λf, hc per block, return GeoDataFrame |
| `initBuildingsBlock(blockDict)` | Clip buildings to a single block extent |
| `_LambdaP()` | Plan area fraction: total building footprint / block area |
| `_LambdaF(windDirection)` | Frontal area density: rotate each building by wind direction, project width × height |
| `_A_f(buildingGeometry, height, windDirection)` | Frontal area of a single building |
| `getHc()` | Area-weighted mean building height |

**Cache strategy:** `LambdaFromBuildingData` stores results as `type="Lambda_Buildings"` documents in the Cache collection, keyed by `(bounds, wind, resolution, crs)`. Uses the project config key `analysis_CacheCounter` to name the GeoJSON files (`LambdaGeoJson_<n>.geojson` in the project `filesDirectory`).

**Presentation layer (`presentation.py`):** `plotBuildings`, `plotBuildingHeights`, `plotMorphology`, `plotBuildingsWithMorphology`, `plotHeightHistogram`, `plotBuildingsOnMap` (overlay on a `TilesToolkit` image).

### DemographyToolkit

**Source:** `vector/demography.py`

Key internals:
- Population type mapping: `{"All": "total_pop", "Children": "age_0_14", ...}`
- `analysis.calculatePopulationInPolygon()` — area-weighted spatial intersection
- `analysis.createNewArea()` — aggregate population within a custom polygon
- `presentation` layer with 6 plotting methods (density, counts, age groups, polygon intersection, area annotation, map overlay)

---

## Cross-references

| What | Where |
|------|-------|
| User guide (per-toolkit usage) | [Measurements > GIS](../../toolkits/measurements/gis/topography_raster.md) |
| API reference (auto-generated) | [API > Measurements](../api/measurements.md) |
| Data source formats | [Measurements > GIS > each toolkit](../../toolkits/measurements/gis/topography_raster.md) (Data source format tables) |
| CLI commands | [CLI Reference > hera-GIS](../../cli/reference.md#hera-gis) |
