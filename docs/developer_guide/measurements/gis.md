# GIS Toolkits Implementation

Implementation details for the GIS toolkit package (`hera/measurements/GIS/`).

For user-facing documentation, see the individual toolkit pages under [User Guide > Toolkits > Measurements > GIS](../../toolkits/measurements/gis/topography_raster.md).

---

## Package structure

```
hera/measurements/GIS/
    __init__.py              # exports WSG84, ITM, convertCRS
    utils.py                 # CRS constants, coordinate conversion, STL factory
    CLI.py                   # hera-GIS CLI entry points
    raster/
        topography.py        # TopographyToolkit — SRTM elevation
        landcover.py         # LandCoverToolkit — MODIS land cover + roughness
        tiles.py             # TilesToolkit — tile server map images
    vector/
        toolkit.py           # VectorToolkit — base class for vector GIS
        topography.py        # TopographyToolkit — contour lines
        demography.py        # DemographyToolkit — census population + presentation
        buildings/
            toolkit.py       # BuildingsToolkit — footprints + 3D STL
            analysis.py      # Building analysis (lambda calculations, block iteration)
```

---

## Inheritance hierarchy

All GIS toolkits inherit from `abstractToolkit`. Vector-based toolkits share a common `VectorToolkit` base class.

| Class | Module | Inherits from | Layers |
|-------|--------|--------------|--------|
| `TopographyToolkit` (raster) | `raster/topography.py` | `abstractToolkit` | analysis |
| `LandCoverToolkit` | `raster/landcover.py` | `abstractToolkit` | analysis, presentation |
| `TilesToolkit` | `raster/tiles.py` | `abstractToolkit` | presentation |
| `VectorToolkit` | `vector/toolkit.py` | `abstractToolkit` | (base class) |
| `TopographyToolkit` (vector) | `vector/topography.py` | `VectorToolkit` | analysis |
| `BuildingsToolkit` | `vector/buildings/toolkit.py` | `VectorToolkit` | analysis |
| `DemographyToolkit` | `vector/demography.py` | `VectorToolkit` | analysis, presentation |

---

## Coordinate utilities (`GIS/utils.py`)

Key constants and functions used across all GIS toolkits:

| Item | Value | Description |
|------|-------|-------------|
| `WSG84` | `4326` | WGS84 CRS identifier (EPSG:4326) |
| `ITM` | `2039` | Israeli Transverse Mercator CRS identifier (EPSG:2039) |
| `convertCRS(points, inputCRS, outputCRS)` | — | Transform coordinates between CRS |

These are exported from `hera.measurements.GIS`:
```python
from hera.measurements.GIS import WSG84, ITM, convertCRS
```

### STL factory (`stlFactory` in `utils.py`)

The `stlFactory` class generates STL binary meshes from elevation grids. Used by:
- `TopographyToolkit.createElevationSTL()` — terrain mesh
- `BuildingsToolkit.regionToSTL()` — building 3D extrusions

Key properties:
- `heightColumnsNames` — column names for height data in GeoDataFrames

---

## Raster toolkits

### TopographyToolkit (raster)

**Source:** `raster/topography.py`

Reads SRTM HGT files (30m resolution binary elevation data). Key internals:
- Uses GDAL (`osgeo.gdal`) to open HGT files and read bands
- Data source is a folder path containing `.hgt` files
- `topographyAnalysis` inner class provides `calculateStatistics()`

### LandCoverToolkit

**Source:** `raster/landcover.py`

Reads MODIS MCD12Q1 GeoTIFF land cover classification. Key internals:
- Uses `rasterio` for raster access
- Maps IGBP class codes (0–16) to roughness lengths via lookup tables
- `presentation` inner class provides rectangle-based plotting methods
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
- `geopandasToGeoJson()` — convert GeoDataFrame to GeoJSON string
- Common patterns for data source management of vector data

### BuildingsToolkit

**Source:** `vector/buildings/toolkit.py` + `vector/buildings/analysis.py`

**Toolkit methods:**

| Method | Purpose |
|--------|---------|
| `getBuildingsFromRectangle(minx, miny, maxx, maxy, dataSourceName, inputCRS, withElevation)` | Query footprints by bounding box, optionally add raster topography elevation |
| `getBuildingHeightFromRasterTopographyToolkit(buildingData, topographyDataSource)` | Add ground elevation column from SRTM at each building centroid |
| `get_buildings_height(gdf)` | Extract height from `BLDG_HT` column; fallback to `building:levels * 3` |
| `filter_buildings_in_area(buildings_data, min_lon, min_lat, max_lon, max_lat)` | Filter GeoJSON features by bounding box (static method) |
| `buildingsGeopandasToSTLRasterTopography(buildingData, ...)` | Convert GeoDataFrame to 3D STL mesh via FreeCAD |
| `regionToSTL(minx, miny, maxx, maxy, dxdy, inputCRS, outputCRS, solidName, fileName)` | High-level STL generation for a bounding box |

**Height resolution logic** (in `getBuildingsFromRectangle` and `_LambdaP`):

1. Use `BLDG_HT` column if present and > 0
2. If `BLDG_HT` is 0 but `HI_PNT_Z` and `HT_LAND` exist, compute `HI_PNT_Z - HT_LAND`
3. If `HT_LAND` is 0, use nearest building's `HT_LAND` as fallback
4. If no height data, `get_buildings_height` uses `building:levels * 3`
5. Buildings with `FTYPE` 14 or 16 are excluded from morphology calculations

**Analysis layer (`analysis.py`):**

| Method | Purpose |
|--------|---------|
| `LambdaFromBuildingData(windDirection, resolution, buildingsData, overwrite, saveCache)` | Compute block-averaged λp, λf, hc for a domain; caches results in DB |
| `ConvexPolygons(regionNameOrData, buffer)` | Group nearby buildings into convex hulls (recursive merge of overlapping hulls) |

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

**Cache strategy:** `LambdaFromBuildingData` stores results as `type="Lambda_Buildings"` documents in the Cache collection, keyed by `(bounds, wind, resolution, crs)`. Uses a project counter (`analysis_CacheCounter`) for unique GeoJSON filenames.

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
