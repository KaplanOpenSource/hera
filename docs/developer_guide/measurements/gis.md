# GIS toolkits

## Inheritance hierarchy

All GIS toolkits inherit from `abstractToolkit`. Vector-based toolkits share a common `VectorToolkit` base class.

| Class | Inherits from | Purpose |
|-------|--------------|---------|
| `TopographyToolkit` (raster) | `abstractToolkit` | SRTM HGT file access, elevation grids, STL generation |
| `LandCoverToolkit` | `abstractToolkit` | MODIS MCD12Q1 land cover classification, roughness |
| `TilesToolkit` | `abstractToolkit` | Tile server images, WGS84/ITM coordinate support |
| `VectorToolkit` | `abstractToolkit` | Base class for vector GIS toolkits |
| `TopographyToolkit` (vector) | `VectorToolkit` | Contour lines, vector DEM |
| `BuildingsToolkit` | `VectorToolkit` | Building footprints, 3D STL with ground elevation |
| `DemographyToolkit` | `VectorToolkit` | Census population data, polygon intersections |

## Coordinate utilities (`GIS/utils.py`)

Key constants and functions used across GIS toolkits:

| Item | Description |
|------|-------------|
| `WSG84` | WGS84 CRS identifier (EPSG:4326) |
| `ITM` | Israeli Transverse Mercator CRS identifier (EPSG:2039) |
| `convertCRS(points, inputCRS, outputCRS)` | Transform coordinates between CRS |

## STL generation (`GIS/raster/hill2stl.py`)

The `hill2stl` module converts elevation grids to STL meshes for CFD simulations. Used by both the raster topography and buildings toolkits.
