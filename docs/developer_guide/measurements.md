# Measurements Implementation

This page covers the internal architecture of the measurement toolkits for developers extending or maintaining them.

---

## Package structure

```
hera/measurements/
    GIS/
        raster/
            topography.py      # TopographyToolkit — SRTM elevation data
            landcover.py       # LandCoverToolkit — MODIS land cover
            tiles.py           # TilesToolkit — tile server map images
            hill2stl.py        # STL mesh generation from elevation
        vector/
            toolkit.py         # VectorToolkit — base class for vector GIS
            topography.py      # TopographyToolkit (vector contours)
            buildings/
                toolkit.py     # BuildingsToolkit — footprints + 3D STL
                analysis.py    # Building analysis layer
            demography.py      # DemographyToolkit — population data
            abstractLocation.py
        shapes.py              # GIS shapes toolkit
        utils.py               # CRS conversion, coordinate utilities
        CLI.py                 # hera-GIS CLI entry points
    meteorology/
        lowfreqdata/
            toolkit.py         # lowFreqToolKit — hourly/daily station data
            analysis.py        # Statistical analysis methods
            presentationLayer.py  # Daily and seasonal plots
        highfreqdata/
            toolkit.py         # HighFreqToolKit — sonic anemometer data
            analysis/
                analysislayer.py       # Analysis dispatcher
                abstractcalculator.py  # Base calculator class
                meandatacalculator.py  # Mean statistics
                turbulencestatistics.py # Turbulence calculations
            parsers/
                CampbellBinary.py  # Campbell Scientific binary parser
                TOA5.py            # TOA5 CSV format parser
        radiosonde.py          # Radiosonde atmospheric data
        analysis.py            # Shared meteorology analysis
        GFS.py                 # GFS weather model data
    experiment/
        experiment.py          # experimentHome — experiment toolkit
        analysis.py            # Experiment analysis layer
        dataEngine.py          # Parquet data engine for experiments
        presentation.py        # Experiment visualizations
        parsers.py             # Data file parsers
        CLI.py                 # hera-experiment CLI entry points
```

---

## GIS toolkits

### Inheritance hierarchy

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

### Coordinate utilities (`GIS/utils.py`)

Key constants and functions used across GIS toolkits:

| Item | Description |
|------|-------------|
| `WSG84` | WGS84 CRS identifier (EPSG:4326) |
| `ITM` | Israeli Transverse Mercator CRS identifier (EPSG:2039) |
| `convertCRS(points, inputCRS, outputCRS)` | Transform coordinates between CRS |

### STL generation (`GIS/raster/hill2stl.py`)

The `hill2stl` module converts elevation grids to STL meshes for CFD simulations. Used by both the raster topography and buildings toolkits.

---

## Meteorology toolkits

### lowFreqToolKit

The toolkit follows the standard analysis + presentation pattern:

| Layer | Class | Key methods |
|-------|-------|------------|
| Toolkit | `lowFreqToolKit` | `getDataSourceData`, standard datasource API |
| Analysis | `lowFreqAnalysis` | `addDatesColumns`, `calcHourlyDist`, `resampleSecondMoments`, `matchDataWithOther` |
| Presentation | `lowFreqPresentation` | `dailyPlots.plotScatter`, `dailyPlots.plotDaily`, `seasonalPlots.plotProbContourf_bySeason`, `seasonalPlots.plotSeasonalHourly` |

### HighFreqToolKit

High-frequency data uses a calculator-based analysis architecture:

| Layer | Class | Purpose |
|-------|-------|---------|
| Toolkit | `HighFreqToolKit` | Data source management, station metadata |
| Analysis | `RawdataAnalysis` | Dispatches to specialized calculators |
| Calculator | `AbstractCalculator` | Base class for analysis calculators |
| Calculator | `MeanDataCalculator` | Mean statistics (wind speed, direction, temperature) |
| Calculator | `TurbulenceStatistics` | Second-order moments, eddy covariance |
| Parsers | `CampbellBinary`, `TOA5` | Raw data format parsers |

---

## Experiment toolkit

The experiment system organizes raw data files into structured experiments:

| Component | Class | Purpose |
|-----------|-------|---------|
| Toolkit | `experimentHome` | Create/list/get experiments |
| Analysis | `experimentAnalysis` | Experiment-specific analysis |
| Data engine | `parquetDataEngineHera` | Parquet-based data storage |
| Presentation | Experiment presentation | Visualization layer |

Experiments are stored as measurement documents and can be exposed as toolkits through `ToolkitHome.getExperimentToolkitDocuments()`.

---

## Adding a new measurement toolkit

1. Create a new module under `hera/measurements/<domain>/`
2. Create a class inheriting from `abstractToolkit` (or `VectorToolkit` for vector GIS)
3. Set `toolkitName` in `__init__`
4. Optionally add `_analysis` and `_presentation` layers
5. Register in `ToolkitHome._toolkits` dict in `hera/toolkit.py`

For the full API, see the [Measurements API Reference](api/measurements.md).
