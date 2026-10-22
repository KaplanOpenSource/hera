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
        vector/
            toolkit.py         # VectorToolkit — base class for vector GIS
            topography.py      # TopographyToolkit (vector contours)
            buildings/
                toolkit.py     # BuildingsToolkit — footprints + 3D STL
                analysis.py    # Building analysis layer
                presentation.py  # Building presentation layer
            demography.py      # DemographyToolkit — population data
            _io_utils.py       # Engine-agnostic GeoJSON -> GeoDataFrame helpers
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

## Toolkit sub-pages

- [GIS toolkits](gis.md) — Inheritance hierarchy, coordinate utilities, STL generation
- [Meteorology toolkits](meteorology.md) — lowFreqToolKit, HighFreqToolKit
- [Experiment toolkit](experiment.md) — Experiment system architecture

---

## Adding a new measurement toolkit

1. Create a new module under `hera/measurements/<domain>/`
2. Create a class inheriting from `abstractToolkit` (or `VectorToolkit` for vector GIS)
3. Pass `toolkitName` to `super().__init__(...)`
4. Optionally add `_analysis` and `_presentation` layers
5. Register it:
    - **External / project toolkits** — register dynamically with `toolkitHome.registerToolkit(toolkit_name=..., toolkit_path=...)` (stored as a `ToolkitDataSource`; `getToolkit` then imports `<toolkit_name>.<toolkit_name>` from that path), or in bulk with `toolkitHome.import_toolkits_from_json(...)`.
    - **Built-in toolkits** — the static registry is the `self._toolkits` dict built in `ToolkitHome.__init__` (`hera/toolkit.py`), together with a name constant on `ToolkitHome` (e.g. `METEOROLOGY_LOWFREQ = "MeteoLowFreq"`). Prefer dynamic registration over editing this dict.

For the full API, see the [Measurements API Reference](../api/measurements.md).
