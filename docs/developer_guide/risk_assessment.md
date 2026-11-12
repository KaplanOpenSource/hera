# Risk Assessment Implementation

This page covers the internal architecture of the risk assessment toolkit for developers extending or maintaining it.

---

## Package structure

```
hera/riskassessment/
    riskToolkit.py             # RiskToolkit + analysis layer class
    agents/
        Agents.py              # Agent + PhysicalPropeties classes
        effects/
            Injury.py          # InjuryFactory + Injury base + concrete injuries
            InjuryLevel.py     # InjuryLevel base + Lognormal10DoseResponse, Threshold, Exponential
            Calculator.py      # Dose calculators (Haber, TenBerge, MaxConcentration)
            thresholdGeoDataFrame.py  # Spatial threshold results
    protectionpolicy/
        ProtectionPolicy.py    # ProtectionPolicy + action pipeline (Indoor, Masks)
    analysis/
        riskAreas.py           # Risk area algorithms (Sweep)
    presentation/
        casualtiesFigs.py      # Casualty plots (roses, bar charts)
    CLI.py                     # hera-riskassessment CLI (agents createRepository)
    docs/                      # Internal documentation and helper scripts
```

---

## Core classes

### RiskToolkit

The main entry point (`toolkitHome.RISKASSESSMENT`, `toolkitName="RiskAssessment"`). Inherits from `abstractToolkit` and manages agents as versioned data sources.

| Component | Type | Purpose |
|-----------|------|---------|
| `RiskToolkit` | Toolkit | `getAgent(nameOrDesc, version=None)`, `listAgentsNames()`, `loadAgent(name, agentDescription, version, saveMode=...)`, `loadData(fileNameOrData, saveMode=..., **kwargs)` |
| `analysis` | `analysis` class in `riskToolkit.py` | `getRiskAreas(tenbergeCoefficient, levels, Q=1*ureg.kg, protectionPolicy=None, LSMfile=None, **LSMparams)`; holds an `LSMToolkit` (`analysis.LSM`) |
| `presentation` | `casualtiesPlot` | `plotCasualtiesRose(...)`, `plotCasualtiesProjection(...)` |
| `ProtectionPolicy` | Class reference | Protection action pipeline builder (the class, not an instance) |

`getAgent` accepts either a datasource name (the descriptor is loaded with `getDataSourceData`; missing `name`/`version` are filled from the datasource document) or a descriptor `dict`, and returns an `Agent`.

### Agent (`agents/Agents.py`)

An `Agent` represents a hazardous material with its effects and physical properties.

| Component | Class | Purpose |
|-----------|-------|---------|
| `Agent` | Container | Holds effects, physical properties, Ten Berge coefficient |
| `PhysicalPropeties` | Properties | Molecular weight, density, volatility, vapor pressure, sorption |
| Effects (dict) | `Injury` instances | Accessed by name: `agent["inhalation"]` |

#### Agent initialization flow

When `Agent(descriptor)` is called:

1. Extract `effectParameters` from the descriptor (e.g., `tenbergeCoefficient`)
2. For each entry in `effects`:
   - Call `InjuryFactory.getInjury(name, cfgJSON, **effectParameters)` to build the `Injury` object
   - Store in `self._effects[name]`
3. Add all effects as attributes on the Agent (`self.__dict__.update(self._effects)`) — this allows `agent.RegularPopulation` as well as `agent["RegularPopulation"]`
4. Initialize `PhysicalPropeties` from the `physicalProperties` section

#### Ten Berge coefficient mutation

When `agent.tenbergeCoefficient = value` is set:
- Updates `self._effectParameters["tenbergeCoefficient"]`
- **Rebuilds all effects** by re-calling `InjuryFactory.getInjury()` for each effect with the new coefficient
- This means changing the coefficient propagates to all calculators immediately

#### PhysicalPropeties class

Stores and computes physical properties using empirical correlations:

| Property/Method | What it computes |
|----------------|-----------------|
| `molecularWeight` | Molecular weight in g/mol (auto-converted from string with units) |
| `sorptionCoefficient` | Sorption coefficient in cm/s |
| `spreadFactor` | Dimensionless spread factor |
| `molecularVolume` | Molecular volume from descriptor |
| `getVolatility(temperature)` | Vapor saturation concentration [g/cm³] using Antoine-type equation |
| `getDensity(temperature)` | Liquid density [g/cm³] using linear temperature model |
| `vaporPressure(temperature)` | Vapor pressure [bar] using extended Antoine equation |

Temperature inputs must carry units (`pint.Quantity`, or `Unum` when installed); the methods convert to Celsius (volatility/density) or Kelvin (vapor pressure) via `unumToPint`. A raw float raises `AttributeError` (`unumToPint` returns it unchanged and `.m_as` fails).

#### Agent descriptor JSON format

```json
{
    "name": "Chlorine",
    "effectParameters": {
        "tenbergeCoefficient": 2.0
    },
    "physicalProperties": {
        "molecularWeight": "70.9*g/mol",
        "sorptionCoefficient": "0.5*cm/s",
        "spreadFactor": 1.0,
        "volatilityConstants": [7.58, 1500, 230, 0],
        "densityConstants": [1.56, 0.002, 20],
        "molecularVolume": 45.5,
        "vaporPressure": {"A": 6.93, "B": 861, "C": -27.01, "units": "mmHg"}
    },
    "effects": {
        "RegularPopulation": {
            "type": "Lognormal10",
            "calculator": {"TenBerge": {"breathingRate": 10}},
            "parameters": {
                "type": "Lognormal10DoseResponse",
                "levels": ["Severe", "Mild"],
                "parameters": {
                    "Severe": {"TL_50": 1000, "sigma": 0.5},
                    "Mild": {"TL_50": 100, "sigma": 0.4}
                }
            }
        }
    }
}
```

---

## Effects system

The effects system follows a layered architecture. Each layer is resolved dynamically from JSON descriptors using `pydoc.locate`.

### Calculators (`agents/effects/Calculator.py`)

Calculators compute the toxic load (dose) from a concentration field. All inherit from `AbstractCalculator`, which stores the reference `injuryBreathingRate` (default `10*ureg.L/ureg.min`). Each implements `calculate(concentrationField, field, breathingRate=10*ureg.L/ureg.min, time="datetime", inUnits=None)`.

| Calculator | Class | Formula | Input |
|-----------|-------|---------|-------|
| Haber | `CalculatorHaber` | `D(T) = integral(C × dt) × breathingRatio` | `pandas.DataFrame` or `xarray.Dataset` |
| Ten Berge | `CalculatorTenBerge` | `D(T) = integral(C^n × dt) × breathingRatio` | Same, plus `tenbergeCoefficient` from agent |
| Max Concentration | `CalculatorMaxConcentration` | Maximum of the rolling-mean concentration over a `sampling` window | Same, plus `sampling` parameter (e.g., `"10min"`) |

#### Implementation details

- **Unit handling**: All calculators convert concentrations to `mg/m³` internally using `inUnits.m_as(ureg.mg / ureg.m**3)`. xarray attributes are checked for field-specific units.
- **Breathing rate ratio**: `breathingRatio = (breathingRate / self.injuryBreathingRate).magnitude` — scales the dose based on actual vs. reference breathing rate.
- **Time integration**: For `pandas.DataFrame`, time steps are computed via `diff().apply(lambda x: x.seconds)/60.` (non-uniform spacing supported). For `xarray.Dataset`, uses `cumsum(dim=time) * concentrationField.dt` (the dataset must expose a `dt` quantity).

### Injury levels (`agents/effects/InjuryLevel.py`)

Injury levels map toxic loads to the fraction of population affected. All inherit from `InjuryLevel`.

| Level type | Class | `getPercent(ToxicLoad)` implementation |
|-----------|-------|---------------------------------------|
| Lognormal10 | `InjuryLevelLognormal10DoseResponse` | `lognorm.cdf(ToxicLoad/a, sigma/a, scale=TL_50/a)` where `a = log10(e)` (`getToxicLoad(Percent)` is the inverse) |
| Threshold | `InjuryLevelThreshold` | `1 if ToxicLoad > threshold else 0` (vectorized via `numpy.atleast_1d`) |
| Exponential | `InjuryLevelExponential` | `1 - numpy.exp(-k * ToxicLoad)` |

#### Contour generation (`calculateContours`)

Each `InjuryLevel` can generate spatial contour polygons from a 2D toxic load field:

1. The `_getGeopandas()` method calls `getContoursFromLevels()`, which runs `matplotlib.pyplot.contour()` on the toxic load grid
2. For **Lognormal10DoseResponse**: generates contours at the toxic loads of 5%, 10%, ..., 95% (plus the higher-severity level's toxic loads to avoid double-counting)
3. For **Threshold**: generates a single contour at the threshold value
4. For **Exponential**: generates a contour at the `k` value
5. The contour collection is converted to a `GeoDataFrame` using `hera.utils.matplotlibCountour.toGeopandas()`

`calculateContours(toxicLoads, time="datetime", x="x", y="y")` iterates over time steps, calling `_getGeopandas()` for each, drops polygons under 1 m², dissolves by level, and returns a `GeoDataFrame` with columns `ToxicLoad`, `TotalPolygon` (geometry), `severity`, `datetime` (or `None` if no contour). The injury's `_postCalculate()` then adds `percentEffected` and the differential `ThresholdPolygon` (via `Injury.calculateThresholdPolygon`).

#### Higher-severity subtraction

When computing `Injury.getPercent(level, ToxicLoad)`, the percentage at each severity level is the **difference** from the level above:

```python
# For the most severe level (index 0):
val = levels[0].getPercent(ToxicLoad)

# For subsequent levels:
val = levels[i].getPercent(ToxicLoad) - levels[i-1].getPercent(ToxicLoad)
```

This ensures that the sum of all severity levels doesn't exceed 100%.

### Injuries (`agents/effects/Injury.py`)

An `Injury` combines a calculator with one or more injury levels:

| Class | Description |
|-------|-------------|
| `Injury` | Base class holding calculator + levels dict |
| `InjuryLognormal10` | Injury using Lognormal10 dose-response (levels: `InjuryLevelLognormal10DoseResponse`) |
| `InjuryThreshold` | Injury using threshold dose-response |
| `InjuryExponential` | Injury using exponential dose-response |

#### Injury initialization

When `Injury(name, cfgJSON, calculator, units=None)` is created:

1. Resolve the `InjuryLevel` class from `cfgJSON["type"]` via `pydoc.locate("hera.riskassessment.agents.effects.InjuryLevel.InjuryLevel<type>")`
2. For each level name in `cfgJSON["levels"]` (ordered from most severe to least):
   - Extract parameters from `cfgJSON["parameters"][levelName]`
   - Inject `higher_severity` reference to the previous level (for subtraction in `getPercent`)
   - Create the `InjuryLevel` instance
3. Store in `self._levels` (ordered list) and `self._levelsmap` (name → level dict)

#### Key methods

| Method | What it does |
|--------|-------------|
| `calculateRegionOfInjured(concentrationField, field, time="datetime", x="x", y="y", breathingRate=..., sel={}, isel={})` | Full pipeline: calculator → toxic loads → contour polygons for each level → `thresholdGeoDataFrame` |
| `calculate(...)` | Deprecated alias of `calculateRegionOfInjured` (emits a warning) |
| `calculateToxicLoads(concentrationField, time="datetime", breathingRate=..., field=None)` | Calculator step only — returns toxic loads (xarray/pandas) |
| `calculatePointWiseFractionInjured(timeConcentration, time="datetime", breathingRate=..., field=None)` | Fraction injured per level at point (device) time series |
| `getPercent(level, ToxicLoad)` | Get fraction affected at a specific severity level (with subtraction) |

### Factory pattern (`InjuryFactory`)

`InjuryFactory.getInjury(name, cfgJSON, **additionalparameters)` dynamically resolves:

1. **Calculator**: `pydoc.locate("hera.riskassessment.agents.effects.Calculator.Calculator<calcType>")` — instantiated with `calcParam` + `additionalparameters` (e.g., `tenbergeCoefficient`)
2. **Injury**: `pydoc.locate("hera.riskassessment.agents.effects.Injury.Injury<injuryType>")` — instantiated with `(name, cfgJSON["parameters"], calculator, units=cfgJSON.get("units"))`

A module-level instance is exported as `hera.riskassessment.agents.effects.injuryfactory` and is what `Agent` uses.

The naming convention is strict — class names must match `Calculator<Name>`, `Injury<Name>` and `InjuryLevel<Name>` exactly (e.g. `"type": "Lognormal10"` for the injury and `"type": "Lognormal10DoseResponse"` for its levels).

---

## Protection policies

`ProtectionPolicy` implements a chainable action pipeline that modifies concentration fields:

```python
policy = ProtectionPolicy()          # or ProtectionPolicy(actionList=[{"name": "indoor", "params": {...}}, ...])
policy.indoor(turnover=2*ureg.hour, enter="5min", stay="2h")
policy.masks(protectionFactor=1000, wear="0min", duration="3h")
result = policy.compute(concentration_data, C="C")
```

Actions can also be added with `addAction(name, params)` or `addActions(jsonStrOrFile)` (a dict/JSON with an `"actions"` list). The action class is resolved as `Action<name.title()>` via `pydoc.locate`.

### Action classes

| Action | Class | Parameters |
|--------|-------|-----------|
| Indoor sheltering | `ActionIndoor` | `turnover` (or `alpha`; see note below), time window (`begin/end` or `enter/stay`) |
| Masking | `ActionMasks` | `protectionFactor`, time window (`begin/end` or `wear/duration`) |

Each action:
1. Reads the current concentration field (`policy.finalname`, always `"C"`)
2. Applies its transformation (indoor air exchange model, or division by protection factor) and stores it as a new field named `<type>_<actionid>` (e.g. `indoor_1`, `masks_2`)
3. Overwrites the `finalname` field with the transformed values inside the action's time window (outside it the previous values are kept), so the next action builds on it
4. Records its parameters in `data.attrs[str(actionid)]` for traceability

!!! warning "Known issue in master"
    `ActionIndoor(alpha=...)` calls `kwargs["alpha"].ureg(1/ureg.h)` (`protectionpolicy/ProtectionPolicy.py:355`), which raises `AttributeError`; use `turnover` instead.

### Indoor model

The indoor concentration `Cin` is computed iteratively:

```
Cin[t] = (Cin[t-1] + alpha * dt * Cout[t]) / (1 + alpha * dt)
```

Where `alpha = 1/turnover` is the air exchange rate.

---

## Risk area calculation

`getRiskAreaAlgorithm(algorithmName, **kwargs)` resolves `riskAreaAlgorithm_<Name>` (only `Sweep` exists). `riskAreaAlgorithm_Sweep(dxdy=150, outlayers=3, parallel=True)` computes casualty counts across a spatial grid via `calculate(effectIsopleths, demog, mathematical_angle=None, meteorological_angle=None, severityColumn="severity", valueColumn="effectedPopulation")`:

1. **Build bounding box** — determine the search area from the demographic region + effect polygon extent
2. **Grid the area** — create a mesh of points at `dxdy` spacing
3. **For each grid point** — project the effect isopleths onto the demographic data at that release location and wind direction
4. **Aggregate** — sum affected population by severity level

Supports parallel execution via `multiprocessing.Pool` (`workerCount` defaults to `cpu_count()`).

!!! warning "Known issue in master"
    `_doCalculation` calls `effectIsopleths.datalayer(...)` (`analysis/riskAreas.py:175`), a method `thresholdGeoDataFrame` does not have (it should be `project`), so `calculate()` fails. The default `valueColumn="effectedPopulation"` also does not match the `effectedtotal_pop` column produced by `project()`.

---

## Spatial results: thresholdGeoDataFrame (`agents/effects/thresholdGeoDataFrame.py`)

`thresholdGeoDataFrame` extends `geopandas.GeoDataFrame` and is the primary data structure for spatial injury results. It holds contour polygons generated by `InjuryLevel.calculateContours()` and provides methods to place them on a real map and project onto population data.

### Expected columns

A `thresholdGeoDataFrame` typically has these columns (produced by `calculateContours`):

| Column | Type | Description |
|--------|------|-------------|
| `datetime` | timestamp | The time step |
| `severity` | str | Injury level name (e.g., `"Severe"`, `"Mild"`) |
| `percentEffected` | float | Fraction of population affected (0–1) |
| `TotalPolygon` | geometry | Cumulative contour polygon |
| `ThresholdPolygon` | geometry | Difference from the next-higher toxic load contour (avoids double-counting); the default geometry for `shiftLocationAndAngle`/`project` |
| `ToxicLoad` | float | The toxic load value for this contour |

### Methods

#### `shiftLocationAndAngle(loc, meteorological_angle=None, mathematical_angle=None, geometry="ThresholdPolygon")`

Returns a **new** `thresholdGeoDataFrame` with polygons rotated and translated to a release location. The original contours are computed relative to the source at `(0, 0)` — this method places them on a real map.

**Implementation:**

1. Set the active geometry to `geometry` column
2. Call `_shiftPolygons()` which:
   - Converts meteorological angle to mathematical angle if needed
   - Rotates all polygons around `(0, 0, 0)` by the wind angle using `self.rotate(angle, origin=(0,0,0))`
   - Translates to the release location using `self.translate(*loc)`
3. Return a copy with the shifted geometry

#### `project(demographic, loc, meteorological_angle=None, mathematical_angle=None, geometry="ThresholdPolygon", population="total_pop", projectName=None)`

Projects the threshold polygons onto a demographic dataset to estimate casualties. This is the main method for going from "injury zones" to "number of people affected".

**Supports three calling patterns:**

1. **Single angle**: `project(demographic, loc, mathematical_angle=45)` — single projection
2. **List of meteorological angles**: `project(demographic, loc, meteorological_angle=[0, 90, 180, 270])` — iterates and concatenates results, adding `meteorological_angle`, `mathematical_angle`, `mathematical_angle_rad` columns
3. **List of mathematical angles**: Same as above with `mathematical_angle=[...]`

In the list forms, `population` is not forwarded, so the default `"total_pop"` is always used. `projectName` is unused (kept for backward compatibility).

**Implementation of `_project()` (single angle):**

1. Convert demographic data to ITM (EPSG:2039) coordinate system
2. Call `_shiftPolygons()` to position the contours at the release location and wind angle
3. Group the contours by `(severity, datetime)`
4. For each contour polygon:
   - Call the static `_calculatePopulationInPolygon(demog_data, polygon, population)` to intersect with the census data (area-weighted)
   - Multiply each population column by `percentEffected` to get `effected<population>` (e.g. `effectedtotal_pop`)
   - Record `ToxicLoad`, `severity`, `datetime`
5. Concatenate all results into a single DataFrame

**Output columns:**

| Column | Description |
|--------|-------------|
| `effectedtotal_pop` | Number of people affected (population × percentEffected) |
| `percentEffected` | Fraction affected at this contour level |
| `ToxicLoad` | The toxic load value |
| `severity` | Injury severity level name |
| `datetime` | Time step |
| `total_pop` | Area-weighted population in the intersected area (one column per requested population type) |
| `areaFraction` | Intersected area / census polygon area |
| `geometry` | Intersected census polygon geometry |

#### Population calculation

The spatial intersection with demographic data is done directly inside `_project()` via a static `_calculatePopulationInPolygon()` method, which performs area-weighted population calculation without depending on `DemographyToolkit` or any project. The `demographic` GeoDataFrame is passed in by the caller, so the data source is fully under user control.

!!! warning "Known issue in master"
    `_calculatePopulationInPolygon` filters with `demography.loc[not demography["geometry"].intersection(poly).is_empty]` (`agents/effects/thresholdGeoDataFrame.py:125`); `not` on a multi-row boolean Series raises `ValueError`, so `project()` fails for demographic layers with more than one polygon.

### Data flow

```
InjuryLevel.calculateContours()
    ↓ produces
thresholdGeoDataFrame (contours at origin)
    ↓ shiftLocationAndAngle(loc, angle)
thresholdGeoDataFrame (positioned on map)
    ↓ project(demographic, loc, angle)
pandas.DataFrame (casualties per severity per time step)
    ↓ groupby("severity").sum()
Total casualties per severity level
```

---

## Adding new effects

1. **New calculator:** Create a class in `agents/effects/Calculator.py` named `Calculator<Name>` inheriting from `AbstractCalculator` and implementing `calculate(concentrationField, field, breathingRate, time, inUnits)`
2. **New injury level:** Create a class in `agents/effects/InjuryLevel.py` named `InjuryLevel<Name>` inheriting from `InjuryLevel`
3. **New injury type:** Create a class in `agents/effects/Injury.py` named `Injury<Name>` inheriting from `Injury`

The factory uses `pydoc.locate` to find classes by name, so naming conventions must be followed exactly.

## Adding new protection actions

1. Create a class in `protectionpolicy/ProtectionPolicy.py` named `Action<Name>` inheriting from `abstractAction`
2. Implement `compute()` to modify `self.policy.data`
3. Add a convenience method on `ProtectionPolicy` (e.g., `def evacuation(self, **kwargs)`)

## CLI

`hera-riskassessment agents createRepository <repository_name> [--path DIR]` (`CLI.createRepository`) scans `DIR` (default: current directory) for JSON files containing an `effectParameters` key and writes a repository JSON with one `RiskAssessment` → `DataSource` entry per agent (`dataFormat: "JSON_dict"`, the descriptor embedded as `resource`, `isRelativePath: "False"`).

For the full API, see the [Risk Assessment API Reference](api/risk.md).
