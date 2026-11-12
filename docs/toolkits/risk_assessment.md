# Risk Assessment Toolkit

**Toolkit name:** `RiskAssessment` | **Constant:** `toolkitHome.RISKASSESSMENT`

The Risk Assessment toolkit provides an agent-based framework for modeling hazards, calculating injury effects, applying protection policies, and estimating casualties.

```python
from hera import toolkitHome

# Tip: if you created the project with `hera-project project create`, you can omit projectName
risk = toolkitHome.getToolkit(toolkitHome.RISKASSESSMENT, projectName="MY_PROJECT")
```

---

## Concepts

The risk assessment framework models a chain of events:

1. **Hazard source** releases a dangerous agent (chemical, thermal, blast)
2. The agent **disperses** through the environment (using simulation toolkits like Gaussian or LSM)
3. **Effects** are calculated on the exposed population based on dose, concentration, or overpressure
4. **Protection policies** (sheltering, evacuation) modify the exposure
5. **Casualties** are estimated based on injury levels

```mermaid
flowchart LR
    Agent["Agent\n(hazardous material)"] --> Effects["Effects\n(injury models)"]
    Effects --> Protection["Protection\n(policies)"]
    Protection --> Casualties["Casualties\n(estimates)"]
```

---

## Agents

An **Agent** represents a hazardous material with its associated effects (injury models). Agents are stored as `JSON_dict` data sources of the toolkit, usually loaded from a repository.

```python
# List available agents (the data source names)
risk.listAgentsNames()
# ['Chlorine', 'Ammonia', 'Propane']

# Load an agent (by data source name, optionally a version; or pass the descriptor dict itself)
agent = risk.getAgent("Chlorine", version=None)

# Access agent properties
print(agent.name)                # "Chlorine"
print(agent.effectNames)         # ['inhalation', 'thermal']
print(agent.tenbergeCoefficient) # Tenberge coefficient for the agent

# Access a specific effect
inhalation_effect = agent["inhalation"]   # effects are also attributes: agent.inhalation
```

An agent descriptor has the form:

```json
{
    "name": "Chlorine",
    "version": [0, 0, 1],
    "effectParameters": {"tenbergeCoefficient": 2.0},
    "effects": {"RegularPopulation": {"type": "Lognormal10", "calculator": {...}, "parameters": {...}}},
    "physicalProperties": {...}
}
```

### Adding agents to the project

```python
# From a JSON file, a JSON string or a dict (the "name" key is required)
risk.loadData("/path/to/Chlorine.json")
# saveMode=risk.TOOLKIT_SAVEMODE_FILEANDDB (default) raises ValueError if the agent/version exists;
# use risk.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE to replace it

# Or with an explicit name and version
risk.loadAgent("Chlorine", agentDescription=descriptor, version=[0, 0, 1])
```

### Creating an agent repository

From the CLI:

```bash
# Create a repository JSON named "myAgents" in --path (default: the current directory)
# from the agent JSON files found there (files with an "effectParameters" key)
hera-riskassessment agents createRepository myAgents --path /path/to/agents/
```

Register and load the resulting repository with `hera-project repository add` / `hera-project repository load` (see the [CLI Reference](../cli/reference.md)).

!!! note
    `createRepository` opens the listed file names relative to the **current directory**, so run it from inside the agents directory (`--path` must point to the same place). The agents are embedded in the repository (`isRelativePath: "False"`).

---

## Effects and injury levels

Each agent has one or more **effects** — named injury models that describe how exposure to the agent harms a population. An effect combines two things: a **calculator** that computes the toxic load from concentration data, and one or more **injury levels** that map the toxic load to a percentage of people affected.

### The chain: Concentration → Toxic Load → Injury

```
Concentration field    →    Calculator    →    Toxic Load    →    Injury Levels    →    % Affected
(from dispersion)          (Haber/TenBerge)    (cumulative)       (Severe/Mild/...)     (per level)
```

1. A dispersion simulation (LSM or Gaussian) produces a **concentration field** — concentration values over space and time
2. A **calculator** integrates the concentration over time to produce a **toxic load** (dose)
3. **Injury levels** convert the toxic load into a percentage of the population affected at each severity

### Calculators

Calculators determine *how* concentration is converted to dose. The choice depends on the toxicological model for the agent:

| Calculator | What it computes | When to use |
|-----------|-----------------|-------------|
| **Haber** (`"Haber"`, `CalculatorHaber`) | `D(T) = integral(C dt)` — simple time integral | Agents where dose is proportional to concentration × time |
| **Ten Berge** (`"TenBerge"`, `CalculatorTenBerge`) | `D(T) = integral(C^n dt)` — concentration raised to power n | Agents where higher concentrations are disproportionately more dangerous (most real chemicals) |
| **Max Concentration** (`"MaxConcentration"`, `CalculatorMaxConcentration`) | Maximum over time of the rolling mean over `sampling` (e.g. `"10min"`) | Agents where the peak matters, not cumulative dose |

The quoted name is the key used under `"calculator"` in the effect JSON; the class lives in `hera.riskassessment.agents.effects.Calculator`.

The Ten Berge exponent `n` (the `tenbergeCoefficient`) is an agent property — it controls how much weight is given to high vs. low concentrations. When `n=1`, Ten Berge reduces to Haber.

All calculators accept concentration data as either `pandas.DataFrame` (time index, one column per point, mg/m³) or `xarray.Dataset`, and scale the result by the ratio of the population's breathing rate to the calculator's `breathingRate` (default 10 L/min). For an `xarray.Dataset` the concentration units are read from `attrs[field]` and the time step from `attrs["dt"]` — the LSM `getConcentration` output provides both. Haber/Ten Berge return the cumulative toxic load at every time step; MaxConcentration returns a single map (no time dimension).

### Injury levels

Injury levels define the **dose-response relationship** — given a toxic load, what percentage of the population is affected? Each effect can have multiple severity levels (e.g., Severe, Mild, Light):

| Injury level type | Dose-response model | Parameters |
|------------------|--------------------|-----------|
| **Lognormal10** | Log-normal CDF (base 10) | `TL_50` (toxic load at 50% effect), `sigma` (spread) |
| **Threshold** | Binary: 0% below, 100% above | `threshold` value |
| **Exponential** | Exponential curve | `k` (rate parameter) |

#### Lognormal10

The most common model. Uses a log-normal cumulative distribution function (base 10) with two parameters per severity level:

- **TL_50** — the toxic load at which 50% of the population is affected
- **sigma** — how spread out the dose-response curve is (larger = more gradual transition)

The percentage affected at a given toxic load `D` is: `P(D) = CDF_lognormal(D; TL_50, sigma)`

**JSON format:**

```json
{
    "effects": {
        "RegularPopulation": {
            "type": "Lognormal10",
            "calculator": {
                "TenBerge": {"breathingRate": 10}
            },
            "parameters": {
                "type": "Lognormal10DoseResponse",
                "levels": ["Severe", "Mild", "Light"],
                "parameters": {
                    "Severe": {"TL_50": 1000, "sigma": 0.5},
                    "Mild":   {"TL_50": 100,  "sigma": 0.4},
                    "Light":  {"TL_50": 10,   "sigma": 0.3}
                }
            }
        }
    }
}
```

This defines an effect called `"RegularPopulation"` that uses the Ten Berge calculator and has three severity levels. At a toxic load of 1000, 50% of the population has Severe injuries.

#### Threshold

A binary model — 0% affected below the threshold, 100% above. Commonly used with the **MaxConcentration** calculator for acute exposure standards (e.g., AEGL, ERPG levels).

- **threshold** — the concentration value (with units) above which the population is affected

**JSON format:**

```json
{
    "effects": {
        "AEGL": {
            "type": "Threshold",
            "calculator": {
                "MaxConcentration": {"sampling": "10min"}
            },
            "parameters": {
                "type": "Threshold",
                "levels": ["AEGL-3", "AEGL-2", "AEGL-1"],
                "parameters": {
                    "AEGL-3": {"threshold": "50.0*mg/m**3"},
                    "AEGL-2": {"threshold": "8.6*mg/m**3"},
                    "AEGL-1": {"threshold": "1.5*mg/m**3"}
                }
            }
        }
    }
}
```

This defines three AEGL levels with concentration thresholds. The MaxConcentration calculator uses the peak concentration over a sampling period instead of a cumulative dose. At any point where the concentration exceeds the threshold, 100% of the population at that location is affected at that level.

#### Exponential

A smooth dose-response curve: `P(D) = 1 - exp(-k × D)`. The percentage affected increases continuously with toxic load, approaching 100% asymptotically.

- **k** — the rate parameter (larger = steeper response)

**JSON format:**

```json
{
    "effects": {
        "ExponentialEffect": {
            "type": "Exponential",
            "calculator": {
                "Haber": {"breathingRate": 10}
            },
            "parameters": {
                "type": "Exponential",
                "levels": ["HighExposure", "LowExposure"],
                "parameters": {
                    "HighExposure": {"k": 0.01},
                    "LowExposure":  {"k": 0.001}
                }
            }
        }
    }
}
```

### Combining calculators and injury levels

Any calculator can be paired with any injury level type. Common combinations:

| Use case | Calculator | Injury level | Example |
|---------|-----------|-------------|---------|
| Toxic gas (cumulative dose) | TenBerge | Lognormal10 | Chlorine inhalation over time |
| Acute exposure standards | MaxConcentration | Threshold | AEGL/ERPG concentration limits |
| Simple dose-response | Haber | Exponential | Linear dose accumulation |
| Sensitivity analysis | TenBerge | Lognormal10 | Vary `tenbergeCoefficient` to test assumptions |

### Working with effects in code

#### Loading an agent and exploring its effects

```python
risk = toolkitHome.getToolkit(toolkitHome.RISKASSESSMENT, projectName="MY_PROJECT")

# Load an agent
agent = risk.getAgent("Chlorine")

# Agent-level properties
print(agent.name)                  # 'Chlorine'
print(agent.tenbergeCoefficient)   # 2.0
print(agent.effectNames)           # ['RegularPopulation']
```

#### Inspecting an effect

```python
# Access an effect by name (dictionary-style)
effect = agent["RegularPopulation"]

# The effect's calculator
print(type(effect.calculator))     # <class 'CalculatorTenBerge'>

# List all severity levels in this effect
print(effect.levelNames)           # ['Severe', 'Mild', 'Light']
```

#### Inspecting injury levels

Each severity level has its own dose-response parameters:

```python
# Access a specific injury level
severe = effect["Severe"]
mild = effect["Mild"]
light = effect["Light"]

# For Lognormal10 levels, inspect the parameters:
# TL_50 — the toxic load at which 50% of the population is affected
# sigma — the spread of the dose-response curve
print(f"Severe: TL_50={severe.TL_50}, sigma={severe.sigma}")
print(f"Mild:   TL_50={mild.TL_50}, sigma={mild.sigma}")
print(f"Light:  TL_50={light.TL_50}, sigma={light.sigma}")

# For Threshold levels:
# threshold — the binary cutoff value
# print(f"Threshold: {level.threshold}")

# For Exponential levels:
# k — the rate parameter
# print(f"k: {level.k}")
```

#### Computing toxic loads and injury percentages

```python
# Given a concentration field from a dispersion simulation:
# concentration_data is an xarray.Dataset with a "C" field

# Step 1: Calculate toxic loads using the effect
toxic_loads = effect.calculateToxicLoads(concentration_data, field="C")

# Step 2: Get the percentage affected at each severity level
# for a given toxic load value
toxic_load_value = 500
percent_severe = effect.getPercent("Severe", toxic_load_value)
percent_mild = effect.getPercent("Mild", toxic_load_value)
print(f"At toxic load {toxic_load_value}: {percent_severe:.1%} severe, {percent_mild:.1%} mild")

# Step 3: Calculate the injury regions on a 2D map (all severity levels)
regions = effect.calculateRegionOfInjured(concentration_data, field="C", time="datetime", x="x", y="y")
# Returns a thresholdGeoDataFrame with the columns:
#   datetime | severity | ToxicLoad | TotalPolygon | ThresholdPolygon | percentEffected
# (ThresholdPolygon = the area of this level minus the area of the level above)

# The contours of a single level from precomputed toxic loads (a plain GeoDataFrame
# with ToxicLoad, TotalPolygon, severity, datetime — no percentEffected/ThresholdPolygon)
contours = severe.calculateContours(toxic_loads, time="datetime", x="x", y="y")
```

For point time-series (a `pandas.DataFrame` with a time index and one column per point), `effect.calculatePointWiseFractionInjured(df)` returns the injured fraction per point, time and level.

#### Modifying the Ten Berge coefficient

The Ten Berge coefficient can be changed after loading — this rebuilds all effects:

```python
# Change the exponent (e.g., for sensitivity analysis)
agent.tenbergeCoefficient = 1.5
# All effects are automatically recalculated with the new coefficient
```

#### Serializing an agent

```python
# Export agent to JSON (for saving or inspection)
import json
print(json.dumps(agent.toJSON(), indent=2))
```

### Physical properties

Agents can also have physical properties (the `physicalProperties` key of the descriptor: `molecularWeight`, `sorptionCoefficient`, `spreadFactor`, `volatilityConstants`, `densityConstants`, `vaporPressure`) used for evaporation and dispersion modeling. Pass temperatures as pint quantities:

```python
# Access physical properties
props = agent.physicalproperties

# Molecular weight, density, vapor pressure at a temperature
mw = props.molecularWeight                                 # pint quantity in g/mol
density = props.getDensity(ureg.Quantity(20, ureg.degC))     # g/cm³ at 20°C
volatility = props.getVolatility(ureg.Quantity(20, ureg.degC))  # vapor saturation, g/cm³
vp = props.vaporPressure(293 * ureg.K)                       # vapor pressure in bar (float)
```

---

## Protection policies

Protection policies model how sheltering, evacuation, or other protective actions reduce exposure. A `ProtectionPolicy` builds a **pipeline of actions** that modify the concentration field before injury calculations.

### Indoor sheltering

The indoor model computes the indoor concentration `Cin` based on the outdoor concentration `Cout` and the building's air exchange rate:

```
Cin[t] = (Cin[t-1] + alpha × dt × Cout[t]) / (1 + alpha × dt)
```

Where `alpha = 1/turnover` is the air exchange rate. A longer turnover time means better protection (slower air exchange).

```python
from hera.riskassessment import ProtectionPolicy   # also available as risk.ProtectionPolicy
from hera.utils import ureg

# Create a policy with indoor sheltering
policy = ProtectionPolicy()
policy.indoor(
    turnover=2*ureg.hour,    # air exchange every 2 hours
    enter="5min",            # population enters buildings 5 min after release
    stay="2h"                # stays indoor for 2 hours
)

# Apply to a concentration field (an xarray.Dataset with x, y, datetime)
protected = policy.compute(concentration_data, C="C")
# Returns an xarray.Dataset: "outdoor_0", one variable per action, and the final
# concentration in "C" (the protected value inside the action's time window)
```

Instead of `enter`/`stay` (timedelta strings relative to the first time step) you can give absolute `begin`/`end` times; `alpha` (1/time) can replace `turnover`. Without a time window the action applies to the whole period. A policy can also be built from a list of actions: `ProtectionPolicy([{"name": "Indoor", "params": {...}}, {"name": "Masks", "params": {...}}])`, or with `addActions(jsonFileOrDict)` (actions under an `"actions"` key).

### Masking

Masks reduce the inhaled concentration by a protection factor:

```python
policy.masks(
    protectionFactor=1000,   # mask reduces concentration by 1000x
    wear="0min",             # masks on immediately (or begin=/end=)
    duration="3h"            # worn for 3 hours
)
```

### Chaining actions

Actions can be chained — each modifies the concentration field for the next:

```python
policy = ProtectionPolicy()
policy.indoor(turnover=2*ureg.hour, enter="5min", stay="2h") \
      .masks(protectionFactor=100, wear="0min", duration="3h")
result = policy.compute(concentration_data, C="C")
```

---

## Analysis

`risk.analysis.getRiskAreas` computes the extent of the area where the Ten Berge toxic load of an LSM simulation (at the lowest `z` level, last time step) exceeds given levels:

```python
areas = risk.analysis.getRiskAreas(
    tenbergeCoefficient=2.0,
    levels=[10, 100],              # toxic-load levels
    Q=1 * ureg.kg,                 # released mass
    protectionPolicy=None,         # optional action list for ProtectionPolicy
    LSMfile=None,                  # a simulation directory, or None to query the LSM toolkit:
    simulationName="run_001",      # ... extra keyword arguments go to lsm.getSimulations
)
# GeoDataFrame: tenBergCoefficient | level | bounds | geometry (convex hull)
```

For casualty estimates over many release points, `hera.riskassessment.getRiskAreaAlgorithm("Sweep", dxdy=150, outlayers=3, parallel=True)` returns an estimator whose `calculate(effectIsopleths, demog, meteorological_angle=...)` sweeps release locations around the populated area.

---

## Presentation

`risk.presentation` (`casualtiesPlot`) draws casualty figures from the result of `calculateRegionOfInjured`:

- `plotCasualtiesRose(results, area, severityList, loc, meteorological_angles=None, mathematical_angles=None, effectedPopulation="effectedtotal_pop", ...)` — radial bars of the affected population for each wind direction
- `plotCasualtiesProjection(results, area, severityList, loc, meteorological_angle=None, mathematical_angle=None, plumSeverity=[], ...)` — the injury polygons projected on the population map

Both rely on `thresholdGeoDataFrame.project`, which is currently broken (see the note in the examples).

---

## Integration with other toolkits

Risk assessment typically uses data from several other toolkits:

| Input | Source toolkit | Purpose |
|-------|--------------|---------|
| Population data | `GIS_Demography` | Number and distribution of people at risk |
| Building footprints | `GIS_Buildings` | Sheltering factors and indoor/outdoor ratios |
| Terrain | `GIS_Raster_Topography` | Elevation effects on dispersion |
| Weather | `MeteoLowFreq` | Wind speed, direction, stability for dispersion |
| Dispersion results | `GaussianDispersion` or `LSM` | Concentration fields from source to receptor |

```python
# A typical risk assessment workflow
topo = toolkitHome.getToolkit(toolkitHome.GIS_RASTER_TOPOGRAPHY, projectName="MY_PROJECT")
demo = toolkitHome.getToolkit(toolkitHome.GIS_DEMOGRAPHY, projectName="MY_PROJECT")
meteo = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName="MY_PROJECT")
risk = toolkitHome.getToolkit(toolkitHome.RISKASSESSMENT, projectName="MY_PROJECT")

# Each toolkit contributes its domain data to the risk analysis
```

For the full API details, see the [Toolkit Catalog](overview.md) and the [API Reference](../developer_guide/api/risk.md).

---

## Complete Examples

### Example 1: Toxic gas release with Lognormal10 dose-response

This example models a chlorine gas release using the LSM dispersion model and Lognormal10 injury levels with the Ten Berge calculator. It covers the full pipeline: dispersion → concentration → toxic load → injury contours → casualty estimation.

!!! bug "Known issue"
    The casualty steps (`thresholdGeoDataFrame.project`, and therefore `risk.presentation.plotCasualtiesRose` / `plotCasualtiesProjection`) currently fail with `ValueError: The truth value of a Series is ambiguous` (`hera/riskassessment/agents/effects/thresholdGeoDataFrame.py:125`). The steps up to the injury regions work.

```python
from hera import toolkitHome
from hera.utils import ureg

# ---------------------------------------------------------------
# Step 1: Set up toolkits
# ---------------------------------------------------------------
PROJECT = "ChlorineRelease"

lsm  = toolkitHome.getToolkit(toolkitHome.LSM, projectName=PROJECT)
risk = toolkitHome.getToolkit(toolkitHome.RISKASSESSMENT, projectName=PROJECT)
demo = toolkitHome.getToolkit(toolkitHome.GIS_DEMOGRAPHY, projectName=PROJECT)

# ---------------------------------------------------------------
# Step 2: Run or retrieve an LSM dispersion simulation
# ---------------------------------------------------------------
# Option A: run a new simulation from a template
template = lsm.getTemplateByName("urban_chlorine")
sim = template.run(                    # returns a SingleSimulation (to_xarray=True)
    topography=topography_str,         # topography file content, or None for flat terrain
    stations=weather_stations_df,
    simulationName="chlorine_run_001",
    windSpeed=3.0,
    windDirection=225,
    releaseRate=5.0
)

# Option B: retrieve an existing simulation
sim = lsm.getSimulations(simulationName="chlorine_run_001")[0]

# ---------------------------------------------------------------
# Step 3: Get concentration field
# ---------------------------------------------------------------
# Q = total released mass (scales the normalized LSM output)
concentration = sim.getConcentration(Q=500 * ureg.kg)
# Result: xarray.Dataset with 'C' field in mg/m³ over (datetime, x, y, z),
# with attrs dt, Q and C (units) used by the calculators

# ---------------------------------------------------------------
# Step 4: Load the agent and compute toxic loads
# ---------------------------------------------------------------
agent = risk.getAgent("Chlorine")
print(f"Agent: {agent.name}")
print(f"Ten Berge coefficient: {agent.tenbergeCoefficient}")
print(f"Effects: {agent.effectNames}")

# The agent has a Lognormal10 effect with Ten Berge calculator
effect = agent["RegularPopulation"]
print(f"Calculator: {type(effect.calculator).__name__}")  # CalculatorTenBerge
print(f"Severity levels: {effect.levelNames}")            # ['Severe', 'Mild', 'Light']

# Use the concentration near the ground (the z levels depend on the template)
ground = concentration.sel(z=10, method="nearest")

# Compute toxic loads: integral(C^n dt) over time (cumulative, per time step)
toxic_loads = effect.calculateToxicLoads(ground, field="C")

# ---------------------------------------------------------------
# Step 5: Calculate injury regions
# ---------------------------------------------------------------
# Polygons of all the severity levels (computes the toxic loads internally)
regions = effect.calculateRegionOfInjured(ground, field="C", time="datetime", x="x", y="y")
# Result: thresholdGeoDataFrame with columns:
#   datetime | severity | ToxicLoad | TotalPolygon | ThresholdPolygon | percentEffected

# ---------------------------------------------------------------
# Step 6: Project onto population data for casualty estimates
# ---------------------------------------------------------------
# Load demographic data
population = demo.getDataSourceData("census_2020")

# Project injury contours onto population
release_location = (178000, 665000)  # ITM coordinates
wind_direction_math = 45  # mathematical angle

casualties = regions.project(
    demographic=population,
    loc=release_location,
    mathematical_angle=wind_direction_math
)

# Result: DataFrame with affected population per severity level
print(casualties.groupby("severity")["effectedtotal_pop"].sum())

# ---------------------------------------------------------------
# Step 7: Apply protection policy (optional)
# ---------------------------------------------------------------
from hera.riskassessment import ProtectionPolicy
from hera.utils import ureg

policy = ProtectionPolicy()
policy.indoor(turnover=2 * ureg.hour, enter="10min", stay="3h")

# Apply protection to the concentration field before computing toxic loads
protected_concentration = policy.compute(ground, C="C")   # keeps the attrs of the input

# Recalculate with protected concentrations
protected_regions = effect.calculateRegionOfInjured(protected_concentration, field="C")

# Compare casualties with and without protection
protected_casualties = protected_regions.project(
    demographic=population,
    loc=release_location,
    mathematical_angle=wind_direction_math
)
print("Without protection:", casualties["effectedtotal_pop"].sum())
print("With indoor shelter:", protected_casualties["effectedtotal_pop"].sum())

# ---------------------------------------------------------------
# Step 8: Visualize
# ---------------------------------------------------------------
risk.presentation.plotCasualtiesRose(
    results=regions,
    area=population,
    severityList=["Severe", "Mild"],
    loc=release_location,
    meteorological_angles=[0, 45, 90, 135, 180, 225, 270, 315]
)
```

---

### Example 2: Industrial accident with Threshold (AEGL) dose-response

This example models an industrial chemical release using AEGL threshold levels with the MaxConcentration calculator. Instead of cumulative dose, AEGL uses peak concentration over a sampling period.

```python
from hera import toolkitHome
from hera.utils import ureg

# ---------------------------------------------------------------
# Step 1: Set up toolkits
# ---------------------------------------------------------------
PROJECT = "IndustrialAccident"

lsm  = toolkitHome.getToolkit(toolkitHome.LSM, projectName=PROJECT)
risk = toolkitHome.getToolkit(toolkitHome.RISKASSESSMENT, projectName=PROJECT)
demo = toolkitHome.getToolkit(toolkitHome.GIS_DEMOGRAPHY, projectName=PROJECT)

# ---------------------------------------------------------------
# Step 2: Get concentration from LSM simulation
# ---------------------------------------------------------------
sim = lsm.getSimulations(simulationName="factory_release")[0]
concentration = sim.getConcentration(Q=100 * ureg.kg)
ground = concentration.sel(z=10, method="nearest")

# ---------------------------------------------------------------
# Step 3: Load an agent with Threshold effects
# ---------------------------------------------------------------
# This agent uses AEGL levels — binary thresholds based on peak
# concentration, not cumulative dose.
#
# Agent JSON looks like:
# {
#     "name": "HydrogenFluoride",
#     "effectParameters": {"tenbergeCoefficient": 1},
#     "effects": {
#         "AEGL10min": {
#             "type": "Threshold",
#             "calculator": {"MaxConcentration": {"sampling": "10min"}},
#             "parameters": {
#                 "type": "Threshold",
#                 "levels": ["AEGL-3", "AEGL-2", "AEGL-1"],
#                 "parameters": {
#                     "AEGL-3": {"threshold": "139*mg/m**3"},
#                     "AEGL-2": {"threshold": "34*mg/m**3"},
#                     "AEGL-1": {"threshold": "1.3*mg/m**3"}
#                 }
#             }
#         }
#     }
# }

agent = risk.getAgent("HydrogenFluoride")
print(f"Agent: {agent.name}")
print(f"Effects: {agent.effectNames}")  # ['AEGL10min']

# Access the AEGL effect
aegl = agent["AEGL10min"]
print(f"Calculator: {type(aegl.calculator).__name__}")  # CalculatorMaxConcentration
print(f"Levels: {aegl.levelNames}")  # ['AEGL-3', 'AEGL-2', 'AEGL-1']

# ---------------------------------------------------------------
# Step 4: Inspect threshold values
# ---------------------------------------------------------------
for level_name in aegl.levelNames:
    level = aegl[level_name]
    print(f"{level_name}: threshold = {level.threshold}")
    # AEGL-3: threshold = 139 mg/m³ (life-threatening)
    # AEGL-2: threshold = 34 mg/m³  (irreversible health effects)
    # AEGL-1: threshold = 1.3 mg/m³ (notable discomfort)

# ---------------------------------------------------------------
# Step 5: Calculate injury regions
# ---------------------------------------------------------------
# The MaxConcentration calculator finds the peak concentration
# at each location over the sampling period.
# The Threshold injury level then marks locations where
# peak concentration exceeds each AEGL threshold.

# The threshold levels must be applied to the toxic loads (the max rolling
# concentration), not to the raw concentration: calculateRegionOfInjured does both.
aegl_regions = aegl.calculateRegionOfInjured(ground, field="C", time="datetime", x="x", y="y")
# One row per level that is exceeded somewhere; percentEffected = 1.
# (None/empty if the concentration is below all thresholds)

aegl3_contours = aegl_regions.query("severity == 'AEGL-3'")
aegl2_contours = aegl_regions.query("severity == 'AEGL-2'")
aegl1_contours = aegl_regions.query("severity == 'AEGL-1'")

# ---------------------------------------------------------------
# Step 6: Estimate affected population
# ---------------------------------------------------------------
population = demo.getDataSourceData("census_2020")
release_location = (180000, 663000)  # ITM coordinates

# Project the AEGL zones onto the population (project() is defined on the full
# thresholdGeoDataFrame; slices made with query() are plain GeoDataFrames)
casualties = aegl_regions.project(
    demographic=population,
    loc=release_location,
    mathematical_angle=45
)
if casualties is None:
    print("no affected population")
else:
    print(casualties.groupby("severity")["effectedtotal_pop"].sum())

# ---------------------------------------------------------------
# Step 7: Compare with protection (masks)
# ---------------------------------------------------------------
from hera.riskassessment import ProtectionPolicy
from hera.utils import ureg

# Gas masks with protection factor 1000
policy = ProtectionPolicy()
policy.masks(protectionFactor=1000, wear="5min", duration="4h")

protected = policy.compute(ground, C="C")

# With masks, the effective concentration is 1000x lower (inside the wear window)
# so fewer locations exceed the AEGL thresholds
aegl_protected = aegl.calculateRegionOfInjured(protected, field="C")
aegl3_protected = aegl_protected.query("severity == 'AEGL-3'")
print(f"AEGL-3 area without masks: {aegl3_contours.area.sum():.0f} m²")
print(f"AEGL-3 area with masks:    {aegl3_protected.area.sum():.0f} m²")
```

### Key differences between the two examples

| Aspect | Example 1 (Lognormal10) | Example 2 (Threshold) |
|--------|------------------------|----------------------|
| **Agent** | Chlorine | Hydrogen Fluoride |
| **Calculator** | TenBerge (`integral(C^n dt)`) | MaxConcentration (peak C) |
| **Injury model** | Lognormal10 — gradual % affected | Threshold — binary 0%/100% |
| **Parameters** | TL_50, sigma per severity | threshold (with units) per AEGL level |
| **Result** | % of population affected at each severity | Zone where concentration exceeds standard |
| **Use case** | Toxic exposure over time | Acute exposure standards (AEGL, ERPG) |
| **Protection** | Indoor sheltering (air exchange) | Gas masks (protection factor) |
