# Risk Assessment Implementation

This page covers the internal architecture of the risk assessment toolkit for developers extending or maintaining it.

---

## Package structure

```
hera/riskassessment/
    riskToolkit.py             # RiskToolkit + analysis inner class
    agents/
        Agents.py              # Agent + PhysicalPropeties classes
        effects/
            Injury.py          # InjuryFactory + Injury base + concrete injuries
            InjuryLevel.py     # InjuryLevel base + Lognormal10, Threshold, Exponential
            Calculator.py      # Dose calculators (Haber, TenBerge, MaxConcentration)
            thresholdGeoDataFrame.py  # Spatial threshold results
    protectionpolicy/
        ProtectionPolicy.py    # ProtectionPolicy + action pipeline (Indoor, Masks)
    analysis/
        riskAreas.py           # Risk area algorithms (Sweep)
    presentation/
        casualtiesFigs.py      # Casualty plots (roses, bar charts)
    CLI.py                     # hera-riskassessment CLI
    docs/                      # Internal documentation and helper scripts
```

---

## Core classes

### RiskToolkit

The main entry point. Inherits from `abstractToolkit` and manages agents as versioned data sources.

| Component | Type | Purpose |
|-----------|------|---------|
| `RiskToolkit` | Toolkit | Agent loading, listing, data source management |
| `analysis` | Inner class | Risk area calculation, LSM integration |
| `presentation` | `casualtiesPlot` | Casualty roses, bar plots |
| `ProtectionPolicy` | Class reference | Protection action pipeline builder |

### Agent

An `Agent` represents a hazardous material with its effects and physical properties.

| Component | Class | Purpose |
|-----------|-------|---------|
| `Agent` | Container | Holds effects, physical properties, Ten Berge coefficient |
| `PhysicalPropeties` | Properties | Molecular weight, density, volatility, vapor pressure, sorption |
| Effects (dict) | `Injury` instances | Accessed by name: `agent["inhalation"]` |

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

The effects system follows a layered architecture:

### Calculators

Calculators compute the toxic load (dose) from a concentration field:

| Calculator | Class | Formula |
|-----------|-------|---------|
| Haber | `CalculatorHaber` | Simple time integral of concentration |
| Ten Berge | `CalculatorTenBerge` | `integral(C^n * dt)` where n is the Ten Berge exponent |
| Max Concentration | `CalculatorMaxConcentration` | Peak concentration value |

### Injury levels

Injury levels map toxic loads to percent affected:

| Level type | Class | Dose-response model |
|-----------|-------|-------------------|
| Lognormal10 | `InjuryLevelLognormal10` | Log-normal CDF (base 10) with `TL_50` and `sigma` |
| Threshold | `InjuryLevelThreshold` | Binary: 0% below threshold, 100% above |
| Exponential | `InjuryLevelExponential` | Exponential dose-response with parameter `k` |

### Injuries

An `Injury` combines a calculator with one or more injury levels:

| Class | Description |
|-------|-------------|
| `Injury` | Base class holding calculator + levels dict |
| `InjuryLognormal10` | Injury using Lognormal10 dose-response |
| `InjuryThreshold` | Injury using threshold dose-response |
| `InjuryExponential` | Injury using exponential dose-response |

### Factory pattern

`InjuryFactory.getInjury(name, cfgJSON)` dynamically resolves the calculator and injury classes from the JSON descriptor using `pydoc.locate`.

---

## Protection policies

`ProtectionPolicy` implements a chainable action pipeline that modifies concentration fields:

```python
policy = ProtectionPolicy()
policy.indoor(turnover=2*ureg.hour, enter="5min", stay="2h")
policy.masks(protectionFactor=1000, wear="0min", duration="3h")
result = policy.compute(concentration_data, C="C")
```

### Action classes

| Action | Class | Parameters |
|--------|-------|-----------|
| Indoor sheltering | `ActionIndoor` | `alpha` or `turnover`, time window (`begin/end` or `enter/stay`) |
| Masking | `ActionMasks` | `protectionFactor`, time window (`begin/end` or `wear/duration`) |

Each action:
1. Reads the current concentration field (`policy.finalname`)
2. Applies its transformation (indoor air exchange model, or division by protection factor)
3. Writes a new field and updates `finalname` for the next action
4. Records its parameters in `data.attrs` for traceability

### Indoor model

The indoor concentration `Cin` is computed iteratively:

```
Cin[t] = (Cin[t-1] + alpha * dt * Cout[t]) / (1 + alpha * dt)
```

Where `alpha = 1/turnover` is the air exchange rate.

---

## Risk area calculation

The `riskAreaAlgorithm_Sweep` class computes casualty counts across a spatial grid:

1. **Build bounding box** — determine the search area from the demographic region + effect polygon extent
2. **Grid the area** — create a mesh of points at `dxdy` spacing
3. **For each grid point** — project the effect isopleths onto the demographic data at that release location and wind direction
4. **Aggregate** — sum affected population by severity level

Supports parallel execution via `multiprocessing.Pool`.

---

## Spatial results: thresholdGeoDataFrame

`thresholdGeoDataFrame` extends `geopandas.GeoDataFrame` with methods for:

| Method | Purpose |
|--------|---------|
| `shiftLocationAndAngle(loc, angle)` | Translate and rotate threshold polygons to a release location |
| `project(demographic, loc, angle)` | Project threshold polygons onto demographic data to estimate casualties |

---

## Adding new effects

1. **New calculator:** Create a class in `agents/effects/Calculator.py` named `Calculator<Name>` inheriting from `AbstractCalculator`
2. **New injury level:** Create a class in `agents/effects/InjuryLevel.py` named `InjuryLevel<Name>` inheriting from `InjuryLevel`
3. **New injury type:** Create a class in `agents/effects/Injury.py` named `Injury<Name>` inheriting from `Injury`

The factory uses `pydoc.locate` to find classes by name, so naming conventions must be followed exactly.

## Adding new protection actions

1. Create a class in `protectionpolicy/ProtectionPolicy.py` named `Action<Name>` inheriting from `abstractAction`
2. Implement `compute()` to modify `self.policy.data`
3. Add a convenience method on `ProtectionPolicy` (e.g., `def evacuation(self, **kwargs)`)

For the full API, see the [Risk Assessment API Reference](api/risk.md).
