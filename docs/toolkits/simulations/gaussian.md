# Gaussian Dispersion

**Toolkit name:** `GaussianDispersion` | **Constant:** `toolkitHome.GAUSSIANDISPERSION`

The Gaussian dispersion toolkit computes the concentration and dosage of a gas cloud released from a single source, using Gaussian puff equations with Briggs (rural) dispersion coefficients. It distinguishes **instantaneous** releases (a mass, e.g. `1 kg`) from **continuous** releases (a rate, e.g. `1 g/s` — currently broken, see below), and handles ground reflections, and an optional depletion factor for dry deposition. Results are `xarray.DataArray` objects on a regular `(time, x, y, z)` grid, with the wind blowing along the `x` axis from a source at `x = 0`, `y = 0`.

The toolkit does not store anything in the database: it is a set of calculation helpers (meteorology, sigmas, gas cloud) plus a plotting layer.

```python
from hera import toolkitHome

gauss = toolkitHome.getToolkit(toolkitHome.GAUSSIANDISPERSION, projectName="MY_PROJECT")
```

All dimensional inputs are pint quantities (`hera.utils.ureg`); plain numbers are interpreted in the default units noted below.

---

## Workflow

A calculation has four steps:

1. **Meteorology** — wind speed at a reference height, stability class, roughness.
2. **Space-time grid** — the extents and resolution of the output grid.
3. **Gas cloud** — the source (amount, height, initial size) bound to the meteorology and the grid.
4. **Concentration / dosage** — computed from the gas cloud.

### 1. Meteorology

```python
from hera.utils import ureg

# Wind speed given at 10 m
met = gauss.getMeteorologyFromU10(
    u10=2 * ureg.m / ureg.s,
    inversion=300 * ureg.m,          # inversion (mixing) height
    verticalProfileType="powerLaw",  # "powerLaw" (default), "log" or "uniformWind"
    stability="D",                   # Pasquill class A-F
    z0=0.1 * ureg.m,                 # roughness length
    temperature=ureg.Quantity(20, ureg.degC),
)

# Wind speed given at another reference height
met = gauss.getMeteorologyFromURefHeight(u=3 * ureg.m / ureg.s, refHeight=20 * ureg.m,
                                         inversion=300 * ureg.m, stability="D")

met.getWindVelocity(height=10 * ureg.m)   # wind speed at a height (profile clipped to 0-300 m)
```

The remaining keyword arguments are `ustar` (default `0.3 m/s`) and `skinSurfaceTemperature` (default `35 °C`).

### 2. Space-time grid

`getSpaceTime` builds the grid dictionary used by the gas cloud. It chooses the horizontal cell size as a multiple of the distance the cloud travels in one time step (`dxdy = dt · u · dxdy_multiplier`), and a time span long enough for the cloud to pass `maxx`:

```python
initialCloudSize = [1 * ureg.m, 1 * ureg.m, 1 * ureg.m]   # sigma_x, sigma_y, sigma_z at the source

spaceTime = gauss.getSpaceTime(
    meteorology=met,
    sourceHeight=10 * ureg.m,
    wind_profile_type="default",   # "default" (the meteorology profile) or "HotSpot"
    maxx=1000 * ureg.m,
    maxz=50 * ureg.m,
    dt=10 * ureg.s,
    dz=1 * ureg.m,
    dxdy_multiplier=1,
    minimal_maxy=200 * ureg.m,     # the y extent is rounded up to a multiple of dxdy
    initialCloudSize=initialCloudSize,
)
# {'minx', 'maxx', 'miny', 'maxy', 'minz', 'maxz', 'dxdy', 'dz', 'timeSpan', 'dt'}
```

You can also build the dictionary yourself with the same keys (all pint quantities).

### 3. Gas cloud

```python
cloud = gauss.getGasCloud(
    sourceQ=1 * ureg.kg,               # [mass] -> instantaneous, [mass/time] -> continuous
    sourceHeight=10 * ureg.m,
    initialCloudSize=initialCloudSize,
    meteorology=met,
    wind_profile_type="default",
    spaceTime=spaceTime,
    deposition_velocity=0 * ureg.m / ureg.s,   # used only when DF=True
    sigmaTypeName="briggsRural",               # see gauss.listSigmaTypes()
)
```

`getGasCloud` returns an `instantaneousReleaseGasCloud` or a `continuousReleaseGasCloud` (from `hera.simulations.gaussian.gasCloud`), depending on the units of `sourceQ`. A `ValueError` is raised if `sourceQ` is neither a mass nor a mass rate.

### 4. Concentration and dosage

Instantaneous release:

```python
C = cloud.getConcentration_inst(numOfReflections=3)   # mg/m³, dims (time, x, y, z)
D = cloud.getDosage_inst(numOfReflections=3)          # mg·min/m³ (time-integrated, erf form)
D2 = cloud.getDosage_inst_NoERF()                     # dosage by numerical integration of C

# Normalised (per unit mass) versions
C1 = cloud.getConcentration_inst_noQ()                # 1/m³
```

Continuous release:

```python
cloud = gauss.getGasCloud(sourceQ=1 * ureg.g / ureg.s, ...)
C = cloud.getConcentration_cont()        # concentration of a constant-rate release
D = cloud.getDosage_cont_NoERF()
```

!!! bug "Known issue"
    In master every `continuousReleaseGasCloud` method (`getConcentration_cont`, `getConcentration_cont_NoERF`, `getDosage_cont_NoERF`, `getDosage_cont_doubleNoERF`) raises `AttributeError`: they call the `getDosage_inst_*_noQ` methods, which are defined on `instantaneousReleaseGasCloud`, not on the common base class (`hera/simulations/gaussian/gasCloud.py:596`). Only instantaneous releases work.

The time coordinate is in minutes from the release; `x`, `y`, `z` are in metres. The units are stored in `C.attrs['Q']`. Pass `DF=True` to the dosage functions to apply the dry-deposition depletion factor computed from `deposition_velocity`.

For radiological releases the instantaneous cloud also provides time-integrated activity concentrations (`getTIAC_inst(specifitActivity, ...)`, `getTIAC_inst_NoERF(...)` and their `_noQ` variants, output in `Bq·s/m³`).

---

## Sigmas

```python
gauss.listSigmaTypes()                 # ['briggsRural']
sigma = gauss.getSigmaType("briggsRural")
sigma.getSigma(x=[100, 500, 1000], stability="D", sigma0=initialCloudSize)
# {'sigmaX': ..., 'sigmaY': ..., 'sigmaZ': ..., 'distance': ...}
```

---

## Presentation

`gauss.presentation` plots the results with matplotlib (each call shows the figure):

| Method | Plot |
|--------|------|
| `plotMaxConcentrationOverTime(C, y, z, x_min=None, x_max=None)` | Maximum concentration along x at (y, z) |
| `plotFixedPointConcentrationOverTime(C, x, y, z, t_min=None, t_max=None)` | Concentration time series at a point |
| `plotDosagePerDistance(D, y, z, time)` | Dosage along x at a given time |
| `plotFixedPointDosageOverTime(D, x, y, z)` | Dosage time series at a point |

The `_noQ` variants (`plotMaxConcentrationOverTime_noQ`, `plotFixedPointConcentrationOverTime_noQ`, `plotTIACPerDistance_noQ`, `plotFixedPointTIACOverTime_noQ`) plot the normalised results.

---

## Complete example

```python
from hera import toolkitHome
from hera.utils import ureg

gauss = toolkitHome.getToolkit(toolkitHome.GAUSSIANDISPERSION, projectName="MY_PROJECT")

met = gauss.getMeteorologyFromU10(u10=3 * ureg.m / ureg.s, inversion=500 * ureg.m, stability="D")
sigma0 = [1 * ureg.m, 1 * ureg.m, 1 * ureg.m]

spaceTime = gauss.getSpaceTime(meteorology=met, sourceHeight=2 * ureg.m, wind_profile_type="default",
                               maxx=2000 * ureg.m, maxz=30 * ureg.m, dt=10 * ureg.s, dz=1 * ureg.m,
                               dxdy_multiplier=1, minimal_maxy=300 * ureg.m, initialCloudSize=sigma0)

cloud = gauss.getGasCloud(sourceQ=10 * ureg.kg, sourceHeight=2 * ureg.m, initialCloudSize=sigma0,
                          meteorology=met, wind_profile_type="default", spaceTime=spaceTime,
                          deposition_velocity=0 * ureg.m / ureg.s)

C = cloud.getConcentration_inst()
D = cloud.getDosage_inst()

gauss.presentation.plotMaxConcentrationOverTime(C, y=0 * ureg.m, z=1.5 * ureg.m)
gauss.presentation.plotDosagePerDistance(D, y=0 * ureg.m, z=1.5 * ureg.m, time=D.time[-1].item() * ureg.min)
```

There is no CLI for this toolkit.

For implementation details, see the [Developer Guide](../../developer_guide/simulations/gaussian.md). For the full API, see the [API Reference](../../developer_guide/api/simulations.md).
