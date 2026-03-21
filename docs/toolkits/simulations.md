# Simulations Toolkits

Simulation toolkits manage the lifecycle of computational models — setting up cases, running simulations, and post-processing results.

All simulation toolkits are accessed via `toolkitHome.getToolkit()` and bound to a project:

```python
from hera import toolkitHome

toolkit = toolkitHome.getToolkit(toolkitHome.SIMULATIONS_OPENFOAM, projectName="MY_PROJECT")
```

---

## OpenFOAM

**Toolkit name:** `OpenFOAM`

Full lifecycle management for OpenFOAM CFD simulations: templates, case setup, mesh generation, running, and post-processing.

### Templates

Templates are saved OpenFOAM case configurations that can be reused across simulations.

```python
of = toolkitHome.getToolkit(toolkitHome.SIMULATIONS_OPENFOAM, projectName="MY_PROJECT")

# Save a case directory as a template
of.saveTemplate(templateName="simpleFoam_base", caseDirectory="/path/to/case")

# Load a template to a new location
of.loadTemplate(
    templateName="simpleFoam_base",
    toDirectory="/data/simulations",
    caseName="run_001"
)
```

### Running simulations

```python
# Run an OpenFOAM case
of.runCase(casePath="/data/simulations/run_001", nproc=4)
```

### Post-processing

```python
# Access the analysis layer for data extraction
of.analysis.extractData(casePath="/data/simulations/run_001")
```

From the CLI:

```bash
# List templates
hera-openFoam simpleFoam templates list --projectName MY_PROJECT

# Save a template
hera-openFoam simpleFoam templates save myTemplate \
    --projectName MY_PROJECT --directory /path/to/case

# Load a template
hera-openFoam simpleFoam templates load myTemplate \
    --projectName MY_PROJECT --toDirectory /path/to/output --caseName run_001
```

---

## LSM (Lagrangian Stochastic Model)

**Toolkit name:** `LSM`

Atmospheric dispersion simulations using Lagrangian particle tracking.

```python
lsm = toolkitHome.getToolkit(toolkitHome.LSM, projectName="MY_PROJECT")

# Load model data
lsm.loadData(
    fileNameOrData="/data/lsm/particles.nc",
    saveMode=lsm.TOOLKIT_SAVEMODE_FILEANDDB
)

# Access analysis
lsm.analysis  # LSM-specific analysis methods

# Access presentation
lsm.presentation  # concentration maps, particle trajectories
```

From the CLI:

```bash
# Load an LSM model
hera-LSM load MY_PROJECT /path/to/resource modelFolder "0,0,1" params.json

# List LSM simulations
hera-LSM list MY_PROJECT
```

---

## Gaussian Dispersion

**Toolkit name:** `GaussianDispersion`

Gaussian puff and plume models for atmospheric pollutant transport.

```python
gauss = toolkitHome.getToolkit(toolkitHome.GAUSSIANDISPERSION, projectName="MY_PROJECT")

# The Gaussian toolkit provides methods for:
# - Cloud dispersion modeling
# - Meteorology integration (wind field, stability class)
# - Downwind concentration calculations
```

---

## Wind Profile

**Toolkit name:** `WindProfile`

Vertical wind profile modeling and analysis.

```python
wp = toolkitHome.getToolkit(toolkitHome.WINDPROFILE, projectName="MY_PROJECT")

# Wind profile analysis and modeling
```

---

## Hermes Workflows

**Toolkit name:** `hermesWorkflows`

Manages simulation groups and workflow pipelines that chain multiple simulation steps together.

```bash
# List workflow groups
hera-workflows list groups --projectName MY_PROJECT
```

For the full API details of each toolkit, see the [Toolkit Catalog](overview.md) and the [API Reference](../developer_guide/api/simulations.md).
