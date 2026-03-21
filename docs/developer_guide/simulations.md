# Simulations Implementation

This page covers the internal architecture of the simulation toolkits for developers extending or maintaining them.

---

## Package structure

```
hera/simulations/
    openFoam/
        toolkit.py             # OFToolkit — main OpenFOAM toolkit
        eulerian/
            abstractEulerianSolver.py  # Base for Eulerian solvers
            simpleFoam.py      # simpleFoam solver interface
            buoyantReactingFoam.py     # buoyantReactingFoam solver
            NavierStokes.old/  # Legacy Navier-Stokes (preprocess, postprocess, GUI)
        lagrangian/
            LSM/
                toolkit.py     # OFLSMToolkit — OpenFOAM + LSM coupling
        preprocessOFObjects/   # Mesh and case preprocessing utilities
        postProcess/
            VTKPipeline.py     # VTK post-processing pipeline
            VTKPipelineExecutionContext.py
    LSM/
        toolkit.py             # LSMToolkit — Lagrangian Stochastic Model
        singleSimulation.py    # Single LSM simulation handler
        template.py            # LSM template management
        hermesWorkflowToolkit.py  # LSM workflow integration
        CLI.py                 # hera-LSM CLI entry points
    gaussian/
        toolkit.py             # gaussianToolkit — Gaussian dispersion
    windProfile/
        toolkit.py             # WindProfileToolkit — vertical wind profiles
    machineLearningDeepLearning/
        toolkit.py             # ML/DL toolkit
        torch/
            modelContainer.py  # PyTorch model container
    hermesWorkflowToolkit.py   # hermesWorkflowToolkit — workflow management
    evaporation/               # Evaporation models
    deposition/
        models.py              # Particle deposition models
    hydrodynamics/
        nearWallFlow.py        # Near-wall flow calculations
    WRF/
        wrfDatalayer.py        # WRF weather model data layer
    utils/
        interpolations.py      # Spatial interpolation
        coordinateHandler.py   # Coordinate transformations
        canopyWindProfile.py   # Canopy wind profile model
        inputForModelsCreation.py  # Model input generation
```

---

## OpenFOAM toolkit

### Architecture

The OpenFOAM toolkit (`OFToolkit`) manages the full CFD simulation lifecycle:

| Component | Module | Purpose |
|-----------|--------|---------|
| Main toolkit | `toolkit.py` | Template management, case setup, running, data source management |
| Eulerian solvers | `eulerian/` | Solver-specific interfaces (simpleFoam, buoyantReactingFoam) |
| Abstract solver | `abstractEulerianSolver.py` | Base class for Eulerian solver interfaces |
| Lagrangian coupling | `lagrangian/LSM/toolkit.py` | `OFLSMToolkit` — OpenFOAM + LSM particle tracking |
| Preprocessing | `preprocessOFObjects/` | Mesh generation, boundary condition setup |
| Post-processing | `postProcess/` | VTK pipeline for field extraction and visualization |

### Solver hierarchy

| Class | Inherits from | Solver |
|-------|--------------|--------|
| `OFToolkit` | `hermesWorkflowToolkit` | Main toolkit (not solver-specific) |
| `abstractEulerianSolver` | — | Base for all Eulerian solvers |
| `simpleFoam` | `abstractEulerianSolver` | Steady-state incompressible turbulent flow |
| `buoyantReactingFoam` | `abstractEulerianSolver` | Buoyancy-driven reacting flow |
| `OFLSMToolkit` | — | Coupled OpenFOAM + Lagrangian particle tracking |

### Template system

Templates are saved OpenFOAM case directories stored as data sources. The toolkit provides:
- `saveTemplate(templateName, caseDirectory)` — save a case as a template
- `loadTemplate(templateName, toDirectory, caseName)` — instantiate a template

---

## LSM toolkit

### Architecture

| Component | Module | Purpose |
|-----------|--------|---------|
| Toolkit | `toolkit.py` | `LSMToolkit` — simulation management, data source API |
| Single simulation | `singleSimulation.py` | Handles one LSM run (concentration, particle trajectories) |
| Template | `template.py` | LSM case template management |
| Workflow | `hermesWorkflowToolkit.py` | Integration with Hermes workflow system |
| CLI | `CLI.py` | `hera-LSM` command-line interface |

### Key methods

| Method | Description |
|--------|-------------|
| `loadData(file, saveMode)` | Load LSM simulation data |
| `getSimulations(**params)` | Query stored simulations |
| `singleSimulation(file)` | Create a handler for one simulation |
| `getConcentration(Q)` | Calculate concentration field |

---

## Gaussian dispersion toolkit

The `gaussianToolkit` provides Gaussian puff/plume models for atmospheric pollutant transport. It implements standard Gaussian equations for:

- Cloud dispersion
- Meteorology integration (wind field, stability classes)
- Downwind concentration estimation

---

## Wind profile toolkit

The `WindProfileToolkit` handles vertical wind profile modeling and analysis using standard atmospheric boundary layer theory.

---

## Hermes workflow toolkit

The `hermesWorkflowToolkit` is the base class for simulation toolkits that support workflow pipelines:

| Feature | Description |
|---------|-------------|
| Workflow groups | Organize simulations into named groups |
| Chained steps | Chain multiple simulation steps (preprocessing → solving → post-processing) |
| Template support | Use templates for reproducible case setup |

`OFToolkit` inherits from `hermesWorkflowToolkit`, gaining workflow support automatically.

---

## Machine Learning / Deep Learning toolkit

The ML/DL toolkit (`machineLearningDeepLearningToolkit`) provides:

| Component | Module | Purpose |
|-----------|--------|---------|
| Toolkit | `toolkit.py` | Model management as data sources |
| PyTorch container | `torch/modelContainer.py` | Model packaging, checkpoint loading, submodule management |

---

## Adding a new simulation toolkit

1. Create a new module under `hera/simulations/<domain>/`
2. Inherit from `abstractToolkit` (or `hermesWorkflowToolkit` for workflow support)
3. Set `toolkitName` in `__init__`
4. Add `_analysis` and `_presentation` layers as needed
5. Register in `ToolkitHome._toolkits` dict in `hera/toolkit.py`

For the full API, see the [Simulations API Reference](api/simulations.md).
