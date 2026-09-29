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

## Toolkit sub-pages

- [OpenFOAM toolkit](openfoam.md) — Overall design, solvers, workflows, VTK pipeline, templates
- [LSM toolkit](lsm.md) — Lagrangian Stochastic Model architecture
- [Gaussian dispersion toolkit](gaussian.md) — Gaussian puff/plume models
- [Wind profile toolkit](windprofile.md) — Vertical wind profile modeling
- [Hermes workflow toolkit](workflows.md) — Workflow pipeline base class
- [Machine Learning / Deep Learning toolkit](ml.md) — ML/DL model management

---

## Adding a new simulation toolkit

1. Create a new module under `hera/simulations/<domain>/`
2. Inherit from `abstractToolkit` (or `hermesWorkflowToolkit` for workflow support)
3. Set `toolkitName` in `__init__`
4. Add `_analysis` and `_presentation` layers as needed
5. Register in `ToolkitHome._toolkits` dict in `hera/toolkit.py`

For the full API, see the [Simulations API Reference](api/simulations.md).
