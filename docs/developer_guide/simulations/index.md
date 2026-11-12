# Simulations Implementation

This page covers the internal architecture of the simulation toolkits for developers extending or maintaining them.

---

## Package structure

```
hera/simulations/
    CLI.py                     # hera-workflows CLI handlers
    openFoam/
        toolkit.py             # OFToolkit — main OpenFOAM toolkit
        OFWorkflow.py          # Hermes workflow wrappers for OpenFOAM cases
        CLI.py                 # hera-openFoam CLI entry points
        eulerian/
            abstractEulerianSolver.py  # Base for Eulerian solvers
            simpleFoam.py      # simpleFoam solver interface
            buoyantReactingFoam.py     # buoyantReactingFoam solver
            NavierStokes.old/  # Legacy Navier-Stokes (preprocess, postprocess, GUI)
        lagrangian/
            abstractLagrangianSolver.py   # Base for Lagrangian solvers
            StochasticLagrangianSolver.py # StochasticLagrangianSolver toolkit extension
            LSM/
                toolkit.py     # OFLSMToolkit — OpenFOAM + LSM coupling
                sourcesFactoryTool.py  # Particle source factory
        preprocessOFObjects/   # OpenFOAM field/object writers (OFObject, OFField, datasetToOF)
        postProcess/
            VTKPipeline.py     # VTK post-processing pipeline
            pvOpenFOAMBase.py  # ParaView OpenFOAM reader base
        toberewritten/         # Legacy WRF -> OpenFOAM helpers (pending rewrite)
    LSM/
        toolkit.py             # LSMToolkit — Lagrangian Stochastic Model
        singleSimulation.py    # Single LSM simulation handler
        template.py            # LSM template management
        CLI.py                 # hera-LSM CLI entry points
    gaussian/
        toolkit.py             # gaussianToolkit — Gaussian dispersion
        gasCloud.py            # Gaussian gas cloud (puff/plume) model
        DropletCloud.py        # Droplet cloud model
        FallingNonEvaporatingDroplets.py  # Settling droplet model
        Meteorology.py         # Meteorology / stability-class wind profiles
        Sigma.py               # Dispersion coefficients (sigma)
        source.py              # Source definitions
        MeshUtils.py           # Spatial mesh helpers
    windProfile/
        toolkit.py             # WindProfileToolkit — vertical wind profiles
    machineLearningDeepLearning/
        toolkit.py             # ML/DL toolkit
        LightingModule.py      # PyTorch Lightning module helpers
        torch/
            modelContainer.py  # PyTorch model container
        dataanalisys/          # Data-analysis helpers (ml.py, stat.py)
    hermesWorkflowToolkit.py   # hermesWorkflowToolkit — workflow management
    analysis/
        errorCalculation.py    # Model-vs-observation error metrics
    evaporation/               # Evaporation models (models.py, monaghan.py)
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
3. Pass `toolkitName` to `super().__init__(...)`
4. Add `_analysis` and `_presentation` layers as needed
5. Register it:
    - **External / project toolkits** — register dynamically with `toolkitHome.registerToolkit(toolkit_name=..., toolkit_path=...)` (or `toolkitHome.import_toolkits_from_json(...)`); `getToolkit` resolves these after the built-ins.
    - **Built-in toolkits** — the static registry is the `self._toolkits` dict built in `ToolkitHome.__init__` (`hera/toolkit.py`), with a matching name constant on `ToolkitHome` (e.g. `SIMULATIONS_OPENFOAM = "OpenFOAM"`). Prefer dynamic registration over editing this dict.

For the full API, see the [Simulations API Reference](../api/simulations.md).
