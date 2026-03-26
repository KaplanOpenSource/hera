# LSM toolkit

## Architecture

| Component | Module | Purpose |
|-----------|--------|---------|
| Toolkit | `toolkit.py` | `LSMToolkit` — simulation management, data source API |
| Single simulation | `singleSimulation.py` | Handles one LSM run (concentration, particle trajectories) |
| Template | `template.py` | LSM case template management |
| Workflow | `hermesWorkflowToolkit.py` | Integration with Hermes workflow system |
| CLI | `CLI.py` | `hera-LSM` command-line interface |

## Key methods

| Method | Description |
|--------|-------------|
| `loadData(file, saveMode)` | Load LSM simulation data |
| `getSimulations(**params)` | Query stored simulations |
| `singleSimulation(file)` | Create a handler for one simulation |
| `getConcentration(Q)` | Calculate concentration field |
