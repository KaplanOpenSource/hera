# Sensitivity Analysis

The `hera.utils.SALibUtils` module provides the `SALibUtils` class, which bridges Hera configurations with the [SALib](https://salib.readthedocs.io/) sensitivity analysis library. It handles building properly formatted problem definitions, supporting typed parameters (float, int, log-scale, and categorical), and transforming raw SALib samples back to their actual values.

## Importing

```python
from hera.utils.SALibUtils import SALibUtils
```

## Defining a Problem

Use `SALibUtils.buildSAProblem` to define the parameters and their types. Plain two-element lists are treated as float ranges. Use the `typeInt`, `typeList`, and `typeLog` helpers for other types.

```python
problem_container = SALibUtils.buildSAProblem(
    # Float parameter: uniform between 0.001 and 0.1
    viscosity=[0.001, 0.1],

    # Integer parameter: 1 through 10
    mesh_layers=SALibUtils.typeInt(1, 10),

    # Log-scale parameter: 10^1 through 10^5
    reynolds=SALibUtils.typeLog(1, 5),

    # Categorical parameter: pick from a list
    solver=SALibUtils.typeList(["simpleFoam", "pisoFoam", "pimpleFoam"]),
)
```

The returned `problem_container` is a dict with two keys:

- `"problem"` -- the SALib problem dict (`num_vars`, `names`, `bounds`), ready to pass to any SALib sampler.
- `"type"` -- metadata mapping each parameter name to its type descriptor, used by `transformSample`.

## Type Helpers

| Helper | Description | Example |
|---|---|---|
| `SALibUtils.typeInt(lower, upper)` | Integer range (inclusive) | `typeInt(1, 10)` |
| `SALibUtils.typeList(items)` | Categorical choice from a list | `typeList(["a", "b", "c"])` |
| `SALibUtils.typeLog(lower, upper)` | Log-scale: sampled value becomes `10**x` | `typeLog(-3, 0)` for 0.001--1.0 |
| Plain `[lower, upper]` list | Continuous float range | `[0.0, 1.0]` |

## Sampling

Pass the `"problem"` dict to any SALib sampling method:

```python
from SALib.sample import saltelli

samples = saltelli.sample(problem_container["problem"], N=256)
print(samples.shape)  # (N * (2D + 2), D) where D is number of parameters
```

## Transforming Samples

Raw SALib samples are always floats in the specified bounds. Use `transformSample` to convert them to their intended types (rounding integers, indexing into lists, applying log transforms):

```python
transformed = SALibUtils.transformSample(samples, problem_container)

for row in transformed[:3]:
    print(row)
# [0.042, 5, 1000.0, 'pisoFoam']
# [0.087, 2, 31.62, 'simpleFoam']
# ...
```

Each row is a list of parameter values in the same order as `problem_container["problem"]["names"]`.

## Full Workflow Example

```python
from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np
from hera.utils.SALibUtils import SALibUtils

# 1. Define the problem
problem_container = SALibUtils.buildSAProblem(
    amplitude=[0.5, 2.0],
    frequency=SALibUtils.typeInt(1, 10),
    phase=SALibUtils.typeList([0, 90, 180, 270]),
)

# 2. Generate samples
raw_samples = saltelli.sample(problem_container["problem"], N=512)

# 3. Transform to actual parameter values
param_sets = SALibUtils.transformSample(raw_samples, problem_container)

# 4. Run the model (replace with your actual model)
Y = np.array([
    row[0] * np.sin(row[1] * 2 * np.pi + np.radians(row[2]))
    for row in param_sets
])

# 5. Analyze sensitivity
Si = sobol.analyze(problem_container["problem"], Y)
print("First-order indices:", Si["S1"])
print("Total-order indices:", Si["ST"])
```
