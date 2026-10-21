# DataFrame Comparison

The `hera.utils.dataframeutils` module provides `compareDataframeConfigurations`, a function that identifies parameters that differ across multiple named datasets. It supports several input formats and can return results in either long or wide (pivoted) format.

## Importing

```python
from hera.utils.dataframeutils import compareDataframeConfigurations
```

## Basic Example -- List of Tuples

The simplest input is a list of `(name, dataframe)` tuples. Each DataFrame must have a parameter-name column and a value column:

```python
import pandas as pd
from hera.utils.dataframeutils import compareDataframeConfigurations

df_a = pd.DataFrame({
    "parameterName": ["solver.dt", "solver.maxIter", "mesh.resolution"],
    "value": [0.1, 500, 100],
})

df_b = pd.DataFrame({
    "parameterName": ["solver.dt", "solver.maxIter", "mesh.resolution"],
    "value": [0.2, 500, 200],
})

diff = compareDataframeConfigurations(
    [("baseline", df_a), ("refined", df_b)],
    parameterName="parameterName",
    valueName="value",
)
print(diff)
```

The result is a wide-format DataFrame showing only the rows where the value differs between `baseline` and `refined` (here `solver.dt` and `mesh.resolution`). `solver.maxIter` is excluded because it is the same in both. A parameter that is present in only some of the datasets is also reported as a difference.

## Wide vs. Long Format

By default the output is wide (pivoted), with one column per dataset:

```
datasetName      baseline  refined
parameterName
mesh.resolution     100.0    200.0
solver.dt             0.1      0.2
```

Pass `longFormat=True` to get long format instead:

```python
diff = compareDataframeConfigurations(
    [("baseline", df_a), ("refined", df_b)],
    parameterName="parameterName",
    valueName="value",
    longFormat=True,
)
```

```
     parameterName  value datasetName
2  mesh.resolution  100.0    baseline
2  mesh.resolution  200.0     refined
0        solver.dt    0.1    baseline
0        solver.dt    0.2     refined
```

In long format the rows keep their original index from the input DataFrames, and the `datasetName` column is appended.

## Input as a List of Dicts

Each item can also be a dict with the dataset name under the `datasetName` key and the DataFrame under `"data"`:

```python
diff = compareDataframeConfigurations(
    [{"datasetName": "baseline", "data": df_a},
     {"datasetName": "refined", "data": df_b}],
)
```

## Input from a Single DataFrame

If your data is already in a single long-format DataFrame with a column distinguishing datasets, pass it directly:

```python
combined = pd.DataFrame({
    "run": ["A", "A", "B", "B"],
    "parameterName": ["dt", "maxIter", "dt", "maxIter"],
    "value": [0.01, 100, 0.02, 100],
})

diff = compareDataframeConfigurations(
    combined,
    datasetName="run",
    parameterName="parameterName",
    valueName="value",
)
```

## Making Columns Compatible with pandas.query

Dot-separated parameter names (e.g., `solver.dt`) cannot be used in `pandas.query()`. Set `changeDotToUnderscore=True` to convert dots to underscores in the parameter-name index of the wide output (it does not work with `longFormat=True`, where the index is the original integer index):

```python
diff = compareDataframeConfigurations(
    [("a", df_a), ("b", df_b)],
    parameterName="parameterName",
    valueName="value",
    changeDotToUnderscore=True,
)
# Index labels like "solver.dt" become "solver_dt"
```

## Parameters Reference

| Parameter | Type | Default | Description |
|---|---|---|---|
| `data` | list or DataFrame | -- | Datasets to compare: list of `(name, df)` tuples, list of `{datasetName: name, "data": df}` dicts, or one long-format DataFrame. |
| `datasetName` | str | `"datasetName"` | Column (or dict key) identifying each dataset. For tuple input it is the name of the column added to hold the tuple names. |
| `parameterName` | str | `"parameterName"` | Column containing parameter paths. |
| `valueName` | str | `"value"` | Column containing parameter values. |
| `indexList` | list | `None` | Additional grouping columns (sub-indices). |
| `longFormat` | bool | `False` | Return long format instead of wide. |
| `changeDotToUnderscore` | bool | `False` | Replace dots with underscores in parameter names. |

If no parameter differs, an empty `pandas.DataFrame()` is returned. Passing an empty list currently raises `AttributeError` (the code calls `pandas.Dataframe`).
