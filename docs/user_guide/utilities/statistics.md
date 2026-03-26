# 2D Distributions

The `hera.utils.statistics` module provides `calcDist2d`, a function that computes a normalized 2D histogram suitable for contour plots. It wraps `matplotlib.pyplot.hist2d` and adds several normalization modes.

## Importing

```python
from hera.utils.statistics import calcDist2d
```

## Basic Usage

```python
import numpy as np
import matplotlib.pyplot as plt
from hera.utils.statistics import calcDist2d

# Generate sample data
rng = np.random.default_rng(42)
x = rng.normal(0, 1, 5000)
y = 0.5 * x + rng.normal(0, 0.5, 5000)

x_mid, y_mid, hist = calcDist2d(x, y, bins=30)

plt.contourf(x_mid, y_mid, hist, levels=20, cmap="viridis")
plt.colorbar(label="Normalized density")
plt.xlabel("X")
plt.ylabel("Y")
plt.title("2D Distribution")
plt.show()
```

## Using with a DataFrame

When your data lives in a pandas DataFrame, pass column names as strings and the DataFrame as `data`:

```python
import pandas as pd
from hera.utils.statistics import calcDist2d

df = pd.DataFrame({"wind_speed": x, "temperature": y})

x_mid, y_mid, hist = calcDist2d("wind_speed", "temperature", data=df, bins=25)
```

## Normalization Modes

The `normalization` parameter controls how the raw bin counts are scaled:

| Mode | Description |
|------|-------------|
| `"max_normalized"` | Divide all bins by the maximum count so the peak equals 1. This is the **default**. |
| `"density"` | Divide each bin by its area, producing a probability density (counts per unit area). |
| `"y_normalized"` | Normalize each column (fixed x-bin) so its values sum to 1. Useful for showing conditional distributions. |

```python
# Density normalization
x_mid, y_mid, hist = calcDist2d(x, y, bins=30, normalization="density")

# Column-normalized (conditional on x)
x_mid, y_mid, hist = calcDist2d(x, y, bins=30, normalization="y_normalized")
```

## Restricting the Axis Range

Use `x_range` and `y_range` to limit the histogram to a specific region:

```python
x_mid, y_mid, hist = calcDist2d(
    x, y,
    bins=40,
    x_range=(-2, 2),
    y_range=(-2, 2),
)
```

Both must be provided together as `(lower, upper)` tuples.

## Return Values

`calcDist2d` returns three arrays:

- `x_mid` -- 1D array of bin-center x-coordinates.
- `y_mid` -- 1D array of bin-center y-coordinates.
- `hist` -- 2D array (transposed) ready for `plt.contourf(x_mid, y_mid, hist)`.
