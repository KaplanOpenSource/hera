# Contour to GIS

The `hera.utils.matplotlibCountour` module provides `toGeopandas`, a function that converts matplotlib contour output into a GeoDataFrame. This makes it straightforward to export contour plots as GeoJSON, shapefiles, or any other format supported by geopandas.

## Requirements

This module needs `shapely`, `geopandas` and `matplotlib` (all installed with Hera) and the `unum` package (pinned in `requirements.txt`): `toGeopandas` imports `unumToPint` from `hera.utils.unitHandler`, which only exists when `unum` is installed.

## Importing

```python
from hera.utils.matplotlibCountour import toGeopandas
# or, via the lazy re-export
from hera.utils import toGeopandas
```

## Basic Usage

Generate a contour plot with matplotlib and pass the result to `toGeopandas`:

```python
import numpy as np
import matplotlib.pyplot as plt
from hera.utils.matplotlibCountour import toGeopandas

# Create sample data on a grid
x = np.linspace(0, 1000, 100)
y = np.linspace(0, 1000, 100)
X, Y = np.meshgrid(x, y)
Z = np.exp(-((X - 500)**2 + (Y - 500)**2) / (2 * 200**2))

# Generate contours
fig, ax = plt.subplots()
cs = ax.contour(X, Y, Z, levels=10)
plt.close(fig)

# Convert to GeoDataFrame
gdf = toGeopandas(cs)
print(gdf.head())
```

The resulting GeoDataFrame has two columns:

- `Level` -- the contour level value.
- `contour` -- the polygon geometry (set as the active geometry column).

!!! warning "Several separate contours at the same level"
    With matplotlib 3.8+ (the pinned version is 3.9.4), `ContourSet.get_paths()` returns one path per level. `toGeopandas` uses the first ring of that path as the polygon shell and every other ring as a hole. If a level has several disjoint contours (for example two separate peaks), you get one invalid polygon (`gdf.is_valid` is `False`) instead of one polygon per contour. Check `gdf.is_valid` before using the result.

## Saving as GeoJSON

```python
gdf.to_file("contours.geojson", driver="GeoJSON")
```

## Saving as Shapefile

```python
gdf.to_file("contours.shp")
```

## Unit Conversion

If your contour coordinates are not in meters, use the `inunits` parameter to specify the input units. The coordinates are multiplied by the conversion factor to meters:

```python
from hera.utils.unitHandler import ureg
from hera.utils.matplotlibCountour import toGeopandas

# Contour data is in kilometers
gdf = toGeopandas(cs, inunits=1 * ureg.km)
```

The default (`inunits=None`) means meters (no conversion). Only the coordinates are scaled; the `Level` values are not changed. The resulting GeoDataFrame has no CRS set; use `gdf.set_crs(...)` before exporting if needed.

## Working with contourf

The same function works with filled contours (`contourf`):

```python
fig, ax = plt.subplots()
cs = ax.contourf(X, Y, Z, levels=10)
plt.close(fig)

gdf = toGeopandas(cs)
```

## Example: Overlay on a Map

Combine with contextily or folium for quick map visualization:

```python
import matplotlib.pyplot as plt

gdf.plot(column="Level", legend=True, cmap="RdYlGn_r")
plt.title("Concentration Contours")
plt.xlabel("Easting (m)")
plt.ylabel("Northing (m)")
plt.show()
```
