# Wind Profile

**Toolkit name:** `WindProfile` | **Constant:** `toolkitHome.WINDPROFILE`

The Wind Profile toolkit builds a spatial wind field over a region and extrapolates it vertically. It:

1. builds a regular lat/lon grid of the region with the land cover and the roughness length `z0` (via the [Land Cover toolkit](../measurements/gis/landcover.md)),
2. downloads the latest wind speed and direction of the IMS (Israel Meteorological Service) stations listed in `hera/simulations/windProfile/wind_stations.json`,
3. interpolates the station winds onto the grid, and
4. computes a logarithmic (or, in urban canopies, exponential) vertical profile of `u`/`v` at every grid point.

```python
from hera import toolkitHome

wp = toolkitHome.getToolkit(toolkitHome.WINDPROFILE, projectName="MY_PROJECT")
```

**Requirements:** a land-cover data source in the project (the `defaultLandCover` config key of the Land Cover toolkit, or pass `landcover_DataSource`), network access to `api.ims.gov.il`, and an IMS API token.

---

## Building the spatial wind field

### 1. Region grid and station data

```python
IMS_TOKEN = {"Authorization": "ApiToken <your-token>"}

landcover, stations = wp.getSpatialWind(
    minlat=32.0, minlon=34.7, maxlat=32.2, maxlon=34.9,
    IMS_TOKEN=IMS_TOKEN,
    dxdy=30,                      # grid spacing [m]
    inputCRS=WSG84,               # from hera.measurements.GIS import WSG84
    landcover_DataSource=None,    # None = the project's default land cover
)
# landcover : xarray Dataset over (i, j) with lat, lon, landcover and z0
# stations  : list of [lat, lon, height, [ws, wd]] for every station that reported
#             within 15 minutes of the most recent report
```

`getSpatialWind` queries the latest record of every station in `wind_stations.json` (up to 20 retries per station) and drops readings with `|ws| >= 20` or `|wd| > 360`.

### 2. Interpolating the stations onto the grid

The interpolation uses the terrain height, so add an `elevation` coordinate first (for example with the raster [Topography toolkit](../measurements/gis/topography_raster.md)):

```python
topo = toolkitHome.getToolkit(toolkitHome.GIS_RASTER_TOPOGRAPHY, projectName="MY_PROJECT")
landcover = topo.getElevationOfXarray(landcover)

landcover = wp.add_interpolated_ws_wd(landcover, stations)
# adds 'ws' (wind speed) and 'wd' (meteorological wind direction) at 10 m above the terrain
```

### 3. Vertical profiles

```python
profile = wp.getWindProfile(landcover, height=100, dz=10)
# pandas DataFrame, one row per grid point and height:
#   lat, lon, height, u, v, U_z, direction
```

At every grid point the friction velocity is derived from `ws` and `z0`, and the speed at height `z` is (κ = 0.41, β = 0.2, from `hera.measurements.GIS.utils`)

- `U(z) = (u*/κ) · ln(z / z0)` outside canopies, and
- inside an urban canopy (when the dataset has `hc` > 2 m, plus the mixing length `ll`): the log law above the canopy height `hc` and `U(hc) · exp(β (z − hc) / ll)` below it.

The direction is converted with `toMathematicalAngle` before splitting into `u` and `v`.

!!! note
    The heights start at `z = 0`, where `ln(z/z0)` is undefined, so the first row of every point is `inf`/`NaN`; filter it out (`profile.query("height > 0")`). Grid points with `hc` ≤ 2 m in a dataset that has an `hc` field are skipped.

---

## Presentation

```python
fig = wp.presentation.plotWindDirections(
    plot,                 # a matplotlib image (e.g. returned by imshow) in ITM coordinates
    landcover,            # the dataset with lat, lon and wd
    figsize=(28, 28), color="black", arrow_length=10,
    width=10, head_width=50, head_length=70,
    show_photo=False,     # True: plt.show(); False: return the figure
)
```

The background image is redrawn with the extent of `plot`, and an arrow is drawn at every grid point (converted to ITM) using `wd` directly as the arrow angle (no meteorological-to-mathematical conversion).

---

There is no CLI for this toolkit.

For implementation details, see the [Developer Guide](../../developer_guide/simulations/windprofile.md). For the full API, see the [API Reference](../../developer_guide/api/simulations.md).
