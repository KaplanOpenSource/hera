# Tiles

**Toolkit name:** `GIS_Tiles`

The Tiles toolkit plots raster map images from a tile server (Google Maps, OpenStreetMap, or a custom server). This is useful for creating satellite or street-map backgrounds for your GIS visualizations.

**Data source format:**

| Property | Value |
|----------|-------|
| Data type | URL template string (`dataFormat: "string"`) |
| Format | XYZ tile URL with `{z}`, `{x}`, `{y}` placeholders |
| Example | `http://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}` |
| Tile format | Standard Web Mercator (XYZ) 256×256 PNG tiles |
| CRS | WGS84 (EPSG:4326) or ITM (EPSG:2039) for input coordinates |
| Config key | `defaultTileServer` |

```python
from hera import toolkitHome
import matplotlib.pyplot as plt

tiles = toolkitHome.getToolkit(toolkitHome.GIS_TILES, projectName="MY_PROJECT")
```

## Getting a map image (WGS84 coordinates)

Specify a bounding box using WGS84 (lat/lon) coordinates:

```python
region = dict(
    minx=34.775, maxx=34.8,
    miny=32.05,  maxy=32.1,
    zoomlevel=17,
    tileServer=None  # use the default server
)
img = tiles.getImageFromCorners(**region)

# Plot the image
fig, ax = plt.subplots(1, 1, figsize=(12, 12))
tiles.presentation.plot(img, ax=ax)
```

## Getting a map image (ITM coordinates)

You can also use Israeli Transverse Mercator (ITM) coordinates by specifying `inputCRS`:

```python
from hera.measurements.GIS.utils import ITM, WSG84

region = dict(
    minx=178898.481, maxx=181280.365,
    miny=661943.482, maxy=667478.9,
    zoomlevel=17,
    tileServer=None,
    inputCRS=ITM  # input coordinates are in ITM
)
img = tiles.getImageFromCorners(**region)

# Plot with ITM axes (default)
fig, ax = plt.subplots(1, 1, figsize=(12, 12))
tiles.presentation.plot(img, ax=ax)

# Or plot with WGS84 axes
fig, ax = plt.subplots(1, 1, figsize=(12, 12))
tiles.presentation.plot(img, outputCRS=WSG84, ax=ax)
```

## Tile servers

The default server is Google Maps satellite imagery. You can check or change it:

```python
# Check the default
tiles.getConfig()['defaultTileServer']
# 'http://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}'

# Set a new default (e.g., OpenStreetMap)
tiles.setDefaultTileServer("https://tile.openstreetmap.org/{z}/{x}/{y}.png")

# Or pass a server for a single call without changing the default
region = dict(
    minx=34.775, maxx=34.8,
    miny=32.05, maxy=32.1,
    zoomlevel=17,
    tileServer="https://tile.openstreetmap.org/{z}/{x}/{y}.png"
)
img = tiles.getImageFromCorners(**region)
```

If you have a tile server registered as a data source in your repository, you can refer to it by name:

```python
tiles.getDataSourceList()
# ['localTileServer']

region = dict(
    minx=178898.481, maxx=181280.365,
    miny=661943.482, maxy=667478.9,
    zoomlevel=17,
    tileServer="localTileServer",  # use the registered data source
    inputCRS=ITM
)
img = tiles.getImageFromCorners(**region)
```

> **Jupyter tutorial:** See the [Tile Toolkit notebook](../../../tutorials/TileToolkitDocumentation.ipynb) for a full interactive walkthrough with map output.

For the full API, see the [API Reference](../../developer_guide/api/measurements.md).
