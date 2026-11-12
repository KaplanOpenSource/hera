# Tiles

**Toolkit name:** `GIS_Tiles`

**Class:** `hera.measurements.GIS.raster.tiles.TilesToolkit` (internal toolkit name `Tiles`)

The Tiles toolkit plots raster map images from a tile server (Google Maps, OpenStreetMap, or a custom server). This is useful for creating satellite or street-map backgrounds for your GIS visualizations.

**Data source format:**

| Property | Value |
|----------|-------|
| Data type | URL template string (`dataFormat: "string"`) |
| Format | XYZ tile URL with `{z}`, `{x}`, `{y}` placeholders |
| Example | `http://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}` |
| Tile format | Standard Web Mercator (XYZ) 256×256 PNG tiles |
| CRS | WGS84 (EPSG:4326) or ITM (EPSG:2039) for input coordinates |
| Config key | `defaultTileServer` (set with `setDefaultTileServer`, or in the repository `Config` section) |

```python
from hera import toolkitHome
import matplotlib.pyplot as plt

tiles = toolkitHome.getToolkit(toolkitHome.GIS_TILES, projectName="MY_PROJECT")
```

## Getting a map image (WGS84 coordinates)

Specify a bounding box using WGS84 coordinates (`minx`/`maxx` = longitude, `miny`/`maxy` = latitude).
`getImageFromCorners` returns a tuple `(image, extent)`, where `extent` is `[left, right, bottom, top]` in `outputCRS`
(default WGS84). `presentation.plot` accepts that tuple directly:

```python
region = dict(
    minx=34.775, maxx=34.8,
    miny=32.05,  maxy=32.1,
    zoomlevel=17,
    tileServer=None  # use the default server
)
img = tiles.getImageFromCorners(**region)   # (PIL image, [left, right, bottom, top])

# Plot the image (axes in ITM by default)
fig, ax = plt.subplots(1, 1, figsize=(12, 12))
tiles.presentation.plot(img, ax=ax)
```

## Getting a map image (ITM coordinates)

You can also use Israeli Transverse Mercator (ITM) coordinates by specifying `inputCRS`:

```python
from hera.measurements.GIS.utils import ITM

region = dict(
    minx=178898.481, maxx=181280.365,
    miny=661943.482, maxy=667478.9,
    zoomlevel=17,
    tileServer=None,
    inputCRS=ITM  # input coordinates are in ITM
)
img = tiles.getImageFromCorners(**region)

# The returned extent is in WGS84 (the default outputCRS), so plot with the default inputCRS
fig, ax = plt.subplots(1, 1, figsize=(12, 12))
tiles.presentation.plot(img, ax=ax)   # ITM axes (default outputCRS)
```

`presentation.plot(imageNameOrData, inputCRS=WSG84, outputCRS=ITM, extents=None, ax=None, display=True)`
takes the `(image, extent)` tuple, or an image plus `extents` (list `[minX, maxX, minY, maxY]` or a dict).
`inputCRS` must match the CRS of the extent. Only `outputCRS=ITM` is reliable in master: for any
other output CRS the code swaps the x/y extents (`tiles.py:356-359`). Passing an image name (string)
is not supported in master (it calls a missing `getImage` method).

## Tile servers

There is no built-in default server: if `tileServer` is `None` and the project has no `defaultTileServer` config,
`getImageFromCorners` raises `ValueError`. Repositories commonly set it to Google satellite imagery. You can check or change it:

```python
# Check the default (None if not set)
tiles.getConfig().get('defaultTileServer')
# e.g. 'http://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}'

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

## Other methods

| Method | Description |
|--------|-------------|
| `tileScaleAtLatLonZoom(latitude, longitude, zoomlevel)` | Size of one tile pixel in meters at the location and zoom |
| `deg2tile(lat_deg, lon_deg, zoom)` | WGS84 lat/lon to XYZ tile indices `(xtile, ytile)` |
| `tile2deg(xtile, ytile, zoom)` | Tile indices to WGS84 `(lat, lon)` of the tile's upper-left corner |
| `setDefaultTileServer(server)` | Store the default tile server URL in the project config |
| `listImages(**filters)` | List stored tile image documents (type `Tiles_PNG`) |

> **Jupyter tutorial:** See the [Tile Toolkit notebook](../../../tutorials/TileToolkitDocumentation.ipynb) for a full interactive walkthrough with map output.

For the full API, see the [API Reference](../../../developer_guide/api/measurements.md).
