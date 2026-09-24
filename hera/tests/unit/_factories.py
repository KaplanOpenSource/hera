"""Deterministic synthetic inputs for the unit layer.

Synthetic input is not a compromise here — it is what gives control over
the edge cases the S3 data set happens not to contain: NaN, a single row,
a constant column, an empty frame.
"""
import numpy as np
import pandas as pd


def timeseries(n=10, freq="1min", columns=("u", "v", "w"), seed=0):
    """A DatetimeIndex frame with reproducible pseudo-random values."""
    rng = np.random.default_rng(seed)
    index = pd.date_range("2020-01-01 00:00:00", periods=n, freq=freq)
    return pd.DataFrame(
        {name: rng.standard_normal(n) for name in columns},
        index=index,
    )


def elevation_grid(nx=3, ny=2, base=100.0, step=50.0):
    """A tabular elevation grid with columns X, Y, Elevation."""
    xs = np.arange(1.0, nx + 1.0)
    ys = np.arange(10.0, 10.0 * (ny + 1), 10.0)
    rows = []
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            rows.append((x, y, base + step * i + 10.0 * j))
    return pd.DataFrame(rows, columns=["X", "Y", "Elevation"])


def points_geodataframe(points, crs=4326):
    """A GeoDataFrame of point geometries from (lon, lat) tuples."""
    import geopandas as gpd
    from shapely.geometry import Point

    return gpd.GeoDataFrame(
        {"id": list(range(len(points)))},
        geometry=[Point(x, y) for x, y in points],
        crs=crs,
    )
