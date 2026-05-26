"""Engine-agnostic helpers for reading vector GIS payloads.

GeoPandas defaults to the ``pyogrio`` engine on Python 3.11+; ``pyogrio`` does
not accept text streams (``io.StringIO``) and reinterprets them as filesystem
paths, producing a misleading ``DataSourceError``. This module centralises the
"parse a GeoJSON string into a GeoDataFrame" path so callers do not have to
care about the active engine and so a non-JSON input fails cleanly with
``ValueError`` rather than an engine-specific exception.
"""

import io

import geopandas


def readGeoJSONString(geojsonText):
    """Parse a GeoJSON string into a ``geopandas.GeoDataFrame``.

    Wraps the payload in ``io.BytesIO`` so it is accepted by both ``pyogrio``
    (default on modern GeoPandas) and the legacy ``fiona`` engine. Validates
    that the input is a non-empty string whose first non-whitespace character
    is ``{`` or ``[`` — anything else (typos, stray identifiers, paths that
    happen not to exist) raises ``ValueError`` before reaching the engine.
    """
    if not isinstance(geojsonText, str):
        raise ValueError(
            f"Expected a GeoJSON string, got {type(geojsonText).__name__}"
        )

    stripped = geojsonText.lstrip()
    if not stripped or stripped[0] not in "{[":
        raise ValueError(
            f"Value is not a path on disk and does not look like a GeoJSON "
            f"string: {geojsonText!r}"
        )

    return geopandas.read_file(io.BytesIO(geojsonText.encode("utf-8")))
