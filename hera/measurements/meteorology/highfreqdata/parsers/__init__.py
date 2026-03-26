"""
Parsers for high-frequency meteorological data formats.

Provides normalisation helpers that ensure consistent output across all
parsers: lowercase column names, proper DatetimeIndex, float dtypes,
and separated metadata.
"""

import pandas as pd

# Standard column sets
SONIC_COLUMNS = ["u", "v", "w", "T"]
TRH_COLUMNS = ["TC_T", "TRH", "RH"]

# Column rename maps
_SONIC_RENAME = {"U": "u", "V": "v", "W": "w"}
_TRH_RENAME = {"TcT": "TC_T"}  # CampbellBinary uses TcT


def _ensure_datetime_index(df):
    """Ensure the DataFrame has a DatetimeIndex.

    Checks, in order: existing DatetimeIndex, ``TIMESTAMP`` column,
    ``Time`` column, ``timestamp`` column.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame.

    Returns
    -------
    pandas.DataFrame
        DataFrame with a DatetimeIndex named ``Time``.
    """
    if isinstance(df.index, pd.DatetimeIndex):
        if df.index.name is None:
            df.index.name = "Time"
        return df

    for col in ("TIMESTAMP", "Time", "timestamp"):
        if col in df.columns:
            df[col] = pd.to_datetime(df[col], errors="coerce")
            df = df.dropna(subset=[col]).set_index(col).sort_index()
            df.index.name = "Time"
            return df

    # Last resort: try converting the existing index
    try:
        df.index = pd.to_datetime(df.index)
        df.index.name = "Time"
    except Exception:
        pass
    return df


def _to_float(df, columns):
    """Convert specified columns to float, coercing errors to NaN."""
    for col in columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def normalise_sonic_df(df):
    """Normalise a sonic anemometer DataFrame to the standard format.

    - Renames uppercase columns (``U → u``, ``V → v``, ``W → w``)
    - Renames ``TcT → TC_T`` if present (CampbellBinary convention)
    - Ensures a ``DatetimeIndex``
    - Converts data columns to float
    - Drops non-data columns (``height``, ``station``, ``instrument``,
      ``RECORD``) and returns them as a metadata dict

    Parameters
    ----------
    df : pandas.DataFrame
        Raw parser output for sonic data.

    Returns
    -------
    tuple of (pandas.DataFrame, dict)
        Normalised DataFrame with columns ``[u, v, w, T]`` and a
        DatetimeIndex, plus a metadata dict with extracted fields.
    """
    df = df.copy()

    # Extract metadata columns before dropping
    metadata = {}
    for col in ("height", "station", "instrument"):
        if col in df.columns:
            vals = df[col].dropna().unique()
            if len(vals) > 0:
                metadata[col] = vals[0] if len(vals) == 1 else list(vals)
            df = df.drop(columns=[col])

    # Drop RECORD if present
    if "RECORD" in df.columns:
        df = df.drop(columns=["RECORD"])

    # Rename columns
    df = df.rename(columns=_SONIC_RENAME)
    df = df.rename(columns=_TRH_RENAME)

    # Ensure DatetimeIndex
    df = _ensure_datetime_index(df)

    # Convert to float
    data_cols = [c for c in SONIC_COLUMNS if c in df.columns]
    df = _to_float(df, data_cols)

    # Keep only known data columns (drop any extras)
    keep = [c for c in df.columns if c in SONIC_COLUMNS]
    df = df[keep]

    metadata["deviceType"] = "sonic"
    return df, metadata


def normalise_trh_df(df):
    """Normalise a TRH sensor DataFrame to the standard format.

    - Renames ``TcT → TC_T`` if present
    - Ensures a ``DatetimeIndex``
    - Converts data columns to float
    - Drops non-data columns and returns them as metadata

    Parameters
    ----------
    df : pandas.DataFrame
        Raw parser output for TRH data.

    Returns
    -------
    tuple of (pandas.DataFrame, dict)
        Normalised DataFrame with columns ``[TC_T, TRH, RH]`` and a
        DatetimeIndex, plus a metadata dict.
    """
    df = df.copy()

    # Extract metadata columns
    metadata = {}
    for col in ("height", "station", "instrument"):
        if col in df.columns:
            vals = df[col].dropna().unique()
            if len(vals) > 0:
                metadata[col] = vals[0] if len(vals) == 1 else list(vals)
            df = df.drop(columns=[col])

    if "RECORD" in df.columns:
        df = df.drop(columns=["RECORD"])

    # Rename
    df = df.rename(columns=_TRH_RENAME)

    # Ensure DatetimeIndex
    df = _ensure_datetime_index(df)

    # Convert to float
    data_cols = [c for c in TRH_COLUMNS if c in df.columns]
    df = _to_float(df, data_cols)

    # Keep only known data columns
    keep = [c for c in df.columns if c in TRH_COLUMNS]
    df = df[keep]

    metadata["deviceType"] = "TRH"
    return df, metadata


def detect_device_type(df):
    """Detect whether a DataFrame contains sonic or TRH data.

    Parameters
    ----------
    df : pandas.DataFrame

    Returns
    -------
    str
        ``'sonic'`` or ``'TRH'`` or ``'unknown'``.
    """
    cols = set(c.lower() for c in df.columns)
    if {"u", "v", "w"}.issubset(cols):
        return "sonic"
    if "tc_t" in cols or "tct" in cols:
        return "TRH"
    return "unknown"
