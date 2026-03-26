# Data Filtering

The `hera.utils.filter_immediate` module provides a `Filter` class for applying threshold and interval filters to pandas DataFrames. Filters are chainable, making it easy to build multi-step data-cleaning pipelines.

## Importing

```python
from hera.utils.filter_immediate import Filter
```

## Basic Usage

```python
import pandas as pd
from hera.utils.filter_immediate import Filter

df = pd.DataFrame({
    "temperature": [15, 22, 35, -5, 28, 100],
    "wind_speed": [3, 7, 12, 2, 50, 5],
}, index=pd.date_range("2024-01-01", periods=6, freq="h"))

# Keep only rows where temperature > 0 and wind_speed < 20
result = (
    Filter(df)
    .threshold("gt", 0, column="temperature")
    .threshold("lt", 20, column="wind_speed")
)

print(result.data)
```

Each `.threshold()` call returns a new `Filter` wrapping the filtered DataFrame. Access the final DataFrame via `.data`.

## Threshold Operators

The `preposition` argument selects the comparison:

| Operator | Meaning |
|---|---|
| `"lt"` | Less than |
| `"lte"` | Less than or equal |
| `"gt"` | Greater than |
| `"gte"` | Greater than or equal |
| `"eq"` | Equal |
| `"neq"` | Not equal |
| `"abs_lt"` | Absolute value less than |
| `"abs_lte"` | Absolute value less than or equal |
| `"abs_gt"` | Absolute value greater than |
| `"abs_gte"` | Absolute value greater than or equal |

## Filtering by Index

Omit the `column` argument to filter on the DataFrame index:

```python
# Keep rows after a certain timestamp
result = Filter(df).threshold("gte", pd.Timestamp("2024-01-01 02:00"), column=None)
print(result.data)
```

## Removing Intervals with outsideInterval

Use `.outsideInterval(lower, upper)` to remove rows that fall within `[lower, upper)`:

```python
# Remove wind speeds between 5 and 10 (exclusive of upper bound)
result = Filter(df).outsideInterval(5, 10, column="wind_speed")
print(result.data)
```

## In-Place Filtering

By default, each filter step creates a new `Filter` object. Pass `inplace=True` to mutate the same instance:

```python
f = Filter(df, inplace=True)
f.threshold("gt", 0, column="temperature")
f.threshold("lt", 20, column="wind_speed")

# f.data has been modified in place
print(f.data)
```

## Practical Data Cleaning Example

```python
import pandas as pd
from hera.utils.filter_immediate import Filter

# Load raw sensor data
raw = pd.read_csv("sensor_data.csv", parse_dates=["timestamp"], index_col="timestamp")

cleaned = (
    Filter(raw)
    # Remove negative humidity (sensor error)
    .threshold("gte", 0, column="humidity")
    # Cap temperature at physical bounds
    .threshold("gt", -40, column="temperature")
    .threshold("lt", 60, column="temperature")
    # Remove wind speed outliers
    .threshold("abs_lt", 50, column="wind_speed")
    # Drop a known bad-data window
    .outsideInterval(
        pd.Timestamp("2024-03-15 12:00"),
        pd.Timestamp("2024-03-15 14:00"),
    )
)

print(f"Kept {len(cleaned.data)} of {len(raw)} rows")
```
