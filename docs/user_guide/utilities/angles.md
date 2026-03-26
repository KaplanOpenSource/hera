# Angle Conversions

The `hera.utils.angle` module provides functions to convert between the three angle conventions commonly used in meteorology, mathematics, and navigation.

## Importing

```python
from hera.utils import toMeteorologicalAngle, toMathematicalAngle, toAzimuthAngle
```

## Angle Convention Reference

| Convention | 0 degrees | Direction of increase | Typical use |
|---|---|---|---|
| **Mathematical** | East | Counter-clockwise | Trigonometry, polar coordinates |
| **Meteorological** | North | Clockwise (wind *from*) | Wind direction reporting |
| **Azimuth** | North | Clockwise | Navigation, compass bearing |

The mathematical convention is the input for all three functions.

## Functions

### toMeteorologicalAngle

Convert a mathematical angle to meteorological convention. The result indicates the direction wind is coming *from*.

```python
from hera.utils import toMeteorologicalAngle

# East wind in math convention (0 deg) -> coming from the West in meteo (270 deg)
toMeteorologicalAngle(0)    # 270
toMeteorologicalAngle(90)   # 180
toMeteorologicalAngle(180)  # 90
toMeteorologicalAngle(270)  # 0
```

### toMathematicalAngle

Convert a meteorological angle back to mathematical convention. This is the same transformation as `toMeteorologicalAngle` (the conversion is its own inverse):

```python
from hera.utils import toMathematicalAngle

toMathematicalAngle(270)  # 0  (East)
toMathematicalAngle(0)    # 270
```

### toAzimuthAngle

Convert a mathematical angle to an azimuth (compass) bearing (0 = North, increasing clockwise):

```python
from hera.utils import toAzimuthAngle

toAzimuthAngle(0)    # 90   (East -> 90 deg azimuth)
toAzimuthAngle(90)   # 0    (North -> 0 deg azimuth)
toAzimuthAngle(180)  # 270  (West -> 270 deg azimuth)
toAzimuthAngle(270)  # 180  (South -> 180 deg azimuth)
```

## Quick Reference Table

| Math angle | Compass direction | Meteorological | Azimuth |
|---|---|---|---|
| 0 | East | 270 | 90 |
| 90 | North | 180 | 0 |
| 180 | West | 90 | 270 |
| 270 | South | 0 | 180 |
