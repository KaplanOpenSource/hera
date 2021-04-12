
toMeteorologicalAngle = lambda mathematical_angle: (270 - mathematical_angle) if ((270 - mathematical_angle) >= 0) else (630 - mathematical_angle)
toMathematicalAngle  = toMeteorologicalAngle
toAzimuthAngle = lambda mathematical_angle: (90 - mathematical_angle) if ((90 - mathematical_angle) >= 0) else (450 - mathematical_angle)