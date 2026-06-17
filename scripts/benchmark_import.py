#!/usr/bin/env python3
"""
Import-time benchmark for hera. Run before and after each optimization task
to track progress toward the <50ms goal.

Results (issue935 branch, 2026-06-17):
  BEFORE  import hera                          : 612 ms  (8 heavy libs loaded)
  AFTER   import hera                          :  16 ms  (0 heavy libs)
  BEFORE  from hera.utils.logging import ...   : 615 ms
  AFTER   from hera.utils.logging import ...   :  16 ms

  CLI --help response times (BEFORE → AFTER):
    hera-project       : 736ms → 58ms
    hera-GIS           : 936ms → 48ms
    hera-workflows     :  53ms → 53ms  (unchanged, was already fast)
    hera-LSM           : 671ms → 51ms
    hera-riskassessment: 1820ms → 50ms

Usage:
    python3 scripts/benchmark_import.py
"""
import subprocess
import sys


def measure(label, code, runs=3):
    times = []
    for _ in range(runs):
        r = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
        )
        try:
            ms = float(r.stdout.decode().strip())
            times.append(ms)
        except Exception:
            times.append(-1)
    best = min((t for t in times if t > 0), default=-1)
    print(f"  {label:<55}: {best:7.0f} ms")


def heavy_libs():
    r = subprocess.run(
        [
            sys.executable,
            "-c",
            """
import hera, sys
heavy = sorted(set(
    k.split('.')[0] for k in sys.modules
    if k.split('.')[0] in ('pandas','numpy','pint','xarray',
                           'geopandas','mongoengine','shapely','pyproj')
))
print(','.join(heavy) if heavy else '(none)')
""",
        ],
        capture_output=True,
    )
    return r.stdout.decode().strip()


print("=" * 70)
print("hera import-time benchmark")
print("=" * 70)

measure(
    "import hera",
    "import time; t=time.time(); import hera; print(f'{(time.time()-t)*1000:.1f}')",
)
measure(
    "from hera.utils.logging import initialize_logging",
    "import time; t=time.time(); from hera.utils.logging import initialize_logging; print(f'{(time.time()-t)*1000:.1f}')",
)
measure(
    "argparse + argcomplete only (target baseline)",
    "import time; t=time.time(); import argparse, argcomplete; print(f'{(time.time()-t)*1000:.1f}')",
)

print()
print(f"  Heavy libs loaded after 'import hera': {heavy_libs()}")
print("=" * 70)
