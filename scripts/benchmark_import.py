#!/usr/bin/env python3
"""
Import-time benchmark for hera. Run before and after each optimization task
to track progress toward the <50ms goal.

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
