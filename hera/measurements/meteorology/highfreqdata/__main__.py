"""Ad-hoc inspection helper for a high-frequency sonic parquet file.

Run it, do not import it:

    python -m hera.measurements.meteorology.highfreqdata <path-to.parquet>

The path used to be a module-level constant pointing at one developer's home
directory, so importing this module raised FileNotFoundError on every other
machine.  CLAUDE.md forbids absolute paths for exactly that reason.
"""
import sys

import pandas as pd


def main(path=None):
    """Print the column names of a parquet file."""
    if path is None:
        if len(sys.argv) < 2:
            raise SystemExit(
                "usage: python -m hera.measurements.meteorology.highfreqdata "
                "<path-to.parquet>"
            )
        path = sys.argv[1]

    print(pd.read_parquet(path).columns)


if __name__ == "__main__":
    main()
