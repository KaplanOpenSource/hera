"""Ad-hoc inspection helper for a high-frequency sonic parquet file.

Run it, do not import it:

    python -m hera.measurements.meteorology.highfreqdata <path-to.parquet>
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
