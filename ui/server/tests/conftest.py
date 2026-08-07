import os
import sys

# Make the server modules (server.py, api_models.py, ...) importable regardless of
# where pytest is run from.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# server.py calls argparse.parse_args() at import time; give it clean args so it
# doesn't try to parse pytest's argv.
sys.argv = ["server"]
