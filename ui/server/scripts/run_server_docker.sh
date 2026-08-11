#!/bin/bash
# Runs the Hera UI server in Docker.
# Any arguments are forwarded to server.py (e.g. --debug).
#
# Handles the runtime bits the image doesn't bake in yet:
#   - PYTHONPATH includes /app/Hermes so `import hermes` works (submodule, no setup.py)
#   - installs hermes runtime deps (jsonpath_rw_ext, luigi) before start
# See ui/KNOWN_INTEGRATION_ISSUES.md for why.
set -e

# Run from the repo root so the bind-mount (-v .:/app) maps the whole project.
cd "$(dirname "$0")/../../.."

docker start hera-mongo

docker run -it --network host \
  -v "$(pwd)":/app \
  -e PYTHONPATH=/app:/app/hera/bin:/app/Hermes \
  --rm --name hera-server-instance hera-server \
  bash -c "
    python ui/server/server.py $*
  "
