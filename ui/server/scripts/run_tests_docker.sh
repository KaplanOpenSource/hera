#!/bin/bash
# Run the UI server tests in the hera-server image.
# Usage: bash run_tests_docker.sh [pytest args]   (bash, not sh)
set -e
cd "$(dirname "$0")/../../.."
docker run --rm \
  -v "$(pwd)":/app \
  -w /app/ui/server \
  hera-server \
  bash -c '
    pip install pytest &&
    python -m pytest tests "$@"
  ' _ "$@"
