#!/bin/bash
# Start the MkDocs development server for local documentation preview.
# The site will be available at http://127.0.0.1:8000
# Changes to docs/ and mkdocs.yml are auto-reloaded in the browser.
#
# Usage: ./serve_docs.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# Check if mkdocs is installed
if ! command -v mkdocs &> /dev/null; then
    echo "mkdocs not found. Installing docs dependencies..."
    pip install -r docs/requirements-docs.txt
    pip install -e . --no-deps
fi

echo ""
echo "=== Hera Documentation Server ==="
echo "  http://127.0.0.1:8000"
echo "  Press Ctrl+C to stop"
echo ""

mkdocs serve
