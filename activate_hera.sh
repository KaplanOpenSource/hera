#!/bin/bash
# Activate the Hera Python virtual environment and set environment variables.
# Usage: source activate_hera.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="${HOME}/.pyhera/environment"

if [ ! -f "${VENV_DIR}/bin/activate" ]; then
    echo "ERROR: Hera virtual environment not found at ${VENV_DIR}"
    echo "Run 'source ${SCRIPT_DIR}/set_hera_environment.sh' first to create it."
    return 1 2>/dev/null || exit 1
fi

# Set Hera environment variables
export HERA_REPO_ROOT="${SCRIPT_DIR}"
export PYHERA_DIR="${HOME}/.pyhera"
export PYTHONPATH="${SCRIPT_DIR}${PYTHONPATH:+:${PYTHONPATH}}"

# Activate the virtual environment
source "${VENV_DIR}/bin/activate"

echo "Hera environment activated (venv: ${VENV_DIR})"
