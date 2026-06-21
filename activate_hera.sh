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

# Activate the virtual environment first — its activate script calls
# `deactivate nondestructive`, which restores PATH from _OLD_VIRTUAL_PATH
# and would drop any PATH edits made before sourcing it.
source "${VENV_DIR}/bin/activate"

export PATH="${SCRIPT_DIR}/hera/bin${PATH:+:${PATH}}"

echo "Hera environment activated (venv: ${VENV_DIR})"
