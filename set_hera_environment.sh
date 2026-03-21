#!/bin/bash
# Set up Hera environment variables.
# Usage: source set_hera_environment.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export HERA_REPO_ROOT="${SCRIPT_DIR}"
export PYHERA_DIR="${HOME}/.pyhera"
export PYTHONPATH="${SCRIPT_DIR}${PYTHONPATH:+:${PYTHONPATH}}"

echo ""
echo "=== Hera Environment ==="
echo "  HERA_REPO_ROOT=${HERA_REPO_ROOT}"
echo "  PYHERA_DIR=${PYHERA_DIR}"
echo "  PYTHONPATH=${PYTHONPATH}"
echo ""
echo "To create a Python virtual environment with all requirements:"
echo ""
echo "  python -m venv .venv"
echo "  source .venv/bin/activate"
echo "  pip install -r ${SCRIPT_DIR}/requirements.txt"
echo "  pip install -e ${SCRIPT_DIR}"
echo ""

# Offer to persist to .bashrc
BASHRC="${HOME}/.bashrc"
SOURCE_LINE="source ${SCRIPT_DIR}/set_hera_environment.sh"

if grep -qF "${SOURCE_LINE}" "${BASHRC}" 2>/dev/null; then
    echo "Hera environment is already configured in ${BASHRC}."
else
    read -p "Add Hera environment to ${BASHRC} so it loads on every new shell? [y/N] " answer
    if [ "$answer" = "y" ] || [ "$answer" = "Y" ]; then
        echo "" >> "${BASHRC}"
        echo "# Hera environment setup" >> "${BASHRC}"
        echo "${SOURCE_LINE}" >> "${BASHRC}"
        echo "Added to ${BASHRC}. It will take effect in new shell sessions."
    else
        echo "Skipped. Run 'source ${SCRIPT_DIR}/set_hera_environment.sh' manually each session."
    fi
fi
