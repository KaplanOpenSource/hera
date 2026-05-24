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
# Offer to create Python virtual environment in ~/.pyhera/environment
VENV_DIR="${PYHERA_DIR}/environment"

if [ -d "${VENV_DIR}" ] && [ -f "${VENV_DIR}/bin/activate" ]; then
    echo "Python virtual environment already exists at ${VENV_DIR}."
else
    read -p "Create a Python virtual environment at ${VENV_DIR}? [y/N] " venv_answer
    if [ "$venv_answer" = "y" ] || [ "$venv_answer" = "Y" ]; then
        echo "Creating virtual environment..."
        mkdir -p "${PYHERA_DIR}"
        python3 -m venv "${VENV_DIR}"
        echo "Installing requirements (this may take a while)..."
        "${VENV_DIR}/bin/pip" install --upgrade pip
        "${VENV_DIR}/bin/pip" install -r "${SCRIPT_DIR}/requirements.txt"
        "${VENV_DIR}/bin/pip" install -e "${SCRIPT_DIR}"

        # Patch PyFoam's vendored six.py (v1.10 from 2016) which breaks on Python 3.12+.
        SITE_SIX="$("${VENV_DIR}/bin/python" -c 'import six, pathlib; print(pathlib.Path(six.__file__))' 2>/dev/null)"
        PYFOAM_SIX="$("${VENV_DIR}/bin/python" -c 'import PyFoam.ThirdParty, pathlib; print(pathlib.Path(PyFoam.ThirdParty.__file__).parent / "six.py")' 2>/dev/null)"
        if [ -f "${SITE_SIX}" ] && [ -f "${PYFOAM_SIX}" ]; then
            cp "${SITE_SIX}" "${PYFOAM_SIX}"
            echo "Patched PyFoam vendored six.py -> ${PYFOAM_SIX}"
        fi

        echo "Virtual environment created and packages installed at ${VENV_DIR}."
    else
        echo "Skipped virtual environment creation."
        echo "You can create one manually:"
        echo "  python3 -m venv ${VENV_DIR}"
        echo "  source ${VENV_DIR}/bin/activate"
        echo "  pip install -r ${SCRIPT_DIR}/requirements.txt"
        echo "  pip install -e ${SCRIPT_DIR}"
    fi
fi
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

# Offer to auto-activate the Hera environment when entering the project directory
ACTIVATE_LINE="source ${SCRIPT_DIR}/activate_hera.sh"
AUTO_ACTIVATE_MARKER="# Hera auto-activate"

if grep -qF "${AUTO_ACTIVATE_MARKER}" "${BASHRC}" 2>/dev/null; then
    echo "Auto-activation on directory change is already configured in ${BASHRC}."
elif [ -d "${VENV_DIR}" ] && [ -f "${VENV_DIR}/bin/activate" ]; then
    echo ""
    read -p "Auto-activate Hera environment when you cd into ${SCRIPT_DIR}? [y/N] " auto_answer
    if [ "$auto_answer" = "y" ] || [ "$auto_answer" = "Y" ]; then
        cat >> "${BASHRC}" <<AUTOEOF

${AUTO_ACTIVATE_MARKER}
cd() {
    builtin cd "\$@" || return
    if [ -f "\${PWD}/activate_hera.sh" ] && [ "\${PWD}" = "${SCRIPT_DIR}" ]; then
        source "\${PWD}/activate_hera.sh"
    fi
}
AUTOEOF
        echo "Added auto-activation to ${BASHRC}."
        echo "The Hera environment will activate automatically when you cd into ${SCRIPT_DIR}."
    else
        echo "Skipped. You can activate manually with: source ${SCRIPT_DIR}/activate_hera.sh"
    fi
fi
