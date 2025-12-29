#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./tests/run_all_json_tests.sh prepare --result-set <NAME>
#   ./tests/run_all_json_tests.sh run     --result-set <NAME>
#
# Notes:
#   - Sources local variables from ./tests/env.template (gitignored).
#   - --result-set is mandatory (or set RESULT_SET in env).
#   - Auto-detects json_definitions directory in root or under tests/.

MODE="${1:-}"
shift || true

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_FILE="${REPO_ROOT}/tests/env.template"

if [[ -z "${MODE}" || ( "${MODE}" != "prepare" && "${MODE}" != "run" ) ]]; then
  echo "Usage: $0 <prepare|run> --result-set <NAME>"
  exit 2
fi

if [[ ! -f "${ENV_FILE}" ]]; then
  echo "Missing tests/env.template at ${ENV_FILE}"
  exit 1
fi

# Ensure PYTHONPATH is defined even under 'set -u'
: "${PYTHONPATH:=}"

# Load local env (no hard-coded exports in this script)
set -a
source "${ENV_FILE}"
set +a

RESULT_SET_ARG=""
EXTRA_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --result-set) shift; RESULT_SET_ARG="${1:-}"; shift || true ;;
    *) EXTRA_ARGS+=("$1"); shift ;;
  esac
done

if [[ -z "${RESULT_SET_ARG}" && -z "${RESULT_SET:-}" ]]; then
  echo "Missing --result-set. Example: --result-set BASELINE (or set RESULT_SET in tests/env.template)"
  exit 2
fi
export RESULT_SET="${RESULT_SET_ARG:-${RESULT_SET}}"

# Expected directory must live under HERA_UNITTEST_DATA
if [[ -z "${HERA_UNITTEST_DATA:-}" ]]; then
  echo "HERA_UNITTEST_DATA is not set in tests/env.template"
  exit 1
fi
EXP_DIR="${HERA_UNITTEST_DATA}/expected/${RESULT_SET}"

if [[ "${MODE}" == "prepare" ]]; then
  mkdir -p "${EXP_DIR}"
  export PREPARE_EXPECTED_OUTPUT=1
else
  if [[ ! -d "${EXP_DIR}" ]]; then
    echo "Expected results set '${RESULT_SET}' not found at ${EXP_DIR}"
    echo "First create it: $0 prepare --result-set ${RESULT_SET}"
    exit 3
  fi
  export PREPARE_EXPECTED_OUTPUT=0
fi

PYTHON_BIN="${PYTHON_BIN:-python3}"
cd "${REPO_ROOT}"

# Resolve runner location (root or tests)
RUNNER_ROOT="${REPO_ROOT}/run_all_definitions.py"
RUNNER_TESTS="${REPO_ROOT}/tests/run_all_definitions.py"
if [[ -f "${RUNNER_ROOT}" ]]; then
  RUNNER="${RUNNER_ROOT}"
elif [[ -f "${RUNNER_TESTS}" ]]; then
  RUNNER="${RUNNER_TESTS}"
else
  echo "Cannot find run_all_definitions.py at:"
  echo "  - ${RUNNER_ROOT}"
  echo "  - ${RUNNER_TESTS}"
  exit 4
fi

# Resolve json_definitions directory (root or tests)
if [[ -d "${REPO_ROOT}/json_definitions" ]]; then
  JSON_DIR="${REPO_ROOT}/json_definitions"
elif [[ -d "${REPO_ROOT}/tests/json_definitions" ]]; then
  JSON_DIR="${REPO_ROOT}/tests/json_definitions"
else
  echo "Cannot find json_definitions at:"
  echo "  - ${REPO_ROOT}/json_definitions"
  echo "  - ${REPO_ROOT}/tests/json_definitions"
  exit 5
fi

echo "MODE=${MODE} | RESULT_SET=${RESULT_SET} | EXP_DIR=${EXP_DIR}"
echo "RUNNER=${RUNNER}"
echo "JSON_DIR=${JSON_DIR}"
echo

# Build absolute paths to JSON definition files
DEF_TOP="${JSON_DIR}/topography_toolkit_definitions_extended.json"
DEF_LOW="${JSON_DIR}/lowfreq_toolkit.json"
DEF_DEM="${JSON_DIR}/demography_test_definitions.json"
DEF_HF="${JSON_DIR}/test_definitions_highfreq.json"
DEF_LC="${JSON_DIR}/test_definitions_landcover.json"
DEF_DYN_TK="${JSON_DIR}/test_dynamic_toolkits.json"
DEF_EXP="${JSON_DIR}/test_experiment_loading.json"

# Assemble args without passing an empty extra-arg (bash-safe)
ARGS=(
  "${RUNNER}"
  "${DEF_TOP}"
  "${DEF_LOW}"
  "${DEF_DEM}"
  "${DEF_HF}"
  "${DEF_LC}"
  "${DEF_DYN_TK}"
  "${DEF_EXP}"
)
if ((${#EXTRA_ARGS[@]})); then
  ARGS+=("${EXTRA_ARGS[@]}")
fi

"${PYTHON_BIN}" "${ARGS[@]}"

echo
echo "Done."
