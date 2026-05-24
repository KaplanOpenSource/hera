#!/usr/bin/env bash
# Bootstrap hera unit-test data from a public S3 bucket.
#
# Usage:
#   ./scripts/s3/bootstrap_unittest_data.sh [--base-url URL] [--target-dir DIR]
#
# Defaults:
#   BASE_URL   https://s3.eu-west-1.amazonaws.com/hera.kaplanopensource.co.il
#   TARGET_DIR ~/hera_unittest_data
#
# Bucket layout (as deployed by Lior):
#   {BASE_URL}/latest.json
#   {BASE_URL}/manifest.json
#   {BASE_URL}/hera_unittest_data/<file-path>
#
# Requires: curl, python3, sha256sum

set -euo pipefail

BASE_URL="https://s3.eu-west-1.amazonaws.com/hera.kaplanopensource.co.il"
TARGET_DIR="${HOME}/hera_unittest_data"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --base-url)   BASE_URL="$2";   shift 2 ;;
        --target-dir) TARGET_DIR="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

WORK_DIR=$(mktemp -d)
trap 'rm -rf "${WORK_DIR}"' EXIT

echo "=========================================="
echo " Hera S3 Bootstrap"
echo "=========================================="
echo "  Base URL  : ${BASE_URL}"
echo "  Target dir: ${TARGET_DIR}"
echo ""

# ── Step 1: Fetch latest.json ─────────────────────────────────────────────────
echo "[1/4] Fetching latest.json ..."
curl -fsSL --max-time 30 "${BASE_URL}/latest.json" -o "${WORK_DIR}/latest.json"
VERSION=$(python3 -c "import json; print(json.load(open('${WORK_DIR}/latest.json'))['version'])")
echo "      Version: ${VERSION}"

# ── Step 2: Fetch manifest.json ───────────────────────────────────────────────
echo "[2/4] Fetching manifest.json ..."
curl -fsSL --max-time 30 "${BASE_URL}/manifest.json" -o "${WORK_DIR}/manifest.json"
FILE_COUNT=$(python3 -c "import json; print(len(json.load(open('${WORK_DIR}/manifest.json'))['files']))")
echo "      Files in manifest: ${FILE_COUNT}"

# ── Step 3: Download data files ───────────────────────────────────────────────
echo "[3/4] Downloading ${FILE_COUNT} files to ${TARGET_DIR} ..."
mkdir -p "${TARGET_DIR}"

python3 -c "
import json
d = json.load(open('${WORK_DIR}/manifest.json'))
for f in d['files']:
    print(f['path'])
" > "${WORK_DIR}/file_paths.txt"

while IFS= read -r FILE_PATH; do
    FILE_URL="${BASE_URL}/hera_unittest_data/${FILE_PATH}"
    DEST="${TARGET_DIR}/${FILE_PATH}"
    mkdir -p "$(dirname "${DEST}")"
    printf "  %-55s" "${FILE_PATH}"
    curl -fsSL --max-time 120 -o "${DEST}" "${FILE_URL}"
    SIZE=$(du -h "${DEST}" | awk '{print $1}')
    echo " [${SIZE}]"
done < "${WORK_DIR}/file_paths.txt"

# ── Step 4: SHA256 integrity check ────────────────────────────────────────────
echo ""
echo "[4/4] Verifying SHA256 hashes ..."
echo ""

python3 -c "
import json
d = json.load(open('${WORK_DIR}/manifest.json'))
for f in d['files']:
    print(f\"{f['path']}\t{f['sha256']}\")
" > "${WORK_DIR}/hashes.tsv"

PASS=0
FAIL=0

while IFS=$'\t' read -r FILE_PATH EXPECTED_HASH; do
    DEST="${TARGET_DIR}/${FILE_PATH}"
    if [[ ! -f "${DEST}" ]]; then
        echo "  MISSING      : ${FILE_PATH}"
        FAIL=$((FAIL + 1))
        continue
    fi
    ACTUAL_HASH=$(sha256sum "${DEST}" | awk '{print $1}')
    if [[ "${ACTUAL_HASH}" == "${EXPECTED_HASH}" ]]; then
        SIZE=$(du -h "${DEST}" | awk '{print $1}')
        echo "  OK [${SIZE}]  : ${FILE_PATH}"
        PASS=$((PASS + 1))
    else
        echo "  HASH MISMATCH: ${FILE_PATH}"
        echo "    expected : ${EXPECTED_HASH}"
        echo "    actual   : ${ACTUAL_HASH}"
        FAIL=$((FAIL + 1))
    fi
done < "${WORK_DIR}/hashes.tsv"

echo ""
echo "=========================================="
echo " Bootstrap summary"
echo "   Version : ${VERSION}"
echo "   Passed  : ${PASS} / ${FILE_COUNT}"
if [[ ${FAIL} -gt 0 ]]; then
    echo "   Failed  : ${FAIL} / ${FILE_COUNT}"
    echo "   Status  : FAIL"
    exit 1
else
    echo "   Status  : OK"
fi
echo "=========================================="
