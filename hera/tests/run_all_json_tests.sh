#!/bin/bash

# 🪄 מצב הרצה: שמירה או השוואה
PREPARE_EXPECTED_OUTPUT=${PREPARE_EXPECTED_OUTPUT:-0}
echo "🚀 Running all JSON-based tests with PREPARE_EXPECTED_OUTPUT=${PREPARE_EXPECTED_OUTPUT} ..."
echo

# 📁 משתני סביבה
export HERA_DATA_PATH=/home/ilay/hera_unittest_data
export HERA_UNITTEST_DATA=/home/ilay/hera_unittest_data
export PYTHONPATH=$(pwd)/../

# 🧪 הרצת הבדיקות
python3 run_all_definitions.py \
  json_definitions/topography_toolkit_definitions_extended.json \
  json_definitions/lowfreq_toolkit.json \
  json_definitions/demography_test_definitions.json \
  json_definitions/test_definitions_highfreq.json \
  json_definitions/test_definitions_landcover.json

echo
echo "✅ All tests finished!"
