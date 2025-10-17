#!/usr/bin/env bash
set -euo pipefail

# 1) ודא שאנחנו ברוט של ה־repo
cd "$(dirname "$0")/.."

# 2) אתר את python מה־venv
if [[ -x "./heraenv/bin/python" ]]; then
  PY=./heraenv/bin/python
elif [[ -x "$HOME/hera/heraenv/bin/python" ]]; then
  PY="$HOME/hera/heraenv/bin/python"
else
  echo "❌ Could not find hera venv python at ./heraenv/bin/python or ~/hera/heraenv/bin/python"
  exit 1
fi
echo "✅ Using: $("$PY" -c 'import sys; print(sys.executable)')"

# 3) ודא שקיים pytest ב־venv
if ! "$PY" -m pytest --version >/dev/null 2>&1; then
  echo "ℹ Installing pytest into venv..."
  "$PY" -m pip install -q 'pytest>=7,<9'
fi

# 4) בידוד מפלאגינים חיצוניים והתנגשויות
export PYTEST_DISABLE_PLUGIN_AUTOLOAD=1
export PYTHONNOUSERSITE=1

# 5) smoke test קצר ל־hera
"$PY" - <<'PY'
import hera
print("✅ hera import OK")
PY

echo
echo "==== Running dynamic toolkits tests ===="
"$PY" -m pytest -q hera/tests/UNIT_TEST_DYNAMIC_TOOLKITS

echo
echo "==== Running repository/save/load focused tests ===="
"$PY" -m pytest -q hera/tests -k "repository or RepositoryConfig or getDatasourceDocument or (save and load)"

# (אופציונלי) כל הטסטים:
# echo
# echo "==== Running entire test suite ===="
# "$PY" -m pytest -q -x hera/tests
