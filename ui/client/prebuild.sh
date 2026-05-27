#!/bin/bash
FILE="src/buildNumber.ts"
TODAY=$(date +%Y%m%d)

if [ -f "$FILE" ]; then
  CURRENT=$(grep -oP "buildNumber = '\K[^']+" "$FILE")
  CURRENT_DATE=${CURRENT%%.*}
  CURRENT_N=${CURRENT##*.}

  if [ "$CURRENT_DATE" = "$TODAY" ]; then
    N=$((CURRENT_N + 1))
  else
    N=1
  fi
else
  N=1
fi

BUILD="${TODAY}.${N}"
echo "export const buildNumber = '${BUILD}';" > "$FILE"
echo "Build number: ${BUILD}"
