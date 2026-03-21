#!/usr/bin/env bash
cd "$(dirname "$0")/../.." || exit 1
docker build -t hera-server -f dockerfile .
