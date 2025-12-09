#!/bin/bash
cd "$(dirname "$0")/../.." # Go to hera repo root
pip install -r ui/server/requirements-server.txt
python ui/server/server.py
