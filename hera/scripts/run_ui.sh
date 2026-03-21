#!/bin/bash
cd "$(dirname "$0")/../.." # Go to hera repo root
pip install -r requirements.txt
python ui/server/server.py
