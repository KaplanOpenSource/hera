#!/bin/bash
set -e

echo "Building base image..."
docker build -t hera-base -f Dockerfile .

echo "Building server image..."
docker build -t hera-server -f server/Dockerfile .

echo "Build complete."
