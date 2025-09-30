#!/bin/bash
set -e

# echo "Building base image..."
# docker build -t hera-base -f Dockerfile .

echo "Building server image..."
docker build -t hera-server -f Dockerfile .

echo "\nBuild complete.\n\nTo run server use:"
echo "docker run -it -p 8000:8000 --rm --name hera-server-instance hera-server"
