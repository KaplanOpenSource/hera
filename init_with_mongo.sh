#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set up environment variables
source "${SCRIPT_DIR}/set_hera_environment.sh"

CONTAINER_NAME="hera-mongo"
MONGO_IMAGE="mongo:5.0"
DATA_DIR="${HOME}/mongo-db-datadir"
INIT_DIR="${SCRIPT_DIR}/mongo-init.d"
SYSTEM_USER="$(whoami)"

echo "=== Hera MongoDB Initialization ==="

# 1. Create required directories
echo "Creating directories..."
mkdir -p "${DATA_DIR}"
mkdir -p "${PYHERA_DIR}/log"

# 2. Pull Docker images
echo "Pulling ${MONGO_IMAGE}..."
docker pull "${MONGO_IMAGE}"

echo "Pulling mermaid-cli (for diagram rendering)..."
docker pull minlag/mermaid-cli

# 3. Remove existing container if present
if docker ps -a --format '{{.Names}}' | grep -q "^${CONTAINER_NAME}$"; then
    echo "Removing existing '${CONTAINER_NAME}' container..."
    docker rm -f "${CONTAINER_NAME}" >/dev/null
fi

# 4. Start MongoDB container with init scripts
echo "Starting MongoDB container..."
docker run --name "${CONTAINER_NAME}" \
    -v "${DATA_DIR}":/data/db \
    -v "${INIT_DIR}":/docker-entrypoint-initdb.d \
    -p 127.0.0.1:27017:27017 \
    -d "${MONGO_IMAGE}"

# 5. Wait for MongoDB to be ready
echo "Waiting for MongoDB to be ready..."
for i in $(seq 1 30); do
    if docker exec "${CONTAINER_NAME}" mongosh --quiet --eval "db.runCommand({ping:1})" >/dev/null 2>&1; then
        echo "MongoDB is ready."
        break
    fi
    if [ "$i" -eq 30 ]; then
        echo "ERROR: MongoDB did not become ready in time."
        exit 1
    fi
    sleep 1
done

# 6. Create config.json
CONFIG_FILE="${PYHERA_DIR}/config.json"
if [ -f "${CONFIG_FILE}" ]; then
    echo "config.json already exists at ${CONFIG_FILE}, skipping creation."
else
    echo "Creating ${CONFIG_FILE}..."
    cat > "${CONFIG_FILE}" <<EOF
{
    "${SYSTEM_USER}": {
        "dbIP": "127.0.0.1",
        "dbName": "olymp",
        "password": "heracles",
        "username": "hera"
    }
}
EOF
fi

echo ""
echo "=== Done ==="
echo "MongoDB is running in container '${CONTAINER_NAME}'."
echo "Config written to ${CONFIG_FILE}."
echo ""
echo "To stop/start MongoDB later:"
echo "  docker stop ${CONTAINER_NAME}"
echo "  docker start ${CONTAINER_NAME}"
