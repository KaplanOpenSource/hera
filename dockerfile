# Use Ubuntu 22.04 as the base image
FROM ubuntu:22.04

# Set environment variables to avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install base tools and the libs some Python packages need at build time
RUN apt-get update
RUN apt-get install -y \
    nano curl git build-essential iputils-ping software-properties-common
RUN apt-get install -y \
    libcairo2-dev \
    pkg-config \
    libgirepository1.0-dev

# Install Python 3.11 from deadsnakes (prebuilt, matches CI).
RUN add-apt-repository -y ppa:deadsnakes/ppa && apt-get update && \
    apt-get install -y python3.11 python3.11-venv python3.11-dev python3.11-distutils

# Install into an isolated venv so pip stays clear of apt-owned packages.
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="/opt/venv/bin:${PATH}"
RUN python3.11 -m venv "$VIRTUAL_ENV" && \
    python -m pip install --no-cache-dir --upgrade pip
    
# Install GDAL and other system dependencies    
RUN apt-get update
RUN apt-get install -y \
    libgdal-dev \
    gdal-bin \
    python3-gdal

RUN curl -fsSL https://deb.nodesource.com/setup_22.x | bash - && \
    apt-get install -y nodejs

RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# Set working directory and copy project files (exclude .venv, .git via .dockerignore)
WORKDIR /app
COPY requirements.txt .

# Install Python dependencies
RUN python -m pip install --no-cache-dir -r requirements.txt

# Disable the jupytext lab extension; it is built for JupyterLab 4.4+ and errors on 4.3.
# RUN mkdir -p "$VIRTUAL_ENV/share/jupyter/lab/settings" && \
#     printf '{"disabledExtensions": {"jupyterlab-jupytext": true}}\n' \
#       > "$VIRTUAL_ENV/share/jupyter/lab/settings/page_config.json"

ENV PATH="/app:/app/hera/bin:${PATH}"
ENV PYTHONPATH="/app:/app/hera/bin"

# Default DB credentials — override at runtime via --env or .env file
ENV MONGO_HERA_USER=hera
ENV MONGO_HERA_PWD=heracles

# Create necessary folders and configuration file
RUN mkdir -p /root/.pyhera/log && \
    mkdir -p /root/mongo-db-datadir && \
    echo "{ \"root\": { \"dbIP\": \"127.0.0.1\", \"dbName\": \"olymp\", \"username\": \"${MONGO_HERA_USER}\", \"password\": \"${MONGO_HERA_PWD}\" } }" \
    > /root/.pyhera/config.json

# RUN echo 'mongod --fork --logpath /var/log/mongodb.log --dbpath /data/db' >> /root/.bashrc

COPY ./ui/client /client
WORKDIR /client
RUN npm install
RUN npm run build

WORKDIR /app

# COPY . /app

EXPOSE 27017
EXPOSE 8000
EXPOSE 5678

CMD ["python", "ui/server/server.py"]

# run with  --add-host=host.docker.internal:host-gateway
