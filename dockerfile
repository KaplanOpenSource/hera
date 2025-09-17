# Use Ubuntu 22.04 as the base image
FROM ubuntu:22.04

# Set environment variables to avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages for pyenv and building Python
RUN apt-get update
RUN apt-get install -y \
    nano curl git build-essential \
    libreadline-dev libssl-dev libbz2-dev libsqlite3-dev \
    libffi-dev zlib1g-dev iputils-ping && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install pyenv
RUN curl https://pyenv.run | bash

# Set up environment variables for pyenv
ENV PATH="/root/.pyenv/bin:/root/.pyenv/shims:${PATH}"

# Install Python 3.9.13 using pyenv
RUN pyenv install 3.9.13 && \
    pyenv global 3.9.13

# Install pip for the installed Python version
RUN python -m ensurepip && \
    python -m pip install --no-cache-dir --upgrade pip

# Install GDAL and other system dependencies    
RUN apt-get update
RUN apt-get install -y \
    libcairo2-dev \
    pkg-config \
    libgirepository1.0-dev \
    libgdal-dev \
    gdal-bin \
    python3-gdal && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install MongoDB 6.0 + mongosh (new shell)
RUN curl -fsSL https://pgp.mongodb.com/server-6.0.asc | gpg --dearmor -o /usr/share/keyrings/mongodb-server-6.0.gpg && \
    echo "deb [ signed-by=/usr/share/keyrings/mongodb-server-6.0.gpg ] https://repo.mongodb.org/apt/ubuntu jammy/mongodb-org/6.0 multiverse" \
      | tee /etc/apt/sources.list.d/mongodb-org-6.0.list
RUN apt-get update
RUN apt-get install -y mongodb-org mongodb-mongosh liblzma-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
    
# Pre-initialize MongoDB users using mongosh
RUN mkdir -p /data/db /var/run/mongodb && \
    mongod --fork --logpath /var/log/mongodb.log --dbpath /data/db && \
    mongosh admin --eval 'db.createUser({ user: "Admin", pwd: "Admin", roles: [ { role: "userAdminAnyDatabase", db: "admin" }, "readWriteAnyDatabase" ] })' && \
    mongosh dbhera --eval 'db.createUser({ user: "user", pwd: "1234", roles: [ { role: "readWrite", db: "dbhera" } ] })'
    # mongod --shutdown
    
# Set working directory and copy project files (exclude .venv, .git via .dockerignore)
WORKDIR /app
COPY requirements.txt .

# Install Python dependencies with pyenv's Python
RUN python -m pip install --no-cache-dir -r requirements.txt

COPY . /app

# Set PYTHONPATH
ENV PYTHONPATH="/app"

# Create necessary folders and configuration file
RUN mkdir -p /root/.pyhera/log && \
    mkdir -p /root/mongo-db-datadir && \
    echo '{ \
        "root": { \
            "dbIP": "host.docker.internal", \
            "dbName": "dbhera", \
            "password": "1234", \
            "username": "user" \
        } \
    }' > /root/.pyhera/config.json

EXPOSE 27017
CMD ["mongod", "--bind_ip_all", "--dbpath", "/data/db"]
