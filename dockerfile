# Use Ubuntu 22.04 as the base image
FROM ubuntu:22.04

# Set environment variables to avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages for pyenv and building Python
RUN apt-get update
RUN apt-get install -y \
    nano curl git build-essential iputils-ping
RUN apt-get install -y \
    libreadline-dev libssl-dev libbz2-dev libsqlite3-dev \
    libcairo2-dev \
    pkg-config \
    libgirepository1.0-dev \
    libffi-dev \
    zlib1g-dev
RUN apt-get install -y \
    xz-utils

# Install MongoDB 6.0 + mongosh (new shell)
RUN curl -fsSL https://pgp.mongodb.com/server-6.0.asc | gpg --dearmor -o /usr/share/keyrings/mongodb-server-6.0.gpg && \
    echo "deb [ signed-by=/usr/share/keyrings/mongodb-server-6.0.gpg ] https://repo.mongodb.org/apt/ubuntu jammy/mongodb-org/6.0 multiverse" \
      | tee /etc/apt/sources.list.d/mongodb-org-6.0.list
RUN apt-get update
RUN apt-get install -y \
    mongodb-org \
    mongodb-mongosh \
    liblzma-dev
    
# Install pyenv
RUN curl https://pyenv.run | bash
    
# Install Python 3.9.13 using pyenv
ENV PATH="/root/.pyenv/bin:/root/.pyenv/shims:${PATH}"
RUN pyenv install 3.9.13 && \
    pyenv global 3.9.13
    
# Install pip for the installed Python version
RUN python -m ensurepip && \
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

# Pre-initialize MongoDB users using mongosh, and close it - will start with container
RUN mkdir -p /data/db /var/run/mongodb && \
    mongod --fork --logpath /var/log/mongodb.log --dbpath /data/db && \
    mongosh admin --eval 'db.createUser({ user: "Admin", pwd: "Admin", roles: [ { role: "userAdminAnyDatabase", db: "admin" }, "readWriteAnyDatabase" ] })' && \
    mongosh admin --eval 'db.createUser({ user: "user", pwd: "1234", roles: [ { role: "readWrite", db: "dbhera" } ] })' && \
    mongod --shutdown
    
# Set working directory and copy project files (exclude .venv, .git via .dockerignore)
WORKDIR /app
COPY requirements.txt .
COPY ui/server/requirements-server.txt .

# Install Python dependencies with pyenv's Python
RUN python -m pip install --no-cache-dir -r requirements.txt
RUN python -m pip install --no-cache-dir -r requirements-server.txt

ENV PATH="/app:/app/hera/bin:${PATH}"
ENV PYTHONPATH="/app:/app/hera/bin"

# Create necessary folders and configuration file
RUN mkdir -p /root/mongo-db-datadir && \
# RUN mkdir -p /root/.pyhera/log && \
#     echo '{ \
#         "root": { \
#             "dbIP": "127.0.0.1", \
#             # "dbIP": "host.docker.internal", \
#             "dbName": "dbhera", \
#             "username": "user", \
#             "password": "1234" \
#         } \
#     }' > /root/.pyhera/config.json

RUN echo 'mongod --fork --logpath /var/log/mongodb.log --dbpath /data/db' >> /root/.bashrc

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
