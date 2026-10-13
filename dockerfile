# Ubuntu 26.04 LTS with its native Python (3.14) -- issue #1125.
FROM ubuntu:26.04

# Set environment variables to avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Base tools, build deps for the few sdists, the system Python, GDAL and Node.
# 26.04 ships Node 22 and Python 3.14 itself, so no PPAs are needed any more.
RUN apt-get update && apt-get install -y \
    nano curl git build-essential iputils-ping \
    libcairo2-dev pkg-config libgirepository1.0-dev \
    python3 python3-venv python3-dev \
    libgdal-dev gdal-bin libglu1-mesa \
    nodejs npm \
 && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install into an isolated venv: PEP 668 forbids pip into the system Python.
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="/opt/venv/bin:${PATH}"
RUN python3 -m venv "$VIRTUAL_ENV" && \
    python -m pip install --no-cache-dir --upgrade pip setuptools wheel

# Set working directory and copy project files (exclude .venv, .git via .dockerignore)
WORKDIR /app
COPY requirements.txt .

# Install Python dependencies, then the GDAL binding built against the system
# libgdal (same recipe as `make install-deps`; numpy must already be present).
RUN python -m pip install --no-cache-dir -r requirements.txt && \
    python -m pip install --no-cache-dir --no-build-isolation GDAL=="$(gdal-config --version)"

# ParaView for the OpenFOAM VTK pipeline: same recipe as `make install-paraview`.
# python3-paraview is built for the system Python against the system VTK, so the
# venv links those in and drops pip's vtk.  Re-running `pip install -r requirements.txt`
# in the image brings pip's vtk back over the links; repeat this step after it.
RUN apt-get update && apt-get install -y --no-install-recommends python3-paraview && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    python -m pip uninstall -y vtk && \
    SP="$(python -c "import sysconfig; print(sysconfig.get_paths()['purelib'])")" && \
    dpkg -L python3-paraview python3-vtk9 | grep -E '^/usr/lib/python3/dist-packages/[^/]+$' | sort -u | xargs -I{} ln -sf {} "$SP/" && \
    python -c "import paraview.simple"

# Disable the jupytext lab extension; it is built for JupyterLab 4.4+ and errors on 4.3.
# RUN mkdir -p "$VIRTUAL_ENV/share/jupyter/lab/settings" && \
#     printf '{"disabledExtensions": {"jupyterlab-jupytext": true}}\n' \
#       > "$VIRTUAL_ENV/share/jupyter/lab/settings/page_config.json"

ENV PATH="/app:/app/hera/bin:${PATH}"
ENV PYTHONPATH="/app:/app/hera/bin"

# The bind-mounted repo is owned by the host user, but the container runs as
# root. Mark it safe so git (used by jupyterlab-git) does not refuse to run.
RUN git config --global --add safe.directory '*'

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
