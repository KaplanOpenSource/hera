# Hera

[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue?logo=github)](https://KaplanOpenSource.github.io/hera/)
[![Live Demo](https://img.shields.io/badge/Live%20Demo-GitHub%20Pages-green?logo=github)](https://kaplanopensource.github.io)
[![Open Source](https://img.shields.io/badge/Open%20Source-Kaplan-orange)](https://kaplanopensource.co.il/)

**Hera** is a Python scientific data management platform by [Kaplan Open Source Consulting](https://kaplanopensource.co.il/). It provides a unified MongoDB-backed data layer and a set of domain-specific **Toolkits** for GIS, meteorology, atmospheric dispersion (CFD + Lagrangian particle tracking), and risk assessment.

**Stack:** Python 3 · MongoDB (via mongoengine) · pandas · dask · geopandas · xarray · pint · OpenFOAM (optional)

---

## Key Features
*   **Unified data layer:** Three MongoDB collections (Measurements, Simulations, Cache) with a consistent document API across all toolkits.
*   **GIS toolkits:** Topography (SRTM), land cover, vector layers, buildings, demography — powered by geopandas and GDAL.
*   **Meteorology toolkits:** Low-frequency and high-frequency station data, turbulence statistics, WRF output ingestion.
*   **Simulation toolkits:** OpenFOAM (Eulerian + Lagrangian), LSM particle tracking, Gaussian dispersion, wind profiles.
*   **Risk assessment:** Configurable injury/effect models, protection policies, casualty estimation with spatial output.
*   **Luigi workflows:** `hermesWorkflowToolkit` for DAG-based simulation pipelines on HPC (Slurm).

---

## Getting Started

### Prerequisites
*   **OS:** Ubuntu 22.04 LTS (verified)
*   **Python:** 3.9.13+
*   **Docker:** Required for quick install (MongoDB runs in a container)

### Quick Install (recommended)

The fastest way to get Hera running. Requires Docker and Python 3.9+.

1. **Clone the repository:**
   ```bash
   git clone https://github.com/KaplanOpenSource/hera
   cd hera
   ```

2. **Run the initialization script:**
   ```bash
   source init_with_mongo.sh
   ```
   This single script will:
   - Set up Hera environment variables (`HERA_REPO_ROOT`, `PYHERA_DIR`, `PYTHONPATH`)
   - Prompt to create a Python virtual environment at `~/.pyhera/environment` with all dependencies installed
   - Prompt to persist the environment to `~/.bashrc`
   - Prompt to auto-activate the Hera environment when you `cd` into the project directory
   - Pull and start a MongoDB 5.0 Docker container with preconfigured credentials
   - Create the `~/.pyhera/config.json` configuration file automatically

   After running, everything is ready to use.

3. **Install system dependencies and third-party tools (optional):**
   ```bash
   make install-deps       # System packages (libcairo, GDAL, etc.)
   make install-paraview   # ParaView 5.11.0
   make install-freecad    # FreeCad Python3 bindings
   make install-openfoam   # OpenFOAM 10
   make install-deps-all   # All of the above
   ```

#### Helper scripts

| Script | Purpose |
|---|---|
| `source set_hera_environment.sh` | Set environment variables, create venv, configure `.bashrc` (called automatically by `init_with_mongo.sh`) |
| `source activate_hera.sh` | Activate the Hera venv and environment in the current shell |
| `init_with_mongo.sh` | Full initialization: environment setup + MongoDB via Docker |

#### Managing MongoDB with Make

After the initial setup, use the Makefile to manage the MongoDB container:

```bash
make mongo-up              # Start MongoDB (data at ~/mongo-db-datadir)
make mongo-down            # Stop MongoDB
make mongo-status          # Check container status
make mongo-logs            # Tail container logs
make mongo-clean           # Remove container and delete all data
make mongo-up MONGO_DATA=/path/to/data   # Override data directory
```

#### Running tests

```bash
make test-setup            # Create test data directory structure at ~/hera_unittest_data
make test                  # Run all tests (requires MongoDB + test data)
```

Run `make help` to see all available targets.

---

### Manual Installation

If you prefer not to use Docker or need more control over each step, follow these instructions.

#### 1. System packages

```bash
sudo apt install libcairo2-dev pkg-config python3-dev libgirepository1.0-dev libgdal-dev gdal-bin python3-gdal
```

#### 2. Clone and set up Python environment

```bash
git clone https://github.com/KaplanOpenSource/hera
cd hera
python3.9 -m venv heraenv
source heraenv/bin/activate
```

#### 3. Install dependencies

```bash
pip install -r requirements.txt
```

If this fails, reinstall `setuptools` first:
```bash
pip install --upgrade --force-reinstall setuptools
pip install -r requirements.txt
```

#### 4. Install GDAL Python binding

Use `gdalinfo --version` to obtain your OS GDAL version, then:
```bash
pip install GDAL==`gdal-config --version`
```

#### 5. Configure environment variables

Add the following to your shell profile or virtual environment activate script:
```bash
export PYTHONPATH=$PYTHONPATH:/path/to/hera
export PATH=$PATH:/path/to/hera/hera/bin
```

Create required directories:
```bash
mkdir -p ~/.pyhera/log/
```

#### 6. MongoDB setup

Install MongoDB 6.0 following the [official instructions](https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-ubuntu/) and ensure it runs on port 27017.

Start MongoDB with `mongosh` and create users:

```javascript
use admin

db.createUser(
  {
    user: "Admin",
    pwd: "Admin",
    roles: [ { role: "userAdminAnyDatabase", db: "admin" } , "readWriteAnyDatabase"]
  }
)

use admin
db.createUser(
  {
    user: "{username}",
    pwd:  "{password}",
    roles: [ { role: "readWrite", db: "{dbName}" } ]
  }
)
```

#### 7. Hera configuration file

Create `~/.pyhera/config.json`:

```json
{
    "{username}": {
        "dbIP": "{host}",
        "dbName": "{database name}",
        "password": "{password}",
        "username": "{username}"
    }
}
```

* `{username}` - should match your Ubuntu system user
* `{host}` - MongoDB address, typically `127.0.0.1` for local
* `{dbName}` - name of database
* `{password}` - choose a password

---

## Documentation

Full documentation is available at **[https://KaplanOpenSource.github.io/hera/](https://KaplanOpenSource.github.io/hera/)**.

The documentation covers architecture, toolkits, testing, CLI reference, examples, and more.

### Viewing Documentation Locally

```bash
# Activate the Hera environment (if not already active)
source activate_hera.sh

# Start the local development server with live reload
mkdocs serve

# The site will be available at http://127.0.0.1:8000
```

To build a static version of the site:
```bash
mkdocs build

# Build with strict mode (catches broken links and warnings)
mkdocs build --strict
```

### Automated Documentation Deployment

The documentation is automatically deployed to GitHub Pages whenever changes are pushed to the `master` branch.

- A GitHub Actions workflow (`.github/workflows/docs.yml`) monitors the `master` branch
- When changes are detected in `docs/`, `mkdocs.yml`, or `hera/` code, it builds and deploys automatically
- The live site updates within 1-2 minutes after merging to `main`

**Manual deployment** (for testing from other branches):
```bash
mkdocs gh-deploy --force
```

> **Note:** All dependencies (including documentation) are consolidated in `requirements.txt`. Documentation uses [MkDocs Material](https://squidfunk.github.io/mkdocs-material/).

---

## Hera UI
Refer to [ui/README.md](https://github.com/KaplanOpenSource/hera/blob/master/ui/README.md)

---

## Additional software for the Hera ecosystem

All instructions are for Ubuntu OS.

### Paraview

Paraview may be used to view results in a convenient GUI. Download from [paraview.org](https://www.paraview.org/download/). To prevent conflicts between your Python version and Paraview's Python version, make sure to use Paraview with your Python.

Add Paraview libs to PYTHONPATH:
```bash
export PYTHONPATH=/raid/software/ParaView-5.11.0-MPI-Linux-Python3.9-x86_64/lib/python3.9/site-packages/:$PYTHONPATH
```

### FreeCad

FreeCad is an open source CAD software that can be embedded in Python.

Install the `freecad-python3` package:
```bash
sudo apt-get install libfreecad-python3-0.19
```

Then add the library path to PYTHONPATH:
```python
FREECADPATH = '/usr/lib/freecad-python3/lib/'
import sys
sys.path.append(FREECADPATH)
```

[More information on embedding FreeCAD](https://wiki.freecad.org/Embedding_FreeCAD)

### OpenFOAM

```bash
sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key > /etc/apt/trusted.gpg.d/openfoam.asc"
sudo add-apt-repository http://dl.openfoam.org/ubuntu
sudo apt-get -y install openfoam10

echo  ". /opt/openfoam10/etc/bashrc" > of10 # use source of10 to setup OpenFOAM environment
```
