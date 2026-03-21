# Hera

[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue?logo=github)](https://KaplanOpenSource.github.io/hera/)
[![Live Demo](https://img.shields.io/badge/Live%20Demo-GitHub%20Pages-green?logo=github)](https://kaplanopensource.github.io)
[![Open Source](https://img.shields.io/badge/Open%20Source-Kaplan-orange)](https://kaplanopensource.co.il/)

**Hera** is an advanced open-source project by [Kaplan Open Source Consulting](https://kaplanopensource.co.il/) focused on web-based GIS systems and data processing. It serves as a framework for managing complex geographical data and interactive map visualizations.

---

## Live Access
You can view the live deployment of this project here:
**[https://kaplanopensource.github.io](https://kaplanopensource.github.io)**

---

## Key Features
*   **GIS Integration:** Built-in support for OpenStreetMap and custom geographic data layers.
*   **Django Framework:** Robust backend architecture designed for scalability.
*   **Data Analysis:** Flexible pipelines for processing and visualizing spatial information.
*   **Interactive Maps:** Lightweight, mobile-friendly interactive map interfaces.

---

## Getting Started

### Prerequisites
*   **OS:** Ubuntu 22.04 LTS (verified)
*   **Python:** 3.9.13+
*   **Database:** MongoDB version 6.0 (running on default port 27017)
*   **System packages:**
    ```bash
    sudo apt install libcairo2-dev pkg-config python3-dev libgirepository1.0-dev libgdal-dev gdal-bin python3-gdal
    ```

### Installation
1. **Clone the repository:**
   ```bash
   git clone https://github.com/KaplanOpenSource/hera
   cd hera
   ```

2. **Set up a virtual environment:**
   ```bash
   python3.9 -m venv heraenv
   source heraenv/bin/activate
   ```

3. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

   If this fails, you may need to reinstall `setuptools`:
   ```bash
   pip install --upgrade --force-reinstall setuptools
   pip install -r requirements.txt
   ```

4. **Install GDAL:**
   Use `gdalinfo --version` to obtain your OS GDAL version, then install the matching Python binding:
   ```bash
   pip install GDAL==`gdal-config --version`
   ```
   If GDAL is not installed:
   ```bash
   sudo apt-get install -y libgdal-dev gdal-bin python3-gdal
   ```

### Setup after installation

Run the environment setup script from the project root:
```bash
source set_hera_environment.sh
```

This sets `HERA_REPO_ROOT`, `PYHERA_DIR`, and `PYTHONPATH` for the current session. The script will also ask if you want to add it to your `~/.bashrc` so that the environment loads automatically on every new shell.

Create the required directories:
```bash
mkdir -p ~/.pyhera/log/
```

### Hera configuration files

Create the following JSON file within `.pyhera` folder. The file contains the address and credentials for MongoDB. If not created manually, it will be created at first import of hera, but without values, so the import will fail.

`.pyhera/config.json`

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

* `{username}` - should match the name of your user in the Ubuntu system
* `{host}` - location of MongoDB, if local it is typically `127.0.0.1`
* `{dbName}` - name of database
* `{password}` - choose a password

### MongoDB Schema

Start MongoDB with `mongosh` and run the following commands to create admin and regular users:

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

Replace `{username}`, `{password}`, and `{dbName}` with the same values from `config.json`.

### Use predefined names (with Makefile)

This is especially convenient if you don't have MongoDB installed, as it uses Docker (you need Docker installed).

A `Makefile` is provided in the project root to manage MongoDB, the Hera server, and tests. Run `make help` to see all available targets.

**Start MongoDB:**
```bash
make mongo-up
```

This starts a MongoDB 5.0 container with data stored at `~/mongo-db-datadir`. It automatically creates a user "hera" with password "heracles" using the init scripts in `mongo-init.d/`.

To override the data directory:
```bash
make mongo-up MONGO_DATA=/path/to/data
```

**Stop / restart MongoDB:**
```bash
make mongo-down
make mongo-up
```

**Check status or view logs:**
```bash
make mongo-status
make mongo-logs
```

**Remove container and delete all data:**
```bash
make mongo-clean
```

**Set up the test data directory:**
```bash
make test-setup
```

This creates the directory structure at `~/hera_unittest_data` needed for running the test suite. Override with `make test-setup TEST_HERA=/path/to/data`.

**Run tests:**
```bash
make test
```

For the predefined setup, create `.pyhera/config.json` with:
```json
{
    "<username>": {
        "dbIP": "127.0.0.1",
        "dbName": "olymp",
        "password": "heracles",
        "username": "hera"
    }
}
```
Replace `<username>` with your system username.

### Installing third-party dependencies

The Makefile also provides targets for installing third-party dependencies:

```bash
# Install system packages (libcairo, GDAL, etc.) and the GDAL Python binding
make install-deps

# Install individual third-party tools
make install-paraview    # Download and set up ParaView 5.11.0
make install-freecad     # Install FreeCad Python3 bindings
make install-openfoam    # Install OpenFOAM 10

# Install everything at once
make install-deps-all
```

Run `make help` to see all available targets.

---

## Documentation

Full documentation is available at **[https://KaplanOpenSource.github.io/hera/](https://KaplanOpenSource.github.io/hera/)**.

The documentation covers architecture, toolkits, testing, CLI reference, examples, and more.

### Viewing Documentation Locally

```bash
# Activate your virtual environment
source heraenv/bin/activate

# Install documentation dependencies (only needed once)
pip install -r requirements.txt

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

The documentation is automatically deployed to GitHub Pages whenever changes are pushed to the `main` branch.

- A GitHub Actions workflow (`.github/workflows/docs.yml`) monitors the `main` branch
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
