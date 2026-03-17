## Hera

[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue?logo=github)](https://KaplanOpenSource.github.io/hera/)

**Hera** is a Python-based platform for managing scientific data across measurements, simulations, and cached results.

### Documentation

Full documentation is available at **[https://KaplanOpenSource.github.io/hera/](https://KaplanOpenSource.github.io/hera/)**.

The documentation covers architecture, toolkits, testing, CLI reference, examples, and more.

#### Viewing Documentation Locally

To preview the documentation site on your local machine:

```bash
# Activate your virtual environment (if not already active)
source heraenv/bin/activate

# Install documentation dependencies (only needed once)
pip install -r docs/requirements-docs.txt

# Start the local development server with live reload
mkdocs serve

# The site will be available at http://127.0.0.1:8000
# Any changes you make to docs/ will automatically refresh in the browser
```

To build a static version of the site:

```bash
# Build static site into site/ directory
mkdocs build

# Build with strict mode (catches broken links and warnings)
mkdocs build --strict
```

#### Automated Documentation Deployment

The documentation is **automatically deployed** to GitHub Pages whenever changes are pushed to the `main` branch.

**How it works:**
- A GitHub Actions workflow (`.github/workflows/docs.yml`) monitors the `main` branch
- When changes are detected in `docs/`, `mkdocs.yml`, or `hera/` code, it:
  1. Builds the documentation site using MkDocs
  2. Validates all links and Mermaid diagrams with `--strict` mode
  3. Deploys the site to GitHub Pages automatically
- The live site updates within 1-2 minutes after merging to `main`

**Manual deployment** (for testing from other branches):
```bash
mkdocs gh-deploy --force
```

> **Note:** The legacy Sphinx-based documentation configuration in `requirements-doc.txt` is no longer maintained. All documentation now uses [MkDocs Material](https://squidfunk.github.io/mkdocs-material/).

---

## 1. Introduction

## 2. Getting started

### 2.1. Prerequisites

1. Linux OS - Ubuntu 22.04 LTS - System verified only with this version

2. Python 3.9.13 - [python3.9 from Ubuntu packages](https://packages.ubuntu.com/search?keywords=python3.9)

3. latest pip  - [follow instructions](https://packaging.python.org/en/latest/guides/installing-using-linux-tools/#debian-ubuntu)

4. MongoDB version 6.0 - [Download & follow instructions](https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-ubuntu/), and make sure it's running on the default port (27017).

5. Several required ubuntu Packages
   
   sudo apt install libcairo2-dev pkg-config python3-dev libgirepository1.0-dev libgdal-dev gdal-bin python3-gdal

### 2.2. Installation method
### User Level:
Proceed to the installation section.

### Clone the Git
Use the following command to clone the git:  
```console
$ git clone https://github.com/KaplanOpenSource/hera
```

### Virtual Environment:
Setup a virtual environment within hera folder and activate it and proceed to the  installation section.  
```console
python3.9 -m venv heraenv
```

### 2.3. Installation

### The minimum requirements for running Hera include:

```python
testresources==2.0.1
numpy==1.21.5
matplotlib==3.5.1
pandas==1.3.0
dask[dataframe]==2021.2.0
xarray==0.16.2
geopandas==0.8.2
rasterio==1.2.6
mongoengine==0.22.1
seaborn==0.11.1
shapely==1.7.1
scipy==1.6.0
unum==4.1.4
vtk==9.1.0
pyfoam==2020.5
jinja2==3.0.1
netcdf4==1.5.5.1
geojson==2.5.0
fastparquet==0.8.1
descartes==1.1.0
pytest==7.2.0
exceptiongroup==1.0.4
iniconfig==1.1.1
pluggy==1.0.0
tomli==2.0.1
```

Proceed with the python requirements installation:

`pip install -r requirements.txt`

If this fails, you may need to reinstall `setuptools` (a library
that is instrumental in package installation):

`pip install --upgrade --force-reinstall setuptools`

and then try the original again:

`pip install -r requirements.txt`

Install GADL
Use the command: "gadlinfo --version" to obtain OS GDAL version and install same version (or as close as possible) via 
pip install GDAL==\`gdal-config --version\` 
if GDAL is not installed, install from repository: 
`sudo apt-get install -y libgdal-dev gdal-bin python3-gdal`


### 2.4. Setup after installation 
In order for the package to work automaticly each time you enter the virtual enviorment, the following steps are required:

Enter the virtual enviorment bin folder (not activate it):
```console
cd HERAENV_PATH/bin
```
Edit the activate script:
```console
nano activate
```
Add two export commands at the end of the file:
```python
export PYTHONPATH=$PYTHONPATH:/home/YOUR_OS_USER_NAME/PATH_TO_HERA_GIT_FOLDER/hera/hera/bin
export PATH=$PATH:/home/YOUR_OS_USER_NAME/PATH_TO_HERA_GIT_FOLDER/hera/hera/bin
```
The paths should be of the bin folder inside the hera folder you cloned before.


We are still showing two alternatives. The first is the manual one,
written up before; the second is the more streamlined one which
uses pre-defined names and passwords.

#### 2.4.0 For both alternatives

Create the following empty folders within the home folder

`.pyhera/`

`.pyhera/log/`

mkdir -p ~/.pyhera/log/


#### 2.4.1 Hera configuration files

Create the following json file within .pyhera folder. The file contains the address and credentails for the MongoDB. If file not created manually, it will be created at first import of hera, but without values, so the import will fails.

`.pyhera/config.json`

```JavaScript
{
    "{username}": {
        "dbIP": "{host}",
        "dbName": "{database name}",
        "password": "{password}",
        "username": "{username}"
    }
}
```

* username - should match the name of your user in the ubuntu system.
* {host} - should be changed to the location of mongoDB, if ran locally it is typically "127.0.0.1"
* {dbName} - name of database 
* {password} - choose a password 
* {username} - user name

Afterwards save and exit.

#### 2.4.2 MongoDB Schema

Startup mongoDB `mongosh`. (More parameters are requird if mongodb is not local or requires authentication).

Type in the following commands to create admin and regular users:

```JavaScript
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

* {username} - should be replaced with the same username written in config.json (above)
* {password} - should be replaced with the same password written in config.json (above)
* {dbName} - should be replaced with the same dbName written in config.json (above)

Exit mongosh

#### 2.4.3 Use predefined names

This is especially convenient if you don't have MongoDB installed, as it
uses docker for it (you need docker installed, though...)

Create an empty directory to be used for MongoDB, e.g. in your home folder:

```console
$ mkdir mongo-db-datadir
```

The next command will start a MongoDB server in a container on your machine,
while setting up the users expected by Hera. You may need to adapt it, see
below.

*TODO* verify username, dbname

```console
$ docker run --name hera-mongo \
  -v ${HOME}/mongo-db-datadir:/data/db \
  -v ${HOME}/hera/mongo-init.d:/docker-entrypoint-initdb.d \
  -p 127.0.0.1:27017:27017 -d mongo:5.0 
```

Note: In the above, `${HOME}` refers to a pre-defined environment
variable. So:

* `${HOME}/mongo-db-datadir` is the data directory assuming you
  created it as above. If you created it somewhere else, adapt
  accordingly.

* `${HOME}/hera/mongo-init.d` assumes you've placed the Hera project
  code in a folder `hera` in your home dir. Again, adapt if this
  is not the case
  
This creates a MongoDB user named "hera" with password "heracles". The
MongoDB server is not accessible from outside your computer, so this
is not a terrible security issue.

Note: Once you've successfully executed the above command, you can
stop and start the mongo server with 

```console
$ docker stop hera-mongo
$ docker start hera-mongo
```

Finally, create the following json file within `.pyhera` folder:

`.pyhera/config.json`

```JavaScript
{
    "<username>": {
        "dbIP": "127.0.0.1",
        "dbName": "olymp",
        "password": "heracles",
        "username": "hera"
    }
}
```
where `<username>` should be replaced by your username on your system.

# 3. Hera UI
Refer to [README.md](https://github.com/KaplanOpenSource/hera/blob/master/ui/README.md)

# 4. Additional software for the  hera ecosystem

All the instructions are for Ubuntu OS.

### 4.1 Paraview

Paraview may be use to view the results in a convenient GUI. Paraview my be downloaded from [paraview.org](https://www.paraview.org/download/) and includes python libraries. To prevent conflicts between your python version and Paraview pythons version. make sure to use Paraview with you python. If specific paraview is required, it is recommended to manually download and compile the same python version and install hera in it.

Add Paraview libs to PYTHONPATH - Example
```
export PYTHONPATH=/raid/software/ParaView-5.11.0-MPI-Linux-Python3.9-x86_64/lib/python3.9/site-packages/:$PYTHONPATH

```

### 4.2 FreeCad

Freecad is an open source CAD software. FreeCad can be embedded in python. 

It is required to install `freecad-python3` pkg (apt). In ubuntu, 

```sudo apt-get install libfreecad-python3-0.19```

Then, add the library path (default:'/usr/lib/freecad-python3/lib/') to PYTHONPATH env or dynamically in the code like:
```python
FREECADPATH = '/usr/lib/freecad-python3/lib/' # path to your FreeCAD.so or FreeCAD.pyd file,

import sys
sys.path.append(FREECADPATH)
```
[more information on embedding freecad in freecad sitep](https://wiki.freecad.org/Embedding_FreeCAD)

### 4.3 OPENFOAM.ORG
```
sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key > /etc/apt/trusted.gpg.d/openfoam.asc"
sudo add-apt-repository http://dl.openfoam.org/ubuntu
sudo apt-get -y install openfoam10

echo  ". /opt/openfoam10/etc/bashrc" > of10 # use source of10 to setup OpenFaom environemnt or add to .bashrc
```
