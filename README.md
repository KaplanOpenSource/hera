## Hera

## 1. Introduction

## 2. Getting started

### 2.1. Prerequisites

1. Linux OS - Ubuntu 22.04 LTS - System verified only with this version

2. Python 3.9.13 - [python3.9 from Ubuntu packages](https://packages.ubuntu.com/search?keywords=python3.9)

3. latest pip  - [follow instructions](https://packaging.python.org/en/latest/guides/installing-using-linux-tools/#debian-ubuntu)

4. MongoDB version 6.0 - [Download & follow instructions](https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-ubuntu/), and make sure it's running on the default port (27017).

5. Several required ubuntu Packages
   
   sudo apt install libcairo2-dev pkg-config python3-dev libgirepository1.0-dev

### 2.2. Installation method
### User Level:
Proceed to the installation section.

### Virtual Environment:
Setup a virtual environment within hera folder and activate it and proceed to the  installation section.
export HERAENVPATH=~/heraenv
mkdir $HERAENVPATH
virtualenv $HERAENVPATH
source $HERAENVPATH/bin/activate

### 2.3. Installation

### The minimum requirements for running Hera include:

```python
setuptools
numpy
matplotlib
pandas
dask[dataframe]
xarray
geopandas
rasterio
mongoengine
seaborn
shapely
scipy
unum
vtk
pyfoam
jinja2
netcdf4
geojson
fastparquet
descartes
```

Proceed with the python requirements installation:

`pip install -r requirements.txt`

If this fails, you may need to reinstall `setuptools` (a library
that is instrumental in package installation):

`pip install --upgrade --force-reinstall setuptools`

and then try the original again:

`pip install -r requirements.txt`

### 2.4. Setup after installation 
In order for the package to work the following steps are required.

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
    "username": {
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

# 3. Additional software for the  hera ecosystem

All the instructions are for Ubuntu OS.

### 3.1 Paraview

Paraview may be use to view the results in a convenient GUI. Paraview my be downloaded from [paraview.org](https://www.paraview.org/download/) and includes python libraries. To prevent conflicts between your python version and Paraview pythons version. make sure to use Paraview with you python. If specific paraview is required, it is recommended to manually download and compile the same python version and install hera in it.

Add Paraview libs to PYTHONPATH - Example
```
export PYTHONPATH=/raid/software/ParaView-5.11.0-MPI-Linux-Python3.9-x86_64/lib/python3.9/site-packages/:$PYTHONPATH

```

### 3.2 FreeCad

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

### 3.3 OPENFOAM.ORG
```
sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key > /etc/apt/trusted.gpg.d/openfoam.asc"
sudo add-apt-repository http://dl.openfoam.org/ubuntu
sudo apt-get -y install openfoam10

echo  ". /opt/openfoam10/etc/bashrc" > of10 # use source of10 to setup OpenFaom environemnt or add to .bashrc
```
