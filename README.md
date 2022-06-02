## Hera

## 1. Introduction

## 2. Getting started

### 2.1. Prerequisites

1. Linux OS (currently checked: Ubuntu 20.04)

2. python 3.8 - [python3.8 from Ubuntu packages](https://packages.ubuntu.com/search?keywords=python3.8)

3. MongoDB version 5.0 - [Download & follow instructions](https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-ubuntu/), and make sure it's running on the default port (27017).


### 2.2. Installation method
### User Level:
Proceed to the installation section.

### Virtual Environment:
Setup a virtual environment within hera folder and activate it and proceed to the  installation section.

### 2.3. Installation

### The minimum requirements for running Hera include:

```python
testresources
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


### 2.4. Setup after installation 


### 2.5. Mongo DB setup

the installation should follow the standard installation:

https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-ubuntu/

if installed in Ubuntu, the preinstalled pacakged in Ubuntu (clean installation) must be uninstalled first

1.First stop the MongoDB Process

`sudo service mongod stop`

2.Completely remove the installed MongoDB packages.

`sudo apt-get purge mongodb-org*`

3.Remove the data directories, MongoDB database(s), and log files.

`sudo rm -r /var/log/mongodb /var/lib/mongodb`

4.To verify that MongoDB has been successfully uninstalled, type the command below.

`service mongod status`

Afterwards Install MongoDB to the instructions in the link above.

Need to add - How to setup MongoDB after installation.

### 2.6. Setting up MongoDB Database

Start mongo in the command line
`mongosh`

Create two databases and change the user name and password to anything in your preference

Follow the the code bellow as an example.

```JavaScript
use admin

db.createUser(
  {
    user: "MathAdmin",
    pwd: "MathAdmin",
    roles: [ { role: "userAdminAnyDatabase", db: "admin" } , "readWriteAnyDatabase"]
  }
)

use admin
db.createUser(
  {
    user: "shai",
    pwd:  "shai",   
    roles: [ { role: "readWrite", db: "shai" } ]

  }
)
```

### 2.7. Setting up config.json (THIS FILE IS GENERATED AFTER RUNNING HERA ONCE)

The config file should be filled with similar parameters filled in the mongoDB database
```JavaScript
{
    "shai": {
        "dbIP": "127.0.0.1",
        "dbName": "shai",
        "password": "shai",
        "username": "shai"
    }
}
```

### 2.8. setting log folder
creating log folder inside .pyhera