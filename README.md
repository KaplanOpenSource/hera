## Hera

## 1. Introduction

## 2. Getting started

### 2.1. Prerequisites

1. Linux OS (currently checked: Ubuntu 20.04)

2. Python 3.8 - [python3.8 from Ubuntu packages](https://packages.ubuntu.com/search?keywords=python3.8)

3. pip 22.1.2 - [follow instructions](https://packaging.python.org/en/latest/guides/installing-using-linux-tools/#debian-ubuntu)

4. MongoDB version 5.0 - [Download & follow instructions](https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-ubuntu/), and make sure it's running on the default port (27017).


### 2.2. Installation method
### User Level:
Proceed to the installation section.

### Virtual Environment:
Setup a virtual environment within hera folder and activate it and proceed to the  installation section.

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


### 2.4. Setup after installation 
In order for the package to work the following steps are required.

Create the following empty folders within the home folder

`.pyhera/`

`.pyhera/log/`

Create the following json file within .pyhera folder

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
TODO This shold be tested to see if any choice would work.

* username - should match the name of your user in the ubuntu system.
* {host} - should be changed to the location of mongoDB, if ran locally it is typically "127.0.0.1"
* {dbName} - name of database 
* {password} - choose a password 
* {username} - user name

Afterwards save and exit.

Startup mongoDB `mongosh`

Type in the following commands to enter both users:

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
    user: "{username}",
    pwd:  "{password}",   
    roles: [ { role: "readWrite", db: "{dbName}" } ]

  }
)
```
TODO: TEST/ASK if both configurations are necessary
* {username} - should be replaced with the same username written in config.json (above)
* {password} - should be replaced with the same password written in config.json (above)
* {dbName} - should be replaced with the same dbName written in config.json (above)

Exit mongosh