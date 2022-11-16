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


#### 2.4.1 Manual

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

#### 2.4.2 Use predefined names

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
