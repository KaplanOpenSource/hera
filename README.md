## Hera

## 1. Introduction

## 2. Install
### 2.1. Installation of virtual environment
`python3 -m venv ./venv`

if the virtual environment package doesn't exist (on clean systems) run the following command to install it

`apt install python3.8-venv`

activating the envrironment:

`source venv/bin/activate`

Test if the pip is in the right path:

`which pip`


Install from the requirements list:

`pip install -r requirements.txt`

### 2.2. Install GDAL

if you are using a clean Ubuntu 20.4 version you shoudl run the following commands:

Install python development tools:
`apt-get install python3-dev`

```
sudo add-apt-repository ppa:ubuntugis/ppa && sudo apt-get update
sudo apt-get update
sudo apt-get install gdal-bin
sudo apt-get install libgdal-dev
export CPLUS_INCLUDE_PATH=/usr/include/gdal
export C_INCLUDE_PATH=/usr/include/gdal
pip install GDAL==3.0.4
```


### 2.3. MiniConda installation

To Be Added

### 2.4. Setup Hera

`python setup.py install`

Additional Installations: (not tested yet)

`conda install --offline {path}`

gdal-3.3.1-py36h77b1db5_3.tar.bz2

icu-68.1-h58526e2_0.tar.bz2

nodejs-15.11.0-h92b4a50_0.tar.bz2

### 2.5. Mongo DB Installation

the installation should follow the standard installation:

https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-ubuntu/

if installed in Ubuntu, the pacakged provided by Ubuntu must be uninstalled first

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

run the following command:

!!!NEED to be tested in a clean environment

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

### 2.7. Setting up config.json

!!! NEED TO SEE IN WHAT CONDITION CONFIG.JSON is CREATED

### 2.8. setting log folder
creating log folder inside .pyhera