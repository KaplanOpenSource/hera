#!/bin/bash

# should run inside docker
mongod --fork --logpath /var/log/mongodb.log --dbpath /data/db
