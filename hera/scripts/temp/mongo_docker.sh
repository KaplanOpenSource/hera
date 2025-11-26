docker run --rm --name mongo-hera \
  -p 27017:27017 \
  -v $(pwd)/projects/mongo-data:/data/db \
  -e MONGO_INITDB_ROOT_USERNAME=Admin \
  -e MONGO_INITDB_ROOT_PASSWORD=Admin \
  mongo:latest \
  bash -c "\
    mongod --fork --logpath /var/log/mongodb.log --dbpath /data/db && \
    until mongosh --quiet --eval 'db.runCommand({ping:1})' >/dev/null 2>&1; do sleep 1; done && \
    mongosh --eval 'db.createUser({user:\"user\",pwd:\"1234\",roles:[{role:\"readWrite\",db:\"dbhera\"}]})' && \
    mongod --shutdown --dbpath /data/db && \
    exec mongod --bind_ip_all"

    # mongosh --eval 'db.createUser({user:\"Admin\",pwd:\"Admin\",roles:[{role:\"userAdminAnyDatabase\",db:\"admin\"},\"readWriteAnyDatabase\"]})' && \