docker run --name hera-mongo \
  -v ./mongo-db-datadir:/data/db \
  -v ./mongo-init.d:/docker-entrypoint-initdb.d \
  -p 127.0.0.1:27017:27017 \
  -d mongo:5.0