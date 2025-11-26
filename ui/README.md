# Hera UI

## Build
```sh
sh hera/scripts/docker_build.sh
```

## Server
#### Run Mongo DB
You can use docker:
```sh
sh hera/scripts/docker_mongo.sh
```
#### Run Hera Server
```sh
sh hera/scripts/docker_run.sh
```
#### Browse
http://localhost:8000

## Debug

#### Server
```sh
docker run -it --network host -v .:/app --rm --name hera-server-instance hera-server python ui/server/server.py --debug
```

#### Client
Needs nodejs 22+
```sh
cd ui/client
npm install
npm run dev
```
