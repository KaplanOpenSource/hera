# Hera UI

## Using without docker
Run this script, be sure to activate your `venv` or `heraenv` before
```sh
sh hera/scripts/run_ui.sh
```

## Using with docker
#### Build
```sh
sh hera/scripts/docker_build.sh
```
#### Run Hera Server
```sh
sh hera/scripts/docker_run.sh
```
#### Run Mongo DB (optional)
If you dont have mongo DB yet, you can use docker:
```sh
sh hera/scripts/docker_mongo.sh
```
#### Browse
The web UI will be in this URL:  
http://localhost:8000

## Debug with docker
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
