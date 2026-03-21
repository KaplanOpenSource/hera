# Hera UI

## Using without docker
#### Install
```
pip install -r requirements.txt
```

#### Run
Run this script, be sure to activate your `venv` or `heraenv` before
```sh
hera-ui [--cors]
```
This is a shell script on `hera/bin/` folder. [here](../hera/bin/hera-ui)  
Args:  
- `--cors`: not limiting serving the web page only to local browsing.  
If the machine IP is `4.3.2.1` browse thru http://4.3.2.1:8000  

- `--debug`: allowing you to run vscode debugging of hera [more info](https://code.visualstudio.com/docs/python/python-quick-start)  

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
