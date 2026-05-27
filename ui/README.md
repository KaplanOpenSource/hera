# Hera UI

## Using without docker
#### Install
```
pip install -r requirements.txt
```

#### Run
Run this script, be sure to activate your `venv` or `heraenv` before
```sh
hera-ui [options]
```
This is a shell script on `hera/bin/` folder. [here](../hera/bin/hera-ui)

Use `--help` to see all available options.

Options:
- `--cors ORIGINS`: Enable CORS for external origins.  
Use `all` to allow all origins, or pass a comma-separated list of IPs to allow specific ones.  
Example: `--cors 192.168.1.10,10.0.0.5`  
Each IP is prefixed with `http://` and the server port automatically.  
When active, a red "CORS" indicator appears in the top-right of the UI.

- `--port PORT`: Port for the API server (default: `8000`).

- `-y, --yes`: Skip confirmation prompts. Useful for non-interactive environments (e.g. VS Code launch).

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
http://localhost:8000 (default port, configurable with `--port`)

## Notebooks

Documents can be created as Jupyter notebooks. When adding a document, check the "Notebook" checkbox — this creates an `.ipynb` file under the project's `filesDirectory/notebooks/` and opens it in an embedded Jupyter editor. If a notebook file already exists at that path, it will be used as-is.

### Using AI in notebook cells

Notebooks include [Jupyter AI](https://jupyter-ai.readthedocs.io/) with Ollama integration. Use the `%%ai` cell magic to query LLMs directly from notebook cells.

**Text response** — the LLM output appears as text below the cell:

```
%%ai ollama:llama3
how's the weather in Israel?
```

**Code response** — use `-f code` to have the LLM output inserted as a new code cell below, which you can then execute:

```
%%ai ollama:llama3 -f code
Write a python function that calculates wind speed from u and v components then run it with some reasonable params
```

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
