# Docker lockdown plan

Goal: make `docker-compose.yml` safer, like the example from the other project.

## What changes

1. No code folders mounted. Code is copied into the image. `docker compose up --watch` copies edits in.
   - Python code (`hera/`, `ui/server/`, `Hermes/`): copy, then restart the server (it has no auto reload).
   - Client code (`ui/client/`): copy only. Vite reloads it.
   - `requirements.txt` or `package-lock.json` change: rebuild the image.
2. Data folders stay mounted from the repo: `mongo-db-datadir`, `ollama-data`, `projects`, `notebooks`.
3. No root. Mongo runs as `999` (its files already belong to 999). Server, client and Ollama run as `1000`, the host user.
4. All services drop all extra Linux rights (`cap_drop: [ALL]`, `no-new-privileges`).
5. Shared settings are written once with YAML anchors.
6. `mongo-init.d` is passed as a config file, not a mounted folder.

## Dockerfile

All changes come after the `pip install` layer, so apt and pip do not run again.
- Add user `hera` (uid 1000) with its own `~/.pyhera/config.json`.
- `npm install` runs as `hera`, so Vite can write its cache.
- Copy the code into `/app` at the end.
- The image still ends as root, so `test:server` and the integration test image work as before.

## One-time step on the host

`projects`, `notebooks` and `ollama-data` belong to root today. Give them to the host user once:

    sudo chown -R 1000:1000 projects notebooks ollama-data

## Not tested here

There is no Docker in the agent's shell. Run `docker compose up --build --watch` to check.
