# Known Integration Issues

Setup problems that break the UI/server and how to work around them.

## Hermes

Running workflows needs the **hermes** package. Missing it fails with the
misleading error: `ValueError: The workflow type None not found`.

This only happens with an outdated Docker image. Running without Docker (a local
env that already has hermes) is not affected.

hermes is a submodule at `Hermes/` with no `setup.py`, so it can't be
pip-installed. The Docker launch script `ui/server/scripts/run_server_docker.sh`
(used by both `.vscode/launch.json` configs) handles it:

- Add `/app/Hermes` to `PYTHONPATH` so `import hermes` works.
- `pip install jsonpath_rw_ext luigi` (its runtime dependencies) before server
  start. Missing `luigi` shows `No module named luigi` and the workflow's tasks
  don't actually run (you still get a dispatch id).

Verify: `docker exec hera-server-instance python -c "import hermes, luigi"`

TODO: permanent fix = add `/app/Hermes` to `dockerfile` PYTHONPATH +
`jsonpath_rw_ext` and `luigi` to `requirements.txt`, then rebuild.

## MongoDB

`ImportError: cannot import name 'toolkitHome' from 'hera'` is misleading. The
name exists, but it's built lazily and creating it opens a MongoDB connection.
If Mongo is unreachable, the import fails as a side effect. The server is a
long-running process, so once it fails the broken state sticks until restart.

Fix: make sure MongoDB is up and reachable, then restart the server.

