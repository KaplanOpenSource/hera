# Known Integration Issues

Setup problems that break the UI/server and how to work around them.

## Hermes

Running workflows needs the **hermes** package. Missing it fails with the
misleading error: `ValueError: The workflow type None not found`.

This only happens with an outdated Docker image. Running without Docker (a local
env that already has hermes) is not affected.

hermes is a submodule at `Hermes/` with no `setup.py`, so it can't be
pip-installed. The two Docker configs in `.vscode/launch.json` handle it:

- Add `/app/Hermes` to `PYTHONPATH` so `import hermes` works.
- `pip install jsonpath_rw_ext` (its runtime dependency) before server start.

Verify: `docker exec hera-server-instance python -c "import hermes"`

TODO: permanent fix = add `/app/Hermes` to `dockerfile` PYTHONPATH +
`jsonpath_rw_ext` to `requirements.txt`, then rebuild.

## MongoDB

`ImportError: cannot import name 'toolkitHome' from 'hera'` is misleading. The
name exists, but it's built lazily and creating it opens a MongoDB connection.
If Mongo is unreachable, the import fails as a side effect. The server is a
long-running process, so once it fails the broken state sticks until restart.

Fix: make sure MongoDB is up and reachable, then restart the server.

