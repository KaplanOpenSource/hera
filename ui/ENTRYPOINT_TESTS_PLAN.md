# Unit tests for every entrypoint

Written 2026-10-29. Status: **done**.

Every route the client calls gets a unit test on both sides.

## Client: one test file per api class

Each file checks every method: the path, the query or body, the label and
`silent`, and what the method answers when the call works or fails. The fetch
is stubbed with `tests/helpers/mockApiCall.ts`.

- `tests/serverApi.test.ts`: `init`, `cors`, `ensureJupyter` (`geoJson` is in `fetchGeoJson.test.ts`)
- `tests/projectApi.test.ts`: `names`, `documents`, `toolkits`, `create`, `detectNotebooks`, `delete`
- `tests/documentApi.test.ts`: `create`, `update`, `remove`
- `tests/repositoryApi.test.ts`: `get`, `registered`, `jsonsInFolder`, `loadIntoProject`
- `tests/nodeApi.test.ts`: `catalog`, `validateParams`
- `tests/workflowApi.test.ts`: `start`, `poll`

## Server: the routes with no test yet

- `tests/test_jupyter.py`: `/jupyter/ensure` with a fake notebook server.
  Starts one, reuses it for the same folder, restarts it for a new one, and
  binds to localhost unless CORS is on.
- `tests/test_nodes.py`: `/node-catalog` reads the real Hermes catalog.
- `tests/test_server.py`: `/cors` answers the origins the server was started with.
- `tests/test_files.py`: `/file/{path}` serves a file, or 404.

The rest already had route tests: `/init`, `/geojson`, `/healthz`,
`/node-params/validate`, `/start_workflow`, `/workflow_status`, and every
`/project/*`, `/document/*` and `/repository/*` route.
