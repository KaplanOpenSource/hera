# Three tasks

Written 2026-10-19. Status: **done** (task 3 on 2026-10-29).

1. Build the minimal requirements from `requirements.txt`, so no pin is written
   twice.
2. One API package, apart from the client code.
3. Retire `/exec`.

They are independent. The file is in work order.

---

## 1. Minimal requirements without a second pin

`requirements-minimal.txt` repeats pins that `requirements.txt` already has.
They drift. Today `matplotlib`, `rasterio` and `Unum` had to be copied across by
hand.

### The shape

`requirements-minimal.txt` becomes names only, no versions:

```
fastapi
uvicorn[standard]
pydantic
mongoengine
...
```

`scripts/minimal_requirements.py` reads both files and prints the pinned lines:

```sh
python scripts/minimal_requirements.py requirements.txt requirements-minimal.txt > /tmp/req.txt
```

Rules:

- A name in the minimal list takes its line from `requirements.txt`, verbatim.
- Matching ignores case and extras, so `uvicorn[standard]` finds
  `uvicorn[standard]>=0.30.0`, and `Deprecated` finds `Deprecated==1.2.18`.
- A name that is not in `requirements.txt` is an error, and the script exits
  non zero saying which. No silent "install whatever pip likes".
- Comments and blank lines in the minimal list are kept out of the answer.
- Order follows the minimal list, so the file stays readable.

`argcomplete` is now in `requirements.txt`, so every name in every chosen list
finds a pin. A test asserts that, per list, so a new name without a pin fails
the build instead of installing whatever pip likes.

### Where it is used

`dockerfile` today copies both files and runs pip on one of them. It becomes:

```dockerfile
COPY requirements.txt requirements-minimal.txt scripts/minimal_requirements.py ./scripts/
ARG REQUIREMENTS=requirements.txt
RUN python scripts/minimal_requirements.py ... > /tmp/req.txt && \
    python -m pip install --no-cache-dir -r /tmp/req.txt
```

With `REQUIREMENTS=requirements.txt` the script is skipped and pip reads the
full file, so the slow build does not change.

`ui/server/requirements-server.txt` and `ui/server/requirements-integ.txt` are
gone. They were two more lists of the same names, so they merged into
`requirements-minimal.txt`, which already had all of them but `pytest`. The two
CI jobs that used them now build the minimal set with the script, like the
dockerfile does. They install a few packages they do not need, which is the
price of one list, `pytest` among them. The client tests install no python, so
there is nothing to do there.

### Tests

`ui/server/tests/` is the only python test home, but this script is not the
server. Put the test next to the script: `scripts/tests/test_minimal_requirements.py`,
run by the same pytest. Cover: a name with extras, a name with a different
case, a `>=` pin, a comment line, and a name missing from `requirements.txt`
failing loudly.

---

## 2. One API package

Today `src/io/` holds four different kinds of thing:

| Kind | Files |
|---|---|
| Server calls | `entrypoint.ts`, `fetch*.ts`, `createDocument.ts`, `deleteDocuments.ts`, `duplicateDocuments.ts`, `loadRepositoryIntoProject.ts`, `runWorkflow.ts`, `validateNodeParams.ts`, `FetchDocument.tsx` |
| React components | `FetchProjects.tsx`, `ProjectAutoReload.tsx` |
| UI feedback | `snackbar.tsx` |
| A parser | `parseToolkits.ts` |

And four calls still go straight to `fetch` from a component:
`CorsIndicator.tsx` (`/cors`), `DetailsViewNotebook.tsx` (`/jupyter/ensure`),
`useNodeCatalog.ts` (`/node-catalog`).

### What the package is

`src/api/` is every call to the server and nothing else. The rule: a file in
`src/api/` may not import React, a store, or the snackbar. One file per
entrypoint group, named after the route group:

```
src/api/call.ts          # getJson / postJson, the one place that fetches
src/api/project.ts       # /project/*
src/api/document.ts      # /document/*
src/api/repository.ts    # /repository/*
src/api/workflow.ts      # /workflow/*
src/api/node.ts          # /node-catalog, /node-params
src/api/server.ts        # /init, /cors, /jupyter/ensure, /geojson
```

Each function answers data or `null`, and nothing else. That is already how
`fetchRepository` works.

### The two knots

**Toasts.** `entrypoint.ts` pushes the running and error toasts itself, so it
imports the snackbar. The API package cannot. Give `call.ts` a reporter it is
handed once at startup:

```ts
// src/api/call.ts
export type Reporter = {
  started: (label: string) => number,
  finished: (ticket: number) => void,
  failed: (label: string, problem: string) => void,
};
export const setReporter = (reporter: Reporter) => { ... };
```

`src/io/reporter.ts` builds one from the snackbar and installs it in `main.tsx`.
With no reporter installed the calls still work and say nothing, which is what
the tests want.

**Store writes.** `fetchProjectDetails`, `loadRepositoryIntoProject` and
friends write to `useProjectStore` after the call. That is client work, not API
work. The API function answers the project; a thin file beside the store puts
it in. So `src/io/` keeps its name and holds the client half: the two React
components, the snackbar, the reporter, and the store-writing wrappers.

### Order

1. `src/api/call.ts` with the reporter, and `src/io/reporter.ts`. Keep
   `entrypoint.ts` as a re-export so nothing breaks yet.
2. Move the pure fetchers over, one route group per commit: repository,
   document, project, workflow, node, server.
3. Pull the store writes out of the moved functions into `src/io/`.
4. Move the four raw `fetch` calls into `src/api/`.
5. Delete `entrypoint.ts` and the empty `io/fetch*` files.
6. `parseToolkits.ts` is a parser, not IO: it goes to `src/utils/`.

### How we know it holds

- A lint rule or a test that reads `src/api/*` and fails on any import of
  `react`, `../stores/`, or `../io/snackbar`.
- `grep -rn "fetch(" src | grep -v src/api` answers nothing.
- The existing tests keep passing. They mock `../src/io/entrypoint` today, so
  each move updates the mock path to `../src/api/call` in that group's tests.

If a path alias is wanted (`@api/*`), note that `@shared/*` is declared in
`tsconfig.app.json` only. Add it to `vite.config.ts` as well, or keep relative
imports.

---

## 3. Retire `/exec`

No UI code calls `/exec`. What keeps it alive:

| What | Where |
|---|---|
| The `WARMING_UP` constant | `src/io/fetchPython.ts`, imported by the one file that fetches |
| Dead mocks | `tests/addProject.test.tsx`, `tests/addProjectNavigation.test.tsx` |
| Test setup over the wire | `tests/integration/addProject.api.test.ts`, `tests/integration/centralRepo.integ.test.tsx` |
| The route | `ui/server/routers/exec_code.py`, registered in `server.py` |
| Route tests | `ui/server/tests/test_server.py`, `ui/server/tests/test_create_app.py` |

### Step 3a: the client stops shipping it

- Move `WARMING_UP` to `src/shared/constants.ts`. It is a server answer, not
  python. After task 2 the importer is `src/api/call.ts`.
- Delete `src/io/fetchPython.ts`, and `ExecRequest` / `ExecResponse` from
  `src/shared/types.ts` once nothing uses them.
- Delete the `fetchPython` mocks in the two tests above. Both already mock
  `entrypoint`, so the mocks do nothing.
- Check: `grep -rn "fetchPython\|/exec" src` answers nothing.

This alone is worth doing. After it, no browser code can reach `/exec`.

### Step 3b: the integration tests stop using it

Two tests use `/exec` as a back door to set up and tear down state:

- `addProject.api.test.ts`: lists projects and deletes one.
- `centralRepo.integ.test.tsx`: writes repository files on the server, then
  removes the folder.

The server runs in the same container as the tests, through
`tests/integration/globalSetup.ts`. So the setup can be done there instead of
over HTTP: write the files with `node:fs`, and delete the project with
`POST /project/delete`, which already exists.

The folder the central repo test writes to is `/tmp/hera-test-central-repos`,
which both sides can see. Nothing needs a new route.

### Step 3c: the route goes

- Delete `ui/server/routers/exec_code.py` and its line in `server.py`.
- Delete `ExecPayload`, `ExecResponse` and `Problem` from `api_models.py`.
- Delete the `/exec` tests in `test_server.py` and the `/exec` entry in
  `test_create_app.py`.
- `ui/client/tests/execPython-coverage.md` becomes a note that the file is
  history, or goes.

Why it should go and not stay as a dev switch: it runs any python the caller
sends, and with `--cors remote` the browser origins of the whole machine are
allowed. A debugger belongs behind a flag at most, and we no longer need one.

If a dev back door is wanted anyway, keep it behind an explicit
`--allow-exec` argument that defaults to off, and say so in `ui/README.md`.

---

## Order of the day

1. Task 1: one script, one test file, one dockerfile line.
2. Task 2: the big one, one route group at a time.
3. Task 3a: minutes, no behaviour change.
4. Task 3b and 3c, once nothing else needs `/exec`.
