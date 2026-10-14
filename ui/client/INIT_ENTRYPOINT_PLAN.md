# One entrypoint for startup: `/init`

Status: **built**. `/ready` is gone; `/healthz` still answers for liveness.
Toolkits come over as `toolkitDocs` and the client still parses them.

Issue #1028 part 3. The first named server entrypoint, replacing the `/exec`
calls the UI makes at startup. If it works out, the other operations can follow
the same shape later.

## Today

Six requests before the UI is usable:

| Request | Where |
|---|---|
| `/ready` poll | `ServerReadyGate` |
| datatypes | `ServerConstantReader` |
| project names | `FetchProjects`, on mount |
| toolkits + project documents | `FetchProjects`, after the names arrive |
| hera version | `VersionShower`, in the header |
| username + inDocker | `loadServerUser`, from the header |

Two round-trips deep, five of them arbitrary python assembled in the browser.

## The entrypoint

`GET /init?project=<name>` on the server. One python function, not `exec`.

While hera is still warming it answers `{"ready": false}`. Once warm it answers
everything the UI needs to start:

```json
{
  "ready": true,
  "server": { "heraVersion": "2.16.1", "username": "eran", "inDocker": true },
  "projects": [{ "name": "..." }],
  "datatypes": { "...": "..." },
  "toolkits": [ ... ],
  "project": { "name": "...", "documents": [ ... ] }
}
```

So `/ready` disappears: the gate polls `/init` and the answer that says "ready"
is also the data. One request, no `/exec` at startup.

### Choosing the project

The server picks it, so the client needs no second round-trip:

- the `project` parameter if that project exists
- else the first project in the list
- else none

A name that is not in the list means a bad URL, so the first project is a
better answer than an empty screen. The client fixes the URL afterwards.

### Partial failure

Each piece is gathered in its own try/except. A toolkit document missing its
`desc.classpath` (the #1011 crash) or one project whose documents fail to load
must not blank the whole UI; that piece comes back missing and the rest still
loads. This is the thing a single `exec` script cannot do, and the reason the
endpoint is worth building.

## Client changes

- `ServerReadyGate` polls `/init` instead of `/ready`, and on success fills the
  stores: project list, datatypes, toolkits, current project, server info.
- Then it mounts the app, with the data already in place.
- Remove the startup fetches: `ServerConstantReader`, the names fetch in
  `FetchProjects`, `VersionShower`'s fetch, `loadServerUser`'s fetch. Their
  stores are filled by init; the stores and components stay.
- `FetchProjects` keeps its URL-sync effect for later project switches only.
- The URL is corrected when the loaded project is not the one in the URL.

## Server changes

- A named operation module, e.g. `ui/server/operations/init.py`, with plain
  python functions. No code from the client.
- The route in `server.py` only parses the query and calls it.
- Keep `/exec` for now; the rest of the UI still uses it.

## Tests

- Python: the init function against the test project — project chosen three
  ways (named, missing, none), and one piece failing without killing the rest.
- Client: the gate shows the spinner while `ready` is false, fills the stores
  and mounts once ready, and corrects the URL.
- Then update `tests/execPython-coverage.md` and run the `TEST_UI.md` checklist.

## Risk

`ServerReadyGate` becomes the only thing standing between a cold server and the
app, so its retry and error paths matter more than before. A failed `/init` must
keep polling, not hang on a spinner forever with no message.
