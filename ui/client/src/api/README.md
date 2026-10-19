# src/api

Every call to the hera server, and nothing else.

## The rule

A file in this folder must not import React, a store, the snackbar, or
anything from `src/io/`. `tests/apiPackageIsolation.test.ts` checks it and
fails the build if you break it.

Why: a server call is not UI work. Keeping them apart means a test can stub one
place, and a call can be reused from anywhere without dragging a store in.

## The shape

`ApiCall.ts` is the base class. It is the only place in the client that calls
`fetch`. It gives subclasses two protected methods:

```ts
this.getJson(path, params, { label, silent })
this.postJson(path, body, { label, silent })
```

`label` names the call in the running toast, and goes in front of the message
if it fails. `silent` means no toast at all, for background calls like the
auto-reload or a probe.

One class per route group, one file per class:

| Class | Routes |
|---|---|
| `ProjectApi` | `/project/*` |
| `DocumentApi` | `/document/*` |
| `RepositoryApi` | `/repository/*` |
| `WorkflowApi` | `/start_workflow`, `/workflow_status/{token}` |
| `NodeApi` | `/node-catalog`, `/node-params/validate` |
| `ServerApi` | `/init`, `/geojson`, `/cors`, `/jupyter/ensure` |

`api.ts` builds one instance of each. Always import from there:

```ts
import { projectApi } from '../../api/api';

const data = await projectApi.names();
if (!data) {
  return;
}
```

## Every method answers data or null

A method returns its answer, or `null` when the call failed. It never throws,
so a caller needs no try/catch. A failed call has already been logged and
reported, so the caller only decides what to show instead.

One exception: `RepositoryApi.registered()` answers `[]` rather than `null`,
because every caller wanted an empty list anyway.

## Toasts

The api package does not know about notistack. It reports through a `Reporter`
(see `Reporter.ts`): `started`, `finished`, `failed`.

`src/io/reporter.ts` builds one from the snackbar, and `main.tsx` installs it
with `installReporter`. With no reporter installed the calls still work and say
nothing, which is what the unit tests want.

## Store writes do not belong here

A method answers the data. Putting it in a store is client work, so it happens
in the caller, or in a thin wrapper in `src/io/` that keeps its old name and
signature. Those wrappers are `fetchInit`, `deleteDocuments`,
`duplicateDocuments`, `loadRepositoryIntoProject`, and the three fetchers
inside `src/io/FetchProjects.tsx`.

## Adding a call

1. Find the class for its route group, or add a new class and a new instance in
   `api.ts` if the group is new.
2. Add a short method. No `fetch` prefix in the name: `names`, not
   `fetchNames`. The class name already says it is a server call.
3. Give it a `label`. Add `silent: true` only if a toast would be noise.
4. Put the answer's type in the same file, next to the class.
5. Call it from the component or hook. Do not write a wrapper in `src/io/`
   unless the call has to write to a store.
6. Add a row to `tests/execPython-coverage.md`.

## Testing

Use `tests/helpers/mockApiCall.ts`. It spies on `ApiCall.prototype`, so one
call stubs every api class at once:

```ts
let mockPostJson: ReturnType<typeof vi.fn>;

beforeEach(() => {
  vi.clearAllMocks();
  mockPostJson = mockApiCall().postJson;
});
```

Then assert on the path, the body and the label. Call it inside `beforeEach`,
after `vi.clearAllMocks()`, because `vi.spyOn` is per test.

To stub one method instead of the transport, spy on the instance:
`vi.spyOn(documentApi, 'update')`.

`ApiCall.test.ts` covers the base class itself against a stubbed `fetch`.

## What is not here

- `src/io/fetchPython.ts` still calls `/exec` directly. No UI code uses it.
  Retiring it is task 3 of `../../../../NEXT_TASKS_PLAN.md`.
- `DetailsViewNotebook.tsx` probes the jupyter server with
  `fetch(url, { mode: 'no-cors' })`. That is not our server and it wants the
  raw answer, so it stays in the component.

The move to this folder is written up in `../../../../API_CLASSES_PLAN.md`.
