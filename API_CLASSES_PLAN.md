# API package as classes

Written 2026-10-19. Status: **proposed**.

This replaces task 2 of `NEXT_TASKS_PLAN.md`. Same goal: one package that
holds every call to the server and nothing else. The difference is that each
route group is a class, not a set of loose functions.

---

## Why classes

- The reporter (toasts) comes in through the constructor. No module level
  `setReporter`, no hidden global.
- The rule "no React, no store, no snackbar in `src/api/`" is easy to see. A
  class has no hidden imports.
- Related calls sit together, so `/project/*` is one file instead of four.

The cost: every call site changes name, and the 22 test files that mock
`src/io/entrypoint` need a new way to stub. Both are mechanical. See below.

---

## The shape

```
src/api/ApiCall.ts        # base class: the one place that fetches
src/api/Reporter.ts       # the reporter type only
src/api/ProjectApi.ts     # /project/*
src/api/DocumentApi.ts    # /document/*
src/api/RepositoryApi.ts  # /repository/*
src/api/WorkflowApi.ts    # /workflow/*
src/api/NodeApi.ts        # /node-catalog, /node-params
src/api/ServerApi.ts      # /init, /cors, /jupyter/ensure, /geojson
src/api/api.ts            # the instances, built once
```

One class per file, named after the class.

### The base class

```ts
export class ApiCall {
  constructor(private reporter?: Reporter) { }

  protected async getJson(path: string, params: QueryParams, options: CallOptions): Promise<any | null> { }

  protected async postJson(path: string, body: object, options: CallOptions): Promise<any | null> { }
}
```

The body is today's `entrypoint.ts`, with `pushRunning` / `pushError` /
`dismiss` replaced by `this.reporter?.started(...)` and friends. With no
reporter the calls still work and say nothing, which is what the unit tests
want.

### A route group

```ts
export class RepositoryApi extends ApiCall {
  async get(path: string): Promise<RepositoryTrees | null> {
    return await this.getJson('/repository/get', { path }, { label: 'repo JSON' });
  }

  async registered(): Promise<Repository[] | null> { }

  async jsonsInFolder(folder: string): Promise<RepoJsonsInFolder | null> { }
}
```

Every method answers data or `null`. Nothing else. No store write, no toast of
its own.

### The instances

```ts
// src/api/api.ts
export const repositoryApi = new RepositoryApi();
export const projectApi = new ProjectApi();
// ...
```

`src/io/reporter.ts` builds a reporter from the snackbar. `main.tsx` installs
it. Simplest way: a setter on each instance, or one `installReporter(reporter)`
in `api.ts` that hands it to all of them. Pick the second, it is one call.

---

## What stays in `src/io/`

`src/io/` keeps its name and holds the client half:

- `FetchProjects.tsx` and `ProjectAutoReload.tsx`, the two React components.
- `snackbar.tsx`.
- `reporter.ts`, new.
- The wrappers that call an api method and then write to a store. Today that
  is `deleteDocuments`, `duplicateDocuments`, `loadRepositoryIntoProject`,
  `fetchInit`, `fetchProjectDetails`, `fetchProjectData`, `fetchProjectsNames`.

`parseToolkits.ts` is a parser, not IO. It moves to `src/utils/`.

---

## Tests

22 test files do this today:

```ts
vi.mock('../src/io/entrypoint', () => ({
  postJson: (...args: any[]) => mockPostJson(...args),
  getJson: vi.fn(),
}));
```

`vi.mock` on a module does not reach a class method. Replace it with a spy on
the base class prototype, in one shared helper:

```ts
// tests/helpers/mockApiCall.ts
export const mockApiCall = () => {
  const getJson = vi.fn();
  const postJson = vi.fn();
  vi.spyOn(ApiCall.prototype as any, 'getJson').mockImplementation(getJson);
  vi.spyOn(ApiCall.prototype as any, 'postJson').mockImplementation(postJson);
  return { getJson, postJson };
};
```

`protected` is a compile time word only, so the spy works at run time. Each
test file then swaps its `vi.mock` block for one `mockApiCall()` call in
`beforeEach`. The assertions on arguments do not change, because the paths and
bodies do not change.

`tests/entrypoint.test.ts` becomes `tests/ApiCall.test.ts` and tests the class
against a stubbed `fetch`, plus one case with a fake reporter and one with
none.

Fallback if the prototype spy turns out to fight the test setup: keep a tiny
`src/api/transport.ts` with the raw `fetch`, have `ApiCall` call it, and mock
that module instead. Only do this if the spy fails. It brings back a module
level seam we are trying to drop.

---

## Order: base, then six agents in parallel, then fuse

The six route groups barely touch the same files, so they can run at the same
time, one agent each. The base has to exist first, and the cleanup has to come
last.

### Phase 0: the base (serial, one person, no agents)

1. `Reporter.ts`, `ApiCall.ts`, `src/io/reporter.ts`, install in `main.tsx`.
2. `src/api/api.ts` with **all six** instances already listed, and all six
   `*Api.ts` files already created as empty classes extending `ApiCall`.
3. `tests/helpers/mockApiCall.ts` and `tests/ApiCall.test.ts`.

Step 2 matters. If `api.ts` is complete up front, no agent has to edit it, so
the one shared file never conflicts. Each agent only fills in the body of its
own class.

`src/io/entrypoint.ts` stays and still works. Nothing calls the classes yet.
Commit this. It is the base every agent starts from.

### Phase 1: six agents, in parallel

Each agent owns one class, its call sites, and its tests. Nothing else. Each
runs in its own git worktree, so the work cannot collide, and each ends with
the one relevant test file passing.

| Agent | Fills | Moves out of `src/io/` | Call sites it owns |
|---|---|---|---|
| repository | `RepositoryApi` | `fetchRepository.ts`, `fetchRegisteredRepositories.ts`, `fetchRepoJsonsInFolder.ts`, the call part of `loadRepositoryIntoProject.ts` | `DetailsViewRepo.tsx`, `useRegisteredRepositories.ts`, `CentralRepoFolder.tsx`, `DetailsViewMergedRepo.tsx`, `UpdateRepositoriesButton.tsx`, `RepoTreeWhole.tsx`, `RepoTreeAddButton.tsx` |
| document | `DocumentApi` | `createDocument.ts`, `FetchDocument.tsx`, the call part of `deleteDocuments.ts` and `duplicateDocuments.ts` | `AddDocumentButton.tsx`, `useShownDoc.ts`, `TreeContextMenu.tsx`, `DeleteSelectedButton.tsx`, `DeleteDocumentButton.tsx` |
| project | `ProjectApi` | the `getJson` calls inside `FetchProjects.tsx` | `FetchProjects.tsx`, `ProjectAutoReload.tsx` |
| workflow | `WorkflowApi` | `runWorkflow.ts` call part | `RunWorkflowButton.tsx`, `WorkflowRunPoller.tsx` |
| node | `NodeApi` | `validateNodeParams.ts` | `useNodeParamsValidation.ts`, `useNodeCatalog.ts` (raw fetch) |
| server | `ServerApi` | `fetchInit.ts` call part, `fetchGeoJson.ts` | `GeoJsonPreview.tsx`, `CorsIndicator.tsx`, `DetailsViewNotebook.tsx` |

Rules every agent gets:

- Touch only the files in your row. If you need a file in another row, stop
  and say so. Do not edit it.
- Do not touch `src/api/api.ts`, `src/api/ApiCall.ts`, `src/io/entrypoint.ts`,
  or `main.tsx`.
- The store writes stay in `src/io/`. Your class answers data or `null`.
- The paths, bodies and labels sent to the server must not change.
- Run only the test files for your group. Never the whole suite, never a
  build.
- Report: files changed, test files run and their result, anything you had to
  leave for the fuse step.

Two files need a word of care:

- `FetchProjects.tsx` also calls `fetchInit`. It belongs to the project agent
  only. The server agent leaves it alone, because the `fetchInit` wrapper
  stays in `src/io/`.
- `duplicateDocuments.ts` calls `createDocument`. Both are the document
  agent's, so there is no split.

### Phase 2: fuse (serial)

1. Merge the six worktrees in the table's order. Conflicts should be none, or
   only in `src/shared/types.ts` if two agents removed a type.
2. Run the import rule test and the greps under "How we know it holds".
3. Pull out any store write an agent left inside a class.
4. Delete `src/io/entrypoint.ts` and the empty `src/io/fetch*.ts` files.
5. Move `parseToolkits.ts` to `src/utils/`.
6. Run the full `TEST_UI.md` checklist once, here, not per agent.

Step 4 is the only one that can break an import we missed, so it comes last.

### If the parallel run is not wanted

The same table works serially, one row per commit, in the listed order. The
rows do not depend on each other.

---

## How we know it holds

- A test that reads every file in `src/api/` and fails on an import of
  `react`, `../stores/`, or `../io/`.
- `grep -rn "fetch(" src | grep -v src/api` answers nothing.
- The full unit suite passes.
- `src/api/*Api.ts` files have no `import` of anything but `ApiCall`, shared
  types, and the shared base url.

---

## Open question

`DetailsViewNotebook.tsx` also probes the jupyter url with
`fetch(url, { mode: 'no-cors' })`. That is not a call to our server, and it
wants the raw response, not JSON. Either give `ServerApi` a small `probe`
method that answers a boolean, or leave this one `fetch` in the component and
write the reason in a comment. Decide at step 7.
