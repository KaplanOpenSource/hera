# Named entrypoints to replace every `/exec` call

Status: **proposed**, except `/init` (see `INIT_ENTRYPOINT_PLAN.md`) and
`file.geoJson`, which are built.

Today the client builds python text and posts it to `/exec`: **22 call sites in
18 files**. That endpoint runs arbitrary code as the server user, and the client
pastes project names, document ids, file paths and whole repo JSON straight
into the code it sends.

The list below replaces all of them with **15 named entrypoints**. Each one is a
plain python function on the server, taking JSON parameters.

## Why 15 and not 22

Two reasons the count drops:

- **Seven write-then-reload pairs become one call**, by having the write return
  the fresh project: delete documents, delete one document, duplicate, load
  repository, update repositories, save document, create project (which also
  detects notebooks today).
- **Two operations are written out more than once**: `repository.get` lives in
  three places, `repository.jsonsInFolder` in two.

## Projects

| Entrypoint | Parameters | Returns | Replaces |
|---|---|---|---|
| `project.names` | - | the project list, default project first | `FetchProjects.fetchProjectsNames` |
| `project.documents` | `projectName` | the project and its documents | `FetchProjects.fetchProjectDetails`, and the second half of `fetchProjectData` |
| `project.toolkits` | - | the toolkit documents | first half of `fetchProjectData` |
| `project.create` | `projectName`, `filesDirectory`, `loadRepositories` | the project, notebooks already detected | `AddProjectButton` (both calls) |
| `project.delete` | `projectName`, `deleteFiles` | the remaining project list and the next project | `DeleteProjectButton` |
| `project.detectNotebooks` | `projectName`, `filesDirectory` | the project, and how many were added | `DetectNotebooksButton` |

`project.documents` is the busiest: reached from 8 places today.

## Documents

| Entrypoint | Parameters | Returns | Replaces |
|---|---|---|---|
| `document.get` | `docId` | the document | `FetchDocument.fetchDocument` |
| `document.add` | `projectName`, `kind`, `desc`, `resource`, `collection` | the project | `AddDocumentButton` |
| `document.update` | `docId`, `fields` | the document and the project | `FetchDocument.updateDocument` + the reload in `useShownDoc` |
| `document.delete` | `docIds` | the project | `DeleteDocumentButton`, `deleteDocuments` |
| `document.duplicate` | `docIds` | the project | `duplicateDocuments` |

`document.update` takes the changed fields as a dict, so the
`stringifyToPython` None/True/False hack goes away.

`document.add` takes a kind, so the four code branches of
`buildAddDocumentCode` (document, agent, notebook, workflow) become four
branches of python, and the notebook JSON template stops being a nested escaped
string in TypeScript.

## Repositories

| Entrypoint | Parameters | Returns | Replaces |
|---|---|---|---|
| `repository.registered` | - | the repository table | `useRegisteredRepositories` |
| `repository.get` | `nameOrPath` | the repository JSON, and the path-resolved copy | `DetailsViewRepo` (both branches), `RepoContents` |
| `repository.jsonsInFolder` | `folder`, `withContents` | the repo JSON files found, names only or with contents | `CentralRepoFolder`, `DetailsViewMergedRepo` |
| `repository.loadIntoProject` | `projectName`, `repositoryName` or `repositoryJSON`, `overwrite`, `baseDir` | the project | `loadRepositoryIntoProject`, `RepoTreeAddButton`, `UpdateRepositoriesButton` |

`repository.loadIntoProject` covers three callers that all load datasources into
a project, differing only in where the JSON comes from.

## Files

| Entrypoint | Parameters | Returns | Replaces | |
|---|---|---|---|---|
| `file.geoJson` | `path` | the file as GeoJSON in WGS84 | `GeoJsonPreview` | built: `GET /geojson?path=` |

## Already named

`/init` (project list, datatypes, toolkits, the chosen project, server info) and
the workflow endpoints (`/start_workflow`, `/workflow_status`, `/node-catalog`,
`/node-params/validate`).

## Notes

- `datatypes` and the server info need no entrypoint of their own: `/init`
  carries them, and nothing else asks.
- Validation moves to the server: every path and project name arrives as a
  parameter and is checked there, instead of being pasted into code.
- `/exec` can then be a dev-only switch, or go.
- 21 `/exec` call sites are left.
- Each entrypoint needs a python test next to it; the client keeps a small test
  that it sent the right name and parameters.
