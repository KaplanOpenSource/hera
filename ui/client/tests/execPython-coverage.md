# fetchPython call sites — test coverage

> **Maintainer note:** When adding or modifying calls, update this table.
>
> Startup no longer goes through `/exec`: one `GET /init` brings the project
> list, datatypes, toolkits, the chosen project and the server info. The python
> lives in `ui/server/init_data.py`.
>
> The project calls no longer go through `/exec` either: six named routes under
> `/project/*`, with the python in `ui/server/project_api/`. The client calls them
> through `src/io/entrypoint.ts`. See `ENTRYPOINTS_PROJECTS_STEPS.md`.
>
> The document calls went the same way: three named routes under `/document/*`,
> with the python in `ui/server/document_api/`. There is no per-document read:
> `/project/documents` already brings them all, so the store is the only source.
> See `ENTRYPOINTS_DOCUMENTS_STEPS.md`.

## Named entrypoints (no exec)

| Entrypoint | Caller | Unit test? | Integration? | Test file |
|---|---|:---:|:---:|---|
| `GET /init` | `src/io/fetchInit.ts` | yes | no | `fetchInit.test.ts`, `serverReadyGate.test.tsx`, `fetchProjectsAfterInit.test.tsx`, `server/tests/test_init_data.py` |
| `GET /geojson` | `src/io/fetchGeoJson.ts` | yes | no | `fetchGeoJson.test.ts`, `server/tests/test_geojson_file.py` |
| `getJson` / `postJson` | `src/io/entrypoint.ts` | yes | no | `entrypoint.test.ts` |
| `GET /project/names` | `fetchProjectsNames` | yes | yes | `FetchProjects.test.ts`, `server/tests/test_projects.py` |
| `GET /project/documents` | `fetchProjectDetails`, `fetchProjectData` | yes | yes | `FetchProjects.test.ts`, `server/tests/test_projects.py` |
| `GET /project/toolkits` | `fetchProjectData` | yes | yes | `FetchProjects.test.ts`, `server/tests/test_projects.py` |
| `POST /project/create` | `AddProjectButton` | yes | yes | `addProject.test.tsx`, `integration/addProject.integ.test.tsx`, `server/tests/test_projects.py` |
| `POST /project/delete` | `DeleteProjectButton` | yes | yes | `deleteProject.test.tsx`, `integration/deleteProject.integ.test.tsx`, `server/tests/test_projects.py` |
| `POST /project/detectNotebooks` | `DetectNotebooksButton` | yes | no | `detectNotebooks.test.tsx`, `server/tests/test_projects.py` |
| `POST /document/create` | `AddDocumentButton`, `duplicateDocuments` | yes | yes | `addDocument.test.tsx`, `documentOperations.test.ts`, `treeInteractions.test.tsx`, `integration/documents.integ.test.tsx`, `server/tests/test_documents.py` |
| `POST /document/update` | `useShownDoc` | yes | yes | `documentOperations.test.ts`, `saveDocument.test.tsx`, `TwoUiSync.test.tsx`, `integration/documents.integ.test.tsx`, `server/tests/test_documents.py` |
| `POST /document/delete` | `deleteDocuments`, `DeleteDocumentButton`, `DeleteSelectedButton` | yes | yes | `deleteDocumentButton.test.tsx`, `documentOperations.test.ts`, `projectActions.test.tsx`, `treeInteractions.test.tsx`, `integration/documents.integ.test.tsx`, `server/tests/test_documents.py` |

## Still on `/exec`

8 call sites in 7 files, all of them repositories.

| # | Function | File | Unit test? | Latency? | Integration? | Test file |
|---|----------|------|:---:|:---:|:---:|---|
| 1 | load repo JSON (useEffect) | `src/components/details/DetailsViewRepo.tsx:22` | yes | yes | no | `detailsViewRepo.test.tsx` |
| 2 | merged repo JSONs | `src/components/details/DetailsViewMergedRepo.tsx:22` | yes | no | no | `detailsViewRepo.test.tsx` |
| 3 | `addRepo` | `src/components/details/RepoTreeAddButton.tsx:41` | yes | yes | no | `repoTreeAddButton.test.tsx` |
| 4 | repo JSONs in folder | `src/components/repo/CentralRepoFolder.tsx:23` | yes | no | yes | `centralRepoFolder.test.tsx` |
| 5 | repo contents | `src/components/repo/RepoContents.tsx:21` | yes | no | no | `repoTreeInteractions.test.tsx` |
| 6 | `useRegisteredRepositories` | `src/components/repo/useRegisteredRepositories.ts:9` | no | no | no | - |
| 7 | update repositories | `src/components/repo/UpdateRepositoriesButton.tsx:14` | no | no | no | - |
| 8 | `loadRepositoryIntoProject` | `src/io/loadRepositoryIntoProject.ts:8` | yes | no | no | `repoTreeInteractions.test.tsx` |

## Gaps

- **#6, #7**: no test at all.
- **#8**: no latency test and no integration test.
- **Integration missing for**: the repo views (#1-#3, #5).
