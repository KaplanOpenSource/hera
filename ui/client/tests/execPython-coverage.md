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

## Still on `/exec`

15 call sites in 13 files.

| # | Function | File | Unit test? | Latency? | Integration? | Test file |
|---|----------|------|:---:|:---:|:---:|---|
| 1 | `fetchDocument` | `src/io/FetchDocument.tsx:8` | yes | yes | no | `documentOperations.test.ts` |
| 2 | `updateDocument` | `src/io/FetchDocument.tsx:48` | yes | yes | yes | `documentOperations.test.ts` |
| 3 | `doAddDoc` | `src/components/project/AddDocumentButton.tsx:85` | yes | yes | yes | `addDocument.test.tsx` |
| 4 | `deleteDocument` | `src/components/details/DeleteDocumentButton.tsx:25` | yes | yes | yes | `deleteDocumentButton.test.tsx` |
| 5 | `deleteDocuments` | `src/io/deleteDocuments.ts:8` | yes | no | no | `projectActions.test.tsx`, `treeInteractions.test.tsx` |
| 6 | `duplicateDocuments` | `src/io/duplicateDocuments.ts:34` | yes | no | no | `treeInteractions.test.tsx` |
| 7 | load repo JSON (useEffect) | `src/components/details/DetailsViewRepo.tsx:22` | yes | yes | no | `detailsViewRepo.test.tsx` |
| 8 | merged repo JSONs | `src/components/details/DetailsViewMergedRepo.tsx:22` | yes | no | no | `detailsViewRepo.test.tsx` |
| 9 | `addRepo` | `src/components/details/RepoTreeAddButton.tsx:41` | yes | yes | no | `repoTreeAddButton.test.tsx` |
| 10 | repo JSONs in folder | `src/components/repo/CentralRepoFolder.tsx:23` | yes | no | yes | `centralRepoFolder.test.tsx` |
| 11 | repo contents | `src/components/repo/RepoContents.tsx:21` | yes | no | no | `repoTreeInteractions.test.tsx` |
| 12 | `useRegisteredRepositories` | `src/components/repo/useRegisteredRepositories.ts:9` | no | no | no | - |
| 13 | update repositories | `src/components/repo/UpdateRepositoriesButton.tsx:14` | no | no | no | - |
| 14 | `loadRepositoryIntoProject` | `src/io/loadRepositoryIntoProject.ts:8` | yes | no | no | `repoTreeInteractions.test.tsx` |

## Gaps

- **#12, #13**: no test at all.
- **#5, #6, #14**: no latency tests and no integration tests.
- **Integration missing for**: `fetchDocument` (#1), the repo views (#7-#9, #11).
