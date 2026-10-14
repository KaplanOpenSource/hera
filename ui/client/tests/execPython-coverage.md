# fetchPython call sites — test coverage

> **Maintainer note:** All calls now use `fetchPython`. When adding or modifying calls, update this table.
>
> Startup no longer goes through `/exec`: one `GET /init` brings the project
> list, datatypes, toolkits, the chosen project and the server info. The python
> lives in `ui/server/init_data.py`.

| # | Function | File | Unit test? | Latency? | Integration? | Test file |
|---|----------|------|:---:|:---:|:---:|---|
| 0 | `fetchInit` (GET /init, no exec) | `src/io/fetchInit.ts:55` | yes | no | no | `fetchInit.test.ts`, `serverReadyGate.test.tsx`, `fetchProjectsAfterInit.test.tsx`, `server/tests/test_init_data.py` |
| 1 | `fetchProjectsNames` | `src/io/FetchProjects.tsx:12` | yes | yes | yes | `FetchProjects.test.ts` |
| 2 | `fetchProjectDetails` | `src/io/FetchProjects.tsx:15` | yes | yes | yes | `FetchProjects.test.ts` |
| 3 | `fetchToolkits` | `src/io/FetchProjects.tsx:33` | yes | yes | yes | `FetchProjects.test.ts` |
| 4 | `fetchDocument` | `src/io/FetchDocument.tsx:5` | yes | yes | no | `documentOperations.test.ts` |
| 5 | `updateDocument` | `src/io/FetchDocument.tsx:26` | yes | yes | yes | `documentOperations.test.ts` |
| 7 | `doAddProject` | `src/components/header/AddProjectButton.tsx:25` | yes (integ) | no | yes | `integration/addProject.integ.test.tsx` |
| 8 | `deleteProject` | `src/components/header/DeleteProjectButton.tsx:12` | yes | yes | yes | `deleteProject.test.tsx` |
| 9 | `doAddDoc` | `src/components/project/AddDocumentButton.tsx:34` | yes | yes | yes | `addDocument.test.tsx` |
| 10 | `deleteDocument` | `src/components/details/DeleteDocumentButton.tsx` | yes | yes | yes | `deleteDocumentButton.test.tsx` |
| 11 | load repo JSON (useEffect) | `src/components/details/DetailsViewRepo.tsx:23` | yes | yes | no | `detailsViewRepo.test.tsx` |
| 12 | `addRepo` | `src/components/details/RepoTreeAddButton.tsx:23` | yes | yes | no | `repoTreeAddButton.test.tsx` |
| 13 | `deleteDocuments` | `src/io/deleteDocuments.ts:6` | yes | no | no | `projectActions.test.tsx`, `treeInteractions.test.tsx` |
| 14 | `duplicateDocuments` | `src/io/duplicateDocuments.ts:32` | yes | no | no | `treeInteractions.test.tsx` |
| 15 | `loadRepositoryIntoProject` | `src/io/loadRepositoryIntoProject.ts:6` | yes | no | no | `repoTreeInteractions.test.tsx` |

## Gaps

- **#7 `doAddProject`**: Only has integration test, no unit test with latency.
- **#13-#15**: No latency tests and no integration tests.
- **Integration missing for**: `fetchDocument` (#4), `loadRepoJson` (#11), `addRepo` (#12), `/init` (#0).
