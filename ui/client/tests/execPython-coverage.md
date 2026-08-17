# fetchPython call sites — test coverage

> **Maintainer note:** All calls now use `fetchPython`. When adding or modifying calls, update this table.

| # | Function | File | Unit test? | Latency? | Integration? | Test file |
|---|----------|------|:---:|:---:|:---:|---|
| 1 | `fetchProjectsNames` | `src/io/FetchProjects.tsx:9` | yes | yes | yes | `FetchProjects.test.ts` |
| 2 | `fetchProjectDetails` | `src/io/FetchProjects.tsx:15` | yes | yes | yes | `FetchProjects.test.ts` |
| 3 | `fetchToolkits` | `src/io/FetchProjects.tsx:33` | yes | yes | yes | `FetchProjects.test.ts` |
| 4 | `fetchDocument` | `src/io/FetchDocument.tsx:5` | yes | yes | no | `documentOperations.test.ts` |
| 5 | `updateDocument` | `src/io/FetchDocument.tsx:26` | yes | yes | yes | `documentOperations.test.ts` |
| 6 | `readAllConstants` | `src/stores/useServerConstants.ts:13` | yes | yes | no | `serverConstants.test.ts` |
| 7 | `doAddProject` | `src/components/header/AddProjectButton.tsx:25` | yes (integ) | no | yes | `integration/addProject.integ.test.tsx` |
| 8 | `deleteProject` | `src/components/header/DeleteProjectButton.tsx:12` | yes | yes | yes | `deleteProject.test.tsx` |
| 9 | `doAddDoc` | `src/components/project/AddDocumentButton.tsx:34` | yes | yes | yes | `addDocument.test.tsx` |
| 10 | `deleteDocument` | `src/components/details/DeleteDocumentButton.tsx` | yes | yes | yes | `deleteDocumentButton.test.tsx` |
| 11 | load repo JSON (useEffect) | `src/components/details/DetailsViewRepo.tsx:23` | yes | yes | no | `detailsViewRepo.test.tsx` |
| 12 | `addRepo` | `src/components/details/RepoTreeAddButton.tsx:23` | yes | yes | no | `repoTreeAddButton.test.tsx` |

## Gaps

- **#7 `doAddProject`**: Only has integration test, no unit test with latency.
- **Integration missing for**: `fetchDocument` (#4), `readAllConstants` (#6), `loadRepoJson` (#11), `addRepo` (#12).
