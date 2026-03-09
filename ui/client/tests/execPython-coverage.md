# execPython call sites — test coverage

> **Maintainer note:** When adding or modifying `execPython`/`fetchPython` calls, update this table.

| # | Function | File | Unit test? | Latency? | Integration? | Test file |
|---|----------|------|:---:|:---:|:---:|---|
| 1 | `fetchProjectDetails` | `src/io/FetchProjects.tsx:17` | yes | yes | yes | `FetchProjects.test.ts` |
| 2 | `fetchToolkits` | `src/io/FetchProjects.tsx:36` | yes | yes | yes | `FetchProjects.test.ts` |
| 3 | `fetchDocument` | `src/io/FetchDocument.tsx:6` | yes | yes | no | `documentOperations.test.ts` |
| 4 | `updateDocument` | `src/io/FetchDocument.tsx:45` | yes | yes | yes | `documentOperations.test.ts` |
| 5 | `readAllConstants` | `src/stores/useServerConstants.ts:14` | yes | yes | no | `serverConstants.test.ts` |
| 6 | `doAddProject` | `src/components/header/AddProjectButton.tsx:31` | yes (integ) | no | yes | `integration/addProject.integ.test.tsx` |
| 7 | `deleteProject` | `src/components/header/DeleteProjectButton.tsx:13` | yes | yes | yes | `deleteProject.test.tsx` |
| 8 | `doAddDoc` | `src/components/project/AddDocumentButton.tsx:48` | yes | yes | yes | `addDocument.test.tsx` |
| 9 | `deleteDocument` | `src/components/project/ProjectDocumentItem.tsx:26` | yes | yes | yes | `projectDocumentItem.test.tsx` |
| 10 | load repo JSON (useEffect) | `src/components/details/DetailsViewRepo.tsx:25` | yes | yes | no | `detailsViewRepo.test.tsx` |
| 11 | `addRepo` | `src/components/details/RepoTreeAddButton.tsx:24` | yes | yes | no | `repoTreeAddButton.test.tsx` |

Already uses `fetchPython`:

| # | Function | File | Unit test? | Latency? | Integration? | Test file |
|---|----------|------|:---:|:---:|:---:|---|
| — | `fetchProjectsNames` | `src/io/FetchProjects.tsx:10` | yes | yes | yes | `FetchProjects.test.ts` |

## Gaps

- **#6 `doAddProject`**: Only has integration test, no unit test with latency. Integration test covers the flow end-to-end.
- **Integration missing for**: `fetchDocument` (#3), `readAllConstants` (#5), `loadRepoJson` (#10), `addRepo` (#11).
