# fetchPython call sites — test coverage

> **Maintainer note:** When adding or modifying calls, update this table.
>
> Startup no longer goes through `/exec`: one `GET /init` brings the project
> list, datatypes, toolkits, the chosen project and the server info. The python
> lives in `ui/server/init_data.py`.
>
> The project calls no longer go through `/exec` either: six named routes under
> `/project/*`, with the python in `ui/server/project_api/`. See
> `ENTRYPOINTS_PROJECTS_STEPS.md`.
>
> Every server call now lives in `src/api/`, one class per route group, all on
> the `ApiCall` base class. `src/io/entrypoint.ts` is gone. See
> `API_CLASSES_PLAN.md`.
>
> The document calls went the same way: three named routes under `/document/*`,
> with the python in `ui/server/document_api/`. There is no per-document read:
> `/project/documents` already brings them all, so the store is the only source.
> See `ENTRYPOINTS_DOCUMENTS_STEPS.md`.
>
> The repositories went last: four named routes under `/repository/*`, with the
> python in `ui/server/repository_api/`. A repository name is only a nickname for
> a file, so no name crosses the wire: every route takes a path. See
> `ENTRYPOINTS_REPOSITORIES_STEPS.md`.

## Named entrypoints (no exec)

Each api class also has its own test file, one test per method:
`serverApi`, `projectApi`, `documentApi`, `repositoryApi`, `nodeApi` and
`workflowApi` `.test.ts`. See `ui/ENTRYPOINT_TESTS_PLAN.md`.

| Entrypoint | Caller | Unit test? | Integration? | Test file |
|---|---|:---:|:---:|---|
| `GET /init` | `ServerApi.init`, `src/io/fetchInit.ts` | yes | no | `serverApi.test.ts`, `fetchInit.test.ts`, `serverReadyGate.test.tsx`, `fetchProjectsAfterInit.test.tsx`, `server/tests/test_init_data.py` |
| `GET /geojson` | `ServerApi.geoJson`, `GeoJsonPreview` | yes | no | `fetchGeoJson.test.ts`, `server/tests/test_geojson_file.py` |
| `GET /cors` | `ServerApi.cors`, `CorsIndicator` | yes | no | `serverApi.test.ts`, `server/tests/test_server.py` |
| `POST /jupyter/ensure` | `ServerApi.ensureJupyter`, `DetailsViewNotebook` | yes | no | `serverApi.test.ts`, `server/tests/test_jupyter.py` |
| `GET /node-catalog` | `NodeApi.catalog`, `useNodeCatalog` | yes | no | `nodeApi.test.ts`, `useNodeCatalog.test.ts`, `server/tests/test_nodes.py` |
| `POST /node-params/validate` | `NodeApi.validateParams`, `useNodeParamsValidation` | yes | no | `nodeApi.test.ts`, `nodeParamsValidation.test.tsx`, `server/tests/test_node_params.py` |
| `POST /start_workflow` | `WorkflowApi.start`, `RunWorkflowButton` | yes | no | `workflowApi.test.ts`, `runWorkflow.test.tsx`, `server/tests/test_server.py` |
| `GET /workflow_status/{token}` | `WorkflowApi.poll`, `WorkflowRunPoller` | yes | no | `workflowApi.test.ts`, `workflowRunPoller.test.tsx`, `server/tests/test_server.py` |
| `getJson` / `postJson` | `src/api/ApiCall.ts` | yes | no | `ApiCall.test.ts`, `apiPackageIsolation.test.ts` |
| `GET /project/names` | `fetchProjectsNames` | yes | yes | `FetchProjects.test.ts`, `server/tests/test_projects.py` |
| `GET /project/documents` | `fetchProjectDetails`, `fetchProjectData` | yes | yes | `FetchProjects.test.ts`, `server/tests/test_projects.py` |
| `GET /project/toolkits` | `fetchProjectData` | yes | yes | `FetchProjects.test.ts`, `server/tests/test_projects.py` |
| `POST /project/create` | `AddProjectButton` | yes | yes | `addProject.test.tsx`, `integration/addProject.integ.test.tsx`, `server/tests/test_projects.py` |
| `POST /project/delete` | `DeleteProjectButton` | yes | yes | `deleteProject.test.tsx`, `integration/deleteProject.integ.test.tsx`, `server/tests/test_projects.py` |
| `POST /project/detectNotebooks` | `DetectNotebooksButton` | yes | no | `detectNotebooks.test.tsx`, `server/tests/test_projects.py` |
| `POST /document/create` | `AddDocumentButton`, `duplicateDocuments` | yes | yes | `addDocument.test.tsx`, `documentOperations.test.ts`, `treeInteractions.test.tsx`, `integration/documents.integ.test.tsx`, `server/tests/test_documents.py` |
| `POST /document/update` | `useShownDoc` | yes | yes | `documentOperations.test.ts`, `saveDocument.test.tsx`, `TwoUiSync.test.tsx`, `integration/documents.integ.test.tsx`, `server/tests/test_documents.py` |
| `POST /document/delete` | `deleteDocuments`, `DeleteDocumentButton`, `DeleteSelectedButton` | yes | yes | `deleteDocumentButton.test.tsx`, `documentOperations.test.ts`, `projectActions.test.tsx`, `treeInteractions.test.tsx`, `integration/documents.integ.test.tsx`, `server/tests/test_documents.py` |
| `GET /repository/registered` | `useRegisteredRepositories` | yes | no | `registeredRepositories.test.tsx`, `server/tests/test_repositories.py` |
| `GET /repository/get` | `DetailsViewRepo` | yes | yes | `detailsViewRepo.test.tsx`, `integration/centralRepo.integ.test.tsx`, `server/tests/test_repositories.py` |
| `GET /repository/jsonsInFolder` | `CentralRepoFolder`, `DetailsViewMergedRepo` | yes | yes | `centralRepoFolder.test.tsx`, `detailsViewRepo.test.tsx`, `integration/centralRepo.integ.test.tsx`, `server/tests/test_repositories.py` |
| `POST /repository/loadIntoProject` | `RepoTreeAddButton`, `UpdateRepositoriesButton`, `loadRepositoryIntoProject` | yes | no | `repoTreeAddButton.test.tsx`, `updateRepositoriesButton.test.tsx`, `repoTreeInteractions.test.tsx`, `server/tests/test_repositories.py` |

A load answers `added` and `problems` as well as the project: hera skips what it
cannot load without raising, so the route counts the documents before and after
and reports the reason. See `ENTRYPOINTS_REPOSITORIES_STEPS.md`.

## Gaps

- No integration test for loading a repository into a project (`POST
  /repository/loadIntoProject`), or for the repository views in the browser.
