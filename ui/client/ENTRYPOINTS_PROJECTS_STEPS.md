# Steps: the project entrypoints

The **Projects** part of `ENTRYPOINTS_PLAN.md`. Six named entrypoints replace
the project `/exec` calls. Status: **built**.

## The routes

| Entrypoint | Route | Body / query | Answer |
|---|---|---|---|
| `project.names` | `GET /project/names` | - | `{projects}` |
| `project.documents` | `GET /project/documents?projectName=` | - | `{project}` |
| `project.toolkits` | `GET /project/toolkits` | - | `{toolkitDocs}` |
| `project.create` | `POST /project/create` | `projectName`, `filesDirectory`, `loadRepositories` | `{projects, project, addedNotebooks}` |
| `project.delete` | `POST /project/delete` | `projectName`, `deleteFiles` | `{projects, project}` |
| `project.detectNotebooks` | `POST /project/detectNotebooks` | `projectName`, `filesDirectory` | `{project, addedNotebooks}` |

While hera warms up every route answers 503 with `WARMING_UP`, as `/geojson`
does. A bad project name or a missing files directory is a 400 with a message.

## Server files

```
ui/server/project_api/project_api.py         # class ProjectApi: one project
ui/server/project_api/project_files.py       # class ProjectFiles: its files on disk
ui/server/project_api/toolkit_documents.py   # the toolkit documents
ui/server/routers/projects.py                # the six routes
```

`ProjectApi` is built from a project name and dropped when the request ends: it
keeps no state between calls, so each one reads the database again. A project
that does not exist yet still has a name, so `create` is a normal method and
answers the files directory it used. Only `names`, which is about no project at
all, is static.

`ProjectFiles` is built from the project's own files directory. It answers
whether a path is inside it, removes what is, and drops the directory once
nothing but its own `caseConfiguration.json` is left. Nothing outside the
directory is ever removed.

The toolkits are a plain function: they are not about a project.

`init_data.py` uses the same `ProjectApi` for startup, so `/init` and
`/project/*` cannot drift apart.

## Client files

- `src/io/entrypoint.ts` - `getJson` / `postJson`: one named call, with the
  running and error snackbars `fetchPython` showed.
- `src/io/FetchProjects.tsx` - `fetchProjectsNames`, `fetchProjectDetails` and
  `fetchProjectData` call the entrypoints. Their signatures do not change, so
  their eight callers stay as they are.
- `AddProjectButton` - one `project.create` call instead of create + detect, and
  one call in total: it puts the fresh list in the store, then hands the new
  project to the watcher as navigation state, so nothing is loaded again. The URL
  alone still says which project is current.
- `DeleteProjectButton` - `project.delete`; handles "no projects left".
- `DetectNotebooksButton` - `project.detectNotebooks`.
- `FetchProjects` - shows a project that arrived with the navigation instead of
  loading it, and compares the router's location, not `window.location`.
  `navigate()` changes the window at once but the route parameter only on the
  next render, so the old comparison saw a mismatch for one render and undid the
  navigation, leaving the previous project on screen.
- Gone: `src/io/ProjectCommands.ts`,
  `src/components/project/buildDetectNotebooksCode.ts`.

## Tests

- `ui/server/tests/test_projects.py` - `ProjectApi`, `ProjectFiles` and the six
  routes, on a fake hera in memory.
- Client tests assert the entrypoint name and parameters sent, not python text.

## Left over

`/exec` call sites after this: 15.
