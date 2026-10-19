# Steps: the repository entrypoints

The **Repositories** part of `ENTRYPOINTS_PLAN.md` (issue #1129). Four named
entrypoints replace the last `/exec` call sites. Status: **proposed**.

Built the same way as the project and document entrypoints, so read
`ENTRYPOINTS_PROJECTS_STEPS.md` and `ENTRYPOINTS_DOCUMENTS_STEPS.md` first: same
router shape, same warmup guard, same `getJson` / `postJson` on the client.

## A repository is a file

A repository name is a nickname for a file. `addRepository` stores a datasource
document in `defaultProject`, toolkit `heradata`, with `datasourceName` the name
and `resource` the absolute path of the JSON file. `getRepository(name)` looks
the name up, takes `resource`, reads that file. `loadAllDatasourcesInRepository
ToProject(name)` does the same and resolves the paths against the file's folder.

So **no name crosses the wire**. Every route that takes a repository takes its
path. The table still answers the names, for showing.

That is what removes the guess in `DetailsViewRepo`, which today asks
`repoPath.includes('/') || repoPath.endsWith('.json')` and sends one of two
python texts.

## The routes

| Entrypoint | Route | Query / body | Answer |
|---|---|---|---|
| `repository.registered` | `GET /repository/registered` | - | `{repositories}` |
| `repository.get` | `GET /repository/get?path=` | - | `{json, resolved}` |
| `repository.jsonsInFolder` | `GET /repository/jsonsInFolder?folder=&withContents=` | - | `{files}`, or `{files, repoJsons}` |
| `repository.loadIntoProject` | `POST /repository/loadIntoProject` | `projectName`, `path`, `repositoryJSON`, `baseDir`, `overwrite` | `{project}` |

Every route reaches hera, so each answers 503 `WARMING_UP` while it warms, like
`/project/*` and `/document/*`.

## Which call site goes where

| # | Call site today | Entrypoint |
|---|---|---|
| 1 | `DetailsViewRepo` | `repository.get` |
| 2 | `DetailsViewMergedRepo` | `repository.jsonsInFolder`, `withContents=true` |
| 3 | `RepoTreeAddButton` | `repository.loadIntoProject`, with `repositoryJSON` |
| 4 | `CentralRepoFolder` | `repository.jsonsInFolder`, names only |
| 5 | `RepoContents` | deleted: nothing renders it |
| 6 | `useRegisteredRepositories` | `repository.registered` |
| 7 | `UpdateRepositoriesButton` | `repository.loadIntoProject`, no source |
| 8 | `loadRepositoryIntoProject` | `repository.loadIntoProject`, with `path` |

## The tree rows carry the path

`RegisteredRepositories` gives its rows an id made from `datasourceName` today.
The row holds the whole table row, so it has `resource` too, and its label
already shows it. The fetchers do not: `DetailsViewRepo` is handed only the tree
item id, and `RepoTreeWhole` reads the id off the DOM for the double click.
Neither holds the table, which is why the name has to be guessed apart from a
path.

So the row id becomes `resource`, and the label keeps showing the name. The id
is a prefix plus a string, and the central folder rows already put paths in it,
so nothing else changes. After that every repository row in the tree is a path.

This also fixes a real break: double clicking a central folder row sends its
file path into `loadRepositoryIntoProject`, which looks it up as a name and does
not find it.

Temp repository rows are not files and never fetch: `DetailsViewRepo` already
skips them on `isTempRepo`, and a temp tree is loaded with the JSON form.

## `repository.get`

Reads the file at `path` and answers two trees:

- `json`: the file as it is on disk.
- `resolved`: the same with the datasource paths made absolute, against the
  file's own directory.

The client keeps showing `resolved` and peeking at `json`, which is what the
`Link` / `LinkOff` button already does. A registered repository is read the same
way, so it gets the resolved paths too. It does not today, even though hera
resolves them when it loads.

A path that is not there, a file that is not JSON, or a JSON file that is not
shaped like a repository is a 404 with a message.

## `repository.jsonsInFolder`

`LOAD_REPO_JSON_PYTHON` is python text pasted into both folder calls today. It
moves to the server, and the two callers differ only in `withContents`:

- `files`: the repository files found under the folder, sorted. Always answered.
- `repoJsons`: path to contents, only when `withContents=true`.

`CentralRepoFolder` reads `files`, `DetailsViewMergedRepo` reads `repoJsons`.
The merge stays on the client in `RepoJsonMerger`: it is display logic and the
merged tree is edited in place.

`~` in the folder is expanded on the server. A folder that is not there answers
empty lists, not an error: the central folder is a user setting and may well
point at nothing yet. `caseConfiguration.json` is skipped, and so is any JSON
whose shape is not a repository, exactly as the python does today.

## `repository.loadIntoProject`

Three callers load datasources into a project. They differ only in where the
tree comes from:

- `path` - a repository file. The server reads it, takes its folder as the base
  directory, and calls `loadAllDatasourcesInRepositoryJSONToProject`. That is
  the three lines hera's own load-by-name runs once it has turned the name into
  a file.
- `repositoryJSON` - the tree the user just edited in the browser, with
  `baseDir` to resolve the relative paths against. The server makes it absolute,
  as the python does with `os.path.abspath`.
- neither - every registered repository
  (`loadAllDatasourcesInAllRepositoriesToProject`), which is the Update button.

Giving both `path` and `repositoryJSON` is a 400. `baseDir` only means something
with `repositoryJSON`.

`overwrite` comes from the caller: the Add dialog asks the user, the other two
send `False`, as today.

The Update button today goes through `hera.utils.data.CLI.update` with a
`SimpleNamespace`. That wrapper only reads `caseConfiguration.json` when no
project name is given, and the UI always has one, so the server calls
`loadAllDatasourcesInAllRepositoriesToProject` straight.

Every answer carries the project, the same `{name, documents}` shape
`/project/documents` gives, so the client puts it in the store. That is what
collapses the two write-then-reload pairs: `loadRepositoryIntoProject` and
`UpdateRepositoriesButton` each do a second round trip today.

## Server files

```
ui/server/repository_api/__init__.py
ui/server/repository_api/repo_json.py        # class RepositoryJson: a tree, loaded into a project
ui/server/repository_api/repo_json_file.py   # class RepoJsonFile: one repository file on disk
ui/server/repository_api/repo_folder.py      # class RepoFolder: the repository files under a folder
ui/server/repository_api/repository_table.py # class RepositoryTable: the registered rows
ui/server/routers/repositories.py            # the four routes
```

Each is built for one request and dropped, like `ProjectApi` and `DocumentApi`:
every call builds a fresh `dataToolkit()`.

- `RepositoryJson` is built from a tree and a base directory. It loads the tree
  into a project. It is the one place that reaches
  `loadAllDatasourcesInRepositoryJSONToProject`.
- `RepoJsonFile` is built from a path. It answers the tree when the file loads
  and is shaped like a repository, `None` otherwise; the resolved copy; and
  loading itself into a project, through `RepositoryJson` with its own folder as
  the base directory. One place for the shape rule, so `repository.get` and the
  folder scan cannot disagree.
- `RepoFolder` is built from a folder. It expands `~`, walks it for `*.json`,
  and answers the files `RepoJsonFile` accepts, with or without their contents.
- `RepositoryTable` answers the registered rows, and loads them all into a
  project.

They raise `NoSuchRepository` and `Refused`, which the router turns into 404 and
400.

The project in every load answer comes from `ProjectApi(projectName).documents()`,
so the repository routes and `/project/documents` cannot drift apart.

New in `api_models.py`: `LoadRepositoryPayload`.

Register the router in `server.py` next to the document one.

## Client files

- `src/io/fetchRegisteredRepositories.ts` - new: one `getJson`.
- `src/io/fetchRepository.ts` - new: one `getJson`, answering `json` and
  `resolved`.
- `src/io/fetchRepoJsonsInFolder.ts` - new: one `getJson`, with a
  `withContents` flag.
- `src/io/loadRepositoryIntoProject.ts` - one `postJson` for all three sources,
  storing the project it answers instead of calling `fetchProjectDetails`.
  `loadRepositoryIntoProject` now takes a path; two small functions next to it
  cover the JSON and the all-repositories forms.
- `src/components/repo/RegisteredRepositories.tsx` - the row id comes from
  `resource`; the label still shows the name.
- `src/components/details/DetailsViewRepo.tsx` - one `fetchRepository` call and
  no branch; the `isFilePath` test goes.
- `src/components/details/DetailsViewMergedRepo.tsx` - `fetchRepoJsonsInFolder`
  with contents.
- `src/components/details/RepoTreeAddButton.tsx` - sends the filtered tree as
  JSON and stores the project it gets back. `filterTree` stays: it is about what
  the user hid in the browser.
- `src/components/repo/CentralRepoFolder.tsx` - `fetchRepoJsonsInFolder`, names
  only.
- `src/components/repo/useRegisteredRepositories.ts` -
  `fetchRegisteredRepositories`.
- `src/components/repo/UpdateRepositoriesButton.tsx` - one call, and no reload
  after it.
- Gone: `src/shared/repoJsonPython.ts`,
  `src/components/repo/RepoContents.tsx`,
  `src/components/repo/LoadRepositoryButton.tsx`. Nothing imports the last two,
  not even a test.

## Tests

- `ui/server/tests/test_repositories.py` - the four classes and the four routes,
  on a fake `dataToolkit` and the fake hera database `test_documents.py` already
  builds. Cover: a path that is not there, a file that is not JSON, a JSON file
  that is not a repository, `caseConfiguration.json` skipped, a folder that is
  not there, `path` and `repositoryJSON` given at once, no source meaning all
  repositories, and that loading a file resolves against its own folder.
- Client tests assert the entrypoint name and the parameters sent, not python
  text. `detailsViewRepo.test.tsx`, `repoTreeAddButton.test.tsx`,
  `centralRepoFolder.test.tsx` and `repoTreeInteractions.test.tsx` already exist
  and need rewriting that way.
- Add the two tests that are missing today: `useRegisteredRepositories` and
  `UpdateRepositoriesButton`.
- Add a client test that a registered row opens by its path, not its name.
- Add a client test that loading a repository makes one call, not two.
- `tests/integration/centralRepo.integ.test.tsx` already covers the folder
  against a real server; point it at the new route.
- Update `ui/client/tests/execPython-coverage.md`: the sites go, and nothing is
  left on `/exec`.

## Order of work

1. `api_models.py` payload, `repository_api/` classes,
   `routers/repositories.py`.
2. `ui/server/tests/test_repositories.py`, green before any client change.
3. Delete the two dead components.
4. The three read calls: `fetchRegisteredRepositories`, `fetchRepository`,
   `fetchRepoJsonsInFolder`, and their four callers, including the row id
   change. Delete `repoJsonPython.ts`.
5. `loadRepositoryIntoProject.ts` and its three callers.
6. Client tests, then the full `TEST_UI.md` checklist.

## Follow-up, not part of this

- `CentralRepoFolder` and `DetailsViewMergedRepo` each keep their own copy of
  the `hera-central-repo-folder` key and its default. One shared place for the
  setting would do, but that is not this change.
- `/exec` has no callers left after this, so it can become a dev-only switch or
  go. Its own issue: the endpoint, `src/io/fetchPython.ts`, and the `WARMING_UP`
  constant that `entrypoint.ts` still imports from it.

## Left over

`/exec` call sites after this: 0.
