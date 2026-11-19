# Project export / load plan (issue #527)

Goal: export a project to a zip, and load it back.

## Already done

- `Project.export(path)` in `hera/datalayer/project.py` writes all project records to a zip (#574).
- `Project.load(proj, path, is_hard_import)` reads that zip into a project (new or existing).
- Both open the zip with `zipfile.ZipFile`, which also takes a file object, so no Python change is needed.

## Scope

- The zip holds records only. Data files on disk are not copied (they can be many GB).
- After load, a record works only if its data file is still on disk at the same path.

## Part 1: Export

1. **Server: export route.** `POST /project/export` calls `Project.export` with an in-memory buffer (`io.BytesIO`), and sends the bytes back as a download. No file on the server.
2. **Client: export button.** In the project actions menu, next to delete. The browser saves the zip to the user's machine.
3. **Tests.** Server test for the export route. Client test for the export button.

## Part 2: Load

1. **Server: load route.** `POST /project/load` takes an uploaded zip and a project name, and passes the upload's file object straight to `Project.load`. Refuses if the name is already taken (a second check, after the client's).
2. **Client: load dialog.** Pick a zip from the user's machine and type the new project's name. The Load button is disabled while the name is empty or already in the project list. Then load it and switch to it.
3. **Tests.** Server test for the load route. Client test for the load dialog.

### Open questions

- The user picks any name, so the config record (`<oldName>__config__`) must be renamed to `<newName>__config__` on load. Suggested: fix it in `Project.load`, next to `updateProjectNameOnDoc`, so Python users get the fix too. This is a Python change.
