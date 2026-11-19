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

## Part 1.8: Let load read zips exported without the progress bar

`Project.export` with `show_progressbar=False` (what the UI uses) pickles each record alone, not in lists. `load` expects lists, so it crashes on these zips.

1. **hera: `_iter_pickled_docs` in `hera/datalayer/project.py`.** After unpickling an entry, if it is a single record (a dict), wrap it in a list of one. Lists stay as they are. This fixes `load` and the server's `name_in_zip`, and makes zips already exported from the UI loadable.
2. **Tests.** hera test that loads a zip exported with `show_progressbar=False`.

## Part 1.9: One settings record after load

Happens on every load into a new project, so it comes before Part 2.

1. **Tests first, to show the bug.** hera test that exports, deletes, loads with the same name, and checks for one settings record, in both import modes. It must fail (2 records) before the fix.
2. **Fix `Project.load` in `hera/datalayer/project.py`.** The project already has its own settings record when `load` starts. So when `load` reaches the zip's settings record, it copies those settings into the project's own record, and does not add the zip's record. Then the tests from step 1 pass.

## Part 2: Load with the original name

Just load an exported zip. The project gets the name it had when exported.

1. **Server: load route.** `POST /project/load` takes the zip, reads the project name from its records, and passes it to `Project.load`. Refuses if a project by that name already exists.
2. **Client: load button.** Pick a zip from the user's machine, load it, and switch to it. Show the server's error if the name is taken.
3. **Tests.** Server test for the load route. Client test for the load button.

## Part 2.5: Edge-case bugs found in hera

1. **Export writes one zip entry per record without the progress bar.** Not a blocker: `load` reads these zips since Part 1.8. But one entry per record is slow and makes big zip files for big projects. Optional fix in `Project.export`: use `_batched_cursor` in the `show_progressbar=False` path too, so all zips hold lists.
2. **Old test marked as expected to fail.** Done: the `xfail` mark on `test_export_and_load` is removed, and the test passes.

## Part 3: Load with a new name

The user picks the name.

1. **hera: match the old settings record by the old name.** Changes `Project.load` in `hera/datalayer/project.py`. The zip's record is `<oldName>__config__`. `load` must treat it as the settings record and copy it into `<newName>__config__`.
2. **Server: name in the load route.** Take an optional project name. Refuse if it is taken.
3. **Client: name field in the load dialog.** The Load button is disabled while the name is empty or already in the project list.
4. **Tests.** hera test for loading under a new name. Server and client tests for the name.
