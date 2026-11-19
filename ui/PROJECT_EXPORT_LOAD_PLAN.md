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

## Part 2.1: One project menu in the header

1. **Header menu.** A three-dots button next to the project chooser opens a menu with Add project, Load project, Export project, and Delete project. The separate Add and Load icons go away.
2. **Tree menu.** The menu on the project name in the tree keeps only the document actions: Add document and Delete selected.
3. **Tests.** Client test for the header menu. Update tests that expect the old places.

## Part 2.8: Load dialog

1. **Client: load dialog.** "Load project" in the header menu opens the file picker first. Once a `.zip` is picked, a small MUI `Dialog` shows its file name, with Choose another zip, Cancel, and Load. Load sends the zip as today (Part 2) and closes the dialog. Closing the picker without a file shows no dialog.
2. **Tests.** Client test for the dialog: pick a file, choose another, Load, Cancel.

## Part 2.85: Show the zip's name and folder

1. **Client: read the zip in the browser.** When a file is picked, open the zip (a zip library, e.g. `fflate`) and read its records with `pickleparser` (npm, reads Python pickle in the browser). Take the project name and the files folder from the settings record.
2. **Client: show them read only** in the load dialog.
3. **Check first:** `pickleparser` must handle the Python-only objects in our records (Mongo ids, dates). Try it on a real exported zip before building on it. Done: it reads both zip kinds; Mongo ids come out as raw bytes, which we don't need.
4. **Tests.** Client test with a real exported zip as a test file.

## Part 2.9: Check the zip against the other projects in the load dialog

1. **Server: `GET /project/folders`.** Every project's name and files folder, read from each project's settings record. No side effects (does not open projects with `Project(...)`).
2. **Client: check in the dialog.** When the dialog opens, get the folders. Once the zip's name and folder are read: if a project with that name exists, show an error and disable Load; if other projects use the same folder, show a warning naming them, and keep Load enabled.
3. **Tests.** Server test for the route. Client tests for the block and the warning.

## Part 3: Load with a new name

To be planned in detail after Part 2.85 works. Draft:

The user picks the name. The loaded project keeps the old project's files folder (copied with its settings in Part 1.9), because its data files are there.

1. **hera: test only.** Part 1.9 already finds the zip's settings record by its old name. Add a hera test that loads under a new name and checks for one settings record with the old settings.
2. **Server: name in the load route.** Take an optional `?projectName=` in the URL (the zip is the raw body). Without it, use the name in the zip. Refuse if the name is taken.
3. **Server: shared folder warning.** After the load, check if any other project has the same files folder. If so, add a warning to the answer, naming those projects. Do not block the load.
4. **Client: editable name in the load dialog.** The Load button is disabled while the name is empty or already in the project list. After the load, show the server's warning if there is one.
5. **Client: check right away.** As soon as a file is picked, check the name and folder against the other projects in the dialog: block a taken name, warn on a shared folder. This needs the projects' folders on the client (a new server route that lists them).
6. **Tests.** hera test for the new name. Server tests for the name and the warning. Client tests for the dialog, the checks, and the warning.

## Part 3.5 (nice to have): Edge-case bugs found in hera

Not a blocker. Can be done any time.

1. **Export writes one zip entry per record without the progress bar.** Not a blocker: `load` reads these zips since Part 1.8. But one entry per record is slow and makes big zip files for big projects. Optional fix in `Project.export`: use `_batched_cursor` in the `show_progressbar=False` path too, so all zips hold lists.
2. **Old test marked as expected to fail.** Done: the `xfail` mark on `test_export_and_load` is removed, and the test passes.
