# Plan: a "Manage projects" dialog

## Goal
- One dialog that lists all projects. Later it can export and delete several at once.
- The project dropdown stays as it is.
- Built in three steps. Each step ships on its own.

## Code
- `ManageProjectsButton.tsx`: only the menu item. It opens the dialog.
- `ManageProjectsDialog.tsx`: the dialog with the list.
- `ManageProjectsRow.tsx`: one line of the list.
- `isFilesFolderEmpty.ts`: moved out of `DeleteProjectButton.tsx`, so both use it.

## Step 1: see the projects
1. A new item in the project menu: "Manage projects". It opens the dialog.
2. The dialog lists all projects except `defaultProject`.
   - A search field on top filters the list.
   - A line shows how many projects there are, e.g. "12 projects".
3. Each row shows the project name and its files folder.
4. Each row notes when the folder is empty (or holds only `caseConfiguration.json`): "empty folder".
   - When the dialog opens, it checks each project's folder with the existing files call.
   - The row shows "checking..." until it answers. A failed check shows no note.
5. No marking, no actions yet. No server change.

Tests:
- The list shows all projects but `defaultProject`.
- Search filters the list.
- Empty-folder rows show the note; rows with files do not; a failed check shows no note.

## Step 2: export several projects
1. Each row gets a checkbox.
   - A line always shows how many are marked, e.g. "3 of 12 marked" (also "0 of 12 marked").
   - Search never clears marks.
   - "Mark all" and "Clear" buttons.
2. An "Export 3 projects" button under the list. Disabled when nothing is marked.
3. It gets each marked project's zip, one by one, with the existing export call.
4. It puts them all in one zip, `projects.zip`, holding `<name>.zip` for each project, and saves it.
   - The client builds it with `fflate`, which is already used. No server change.
   - If any export fails, stop at once. Nothing is saved.
   - An error message names the project that failed, e.g. "Export stopped: Beta failed."
5. When all succeed, one message says so, e.g. "Exported 3 projects."

Tests:
- Marking, search keeps marks, mark all, clear, and the count stays right.
- One zip is saved, holding one zip per marked project.
- A failed export stops the rest, saves nothing, and names the failed project.

## Step 3: delete several projects
1. A "Delete 3 projects" button next to export. Disabled when nothing is marked.
2. First it exports all marked projects, as in step 2.
   - If any export fails, stop at once. Nothing is saved and nothing is deleted.
   - An error message names the project that failed, e.g. "Delete stopped: Beta export failed."
3. When all exports succeed, it saves `projects.zip`, then deletes the projects one by one.
   - If a row was noted as an empty folder, delete its files too.
4. When done, one message says what happened, e.g. "Deleted 3 projects."
5. Update the project list. If the current project was deleted, switch to another, as today.

Open question for step 3: is the confirm and the "delete files" choice the same as in the single delete?

Tests:
- All projects are exported and saved before any is deleted.
- A failed export stops it all: nothing saved, nothing deleted, the failed project is named.
- Empty-folder projects are sent with "delete files" on.
- The message says how many were deleted.
- The current project is replaced when it is deleted.

## Step 4: export as JSON instead of pickle (to solve later)
Today each project zip holds `chunk_0`, `chunk_1`, ... These are Python pickle files: binary, no extension, and only Python can read them.

Issue #574 (closed as a duplicate of #932):
- erasta asked for JSON chunks, with load still taking the old pickle zips.
- lior-antonov answered: do not fix the old export. Use Ilay's `exportDocumentsToRepository` (in `hera/utils/data/toolkit.py`), which writes JSON. Later remove the other export code.

What `exportDocumentsToRepository` does, and how it differs from today's zip:
- By default it writes all the project's records into one JSON file.
- It leaves out the project's settings record (it holds the files folder).
- Like the zip, it has no data files from disk, only the records that point to them.
- It saves the JSON on the server and registers it as a repository. It gives no file to download.
- If the JSON file already exists, it merges the records into it, not replaces it.
- It is brought back by loading a repository, not by "Load project".

Open questions:
- Can the UI export and load use it, or does it need changes in hera?
- How to keep the files folder, since the settings record is left out?
- Must load still take the old pickle zips?
- Should there be a new issue for this?

## Later (own plan)
- "Load project" accepts a zip of zips, and lets the user pick which project or projects to load.
