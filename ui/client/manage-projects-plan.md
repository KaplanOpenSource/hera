# Plan: a "Manage projects" dialog

## Goal
- One dialog that lists all projects. Later it can export and delete several at once.
- The project dropdown stays as it is.
- Built in three steps. Each step ships on its own.

## Code
- `ManageProjectsButton.tsx`: only the menu item. It opens the dialog.
- `ManageProjectsDialog.tsx`: the dialog with the list.
- `ManageProjectsRow.tsx`: one line of the list.
- `FolderStateNote.tsx`: the note for a folder state, in a row and in the delete confirm.
- `isFilesFolderEmpty.ts`: moved out of `DeleteProjectButton.tsx`, so both use it.
- `exportProjectsToZip.ts`: exports names one by one into one saved `projects.zip`. Used by export and delete.
- `DeleteProjectsButton.tsx`: the delete button and the delete work.
- `DeleteProjectsConfirm.tsx`: the confirm dialog for deleting several projects.
- `ExportProgressDialog.tsx`: covers the manage dialog with the progress while Export runs. Delete shows it inside its confirm.
- `projectsText.ts`: "1 project" or "3 projects".

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
1. A "Delete 3 projects..." button next to export. Disabled when nothing is marked.
   - The "..." shows that it does not delete yet. It only opens a confirm dialog.
2. The confirm dialog:
   - Title: "Delete 3 projects?"
   - Lists the marked projects in groups:
     - Empty folder: "their files are deleted too".
     - Folder with files: shown on their own, with the folder path, and "their files are kept".
     - No folder or missing folder: nothing on disk to delete.
   - Says that `projects.zip` is saved first, and that nothing is deleted if an export fails.
   - A field asks the user to type the number of projects to delete, e.g. "3".
   - Buttons: "Cancel" and "Export and delete 3 projects".
     - The delete button stays off until the typed number matches.
     - Only this button starts the work.
3. On confirm, it first exports all marked projects, as in step 2.
   - If any export fails, stop at once. Nothing is saved and nothing is deleted.
   - An error message names the project that failed, e.g. "Delete stopped: Beta export failed."
4. When all exports succeed, it saves `projects.zip`, then deletes the projects one by one.
   - If a row was noted as an empty folder, delete its files too.
   - A folder with files is never deleted. There is no choice for it.
5. When done, one message says what happened, e.g. "Deleted 3 projects."
6. Update the project list. If the current project was deleted, switch to another, as today.


Tests:
- The button only opens the confirm; Cancel deletes and exports nothing.
- The confirm lists the marked projects in their groups (empty folder, folder with files, no folder).
- The delete button stays off until the typed number matches.
- All projects are exported and saved before any is deleted.
- A failed export stops it all: nothing saved, nothing deleted, the failed project is named.
- Empty-folder projects are sent with "delete files" on; all others with it off.
- The message says how many were deleted.
- The current project is replaced when it is deleted.

## Step 3.5: load several projects from one zip
1. "Load project" also takes a `projects.zip` made by step 2 (a zip holding `<name>.zip` files).
   - The client opens it with `fflate`. A zip of normal hera export files is a single project, as today.
2. For a zip of zips, the load dialog shows a table of the projects inside, all marked.
   - Each row: a checkbox, the name, and the files folder from that project's zip.
   - Each row gets the same checks as the single load: "A project named X already exists" and "This folder is also used by: ...".
   - Two rows in the zip with the same name or folder are flagged too.
   - Shift-click marks a range, as in the manage table.
3. A "Load 3 projects" button. Off when nothing is marked or a marked row has a name that already exists.
4. It sends each marked inner zip one by one with the existing load call. No server change.
   - If any load fails, stop at once. An error names the project that failed, and says which ones were already loaded.
5. When done, one message says so, e.g. "Loaded 3 projects.", and the project list updates.

Open questions for step 3.5:
- Can the user edit each row's name and folder, as in the single load, or only mark rows?
- When a load fails half way, should the ones already loaded stay?

Tests:
- A zip of zips shows one row per inner project, all marked.
- A single export zip still loads as one project, as today.
- Name and folder warnings show per row; a taken name blocks the load button.
- Only marked projects are sent, one by one.
- A failed load stops the rest and names the failed project and the loaded ones.

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
