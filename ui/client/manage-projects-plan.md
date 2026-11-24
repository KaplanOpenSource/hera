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
   - A project whose export failed is left out.
   - If all failed, nothing is saved.
5. When done, one message says what happened, e.g. "Exported 2 of 3. Failed: Beta."

Tests:
- Marking, search keeps marks, mark all, clear, and the count stays right.
- One zip is saved, holding one zip per marked project.
- A failed export is left out and listed in the message.
- When all fail, nothing is saved.

## Step 3: delete several projects
1. A "Delete 3 projects" button next to export. Disabled when nothing is marked.
2. It goes through the marked projects one by one:
   - Export the zip first. If it fails, skip this project and go on.
   - Delete the project. If its row was noted as an empty folder, delete its files too.
3. When done, one message says what happened, e.g. "Deleted 2 of 3. Not deleted: Beta (export failed)."
4. Update the project list. If the current project was deleted, switch to another, as today.

Open question for step 3: is the confirm and the "delete files" choice the same as in the single delete?

Tests:
- Each project is exported before it is deleted; a failed export skips only that project.
- Empty-folder projects are sent with "delete files" on.
- The summary message lists what failed.
- The current project is replaced when it is deleted.

## Later (own plan)
- "Load project" accepts a zip of zips, and lets the user pick which project or projects to load.
