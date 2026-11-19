# Plan: export a project from the delete dialog

## Goal
- The delete dialog can export the project as a zip before it deletes it.
- If the export works and the files stay on disk, the user does not need to type the project name.

## Steps
1. Move `saveFile` out of `ExportProjectButton.tsx` into its own small file, so both buttons can use it.
2. In `DeleteProjectButton.tsx`, add a checkbox: "Export the project as a zip first".
3. The name field is needed only when export is off, or when "delete files" is on.
   Hide it (or mark it optional) in the other case.
4. On confirm:
   - If export is on, call `projectApi.export` and save the zip.
   - If the export fails, stop. Do not delete.
     Show an error with `pushError`: "Export failed, so the project was not deleted. Try again."
   - Then call `projectApi.delete` as today.
5. Under the export checkbox, add a short help text on how to restore:
   "You can restore the project from this zip later. Open the project menu, click 'Load project', and pick the zip."
6. No server change. The export and delete routes already exist.

## Tests
- `tests/deleteProject.test.tsx`:
  - Export on, files kept: deletes with no name typed, and saves the zip.
  - Export fails: nothing is deleted, and the "not deleted" error shows.
  - Export on, files deleted: still needs the name.
  - Export off: same as today.
