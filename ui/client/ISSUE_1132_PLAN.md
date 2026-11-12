# Issue 1132 - Remove the detect notebooks button

Notebook detection should run on its own, not from a button.

## When it runs

- Project create: the server already detects in `/project/create`. No change.
- Project switch, and startup on a project: the client calls
  `/project/detectNotebooks` after the project is loaded.
- It does not block the UI. The project shows first; the call runs after it.
- If notebooks were added and the user is still on that project, the store
  gets the new project and a toast says how many. Otherwise nothing shows.

## Changes

1. New `src/io/detectNotebooks.ts`: `detectNotebooks(projectName)`, the
   logic taken from the button.
2. Call it (not awaited) from `fetchProjectData` in `src/io/FetchProjects.tsx`
   and from `fetchInit` in `src/io/fetchInit.ts`.
3. Delete `src/components/project/DetectNotebooksButton.tsx` and remove it
   from `ProjectActionsButton.tsx`.
4. Tests: replace `tests/detectNotebooks.test.tsx` with tests of the new
   function; add checks that switch and startup start detection.
5. Fix the comment in `tests/integration/integHelpers.tsx`.

The server route stays. No Python changes.
