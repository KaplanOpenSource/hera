# Run history plan (#1114)

Keep earlier workflow runs, so the user can look back at their output.

## Today

- Server: `ActiveWorkflowRunList.add()` drops the old run of the same doc. Memory only.
- Client: `useWorkflowRunStore` keeps one run per workflow name. Memory only.

## What a run saves

One JSON file per run:

```json
{
  "startTime": "2026-11-03T10:15:00",
  "endTime": "2026-11-03T10:17:42",
  "status": "done",
  "error": "",
  "chunks": [{ "name": "task1", "text": "..." }]
}
```

- Chunks are kept as they are. No flat log, no parsing.
- Status and error show how the run ended.
- Start and end time tell the runs apart.

## Where

- `<FilesDirectory>/workflow_runs/<workflow name>/<YYYY-MM-DD_HH-MM-SS>.json`.
- Not inside `<name>_targetFiles`: that folder is deleted at the start of every run.
- Every project has a files folder (`Project` falls back to `~/.hera/<projectName>`).
- A run with no files folder fails before it starts, so there is nothing to save.

## When

- Write the file once, when the run ends (done, error or stopped).
- Write it before the status changes, so a client that sees the end finds the file.
- The one place every run ends is `ActiveWorkflowRunOne._background`.
- While the run goes, the live output stays in memory, as today.
- If the server crashes mid-run, that run is lost. Accepted.
- Saved runs are never deleted.

## Steps

All steps are done.

1. Server: `WorkflowRunHistory` (`ui/server/workflow_runner/workflow_run_history.py`) owns the folder, save, list and read.
2. Server: folder methods on `ActiveWorkflowRunOne`: `files_directory()`, `target_files_folder()`, `run_history_folder()`.
3. Server: `ActiveWorkflowRunOne` records start and end time, and saves the run at the end of `_background`.
4. Server: `GET /workflow_runs?projectName&workflowName` lists saved runs, newest first, without output.
5. Server: `GET /workflow_run?projectName&workflowName&runId` reads one saved run.
6. Client: `WorkflowRunSelect`, a menu button at the right end of the output tab's bottom bar. Its menu opens upward and lists "Latest run" plus older saved runs. With no run in this session, the newest saved run shows.
7. Client: the picked old run lives in `useWorkflowRunStore.oldRuns`; `shownRun()` feeds both the output tab and the canvas frames. A new run drops the pick.
8. Tests: `ui/server/tests/test_workflow_run_history.py`, `ui/client/tests/workflowRunSelect.test.tsx`, `ui/client/tests/workflowApi.test.ts`.
