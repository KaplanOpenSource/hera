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

- `<FilesDirectory>/workflow_runs/<workflow name>/<start time>.json`.
- Not inside `<name>_targetFiles`: that folder is deleted at the start of every run.
- Every project has a files folder (`Project` falls back to `~/.hera/<projectName>`).
- A run with no files folder fails before it starts, so there is nothing to save.

## When

- Write the file once, when the run ends (done, error or stopped).
- The one place every run ends is `ActiveWorkflowRunOne._background`.
- While the run goes, the live output stays in memory, as today.
- If the server crashes mid-run, that run is lost. Accepted.
- Saved runs are never deleted.

## Steps

1. Server: folder methods on `ActiveWorkflowRunOne`:
   - `files_directory()`: the project's `FilesDirectory`, read through the workflow toolkit.
   - `target_files_folder()`: `<FilesDirectory>/<workflow name>_targetFiles`.
   - `run_history_folder()`: `<FilesDirectory>/workflow_runs/<workflow name>`.
   - Store `workflow_name` from `doc["desc"]["workflowName"]` in `__init__`.
2. Server: record start and end time on `ActiveWorkflowRunOne`.
3. Server: write the run JSON file at the end of `_background`.
4. Server: an endpoint to list the saved runs of a workflow (time, status).
5. Server: an endpoint to read one saved run.
6. Client: a simple select box in the output tab to pick a run; latest run by default.
7. Client: when an old run is picked, the canvas frames show that run.
8. Tests for each step.
