# Workflow run states (#1113)

Written 2026-11-03. Status: **done** (2026-11-03).

Goal: the run button shows the workflow's state, with one action:

- Idle: Run.
- Running: Stop.
- Stopped / finished / error: Run again.

Clean up of target files is out of scope. It is a separate task (#1115).

## What the code already does

- The run is a forked child process (`workflow_runner.py`). Stop can just end it.
- One run at a time: the runner's lock plus its RUNNING status.
- Each run has a random token. The client polls with it.
- There is no stop today: no button, no server route.
- `prepareWorkflowRunFromDoc` deletes `<name>_targetFiles` at the start of every run.
  So a stopped or finished run can rerun without any clean up.

## Server

Three classes, one file each, in `ui/server/workflow_runner/`. All in memory only.

- `ActiveWorkflowRunOne` (`active_workflow_run_one.py`): one run.
  Holds doc id, project name, token, child process, pid, status and log.
  Has `start()` and `stop()`. The run logic moves here from `WorkflowRunner`.
- `ActiveWorkflowRunList` (`active_workflow_run_list.py`): holds the active runs,
  keyed by token. Finds by token, or by project + doc id when the client forgot the token.
- `WorkflowRunner` (`workflow_runner.py`): what the routes call: start, poll, stop.
  Keeps the one-run-at-a-time rule. Makes a run and puts it in the list.

The doc id is `doc["_id"]["$oid"]` from the doc the client already sends.
Projects have only a name, so the run keeps the project name with the doc id.

Steps:

1. Add `RunStatus.STOPPED`.
2. Add the two new classes and move the run logic out of `WorkflowRunner`.
3. Set the pid right after the child process starts.
4. `POST /stop_workflow/{token}`: terminate the child process.
5. The read loop waits on the queue forever. Make it wake when the child dies,
   so the lock is freed and the status becomes STOPPED.
6. Tests in `ui/server/tests/test_workflow_runner.py`: stop a slow run, check STOPPED,
   check a new run can start.

## Client

1. Add `Stopped` to `WorkflowRunStatus` in `useWorkflowRunStore.ts`.
2. `workflowApi.stop(token)`.
3. `RunWorkflowButton`: Run when idle or ended, Stop while running.
4. The poller treats STOPPED as an end state.
5. Unit tests for the button and the store.

## Later

- Recover a run after a browser reload: a route that finds the run by doc id.
- Recover a run after a server restart.
