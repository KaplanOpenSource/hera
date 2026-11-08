# Run graph plan (#1114 follow-up)

Keep the workflow graph each run used, next to its log.
Builds on `ui/RUN_HISTORY_PLAN.md`.

## Why

- A workflow can change after a run.
- An old log shown on today's graph can miss or mix up nodes.
- With the graph saved, an old run always shows the graph it ran.

## Part 1: save the doc with the run

- Save the whole workflow doc the run used, as one more field (`doc`) in the run JSON file.
- All fields, not only `desc.workflow`: `_id`, `desc`, `resource`, type and the rest.
- It is the doc the client sent at start, so it holds unsaved edits too, the real thing that ran.
- The server already holds it (`ActiveWorkflowRunOne.doc`), so this is one more field in `save_history`.
- The run list stays small: the list route still leaves out `doc`, the same as `chunks`.
- The read route answers `doc` with the rest of the run.
- Old run files without `doc` still load; they just have no graph.
- Tests: a saved run holds the doc as it was sent; the list leaves it out.

## Part 2: show the saved graph

- When an old run is picked, the workflow's canvas shows the saved graph in place of today's one.
- The saved graph is read-only. It is a record, not something to edit.
- The node frames come from the old run's log, on its own graph, so names always match.
- A clear mark says this is a past run, with its start time, and a way back to the latest.
- The client keeps the picked run's doc in the run store, next to its chunks.
- Old runs with no saved doc fall back to today's graph, as now.

## Decisions

- The old graph shows on the canvas itself, read-only. No separate tab.
- No marks for how the old graph differs from today's.
- Restoring an old graph as the current workflow comes later, not in this plan.

## Where the code is

- Server save, list, read: `ui/server/workflow_runner/workflow_run_history.py`.
- Server run and its sent doc: `ui/server/workflow_runner/active_workflow_run_one.py` (`self.doc`, `save_history`).
- Routes and models: `ui/server/routers/workflow.py`, `ui/server/api_models.py` (`WorkflowRunRecord`).
- Client calls: `ui/client/src/api/WorkflowApi.ts` (`runs`, `run`).
- Client picked run: `ui/client/src/stores/useWorkflowRunStore.ts` (`oldRuns`, `shownRun`).
- Client picker: `ui/client/src/components/workflow/log/WorkflowRunSelect.tsx`.
- Canvas and node frames: `ui/client/src/components/workflow/WorkflowEditor.tsx`.
- Tests: `ui/server/tests/test_workflow_run_history.py`, `ui/client/tests/workflowRunSelect.test.tsx`.
