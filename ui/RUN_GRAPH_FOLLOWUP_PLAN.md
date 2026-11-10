# Run graph follow-up plan (#1114)

One small change on top of `ui/RUN_GRAPH_PLAN.md`.

## A Restore button for a past run's graph

Size: small.

- The read-only bar gets a "Restore" button, next to "Back to latest".
- Click: put the past run's `desc.workflow` into the shown doc, then go back to the latest.
- Only `desc.workflow` is restored, not the whole doc. The `_id`, name and other fields
  stay as today. `parameters` gets rebuilt from the workflow by the normal save step.
- It is not saved yet. The user sees the graph, then saves with "Update Document" as usual.
- If the shown doc has unsaved edits, Restore throws them away. Advised: ask first in that case.
- Note: `ui/RUN_GRAPH_PLAN.md` said restore comes later. This plan is that later.
- Tests: Restore puts the old workflow in the shown doc and leaves read-only mode.

## Where the code is

- Store: `ui/client/src/stores/useWorkflowRunStore.ts`.
- Bar: `ui/client/src/components/workflow/WorkflowPastRunBanner.tsx`.
- Canvas: `ui/client/src/components/workflow/WorkflowEditor.tsx`, `WorkflowGraph.tsx`.
- Shown doc (where Restore writes): `ui/client/src/components/workflow/WorkflowCanvasPanel.tsx`.
