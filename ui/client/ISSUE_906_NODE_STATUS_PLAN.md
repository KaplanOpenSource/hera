# Node run status outline on the workflow canvas (#906)

Show each canvas node's state in the current run.

- Green outline: the node finished well.
- Red outline: the node failed.
- No outline: the node has not run yet.
- Blue moving dashed outline: the node is running now.

## How

1. `src/components/workflow/nodeRunStatus.ts` reads the run output chunks and
   returns a status per node. It parses the `[luigi-event]` lines
   (START / SUCCESS / FAILURE / BROKEN) and maps the Luigi task name back to the
   node name. Last event wins. A task still running when the run ended counts as
   failed.
2. `WorkflowEditor` reads the run from `useWorkflowRunStore` by workflow name and
   passes the statuses to `WorkflowGraph`.
3. `WorkflowGraph` puts the status in each node's data.
4. `WorkflowFlowNode` draws the outline. The status outline is 10px, far thicker
   than the 1px selection border. The running state adds a moving dashed outline
   drawn with four gradient strips on a `::before` overlay.

## Tests

- `tests/nodeRunStatus.test.ts` for the parsing.
- `tests/WorkflowFlowNode.test.tsx` for the outline colors, width, and animation.
