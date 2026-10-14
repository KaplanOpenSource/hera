# Explicit node dependency (issue #1002) - plan

Three open points. Each is independent.

## 1. Requires label on the node - DONE

Show the nodes a node depends on, in small faint text on the node card.

- Read `requires` with `normalizeRequires`.
- Draw it in `WorkflowNodeSummary` under the node name, e.g. "after: nodeA, nodeB".
- Hide the line when `requires` is empty.
- Use `text.secondary` and a small font so it stays quiet.

## 2. Direction labels on the edge - DONE

Say which end is which, so the arrow is not the only hint.

- In `WorkflowRequiresEdge`, add two small labels with `EdgeLabelRenderer`:
  "depends on" near the target end, "required by" near the source end.
- Show them only on hover, like the delete X, so the canvas stays clean.
- Take the end points from the props the edge already gets.

## 3. A reference also fills `requires` - open

A field that references another node's output should add that node to `requires`.

- In `WorkflowCanvasEdits.connect`, after writing the reference, call
  `onAddRequire(source, target)`.
- Do the same for the right-click "reference output" command.
- Skip it when `requires` already has the node.
- On removing the last reference to that node, drop it from `requires` only if
  no other field still references it.
- Open question: should a hand-typed reference do this too? It would need a
  check when the field value changes.

## Tests

- `WorkflowNodeSummary` / `WorkflowFlowNode`: the label shows and hides.
- `WorkflowDisplayEdges` or an edge test: labels appear on hover.
- `WorkflowCanvasEdits`: a dataflow connect adds `requires`; removing the last
  reference drops it.
