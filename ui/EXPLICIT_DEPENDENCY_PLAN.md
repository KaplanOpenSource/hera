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

## 3. A reference also fills `requires` - DONE

A field that references another node's output should add that node to `requires`.

- `nodeWithRequire` / `nodeWithRequiresFromText` in `workflowNodeEdits.ts` add
  the link. The reference and the link go in one write, or the second write
  would read the node as it was and undo the first.
- Covered: the output-to-input drag, the right-click "reference output", and a
  typed or picked `{other.output}` in a field.
- Add only. Removing a reference leaves `requires` alone, since the link may
  also have been drawn by hand. Delete it with the edge's X.

## Tests

- `WorkflowNodeSummary` / `WorkflowFlowNode`: the label shows and hides.
- `WorkflowDisplayEdges` or an edge test: labels appear on hover.
- `WorkflowCanvasEdits`: a dataflow connect adds `requires`; removing the last
  reference drops it.
