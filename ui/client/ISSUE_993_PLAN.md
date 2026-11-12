# Issue 993 - friendlier hermes editor titles (points 1 and 2)

Goal: hide hermes internals in the node editor. Two points only.

## Point 1 - "parameters" label

Today the node's parameter tree has a top row named `input_parameters` with a
type chip ("object" dropdown) next to it.

Plan:
- Show the text `parameters` on that row instead of `input_parameters`.
- Hide the type chip on that row. Its type is not the user's choice.
- Keep the chevron, add button, and JSON button as they are.

Touches: `WorkflowNodeInputs.tsx`, `DetailsViewItem.tsx` (new flag to hide the
type chip).

## Point 2 - friendly parameter names

Each parameter row shows its hermes name. Show a readable name instead, and the
hermes name only while the name is being edited.

Rules: `ProjectName` -> `Project Name`, `foo_bar` -> `Foo Bar`, runs of capitals
stay together (`URL` -> `URL`).

Plan:
- New helper `friendlyParamName.ts` in `src/components/workflow/`.
- Pass a "name to show" callback into the details tree. `RenameField` already
  supports a view-only label, so editing keeps showing the real name.
- Apply it only to top-level parameter rows of a workflow node.

Touches: new `friendlyParamName.ts`, `WorkflowNodeInputs.tsx`,
`DetailsViewItem.tsx` and the three components that forward tree props
(`DetailsViewItemsInObject.tsx`, `DetailsViewItemsInArray.tsx`,
`DetailsViewSortableItem.tsx`).

## Not in this change

- Point 3 (click to edit a field name, delete button on hover).
- Point 4 (drop the type chip, infer the type from what is typed).

## Tests

- `tests/friendlyParamName.test.ts` - name conversion cases.
- `tests/WorkflowFlowNode.test.tsx` - row says `parameters`, a parameter shows
  the friendly name, clicking it shows the hermes name.
- Then the full `TEST_UI.md` checklist.
