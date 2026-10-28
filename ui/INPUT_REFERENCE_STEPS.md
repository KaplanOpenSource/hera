# Low-level steps: reference another node's input parameter

Scope: the six steps in `INPUT_REFERENCE_PLAN.md`. Each step builds and passes tests on its own,
except that step 1 must not ship without step 2 (see the warning there).

## Step 1 - `references/InputReferenceKind.ts`

New file, one class, nothing imports it yet.

- `section = 'Execution.input_parameters'`, `label = 'Input'`, `handleMark = 'param'`.
- `handleMark` may not be `in`: an input row's left dot already uses `<node>:in:<param>`, and
  `ReferenceKindRegistry.ofHandle` would read those as sources.
- `namesOf` returns the node's own top-level parameter keys:
  `Object.keys(node.Execution?.input_parameters ?? {}).filter(isParameterKey)` from `nodeCatalog`.
  Not the catalog's parameter list - the dot only exists on rows the node actually has.
- No `sectionMatch` override. The default escapes the dot, and `PARSE` already takes a dotted
  section, so `{A.Execution.input_parameters.k}` parses with key `k`.

New test file `tests/InputReferenceKind.test.ts`: the written token, the handle id, the keys of a
node with two params, and nothing for a node with none.

**Do not add it to `knownKinds` yet.** The moment it is registered, `buildDataflowEdges` starts
drawing lines whose source handle is `<node>:param:<key>`, and no such handle is rendered until
step 2. Register it at the end of step 2.

## Step 2 - the right dot on an input row

`DetailsViewItem` has `renderBeforeName` but no slot after the value editor.

- Add `renderAfterValue?: (itemKey: string, parentKey: string | undefined, def?: FieldDef) => ReactNode`
  next to `renderBeforeName`, rendered as the last child of the row `Stack`. Thread it through
  `DetailsViewItemsInObject` and `DetailsViewItemsInArray` exactly as `renderBeforeName` is threaded.
- `WorkflowNodeInputs` passes one. On a top-level row
  (`parentKey === keyForDetailsViewItem(INPUT_PARAMETERS_KEY)`) it renders
  `<Handle type="source" id={INPUT.handleId(nodeName, itemKey)} position={Position.Right} />`,
  styled like the dot in `WorkflowNodeOutputRow` (8x8, `primary.main`, `right: -14`,
  `position: relative`). Other rows get nothing.
- Then `references/knownKinds.ts`: `export const INPUT = new InputReferenceKind();` and
  `new ReferenceKindRegistry([OUTPUT, INPUT])`.

Tests: `tests/WorkflowNodeInputs.test.tsx` - a source handle on a top-level param row, none on a
nested one. Extend `tests/ReferenceKindRegistry.test.ts` for the two registered kinds:
`bySection('Execution.input_parameters')`, `couldBe('Exec')`, and `keysOf` returning both kinds.

## Step 3 - a node may not point at itself

- `WorkflowCanvasEdits.canConnect`: a dataflow connection returns
  `connection.source !== connection.target` instead of a bare `true`.
- `buildDataflowEdges`: skip a reference when `reference.node === target`.

Tests in `tests/WorkflowCanvasEdits.test.ts` and `tests/workflowDataflow.test.ts`.

## Step 4 - the typing path learns the section

This is the real work. Today `tokenAtCaret` splits on the first dot and calls everything after it the
output key, so a dotted section cannot be typed at all.

In `workflowDataflow.ts`:

- `ReferenceTokenStage` becomes `Node`, `Section`, `Key`.
- `ReferenceTokenAtCaret` gains `sectionPart: string` and `kind: ReferenceKind | null`.
- `tokenAtCaret`: with no dot it is the `Node` stage, as now. Otherwise take the text after the first
  dot, split it at its *last* dot, and ask `knownKinds.bySection` about the left half. A hit is the
  `Key` stage with that kind and the right half as the seed. A miss is the `Section` stage with the
  whole of it as `sectionPart`.

In `WorkflowReferences.inlineOptions`: `Node` returns node names, as now. `Section` returns
`knownKinds.couldBe(sectionPart).map(kind => kind.label)`. `Key` returns that one kind's keys for
`nodePart`, filtered by the seed.

In `inlineReferencePick.applyInlinePick`: `Node` writes `{option.` and stays open - it no longer
writes a section, because there is more than one. `Section` looks the label up to a kind and writes
`{node.section.`. `Key` builds a `Reference` and completes.

In `useInlineReference.pickInline`: after an incomplete pick, recompute with
`references.inlineOptions(node, picked.value, picked.caret)` rather than assuming the next options
are that node's keys.

Tests: `tests/workflowDataflow.test.ts` (the three stages, including a half-typed dotted section),
`tests/inlineReferencePick.test.ts` (node -> kind -> key writes the right token),
`tests/WorkflowReferences.test.ts` (the options at each stage).

## Step 5 - the right-click menu shows both kinds

`optionsFor` already returns every kind once step 2 registers `INPUT`, so only the display changes.

- The second fly-out's `Autocomplete` gets `groupBy={reference => reference.kind.label}`.
- Its `TextField` label becomes `Parameter`; the root menu item becomes `Reference another node`.

Update `tests/WorkflowContextMenu.test.tsx` for the new wording and for picking an input option.

## Step 6 - validation

- `npx tsc --noEmit` clean, `npm run test` green.
- Then `ui/client/TEST_UI.md` in order.
- `grep -rn "Execution.input_parameters" src/components/workflow` shows only
  `references/InputReferenceKind.ts`.
- By hand: drag an input's right dot onto another node's param, see the token appear, delete the
  line, see the token go. Then type `{A.Exec` in a param and finish it from the menu.
