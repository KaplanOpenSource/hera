# Issue #1120 - connect a node to a key inside a dict parameter

## The bug

A parameter whose value is a dict gets no connection line. The same reference
written as a plain string works.

Reproduce with the template `Dict Parameter Bug` (templates menu).

Two causes:

1. `workflowDataflow.ts` / `buildDataflowEdges` reads only top-level string
   values. `if (typeof value !== 'string') return;` skips the dict.
2. `WorkflowNodeInputs.tsx` puts a handle only on top-level rows
   (`parentKey === input_parameters`). A nested row has no dot to land on.

Values inside a list have the same problem. A reference in a `Command` list
item draws no line either. The fix covers both.

## The idea

Replace the "parameter name" in every dataflow path with a "parameter path".
A path is the dotted keys from `input_parameters` down to the leaf, e.g.
`Parameters.project_name` or `Command.0`. A top-level parameter is a one
segment path, so today's behaviour is the simple case of the new one.

## Steps

### 1. New file `paramPath.ts`

Pure helpers on a node's `input_parameters`, no React.

- `paramPathOf(segments: string[]): string` - joins with `.`.
- `segmentsOf(path: string): string[]` - splits on `.`.
- `valueAtPath(params, path): any`.
- `withValueAtPath(params, path, value): params` - immutable, copies each
  level, keeps lists as lists.
- `withoutPath(params, path): params` - deletes a key, or splices a list item.
- `leafStringsOf(value): { path, text }[]` - walks dicts and lists, returns
  every string leaf with its path.

New test file `tests/paramPath.test.ts`.

### 2. Edge ids and handle ids accept a path

- `inputHandleId(node, path)` - unchanged shape, `<node>:in:<path>`. The
  matcher `/:in:([^:]+)$/` already allows dots, so nothing breaks.
- `ReferenceKindRegistry.EDGE_ID`: change `->(\w+)\.(\w+)$` to
  `->(\w+)\.(.+)$`. The target node is still one word; the rest is the path.
- Check `Reference.edgeIdTo` needs no change.

Extend `tests/Reference.test.ts` and the registry tests with a dotted path.

### 3. Build edges from nested values

In `buildDataflowEdges`, replace the top-level string scan with
`leafStringsOf(params)`. For each `{ path, text }` parse the references as
today and use `path` where `param` was used.

One line per reference per leaf. The `seen` set already drops duplicates.

### 4. Write and clear a reference at a path

- `setInputReference(node, path, reference)` - use `withValueAtPath`.
- `clearInputReference(node, path, refNode, key)` - read the leaf with
  `valueAtPath`, strip the token, write it back.
- `nodeWithParamValue`, `nodeWithReferenceAt`, `nodeWithoutParam` in
  `workflowNodeEdits.ts` - same change, path instead of name.

`WorkflowCanvasEdits.ts` needs no logic change. Only the names of the
variables it passes through (`param` -> `path`).

Update `tests/workflowDataflow.test.ts` and `tests/workflowNodeEdits.test.ts`.

### 5. A handle on every nested row

In `WorkflowNodeInputs.tsx`, `renderBeforeName` currently tests
`parentKey === keyForDetailsViewItem(INPUT_PARAMETERS_KEY)`. Change it to
accept any row under `input_parameters`:

- A row belongs to the inputs tree when `parentKey` is the inputs key or
  starts with `inputs key + '/'`.
- The path is `parentKey` with the `input_parameters/` prefix cut, split on
  `/`, plus `itemKey`, joined with `.`.

Give the same treatment to the right side source dot (`renderAfterValue`), so
another node can read a nested value too.

Draw the left handle on leaf rows only. A dict row holds no value, so a line
to it means nothing. The row does not know if it is a leaf, so pass the value
down or test it in the caller. Decide this when writing the code.

Check what `itemKey` a list element gets (`DetailsViewSortableItem`). The plan
assumes the index as a string.

Update `tests/WorkflowNodeInputs.test.tsx`.

### 6. Right-click and inline autocomplete on nested rows

`onRowContextMenu` and `onValueCaret` have the same top-level-only test. Apply
the same path rule, so the field menu and the `{` autocomplete work inside a
dict. `onFieldContextMenu` and `onFieldInlineEdit` then carry a path, which
flows into `referenceOutput` from step 4.

This is a separate commit. Step 5 fixes the reported bug on its own.

## Out of scope

- Hermes itself. Only the Web UI draws these lines.
- Changing the reference text format. It stays `{node.section.key}`.

## Validation

Follow `TEST_UI.md` in order: `npx tsc --noEmit`, `npm run test`, then the
build. Then load the `Dict Parameter Bug` template and check that both keys of
the dict get a line, and that deleting a line clears only that key.
