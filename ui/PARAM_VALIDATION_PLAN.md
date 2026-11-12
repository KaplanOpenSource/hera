# Node parameter validation in the UI (issue #996)

Goal: when the user edits a Hermes node parameter in the workflow editor, show
Hermes's own error message on the node.

Autocomplete is not in this plan. `getValuesForParam` is dead code - one stub in
`abstractExecuter.py`, no overrides, no callers. Tracked in a comment on #996.

Example node: `RiskAssessment.calculateThresholds`. `Agent` must be one of the
project's agents, `Calculator` one of that agent's effect names.

## Current state

- `testParamValues(params)` exists on ~14 node executers under
  `Hermes/hermes/Resources/`. Returns `(ok, message)`. Nothing calls it.
- The UI fetches `/node-catalog` once. The catalog is a static source parse
  (`Hermes/hermes/utils/node_lookup.py`) - names only, no values, no hera import.
- The UI never validates a parameter value.

## Step 1 - one Hermes fix

`checkParamAgainstList` falls off the end and returns `None` when the value is
valid, so the caller's two-value unpack raises TypeError on every good value.
Return `(True, "")`. This is the only Hermes change needed.

Two other known problems are worked around, not fixed:

- The base `testParamValues` stub is `@staticmethod` but takes `self`. The server
  only calls the method on nodes that override it (`'testParamValues' in
  cls.__dict__`); the rest report "no validation".
- The result carries no parameter name, only a message for the first failure. So
  the message is shown on the node, not on a field. Per-field errors come later.

## Step 2 - server endpoint

One new endpoint in `ui/server/server.py`, next to `/node-catalog`:

`POST /node-params/validate` - body `{ type, params }`, returns
`{ ok, message }`, or `{ ok: true }` when the node has no `testParamValues`.

It imports the node's executer class by its dotted type (the type is the
Resources directory path, and the class name is the last segment). Unlike the
catalog this needs hera and a live project, so it runs in the server process that
already warms hera up. Not cacheable - the answer depends on project data.

## Step 3 - client

- A store like `useNodeCatalog`, but per node and keyed on the node's current
  params. Debounce; the call hits MongoDB.
- Show the message on the node in `WorkflowFlowNode.tsx`, as a soft warning next
  to the existing `nodeTypeIssue` warning. Never block editing.

## Out of scope

- Autocomplete / `getValuesForParam`.
- Per-field error highlighting.
- References (`{otherNode.output}`), which `isParamTestable` skips. Presence of a
  required reference is still checked.

## Order of work

Fix `checkParamAgainstList`, then `calculateThresholds` end to end through steps
2 and 3, then confirm the other nodes report through the same path.
