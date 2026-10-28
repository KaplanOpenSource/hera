# Issue 1064 - References into a node output sub-field

A node like `RunPythonCode` can return a dict or a list. A later node should be
able to point at one value inside it, not only at the whole output.

Today `{A.output.result.station}` is dropped by the UI: the parser takes the
last dotted part as the key and `output.result` as the section, which matches no
kind. So no reference and no line on the canvas.

## What the server already does

Hermes resolves these tokens itself, so the UI only has to stop dropping them.

- `Hermes/hermes/taskwrapper/wrapper.py` `parsePath` treats everything inside
  `{}` as one path, and only splits on `.` to read the node name.
- `Hermes/hermes/engines/luigi/taskUtils.py` `_handle_output` passes everything
  after `output.` straight to `jsonpath_rw_ext.match`.

So the sub-path is JSONPath, not our dotted param path:

```
{A.output.result}             # whole output, works today
{A.output.result.station}     # dict field
{A.output.items[0]}           # list element
{A.output.items[0].name}      # field inside a list element
{A.output.items[*].name}      # every element - match returns a list
```

Note `items[0]`, not `items.0`. Our input side writes `Command.0`, but that is a
different thing (a path inside this node's own params) and does not change.

## 1. Parse a JSONPath key

`references/ReferenceKindRegistry.ts`

The section can itself hold dots (`InputReferenceKind.section` is
`Execution.input_parameters`), so the split cannot stay "last dot wins".

- Change `PARSE` to capture the node name and the whole rest:
  `\{\s*(\w+)\.([\w.\[\]*'"-]+)\s*\}`. The rest must allow brackets, `*` and
  quotes, since JSONPath also writes `items[?(@.x)]`-style filters. Keep the
  class tight enough that plain prose in a value is not read as a reference.
- Add a helper that finds the kind whose `sectionMatch()` matches a left-anchored
  prefix of the rest, and returns that kind plus the remainder as the key.
  Longest matching section wins, so `Execution.input_parameters.x` is not read as
  a shorter section.
- `bySection` stays as it is for the inline autocomplete, which asks about the
  section alone.
- Widen the key group in `HANDLE` and `EDGE_ID` the same way. `HANDLE` already
  ends in `(.+)`; `EDGE_ID` needs the change.

Edge id risk: `df:<node>:<mark>:<key>-><target>.<paramPath>` - the key now holds
dots and brackets. `->` stays the only separator the regex splits on.

## 2. Validate by the root segment

`workflowDataflow.ts`, `buildDataflowEdges`

`isReal` today asks `namesOf(...).includes(reference.key)`, which a sub-path key
never satisfies.

- Add `Reference.rootKey()`: the key up to the first `.` or `[`.
- Compare `rootKey()` against `namesOf`, and keep the full key in the handle id,
  so there is one line per distinct sub-field.

`WorkflowNodeOutputRow.tsx` renders one dot per catalog output name, so lines for
several sub-fields of the same output share that one dot. Fine for a first step.

## 3. Clearing a reference

`workflowDataflow.ts`, `clearInputReference`

- `ReferenceKind.clearToken` escapes the key for the regex, and the escape list
  already covers `[`, `]` and `*`. No change beyond passing the full key.
- Confirm the delete path uses the key from `parseDataflowEdgeId`, which after
  step 1 gives the full key.

## 4. The menus stay as they are

No help is possible after the output name, and none is needed.

`Hermes/hermes/utils/node_lookup.py` finds outputs by reading the executer
source with `ast` and taking the top-level keys of the dict `run()` returns.
There is nothing deeper, and for `RunPythonCode` the shape comes from the user's
own code at run time, so it cannot be known statically at all. The run result
does not help either: the server sends back captured console text per task
(`ui/server/api_models.py:35`), not the returned values.

- The right-click menu and the inline menu keep listing top-level output names.
- Once a name is picked, everything after it is the user's typing.
- `WorkflowReferences.inlineOptions` key stage: filter by the last segment of the
  seed, not the whole seed, so a seed like `result.st` does not make the menu
  come back empty. Nothing new is offered - this only stops the menu from
  fighting the user while they type a sub-path.

## 5. Tests

`ui/client/tests/`

- Parse: dict field, `items[0]`, `items[0].name`, `items[*].name`, and an
  `Execution.input_parameters.x` input reference still read as before.
- Parse rejects: a value with braces that is not a reference stays untouched.
- `buildDataflowEdges`: a line is built for `{A.output.result.station}` when the
  catalog knows only `result`.
- Edge id round trip: build an id with a bracketed key and parse it back.
- `clearInputReference` with a bracketed key removes only that token.
- `inlineOptions` filtering at the key stage with a dotted seed.

## Out of scope

- Any sub-field menu. The shape is not known anywhere.
- Validating the JSONPath itself. A wrong path fails at run time, as today.
- Server side. No Python changes - Hermes already resolves these.

## Order

1. Parse and ids (step 1), with tests.
2. Validation by root segment (step 2), with tests.
3. Clear and edge round trip (step 3).
4. Inline filter tweak (step 4).
5. Full `TEST_UI.md` checklist before reporting done.
