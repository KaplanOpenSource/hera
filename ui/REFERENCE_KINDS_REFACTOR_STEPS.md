# Low-level steps: the reference-kind refactor only

Scope: steps 1-3 of `REFERENCE_KINDS_PLAN.md`. One class owns what a reference is.
No new kind, no new dot, no UI change. Behaviour stays the same except the
dataflow edge id format.

Status: steps 1-3 are done. Next is step 4 of `REFERENCE_KINDS_PLAN.md`
(`InputReferenceKind`).

## Step 1 - new folder `src/components/workflow/references/`

Nothing else changes in this step. Nothing imports it yet.

### `ReferenceKind` (abstract)

```ts
export abstract class ReferenceKind {
  abstract readonly section: string;      // written form, e.g. 'output'
  abstract readonly label: string;        // menu label, e.g. 'Output'
  abstract readonly handleMark: string;   // handle id fragment, e.g. 'out'
  abstract namesOf(node: WorkflowNode, catalog: NodeCatalogEntry[]): string[];

  format(node: string, key: string): string;      // `{${node}.${this.section}.${key}}`
  handleId(node: string, key: string): string;     // `${node}:${this.handleMark}:${key}`
  clearToken(node: string, key: string): RegExp;   // `\{\s*node\.(sectionMatch)\.key\s*\}` g
  sectionMatch(): string;                          // regex source, default: escaped section
  couldBe(typedSection: string): boolean;          // default: section.startsWith(typedSection)
}
```

`handleMark` is one word, no dots and no colons - handle ids and edge ids split on
those.

### `OutputReferenceKind`

Only subclass in this step.

- `section = 'output'`, `label = 'Output'`, `handleMark = 'out'`
- `namesOf` = `nodeOutputNames(node, catalog)`
- `sectionMatch()` = `'parameters?|outputs?'` (keeps reading the older spellings)

### `Reference` (value)

```ts
export class Reference {
  constructor(readonly node: string, readonly kind: ReferenceKind, readonly key: string) {}
  toString(): string            // kind.format(node, key)
  handleId(): string            // kind.handleId(node, key)
  edgeIdTo(target: string, param: string): string
  clearToken(): RegExp          // kind.clearToken(node, key)
}
```

Edge id: `df:<node>:<handleMark>:<key>-><target>.<param>`. Colons, not dots, so a
future section with a dot in it still parses. This is the one visible format
change; edge ids are not persisted anywhere.

### `ReferenceKindRegistry`

Built once in the constructor: `PARSE = /\{\s*(\w+)\.([\w.]+)\.(\w+)\s*\}/g` plus a
section -> kind lookup that tries each kind's `sectionMatch()` anchored. An
unknown section matches no kind and is skipped, same as today.

```ts
all(): ReferenceKind[]
bySection(section: string): ReferenceKind | null
parseAll(value: string): Reference[]                 // every reference in a value
ofHandle(handleId: string): Reference | null         // `${node}:${mark}:${key}`
ofEdgeId(id: string): { reference: Reference, target: string, param: string } | null
couldBe(typedSection: string): ReferenceKind[]
keysOf(node: string, workflowNode: WorkflowNode, catalog: NodeCatalogEntry[]): Reference[]
```

`export const knownKinds = new ReferenceKindRegistry([new OutputReferenceKind()]);`

New test files `tests/Reference.test.ts` and `tests/ReferenceKindRegistry.test.ts`: format/parse round trip,
handle id round trip, edge id round trip, `parseAll` with two refs in one value,
old `parameters` spelling still parsed, unknown section ignored.

## Step 2 - `workflowDataflow.ts` delegates

Keep every exported name and signature it has now, so `WorkflowCanvasEdits`,
`workflowNodeEdits`, `WorkflowNodeOutputRow` and the existing tests do not move
in this step. Bodies only:

| Export | New body |
|---|---|
| `outputHandleId(node, name)` | `new Reference(node, OUTPUT, name).handleId()` |
| `dataflowReference(node, name)` | `new Reference(node, OUTPUT, name).toString()` |
| `parseDataflowConnection` | source via `knownKinds.ofHandle`, target via `INPUT_HANDLE_MATCH`; add `kind` to the returned object |
| `parseDataflowEdgeId(id)` | `knownKinds.ofEdgeId(id)`, mapped to the existing `{ refNode, key, target, param }` shape |
| `clearInputReference` | `new Reference(refNode, OUTPUT, key).clearToken()` |
| `buildDataflowEdges` | `knownKinds.parseAll(value)` per param; keep the "node in graph and key is a real output" filter via `ref.kind.namesOf(...)`; edge id from `ref.edgeIdTo(target, param)` |

Delete the local `REFERENCE`, `OUTPUT_HANDLE_MATCH` and the inline edge-id regex.
`inputHandleId` / `INPUT_HANDLE_MATCH` / `nodeOutputHandleId` / `nodeInputHandleId`
stay untouched - the target dot is not per-kind.

`tokenAtCaret` stays as is here. It already returns the raw section text; it gains
nothing until a second kind exists. Change only `ReferenceTokenStage.Output` ->
keep the name (renaming it is churn with no reader).

Tests: `tests/workflowDataflow.test.ts` and `workflowReferenceRoundTrip.test.ts`
should pass unchanged except any assertion on the literal edge id string. Grep
`'df:'` in tests and fix those.

## Step 3 - the menus take a `Reference`

- `WorkflowReferences.outputsOf(name)` -> `referencesOf(name): Reference[]` using
  `knownKinds.keysOf`. `optionsFor(name)` returns `{ node: string, references: Reference[] }[]`,
  still dropping self and empty nodes.
- `inlineOptions` keeps returning `string[]` (node names, then keys) - it feeds an
  autocomplete of plain strings. Its Output branch reads keys from `referencesOf`.
- `WorkflowContextMenu`: rename `NodeOutputOption` to `ReferenceOption` with
  `{ node: string, references: Reference[] }`. Second submenu options are the node's
  `Reference`s, labelled by their `key`. `onReferenceOutput(node, param, reference, caret)`
  - one `Reference` argument instead of `(sourceNode, output)`.
- `WorkflowCanvasEdits.referenceOutput(nodeName, param, reference, caret?)` and
  `nodeWithReferenceAt(node, param, reference, caret?)`; `insertReferenceAt` and
  `replaceReferenceAt` take a `Reference` too.
- `setInputReference(node, param, reference)` likewise; `WorkflowCanvasEdits.connect`
  builds `new Reference(connection.source, dataflow.kind, dataflow.key)`.
- `applyInlinePick`: the Node branch writes the scaffold from the kind
  (`` `{${option}.${OUTPUT.section}.` ``, no literal `.output.`); the Output branch
  builds a `Reference`.
- `WorkflowGraph.tsx`: `referenceOptions` type follows `optionsFor`; pass the
  `Reference` straight through to `edits.referenceOutput`.
- `useInlineReference`: `references.referencesOf(option).map(r => r.key)`.

Tests to update: `WorkflowReferences.test.ts`, `WorkflowContextMenu.test.tsx`,
`inlineReferencePick.test.ts`, `workflowReferenceRoundTrip.test.ts` - argument
shapes only, no new expectations.

## Done when

- No file outside `references/` contains the literal `output` as part
  of a reference or handle format. Check with
  `grep -rn "\.output\.\|:out:" src/components/workflow`.
- `npx tsc --noEmit` clean, `npm run test` green.
- Then `ui/client/TEST_UI.md` in order.
- Manual check: drag an output dot to a param, see the token appear, delete the
  line, see the token go.
